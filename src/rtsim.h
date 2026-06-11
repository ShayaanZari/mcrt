#ifndef RTSIM_H
#define RTSIM_H

#include <iostream>
#include <random>
#include <chrono>
#include <cmath> 
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>

constexpr double PI = 3.14159265358979323846;

// Defines a three-dimensional vector structure for spatial coordinates and direction.
// Includes a dot product method utilized for geometric projections and boundary calculations.
struct Vec3 {
    double x, y, z;
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};


// Tracks the final/escaped weights of photons in a simulation batch/wave.
// Computes metrics to monitor convergence.
struct Accumulator {
    std::vector<double> power_sums;
    int sample_count = 0; // Number of samples accumulated thus far

    Accumulator(
        int num_moments = 4 // Tracks Σw, Σw², Σw³, Σw⁴ by default — sufficient for mean, variance, skewness, and excess kurtosis
    ) : power_sums(num_moments, 0.0) {}

    // Accumulates a weight sample, updating raw power sums
    void accumulate(double weight) {
        double weight_power = weight; // w^1, w^2, w^3..
        for (double& ps : power_sums) {
            ps += weight_power;
            weight_power *= weight;
        }
        sample_count++;
    }

    void reset() {
        std::fill(power_sums.begin(), power_sums.end(), 0.0);
        sample_count = 0;
    }

    // Merge another accumulator's data into this one (for cumulative tracking).
    void merge(const Accumulator& other) {
        for (size_t i = 0; i < power_sums.size(); ++i) {
            power_sums[i] += other.power_sums[i];
        }
        sample_count += other.sample_count;
    }

    double mean() const { 
        return sample_count > 0 ? power_sums[0] / sample_count : 0.0;
    }

    double variance() const {
        if (sample_count < 2) return 0.0; // need 2nd raw moment to compute variance
        double m = mean();
        return (power_sums[1] - sample_count * m * m) / (sample_count - 1); // Computational trick for computing variance
    }

    // Kish's Effective Sample Size
    double n_eff() const {
        if (power_sums[1] == 0) return 0.0;
        else return (power_sums[0] * power_sums[0]) / power_sums[1];
    }

    // --- More tentative metrics ---

    double variance_penalty() const {
        if (sample_count == 0 || power_sums[0] == 0.0) return 0.0;
        return (power_sums[1] / (power_sums[0] * power_sums[0])) - (1.0 / sample_count);
    }

    double cv() const {
        if (sample_count < 2 || power_sums[0] == 0) return 0.0;
        return std::sqrt(variance()) / mean();
    }

    // Standardized 3rd central moment.
    // Positive = right-skewed (a few high-weight photons dominate the sum).
    // Formula: (s3/n - 3*m*s2/n + 2*m^3) / sigma^3
    // Uses simple /n denominator; bias correction is negligible at n >= 10k.
    double skewness() const {
        if (sample_count < 3) return 0.0;
        double n  = sample_count;
        double m  = power_sums[0] / n;
        double mu3 = power_sums[2] / n
                   - 3.0 * m * (power_sums[1] / n)
                   + 2.0 * m * m * m;
        double sigma = std::sqrt(variance());
        if (sigma == 0.0) return 0.0;
        return mu3 / (sigma * sigma * sigma);
    }

    // Excess kurtosis (raw kurtosis minus 3, so a Gaussian scores 0).
    // Large positive values indicate rare high-weight photons are dominating the estimator —
    // a direct signal that the cell needs better importance sampling coverage.
    // Formula: (s4/n - 4*m*s3/n + 6*m^2*s2/n - 3*m^4) / sigma^4 - 3
    double excess_kurtosis() const {
        if (sample_count < 4) return 0.0;
        double n  = sample_count;
        double m  = power_sums[0] / n;
        double m2 = power_sums[1] / n;
        double m3 = power_sums[2] / n;
        double m4 = power_sums[3] / n;
        double mu4 = m4
                   - 4.0 * m * m3
                   + 6.0 * m * m * m2
                   - 3.0 * m * m * m * m;
        double sigma2 = variance();
        if (sigma2 == 0.0) return 0.0;
        return mu4 / (sigma2 * sigma2) - 3.0;
    }
};
// "Sum of powers" accumulator: power_sums[k] = Σw^(k+1)
// Raw moments: power_sums[k] / sample_count
// Requires weights in a bounded range (e.g. [0,1]) for numerical stability
// of the higher-order computational formulas.


// Represents a discrete radial shell within the spherical domain.
struct Cell {
    double p_phys; // Physical Probability
    double q_bias; // Biased Importance sampling weight
    Accumulator accumulator; // Localized statistical data
    double r_min, r_max;

    // Returns the initial weight assigned to a photon spawned in this cell.
    // Encapsulates the importance sampling ratio: p_phys / (q_bias * n_photons).
    double init_weight(int n_photons) const {
        return (q_bias > 0.0) ? p_phys / (q_bias * n_photons) : 0.0;
    }
};

// Holds invariant simulation parameters.
struct Config {
    int n_photons;  // Photon count
    double kappa, radius, albedo; // properties of medium: extinction coefficient, spherical radius, and scattering albedo.
    int max_scatters = 10000;
    double mixing_ratio = 0.0; // Defines how much of the target proposal is mixed with physical probability
    int max_waves = 200; // The number of photon waves to simulate; currently photons per wave is fixed.
};

// Manages the collection of radial cells, RNG, and global counting.
struct Environment {
    Config config;
    std::mt19937 rng;
    std::vector<Cell> cells;
    Accumulator accumulator;

    // Constructs the environment by partitioning the sphere into shells of equal radial thickness.
    // Initializes the physical and biased probabilities based on the volume fraction of each shell.
    Environment(Config cfg, int num_cells, unsigned int seed);

    // Samples a unit vector uniformly from a $4\pi$ steradian solid angle.
    // Utilizes inverse transform sampling on azimuthal and polar coordinates.
    Vec3 isotropic_direction();
    
    // Generates a random coordinate within a given cell's radial boundaries, accounting for r^2 scaling of spherical volume.
    Vec3 random_point_in_cell(int idx);
    
    // Simulates the transport of a single photon through the participating medium.
    // Handles Beer-Lambert law path sampling, boundary intersections, and isotropic scattering.
    double simulate_photon(Vec3 pos, double weight);

    // Executes a discrete simulation wave using the current biased probability distribution.
    // Updates CDF and tallies outcomes into cell and global stats.
    void run_wave();

    // Adjusts the importance sampling weights based on recorded statistical moments.
    // Intended to minimize variance by shifting sample density toward high-contribution cells.
    void adapt_distribution(double a1, double a2);

    // Prints summary tables of radial grid values
    void print_local_stats();
};

double compute_escape_fraction(
    Config config,
    int num_cells,
    unsigned int seed,
    double nsr_tolerance,
    bool verbose = false
);

#endif // RTSIM_H

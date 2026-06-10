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
// Computes the mean, variance, and effective sample size (N_eff) to monitor MC convergence.
struct WaveStats {
    double weight_sum = 0.0;    // 1st moment: sum of escaped weights
    double weight_sq_sum = 0.0; // 2nd moment: sum of squared escaped weights
    int n = 0;

    void add(double weight) {
        weight_sum += weight;
        weight_sq_sum += weight * weight;
        n++;
    }

    void reset() { weight_sum = 0.0; weight_sq_sum = 0.0; n = 0; }
    double mean() const { return n > 0 ? weight_sum / n : 0.0; }
    double variance() const {
        if (n < 2) return 0.0;
        double m = mean();
        return (weight_sq_sum - n * m * m) / (n - 1); // Bessel's correction
    }
    double n_eff() const {
        return (weight_sq_sum == 0) ? 0.0 : (weight_sum * weight_sum) / weight_sq_sum;
    }
};

// Represents a discrete radial shell within the spherical domain.
struct Cell {
    double p_phys; // Physical Probability
    double q_bias; // Biased Importance sampling weight
    WaveStats stats; // Localized statistical data
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
};

// Manages the collection of radial cells, RNG, and global counting.
struct Environment {
    Config config;
    std::mt19937 rng;
    std::vector<Cell> cells;
    WaveStats global_stats;

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
    int max_waves = 200,
    bool verbose = false
);

#endif // RTSIM_H

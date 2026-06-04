#include <iostream>
#include <random>
#include <chrono>
#include <cmath> 
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>

// g++ sim.cc -o sim
// ./sim 1000 4 0.4 0.05

constexpr double PI = 3.14159265358979323846;
double mixing_ratio;

// Defines a three-dimensional vector structure for spatial coordinates and direction.
// Includes a dot product method utilized for geometric projections and boundary calculations.
struct Vec3 {
    double x, y, z;
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};

// Accumulates statistical moments for a simulation batch/wave.
// Computes the mean, variance, and effective sample size (N_eff) to monitor MC convergence.
struct WaveStats {
    double sum_w = 0.0; // 1st Moment
    double sum_w2 = 0.0; // 2nd Moment
    int n = 0;

    void add(double w) {
        sum_w += w;
        sum_w2 += w * w;
        n++;
    }

    void reset() { sum_w = 0.0; sum_w2 = 0.0; n = 0; }
    double mean() const { return n > 0 ? sum_w / n : 0.0; }
    double variance() const {
        if (n < 2) return 0.0;
        double m = mean();
        return (sum_w2 - n * m * m) / (n - 1); // Bessel's correction
    }
    double n_eff() const {
        return (sum_w2 == 0) ? 0.0 : (sum_w * sum_w) / sum_w2;
    }
};

// Represents a discrete radial shell within the spherical domain.
struct Cell {
    double p_phys; // Physical Probability
    double q_bias; // Biased Importance sampling weight
    WaveStats stats; // Localized statistical data
    double r_min, r_max;
};

// Holds invariant simulation parameters.
struct Config {
    int n_photons;  // Photon count
    double kappa, radius, albedo; // properties of medium: extinction coefficient, spherical radius, and scattering albedo.
    int max_scatters = 10000;
};

// Manages the collection of radial cells, RNG, and global counting.
struct Environment {
    Config config;
    std::mt19937 rng;
    std::vector<Cell> cells;
    WaveStats global_stats;

    // Constructs the environment by partitioning the sphere into shells of equal radial thickness.
    // Initializes the physical and biased probabilities based on the volume fraction of each shell.
    Environment(Config cfg, int num_cells, unsigned int seed) 
        : config(cfg), rng(seed) {
        
        double dr = config.radius / num_cells;
        double total_vol = std::pow(config.radius, 3);

        for (int i = 0; i < num_cells; ++i) {
            double r0 = i * dr;
            double r1 = (i + 1) * dr;
            double vol = std::pow(r1, 3) - std::pow(r0, 3);
            
            // Initializing p and q as unbiased (p = q = vol / total_vol)
            cells.push_back({vol / total_vol, vol / total_vol, {}, r0, r1});
        }
    }

    // Samples a unit vector uniformly from a $4\pi$ steradian solid angle.
    Vec3 isotropic_direction() { // Utilizes inverse transform sampling on azimuthal and polar coordinates.
        std::uniform_real_distribution<double> d01(0, 1), dmu(-1, 1);
        double theta = d01(rng) * 2.0 * PI; 
        double mu = dmu(rng); 
        double s = std::sqrt(1.0 - mu * mu);
        return { s * std::cos(theta), s * std::sin(theta), mu };
    }

    // Generates a random coordinate within a given cell's radial boundaries, accounting for r^2 scaling of spherical volume.
    Vec3 random_point_in_cell(int idx) {
        std::uniform_real_distribution<double> d01(0, 1);
        double r = std::pow(std::pow(cells[idx].r_min, 3) + 
                   (std::pow(cells[idx].r_max, 3) - std::pow(cells[idx].r_min, 3)) * d01(rng), 1./3.);
        Vec3 dir = isotropic_direction();
        return { dir.x * r, dir.y * r, dir.z * r };
    }

    // Simulates the transport of a single photon through the participating medium.
    // Handles Beer-Lambert law path sampling, boundary intersections, and isotropic scattering.
    double simulate_photon(Vec3 pos, double weight) {
        int scatters = 0;
        Vec3 dir = isotropic_direction();
        std::uniform_real_distribution<double> d01(0, 1);

        while (scatters < config.max_scatters) {
            // Distance to boundary calculation
            double mu_r = dir.dot(pos);
            double s_bound = std::sqrt(mu_r * mu_r - pos.dot(pos) + config.radius * config.radius) - mu_r;
            
            // Sample delta to next interaction
            double delta_l = -std::log(std::max(d01(rng), 1e-12)) / config.kappa;

            // Escape condition: Photons reaches the surface
            if (delta_l >= s_bound) return weight;

            // Interaction: Move to collision site
            pos.x += dir.x * delta_l; pos.y += dir.y * delta_l; pos.z += dir.z * delta_l;

            // Discrete absorption: Photon survives with probability 'albedo'
            if (d01(rng) <= config.albedo) {
                dir = isotropic_direction();
                scatters++;
            } else {
                return 0.0; // Absorbed
            }
        }
        return 0.0; // Only reached if max_scatters
    }

    // Executes a discrete simulation wave using the current biased probability distribution.
    // Updates CDF and tallies outcomes into cell and global stats.
    void run_wave() {
        global_stats.reset();
        for (auto& c : cells) c.stats.reset();

        // Build CDF for current biasing
        std::vector<double> cdf(cells.size());
        cdf[0] = cells[0].q_bias;
        for (size_t i = 1; i < cells.size(); ++i) cdf[i] = cdf[i-1] + cells[i].q_bias;

        std::uniform_real_distribution<double> d01(0, 1);
        for (int k = 0; k < config.n_photons; ++k) {
            double xi = d01(rng);
            int idx = std::distance(cdf.begin(), std::upper_bound(cdf.begin(), cdf.end(), xi));
            idx = std::min(idx, (int)cells.size() - 1);

            double w_init = (1.0 / config.n_photons) * (cells[idx].p_phys / cells[idx].q_bias);
            double w_out = simulate_photon(random_point_in_cell(idx), w_init);

            cells[idx].stats.add(w_out);
            global_stats.add(w_out);
        }
    }

    // Adjusts the importance sampling weights based on recorded statistical moments.
    // Intended to minimize variance by shifting sample density toward high-contribution cells.
    void adapt_distribution(double a1, double a2) {
        double sum_qt = 0.0;
        std::vector<double> q_target(cells.size());
        
        for (size_t i = 0; i < cells.size(); ++i) {
            // extract (stable) physical escape fraction (\overline{f_i})
            double w0 = (cells[i].q_bias > 0) ? cells[i].p_phys / (cells[i].q_bias * config.n_photons) : 0.0;
            double f_i = (w0 > 0) ? cells[i].stats.mean() / w0 : 0.0;
            
            // apply polynomial guess to stable metric
            double b_i = 1.0 + (a1 * f_i) + (a2 * f_i * f_i); 
            
            q_target[i] = cells[i].p_phys * b_i;
            sum_qt += q_target[i];
        }

        // 3. Apply the Mixing Ratio for absolute stability (c = 0.5)
        const double c = mixing_ratio; 
        for (size_t i = 0; i < cells.size(); ++i) {
            double q_ideal = (sum_qt > 0) ? (q_target[i] / sum_qt) : cells[i].p_phys;
            cells[i].q_bias = c * q_ideal + (1.0 - c) * cells[i].p_phys;
        }
    }

void print_local_stats() {
    std::cout << "\n--- Local Cell Statistics (Wave by Wave) ---\n";
    std::cout << std::left << std::setw(6) << "Cell" 
              << std::setw(10) << "p_phys" 
              << std::setw(10) << "q_bias" 
              << std::setw(12) << "Exp_Esc_W" // Expected Escaping Weight
              << std::setw(12) << "Local_NSR" << "\n";
    
    for (size_t i = 0; i < cells.size(); ++i) {
        const auto& c = cells[i];
        double mu = c.stats.mean();
        double var = c.stats.variance();
        int n = c.stats.n;

        // The term 'esc_prob' is now 'Expected Escaping Weight' per unit initial weight.
        // w_0 = p / (q * N)
        double w0 = (c.q_bias > 0) ? c.p_phys / (c.q_bias * config.n_photons) : 0.0;
        double exp_esc_w = (w0 > 0) ? mu / w0 : 0.0;
        
        double local_nsr = (mu > 0 && n > 0) ? (std::sqrt(var / n) / mu) : 0.0;

        std::cout << std::left << std::setw(6) << i 
                  << std::setw(10) << c.p_phys 
                  << std::setw(10) << c.q_bias 
                  << std::setw(12) << exp_esc_w
                  << std::setw(12);
        
        if (n == 0) {
            std::cout << "NO_SAMPLES";
        } else if (mu <= 0) {
            // This indicates total weight evaporation before boundary intersection
            std::cout << "W_EVAPORATED"; 
        } else {
            std::cout << std::to_string(local_nsr * 100).substr(0, 6) + "%";
        }
        std::cout << "\n";
    }
    std::cout << "--------------------------------------------------\n\n";
}
};

// Handles command-line argument parsing and manages the iterative optimization loop.
// Moniters the escape fraction and NSR to determine when the simulation target is achieved.
int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "usage: " << argv[0] << " <N_per_wave> <kappa> <albedo> <nsr_tolerance> <mixing_ratio>\n";
        return 1;
    }
    int num_cells = 10;

    Config config = { std::stoi(argv[1]), std::stod(argv[2]), 1.0, std::stod(argv[3]) };
    double nsr_tolerance = std::stod(argv[4]);
    mixing_ratio = std::stod(argv[5]);
    
    Environment env(config, num_cells, std::random_device{}());
    
    // Cumulative Statistics (Long-Memory for Stopping Criterion)
    double f1_cum = 0.0; 
    double f2_cum = 0.0; 
    int N_cum = 0;       

    // Bias coefficients 
    // (Set to 0.0, 0.0 for the Unbiased Baseline test)
    double a1 = 2., a2 = 0.0;

    std::cout << std::fixed << std::setprecision(6);

    for (int wave = 0; wave < 200; ++wave) {
        env.run_wave();
        // env.print_local_stats(); // Keep uncommented for local diagnostics if needed

        // 1. Per-Wave Statistics (Short-Memory for Diagnostics)
        double f1_wave = env.global_stats.sum_w;
        double f2_wave = env.global_stats.sum_w2;
        
        // Calculate Wave Efficiency based on COLT's calc_n_eff logic
        double Neff_wave = (f2_wave > 0.0) ? (f1_wave * f1_wave / f2_wave) : 0.0;
        double eff_wave_percent = (Neff_wave / config.n_photons) * 100.0;

        // 2. Cumulative Statistics Update
        f1_cum += f1_wave;
        f2_cum += f2_wave;
        N_cum += config.n_photons;

        // 3. Convergence Evaluation (NSR)
        double cum_nsr = 1.0; // Default to 100% relative error until populated
        if (f1_cum > 0.0 && N_cum > 0) {
            // Exact finite-N relative standard error: sqrt(f2/f1^2 - 1/N)
            double var_term = (f2_cum / (f1_cum * f1_cum)) - (1.0 / N_cum);
            cum_nsr = (var_term > 0.0) ? std::sqrt(var_term) : 0.0;
        }

        double cum_f_esc = f1_cum / (wave + 1);

        std::cout << "Wave " << std::setw(3) << wave 
                  << " | (Cumulative) Escape Fraction: " << cum_f_esc 
                  << " | (Cumulative) NSR: " << cum_nsr 
                  << " | (Wave) Number of Effective Photons: " << eff_wave_percent << "%\n";

        // 4. Stopping Condition
        if (cum_nsr < nsr_tolerance && N_cum > config.n_photons) {
            std::cout << "Target Tolerance Reached in " << wave + 1 << " waves.\n";
            break;
        }

        // 5. Adapt Distribution for Next Wave
        // (Ensure adapt_distribution() is using wave stats, not cumulative stats)
        env.adapt_distribution(a1, a2);
    }

    return 0;
}
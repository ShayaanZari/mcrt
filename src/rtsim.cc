#include "rtsim.h"
Environment::Environment(Config cfg, int num_cells, unsigned int seed)
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

    Vec3 Environment::isotropic_direction() { 
        std::uniform_real_distribution<double> d01(0, 1), dmu(-1, 1);
        double theta = d01(rng) * 2.0 * PI; 
        double mu = dmu(rng); 
        double s = std::sqrt(1.0 - mu * mu);
        return { s * std::cos(theta), s * std::sin(theta), mu };
    }

    Vec3 Environment::random_point_in_cell(int idx) {
        std::uniform_real_distribution<double> d01(0, 1);
        double r = std::pow(std::pow(cells[idx].r_min, 3) + 
                   (std::pow(cells[idx].r_max, 3) - std::pow(cells[idx].r_min, 3)) * d01(rng), 1./3.);
        Vec3 dir = isotropic_direction();
        return { dir.x * r, dir.y * r, dir.z * r };
    }

    double Environment::simulate_photon(Vec3 pos, double weight) {
        int scatters = 0;
        Vec3 dir = isotropic_direction();
        std::uniform_real_distribution<double> d01(0, 1);

        // If weight < cutoff, the photon is "absorbed" by the medium.
        const double weight_cutoff = 1e-9;

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

            // Continuous absorption: Reduce weight by albedo.
            weight *= config.albedo;
            if (weight < weight_cutoff) {
                return 0.0; // Terminate because photon is too weak to matter.
            }
            // Scatter
            dir = isotropic_direction();
            scatters++;
        }
        return 0.0; // Only reached if max_scatters
    }

    void Environment::run_wave() {
        wave_stats.reset();
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

            double weight_out = simulate_photon(
                random_point_in_cell(idx), cells[idx].init_weight(config.n_photons));

            cells[idx].stats.accumulate(weight_out);
            wave_stats.accumulate(weight_out);
        }
    }

    void Environment::adapt_distribution(double a1, double a2) {
        double sum_qt = 0.0;
        std::vector<double> q_target(cells.size());

        for (size_t i = 0; i < cells.size(); ++i) {
            double h_m  = cells[i].h_mean(config.n_photons);
            double h2_m = cells[i].h2_mean(config.n_photons);

            double b_i = 1.0;
            switch (config.biasing_mode) {
                case BiasingMode::NEYMAN_SECOND_MOMENT:
                    // Exact Neyman allocation: q* ∝ p_phys * sqrt(⟨H²⟩)
                    b_i = std::sqrt(h2_m + config.floor_epsilon);
                    break;
                case BiasingMode::FIRST_MOMENT:
                    // Mean-proportional allocation: q* ∝ p_phys * ⟨H⟩
                    b_i = h_m + config.floor_epsilon;
                    break;
                case BiasingMode::POLYNOMIAL:
                    // Legacy polynomial: 1 + a1*⟨H⟩ + a2*⟨H⟩²
                    b_i = 1.0 + (a1 * h_m) + (a2 * h_m * h_m);
                    break;
                case BiasingMode::UNBIASED:
                    // No biasing; q = p_phys
                    b_i = 1.0;
                    break;
            }

            q_target[i] = cells[i].p_phys * b_i;
            sum_qt += q_target[i];
        }

        // Apply the Mixing Ratio for absolute stability
        const double c = config.mixing_ratio;
        for (size_t i = 0; i < cells.size(); ++i) {
            double q_ideal = (sum_qt > 0.0) ? (q_target[i] / sum_qt) : cells[i].p_phys;
            cells[i].q_bias = c * q_ideal + (1.0 - c) * cells[i].p_phys;
        }
    }

    void Environment::print_local_stats() {
        const int N = config.n_photons;
        std::cout << "\n--- Local Cell Statistics (Wave by Wave) ---\n";
        std::cout << std::left
                  << std::setw(6)  << "Cell"
                  << std::setw(10) << "p_phys"
                  << std::setw(10) << "q_bias"
                  << std::setw(12) << "⟨H⟩"        // 1st moment: mean escape response
                  << std::setw(12) << "⟨H²⟩"       // 2nd raw moment
                  << std::setw(12) << "H_rms"       // sqrt(⟨H²⟩)
                  << std::setw(12) << "Skewness"    // standardised 3rd central moment
                  << std::setw(14) << "ExKurtosis"  // excess kurtosis (0 = Gaussian)
                  << std::setw(12) << "Local_NSR"
                  << "\n";

        for (size_t i = 0; i < cells.size(); ++i) {
            const auto& c = cells[i];
            int n = c.stats.sample_count;

            double h1   = c.h_mean(N);
            double h2   = c.h2_mean(N);
            double hrms = c.h_rms(N);
            double skew = c.h_skewness(N);
            double kurt = c.h_excess_kurtosis(N);

            double weight_mean = c.stats.mean();
            double weight_var  = c.stats.variance();
            double local_nsr   = (weight_mean > 0 && n > 0)
                                 ? (std::sqrt(weight_var / n) / weight_mean) : 0.0;

            std::cout << std::left
                      << std::setw(6)  << i
                      << std::setw(10) << c.p_phys
                      << std::setw(10) << c.q_bias
                      << std::setw(12) << h1
                      << std::setw(12) << h2
                      << std::setw(12) << hrms
                      << std::setw(12) << skew
                      << std::setw(14) << kurt
                      << std::setw(12);

            if (n == 0) {
                std::cout << "NO_SAMPLES";
            } else if (weight_mean <= 0) {
                // This indicates total weight evaporation before boundary intersection
                std::cout << "W_EVAPORATED";
            } else {
                std::cout << std::to_string(local_nsr * 100).substr(0, 6) + "%";
            }
            std::cout << "\n";
        }
        std::cout << std::string(100, '-') << "\n\n";
    }

// Manages iterative optimization loop and monitors escape fraction and NSR to determine when the simulation target is achieved.
double compute_escape_fraction(Config config, int num_cells, unsigned int seed, double nsr_tol, bool verbose) {
    Environment env(config, num_cells, seed);

    Accumulator cumulative; // default 4 moments. .skewness() and .excess_kurtosis() are available if needed.
    // no need to define and set the variables to 0.0. Cumulative sum of weights, Cumulative sum of squares of weights, Number of cumulative photons

    // Bias coefficients (will be more flexible in next commit)
    double a1 = 2., a2 = 0.0; 

    // new! store the escape fraction instead of just printing it
    double escape_frac_cum = 0.0; 

    if (verbose) {                                                       
        std::cout << std::fixed << std::setprecision(6);                 
    }     

    for (int wave = 0; wave < config.max_waves; ++wave) {
        env.run_wave();

        // 1. Per-Wave Statistics
        double weight_sum_wave    = env.wave_stats.power_sums[0];
        double weight_sq_sum_wave = env.wave_stats.power_sums[1];

        // 2. Cumulative Statistics Update
        cumulative.merge(env.wave_stats);

        // 3. Convergence Evaluation (NSR)
        double nsr_cum = 1.0; // Default to 100% relative error until populated
        if (cumulative.power_sums[0] > 0.0 && cumulative.sample_count > 0) {
            nsr_cum = std::sqrt(std::max(0.0, cumulative.variance_penalty()));
        }

        escape_frac_cum = cumulative.power_sums[0] / (wave + 1);

        if (verbose) {
            // Calculate Wave Efficiency based on COLT's calc_n_eff logic
            double n_eff_pct = (env.wave_stats.n_eff() / config.n_photons) * 100.0;

            std::cout << "Wave " << std::setw(3) << wave
                      << " | (Cumulative) Escape Fraction: " << escape_frac_cum
                      << " | (Cumulative) NSR: " << nsr_cum
                      << " | (Wave) N_eff/N: " << n_eff_pct << "%\n";
        }

        // 4. Stopping Condition
        if (nsr_cum < nsr_tol && cumulative.sample_count > config.n_photons) {
            if (verbose) {
                std::cout << "Target Tolerance Reached in " << wave + 1 << " waves.\n";
            }
            break;
        }

        // 5. Adapt Distribution for Next Wave
        env.adapt_distribution(a1, a2);
    }

    if (verbose) {
        std::cout << "\n=== Final Simulation State ===\n";
        env.print_local_stats();
    }

    return escape_frac_cum;
}
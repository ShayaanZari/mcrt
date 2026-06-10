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

            double weight_out = simulate_photon(
                random_point_in_cell(idx), cells[idx].init_weight(config.n_photons));

            cells[idx].stats.add(weight_out);
            global_stats.add(weight_out);
        }
    }

    void Environment::adapt_distribution(double a1, double a2) {
        double sum_qt = 0.0;
        std::vector<double> q_target(cells.size());
        
        for (size_t i = 0; i < cells.size(); ++i) {
            // h_mean: unit-weight escape response, stripped of importance sampling bias.
            double iw = cells[i].init_weight(config.n_photons);
            double h_mean = (iw > 0) ? cells[i].stats.mean() / iw : 0.0;
            
            // apply polynomial guess to stable metric
            double b_i = 1.0 + (a1 * h_mean) + (a2 * h_mean * h_mean); 
            
            q_target[i] = cells[i].p_phys * b_i;
            sum_qt += q_target[i];
        }

        // Apply the Mixing Ratio for absolute stability
        const double c = config.mixing_ratio; 
        for (size_t i = 0; i < cells.size(); ++i) {
            double q_ideal = (sum_qt > 0) ? (q_target[i] / sum_qt) : cells[i].p_phys;
            cells[i].q_bias = c * q_ideal + (1.0 - c) * cells[i].p_phys;
        }
    }

    void Environment::print_local_stats() {
        std::cout << "\n--- Local Cell Statistics (Wave by Wave) ---\n";
        std::cout << std::left << std::setw(6) << "Cell" 
                  << std::setw(10) << "p_phys" 
                  << std::setw(10) << "q_bias" 
                  << std::setw(12) << "h_mean" // Unit-weight escape response
                  << std::setw(12) << "Local_NSR" << "\n";
        
        for (size_t i = 0; i < cells.size(); ++i) {
            const auto& c = cells[i];
            double weight_mean = c.stats.mean();
            double weight_var = c.stats.variance();
            int n = c.stats.n;

            // h_mean: unit-weight escape response, stripped of importance sampling bias.
            double iw = c.init_weight(config.n_photons);
            double h_mean = (iw > 0) ? weight_mean / iw : 0.0;
            
            double local_nsr = (weight_mean > 0 && n > 0) ? (std::sqrt(weight_var / n) / weight_mean) : 0.0;

            std::cout << std::left << std::setw(6) << i 
                      << std::setw(10) << c.p_phys 
                      << std::setw(10) << c.q_bias 
                      << std::setw(12) << h_mean
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
        std::cout << "--------------------------------------------------\n\n";
    }

// Manages iterative optimization loop and monitors escape fraction and NSR to determine when the simulation target is achieved.
double compute_escape_fraction(Config config, int num_cells, unsigned int seed, double nsr_tol, bool verbose) {
    Environment env(config, num_cells, seed);
    
    // Cumulative Statistics
    double weight_sum_cum = 0.0; // Cumulative sum of weights 
    double weight_sq_sum_cum = 0.0;  // Cumulative sum of squares of weights
    int n_photons_cum = 0; // Number of cumulative photons

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
        double weight_sum_wave = env.global_stats.weight_sum;
        double weight_sq_sum_wave = env.global_stats.weight_sq_sum;

        // 2. Cumulative Statistics Update
        weight_sum_cum += weight_sum_wave;
        weight_sq_sum_cum += weight_sq_sum_wave;
        n_photons_cum += config.n_photons;

        // 3. Convergence Evaluation (NSR)
        double nsr_cum = 1.0; // Default to 100% relative error until populated
        if (weight_sum_cum > 0.0 && n_photons_cum > 0) {
            // Exact finite-N relative standard error: sqrt(f2/f1^2 - 1/N)
            double var_term = (weight_sq_sum_cum / (weight_sum_cum * weight_sum_cum))
                            - (1.0 / n_photons_cum);
            nsr_cum = (var_term > 0.0) ? std::sqrt(var_term) : 0.0;
        }

        escape_frac_cum = weight_sum_cum / (wave + 1);

        if (verbose) {
            // Calculate Wave Efficiency based on COLT's calc_n_eff logic
            double n_eff_wave = (weight_sq_sum_wave > 0.0)
                                    ? (weight_sum_wave * weight_sum_wave / weight_sq_sum_wave)
                                    : 0.0;
            double n_eff_pct = (n_eff_wave / config.n_photons) * 100.0;

            std::cout << "Wave " << std::setw(3) << wave
                      << " | (Cumulative) Escape Fraction: " << escape_frac_cum
                      << " | (Cumulative) NSR: " << nsr_cum
                      << " | (Wave) N_eff/N: " << n_eff_pct << "%\n";
        }

        // 4. Stopping Condition
        if (nsr_cum < nsr_tol && n_photons_cum > config.n_photons) {
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
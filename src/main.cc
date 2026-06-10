#include "rtsim.h"
#include <iostream>

// Handles command-line input, builds Config, and starts the simulation
int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "usage: " << argv[0] << " <N_per_wave> <kappa> <albedo> <nsr_tolerance> <mixing_ratio> [seed]\n";
        return 1;
    }
    int num_cells = 10;

    Config config = { 
        std::stoi(argv[1]), 
        std::stod(argv[2]), 
        1.0, 
        std::stod(argv[3]),
        10000,
        std::stod(argv[5]) // mixing_ratio
    };
    
    double nsr_tol = std::stod(argv[4]);
    
    unsigned int seed = (argc > 6) ? std::stoul(argv[6]) : std::random_device{}();
    
    return compute_escape_fraction(config, num_cells, seed, nsr_tol, true) >= 0 ? 0 : 1;
}
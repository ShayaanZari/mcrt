#include "rtsim.h"
#include <cmath>

// Compares against a converged value with medium configuration: kappa=4 albedo=0.4
int main() { 
    Config config = { 10000, 4.0, 1.0, 0.4, 10000, 0.7 };

    // truth computed via 10M photon unbiased run on above params w/ seed 42
    const double truth = 0.257662;
    const double tolerance = 0.005;
    std::cout << "Stored Truth Escape Fraction: " << truth << "\n";
    
    // run biased sim
    std::cout << "Running adaptive sim (10k photons/wave, mix=0.7)...\n";
    double result = compute_escape_fraction(config, 10, 42, 0.01);
    std::cout << "Biased Escape Fraction: " << result << "\n";
    
    double diff = std::abs(result - truth);
    std::cout << "Difference: " << diff << "\n";
    
    if (diff <= tolerance) {
        std::cout << "SUCCESS: Difference is within tolerance (" << tolerance << ").\n";
        return 0; // Success code for CTest
    } else {
        std::cerr << "FAIL: Difference exceeds tolerance (" << tolerance << ").\n";
        return 1; // Failure code for CTest
    }
}
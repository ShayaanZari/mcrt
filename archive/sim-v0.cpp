#include <iostream>
#include <random>
#include <chrono>
#include <cmath> // For sqrt, cos, sin, and M_PI

// for output
#include <fstream>
#include <vector>
#include <iomanip>

struct Vector3D {
    double x, y, z;

    // Compute Euclidean length
    double mag() const {
        return std::sqrt(x*x + y*y + z*z);
    }
};

Vector3D generate_random_direction(std:: mt19937& rng) {
    const double pi = std::acos(-1.0);
    
    // theta in [0, 2*pi)
    std::uniform_real_distribution<double> dist_theta(0.0, 2.0 * pi);
    
    // mu (cos(phi)) in [-1, 1]
    std::uniform_real_distribution<double> dist_mu(-1.0, 1.0);

    double theta = dist_theta(rng);
    double mu = dist_mu(rng);
    double zenith_sine = std::sqrt(1.0 - mu * mu);

    double nx = zenith_sine * std::cos(theta);
    double ny = zenith_sine * std::sin(theta);
    double nz = mu;

    return Vector3D {nx, ny, nz};
}

struct Photon {
    Vector3D pos = Vector3D{0, 0, 0}; // 3D position vector. Initial: origin
    Vector3D dir; // 3D direction vector. Initial: randomly generated

    void move(double step_size) {
        pos.x += dir.x * step_size;
        pos.y += dir.y * step_size;
        pos.z += dir.z * step_size;
    }
};

struct Medium {
    double kappa = 0.5; // absorption pacity (total)
    double R = 1.; // radius of sphere
};

bool simulate(const Medium& m, std::mt19937& rng, std::uniform_real_distribution<double>& dist_uniform) {
    Photon p;
    p.dir = generate_random_direction(rng);
    
    double xi = dist_uniform(rng);
    
    double step = -std::log(xi) / m.kappa; 
    p.move(step);
    
    // Return true if the photon is no longer inside the medium
    return p.pos.mag() >= m.R;
}

struct SimulationResult {
    double kappa;
    double escape_fraction;
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <num_packets> <kappa>\n";
        return 1;
    }

    int num_packets = std::stoi(argv[1]);
    Medium m {std::stod(argv[2]), 1.};
    
    std::mt19937 rng(std::random_device{}()); 
    std::uniform_real_distribution<double> dist_uniform(0.0, 1.0);
    
    //run_sim(num_packets);
    
    int escaped_packets = 0;
    
    for (int i=0; i<num_packets; i++) {
        if (simulate(m, rng, dist_uniform)) {
            escaped_packets++;
        }
        // else: absorbed
    }

    std::cout << "Escaped packets: " << escaped_packets << std::endl;
    std::cout << "Escape fraction: " << (double)escaped_packets/num_packets << std::endl;
    
    return 0;
}

int run_sim(int num_packets) {
    // Define the range of kappa values
    std::vector<double> kappa_values = {0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 10.0};
    std::vector<SimulationResult> results;
    std::mt19937 rng(std::random_device{}()); 
    std::uniform_real_distribution<double> dist_uniform(0.0, 1.0);
    
    // Perform the sweep
    for (double k : kappa_values) {
        Medium m{k, 1.0};
        int escaped_packets = 0;

        for (int i = 0; i < num_packets; i++) {
            if (simulate(m, rng, dist_uniform)) {
                escaped_packets++;
            }
        }

        double fraction = static_cast<double>(escaped_packets) / num_packets;
        results.push_back({k, fraction});
        
        std::cout << "Kappa: " << k << " | Escape Fraction: " << fraction << std::endl;
    }

    // Save to CSV
    std::ofstream outFile("mcrt_results.csv");
    outFile << "kappa,escape_fraction\n"; // Header
    for (const auto& res : results) {
        outFile << res.kappa << "," << res.escape_fraction << "\n";
    }
    outFile.close();

    return 0;
}
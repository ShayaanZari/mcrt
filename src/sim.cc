#include <iostream>
#include <random>
#include <chrono>
#include <cmath> 

constexpr double PI = 3.14159265358979323846; 

struct Vec3 {
    double x, y, z;
    double mag() const {return std::sqrt(x*x + y*y + z*z);}
};

struct Sampler {
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist_01{0.0, 1.0};
    std::uniform_real_distribution<double> dist_mu{-1.0, 1.0};
    Sampler() : rng(std::random_device{}()) {}

    Vec3 isotropic_direction() {
        double theta = dist_01(rng)*2.0*PI;
        double mu = dist_mu(rng);
        double zenith_sine = std::sqrt(1.0 - mu * mu);

        return Vec3 {
            zenith_sine * std::cos(theta),
            zenith_sine * std::sin(theta),
            mu
        };
    }
};

struct Config {
    int num_packets;
    double kappa; // absorption or opacity (total)
    double radius; // radius of sphere

    Config(int n, double r, double k) 
        : num_packets(n), radius(r), kappa(k) {}
};

struct Photon {
    Vec3 pos = Vec3{0, 0, 0}; // 3D position vector. Initial: origin
    Vec3 dir; // 3D direction vector. Initial: randomly generated
    double weight; // packet weight

    Photon(const Config& config, Sampler& smp) {
        dir = smp.isotropic_direction();
        weight = 1./config.num_packets; // 1/N
    }

    void move(double step_size) {
        pos.x += dir.x * step_size;
        pos.y += dir.y * step_size;
        pos.z += dir.z * step_size;
    }
};

double simulate(const Config& config, Sampler& smp) {
    Photon p(config, smp);
    double tau = config.kappa*config.radius;
    p.weight = p.weight*std::exp(-tau);
    return p.weight;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <num_packets> <kappa>\n";
        return 1;
    }

    int num_packets = std::stoi(argv[1]);
    double kappa = std::stod(argv[2]);
    Config config(num_packets, 1., kappa);
    Sampler smp;
    
    double fraction = 0;
    
    for (int i=0; i<num_packets; i++) {
        fraction += simulate(config, smp);
    }

    std::cout << "Escape fraction: " << fraction << std::endl;
    
    return 0;
}
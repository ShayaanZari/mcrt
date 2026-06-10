### Intro
This repo serves as a proof of concept for, or a minimal reproducible example of integrating importance sampling into the COLT (Cosmic Lyman-alpha Transfer) solver:
- [Source code](https://bitbucket.org/aaron_smith/colt/src/master/)
- [Smith et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.4336S/abstract)
- [Documentation](https://colt.readthedocs.io/en/latest/)

Why implement in this sandbox rather than opening a PR? Primarily for the learning process, but also to reduce regression risk by:
1. Verifying mathematical foundations under simple constraints
2. Validating against architectural and performance requirements for the solver (pre-MPI and HDF5, etc)

Development milestones and iteration history are tracked in [plan.md](plan.md).

## Repository Structure

```text
.
├── CMakeLists.txt         # Build configuration for CMake & CTest
├── README.md              # Project overview and file structure
├── plan.md                # Development roadmap, milestones, and version log
└── src/                   # C++ source code directory
    ├── main.cc            # CLI entry point and simulation driver
    ├── regression_test.cc # Validation test checking results against reference truth
    ├── rtsim.cc           # Monte Carlo transport & adaptive distribution updates
    └── rtsim.h            # Grid cells, configuration, and trackers
```

## Core Components
- **Simulation Interface ([src/rtsim.h](src/rtsim.h))**:
    - Defines 3D vector geometry ([Vec3](src/rtsim.h#L17-L20)) and statistical tracking ([WaveStats](src/rtsim.h#L24-L45)).
    - Implements the [Cell](src/rtsim.h#L48-L59) structure for radial shell discretization.
    - Declares the [Environment](src/rtsim.h#L70-L101) class for managing Monte Carlo transport and the [compute_escape_fraction](src/rtsim.h#L103-L110) interface.
- **Simulation Engine ([src/rtsim.cc](src/rtsim.cc))**:
    - Implements Beer-Lambert path sampling, isotropic direction generation, and boundary intersections in [Environment::simulate_photon](src/rtsim.cc#L33-L65).
    - Implements the Monte Carlo wave-by-wave loop in [Environment::run_wave](src/rtsim.cc#L67-L88).
- **Adaptive Biasing**:
    - Located in [Environment::adapt_distribution](src/rtsim.cc#L90-L112).
    - Dynamically adjusts photon spawn distributions (`q_bias`) based on recorded escape responses to minimize variance.
    - Uses a `mixing_ratio` to balance the ideal proposal with the physical probability for numerical stability.
- **Convergence Tracking**:
    - Calculated in [compute_escape_fraction](src/rtsim.cc#L154-L226).
    - Computes the **Normalized Standard Error (NSR)** cumulatively.
    - Target precision is achieved when the cumulative NSR falls below the user-defined tolerance.

## Build Instructions
Standard CMake build process:
```bash
cmake -B build -S .
cmake --build build --config Release
```

## Running the Simulation
Execute the compiled binary with the following parameters:
```bash
./build/sim <N_per_wave> <kappa> <albedo> <nsr_tolerance> <mixing_ratio> [seed]
```
- `N_per_wave`: Number of photons to simulate per iteration (wave).
- `kappa`: Extinction coefficient (opacity).
- `albedo`: Scattering albedo (ratio of scattering to total extinction, 0.0 to 1.0).
- `nsr_tolerance`: Target relative error (e.g., `0.01` for 1%).
- `mixing_ratio`: Balance between physical and biased sampling (e.g., `0.7`).
- `seed` (optional): Random seed for reproducibility.

### Example Run
Running the simulation with 10k photons/wave, $\kappa = 4.0$, albedo = 0.4, NSR tolerance of 1%, and a mixing ratio of 0.7:
```bash
$ ./build/sim 10000 4.0 0.4 0.01 0.7 42
Wave   0 | (Cumulative) Escape Fraction: 0.255754 | (Cumulative) NSR: 0.014531 | (Wave) N_eff/N: 32.139009%
Wave   1 | (Cumulative) Escape Fraction: 0.255455 | (Cumulative) NSR: 0.009989 | (Wave) N_eff/N: 34.730660%
Target Tolerance Reached in 2 waves.

=== Final Simulation State ===

--- Local Cell Statistics (Wave by Wave) ---
Cell  p_phys    q_bias    h_mean      Local_NSR   
0     0.001000  0.000791  0.017715    66.232%     
1     0.007000  0.005573  0.010288    33.232%     
2     0.019000  0.015045  0.053780    24.032%     
3     0.037000  0.030354  0.062371    17.056%     
4     0.061000  0.050862  0.064726    12.420%     
5     0.091000  0.078082  0.121080    7.7254%     
6     0.127000  0.113715  0.153692    5.4490%     
7     0.169000  0.158570  0.187732    4.1640%     
8     0.217000  0.223674  0.282343    2.8221%     
9     0.271000  0.323334  0.455923    1.6362%     
--------------------------------------------------
```

## Testing
Run the regression test using CTest:
```bash
ctest --test-dir build --output-on-failure
```
The test verifies that the adaptive simulation converges to the expected escape fraction within a tolerance (currently set to 0.005) under configuration defined in [src/regression_test.cc](src/regression_test.cc).


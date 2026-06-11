### `rtsim.h` and `rtsim.cc`: stateless metrics and stateful simulation

Code is separated by purpose rather than function length. 
- `rtsim.h` contains type declarations and stateless/const query metrics. 
    - E.g. the lengthy mathematical methods of `Accumulator` such as `skewness()` are defined inline within the header because they operate on frozen snapshot data.
- `rtsim.cc` contains the stateful execution logic, as well as the stochastic sampling methods, and environment mutations. 
    - E.g. `random_point_in_cell()`, while relatively short, belongs to `Environment` and mutates the stateful random number generator, so it is placed in the source file to keep compilation dependencies clean and enforce interface separation.

### Accumulator Struct

The `Accumulator` struct computes the raw moments of escaped photon weights by accumulating the sums of powers ($w^k$), and provides operations such as NSR and Kish's ESS ($N_eff$). It is utilized in three distinct areas of the simulation:
1. **Per-Cell Statistics (`Cell::stats`)**
    - Every `Cell` stores its own `Accumulator` instance, which is reset at the start of each wave.
    - Tracks the escaped weights of photons spawned within a specific cell. 
    - Having cell-wise `Accumulator`s gives the simulation ability to access to local statistical metrics (e.g. local mean escape response, sample count, variance), which may be used in the adjustment of the sampling distribution for subsequent waves.
2. **Per-Wave Statistics (`Environment::wave_stats`)**
    - A member of the `Environment` struct; reset at the beginning of each wave.
    - Tracks the escaped weights of all photons simulated within the current wave.
    - Compute wave-level statistics, such as the wave's ($N_{eff}$).
3. **Cumulative Statistics (`cumulative` in `compute_escape_fraction`)**:
    - Declared locally within the main simulation loop.
    - Merges the stats of each completed wave using `Accumulator::merge()`. 
    - Tracks the overall progress of the simulation, computing the cumulative escape fraction and the cumulative relative error (NSR) to determine when the convergence tolerance threshold has been reached.
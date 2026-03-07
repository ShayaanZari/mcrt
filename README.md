### Intro
This repo serves as a proof of concept for, or a minimal reproducible example of integrating importance sampling into the COLT (Cosmic Lyman-alpha Transfer) solver:
- [Source code](https://bitbucket.org/aaron_smith/colt/src/master/)
- [Smith et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.4336S/abstract)
- [Documentation](https://colt.readthedocs.io/en/latest/)

Why implement in this sandbox rather than opening a PR? Primarily for the learning process, but also to reduce regression risk by
1. Verifying mathematical foundations under simple constraints
2. Validating against architectural and performance requirements for the solver (pre-MPI and HDF5, etc)

Design principles:
1. Where obvious, match structure names (e.g. `Vec3`) and function signatures (e.g. `isotropic_direction`) of the COLT source code to reduce PR friction, also acting as an on-ramp to the solver.
2. AI as an accelerator, never the architect. Given a problem $(X)$, I decide that I need to solve it by doing $(Y)$, so I use AI to rapidly understand various methods for accomplishing $(Y)$, and as a fast syntax search for the implementation of $(Y)$.
3. Minimal but efficient code: OOP principles to decouple simulation, configuration, RNG, computation methods, and I/O.  

I track my completed iterations as well as future iterations in `plan.md`.

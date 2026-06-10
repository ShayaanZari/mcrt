### Intro
This repo serves as a proof of concept for, or a minimal reproducible example of integrating importance sampling into the COLT (Cosmic Lyman-alpha Transfer) solver:
- [Source code](https://bitbucket.org/aaron_smith/colt/src/master/)
- [Smith et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.4336S/abstract)
- [Documentation](https://colt.readthedocs.io/en/latest/)

Why implement in this sandbox rather than opening a PR? Primarily for the learning process, but also to reduce regression risk by:
1. Verifying mathematical foundations under simple constraints
2. Validating against architectural and performance requirements for the solver (pre-MPI and HDF5, etc)

Thought process during designing is explained in `design.md`.

I track my completed iterations as well as future iterations in `plan.md`.

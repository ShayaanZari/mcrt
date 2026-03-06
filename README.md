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

---

### Version log and plan of implementation

#### v0 - event-based pure absorption (complete: `archive/sim-v0.cc`)
Assumptions:
1. Optical depth is homogeneous $\tau=\int_0^R\kappa ds=\kappa R=\kappa$
2. Scattering albedo $\omega=0$

`Photon` packet with `position` $\vec r_0$ initialized to origin and `direction` $\hat n$ initialized to random isotropic direction (using the $\phi$ and $\mu$ sampling scheme outlined in section 4.2.1 of the MCRT paper)

Compute distance to interaction: 
$$\triangle l=-\frac{\ln\xi}{\kappa}$$.
Move the packet: $\vec r=\vec r_0+\hat n\triangle l$.
If $|\vec r|\geq R$, where $R$ is radius of the sphere, photon escapes and count is incremented.
Otherwise, photon is absorbed (do nothing).

Graph:
![alt2](plots/exponential_curve.png?raw=true "Title")

#### v1 - continuous attenuation (complete: `src/sim.cc`)
Start `Photon.weight` with $1/N$. On each step, multiply by $e^{-\tau}$. As there is only one step currently, the sum of the final weights is trivially:

$$
\sum_{i=1}^N(\frac{1}{N}e^{-\tau})=e^{-\tau}
$$

#### v2 - off-center point source, multiple sources
Vary the location of the point source, so distance $s$ from source position to edge of sphere is not equal to $R$. Calculate $s$ using `face_distance` function of `spherical.cc`. Compute $e^{-\kappa s}$ instead of $e^{-kappa}$. 

Escape fraction: Let $r$ be the distance of the source from the origin. $f_\text{esc}$ should be higher at $r=R/2$ than $r=0$. For $r=R$ (i.e. source at edge of sphere), $f_\text{esc}$ should approach $0.5$ as $\tau\to\infty$, as exactly half protons should have a direction pointing outwards.

#### v3 - multiple point sources
Create `Source` struct with position and luminosity fields. Luminosity is a simple weight. Create vector of `Source`s. Run $N$ packets for each `Source`. Compute total escape fraction:

$$
f_\text{esc,tot}=\dfrac{\sum L_i f_\text{esc,i}}{\sum L_i}
$$

Escape fraction should be skewed towards the escape fraction of sources with higher luminosity.

#### v4 - source importance sampling
Compute total luminosity of all sources $L_\text{tot}$. Initialize CDF array $C$, where the $i$-element of the array is

$$
C[i]=C[i-1] + \dfrac{L_i}{L_\text{tot}}
$$

where $C[0]=\frac{L_1}{L_\text{tot}}$.

Make two changes to code:
- Sourcing: for each proton, draw a random number $\xi=[0,1)$. Find first index $j$ where $\xi<C[j]$. The proton's source is determined by the source which corresponds to $C[j]$.
- Weighting for conservation of energy: initial weight $W_i$ of each photon is inversely proportional to the probability of its source being sampled from:

$$
W_i\propto L_\text{tot}/N
$$

#### v5 - radial shells
Reference: `random_point_in_cell` function from `spherical.cc` for uniform volumetric sampling. $\vec r =R\cdot \sqrt[3]{\xi}$ for random number $\xi\in[0,1)$.

Compute noise of each shell: 

$$
\dfrac{\sum f_\text{esc}^2}{(\sum f_\text{esc})^2}
$$


#### v6 - scattering
Scattering albedo:

$$
\omega=\dfrac{\kappa_\text{scat}}{\kappa}
$$

Random walk: For each photon, after initial movement step, if it has not escaped, draw $\xi\in[0,1)$. If $\xi\leq\omega$, it scatters: draw a new $\va n$, take a step in the corresponding direction, and repeat until absorption $\xi>\omega$ or escape $\va r\geq R$. 
Add a limit on number of scatters.

At $\omega=1.0$, $f_\text{esc}$ should asymptote to $1.0$.

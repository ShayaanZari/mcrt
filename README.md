### Version log and plan of implementation

#### v0 - event-based pure absorption (done, `archive/sim-v1.cc`)
Assumptions:
1. Optical depth is homogeneous $\tau=\int_0^R\kappa ds=\kappa R=\kappa$
2. Scattering albedo $\omega=0$

`Photon` packet with `position` $\vec r_0$ initialized to origin and `direction` $\hat n$ initialized to random isotropic direction.

Compute distance to interaction: 
$$\triangle l=-\frac{\ln\xi}{\kappa}$$.
Move the packet: $\vec r=\vec r_0+\hat n\triangle l$.
If $|\vec r\geq R$, where $R$ is radius of the sphere, photon escapes and count is incremented.
Otherwise, photon is absorbed (do nothing).

Graph:
![alt2](plots/exponential_curve.png?raw=true "Title")

#### v1 - continuous attenuation (done, `src/sim.cc`)
Start `Photon.weight` with $1/N$. On each step, multiply by $e^{-\tau}$. As there is only one step currently, the sum of the final weights is trivially:

$$
\sum_{i=1}^N(\frac{1}{N}e^{-\tau})=e^{-\tau}
$$

#### v2 - off-center point source
Reference: `face_distance` function of `spherical.cc`. Return $e^{-\kappa d}$.

#### v3 - multiple point sources
Create `Source` struct with position and luminosity fields. Create vector of `Source`s. Compute total escape fraction:

$$
f_\text{esc,tot}=\dfrac{\sum L_i f_\text{esc,i}}{\sum L_i}
$$

#### v4 - importance sampling
Compute total luminosity. Create a CDF array...

#### v5 - radial shells
Reference: `random_point_in_cell` function from `spherical.cc`.

implementation plan

### v0 - pure absorption
Photon packet with fields Vec3 position initialized to $\vec r_0=(0,0,0)$ and Vec3 direction initialized to random isotropic direction $\hat n$.

Compute distance to interaction: $\triangle l=-\frac{\ln\xi}{\kappa}$.
Move the packet: $\vec r=\vec r_0+\hat n\triangle l$.
If $|\vec r\geq R$, where $R$ is radius of the sphere, photon escapes and increment the count.
Otherwise, photon is absorbed (do nothing).

Graph:
![alt2](plots/exponential_curve.png?raw=true "Title")
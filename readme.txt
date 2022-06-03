Calculating the scattering reflected and transmitted via a slab of Omega-medium.

This code applies the general theory of propagation in multilayered structures [1] to reproduce the results shown in [2], Figure 4.

[1] EI3302 Propagation and scattering in multilayered structures, Martin Norgren, KTH, May 10, 2022. Handouts for the course EI3302 Electromagnetic Waves in Complex Media.

[2] Norgren M. and He S., “Electromagnetic reflection and transmission for a dielectric-Ω interface and an Ω slab”, International Journal of Infrared and Millimeter Waves, 15(9), 1537-1554, 1994.

-----------------------------------------------------------------------

Important differences between [1] and [2], and other notes:
- in [1] material tensors are defined wrt the dimensionless relative parameters (e.g. relative permittivity) while in [2] the vacuum permittivity and permeability are infused in the tensor matrices.
- Vacuum wavenumber k_0 is not needed, since in calculating the W matrix in [1] all occurences of k_t vector is accompanied with division by k_0.
- In vacuum the exponential matrix M and thus also propagator matrix P is the identity matrix I. It is therefore sufficient to calculate the propagator for the slab only. 
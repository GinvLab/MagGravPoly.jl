[![Docs Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliageoph.gitlab.io/MagGrav2Dpoly.jl/stable)
[![Docs Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliageoph.gitlab.io/MagGrav2Dpoly.jl/dev)


# MagGrav2Dpoly

**MagGrav2Dpoly** is a Julia package to perform magnetic and gravity anomaly calculations using a 2D or 2.75D parameterization in terms of polygons with uniform arbitrary magnetizations and density contrasts. It provides functions to both solve the forward problem and to calculate the gradient of a given misfit function. Such functions can be used to solve inverse problems both in the deterministic and probabilistic approach. In particular, this package provides some functions to solve inverse problems using the Hamiltonian Monte Carlo (HMC) method, as part of the `HMCLab` project (see the  `HMCSampler.jl` package). Gradients are calculated using the technique of automatic differentiation.
With this package it is also possible to perform joint magnetic and gravity forward and gradient calculations and hence solve joint inverse problems, see the tutorials below.

The forward problem formulations for the magnetic case implemented in this package are the following:
* 2D case: Talwani & Heirtzler (1962, 1964), Won & Bevis (1987) and revised Kravchinsky, Hnatyshin, Lysak, & Alemie (2019);
* 2.75D case: revised Rasmussen & Pedersen (1979) and Campbell (1983).

The forward problem formulations for the gravity case implemented in this package are the following:
* 2D case: Talwani, Worzel, & Landisman (1959), with background theory derived from the paper of Hubbert (1948);
* 2.75D case:  Rasmussen & Pedersen (1979).

  
If you use this code for research or else, please cite the related papers:

* Ghirotto, Zunino, Armadillo &, Mosegaard (2021). **Magnetic Anomalies Caused by 2D Polygonal Structures with Uniform Arbitrary Polarization: new insights from analytical/numerical comparison among available algorithm formulations**. *Geophysical Research Letters, 48*(7), e2020GL091732, https://doi.org/10.1029/2020GL091732.

* Zunino, Ghirotto, Armadillo, & Fichtner (2022). **Hamiltonian Monte Carlo probabilistic joint inversion of 2D (2.75D) gravity and magnetic data**. *Geophysical Research Letters*,  49, e2022GL099789. https://doi.org/10.1029/2022GL099789.

Regarding solving the inverse problem with the HMC method, please see the following paper and check out the package `HMCSampler.jl`:

* Zunino, Gebraad, Ghirotto, & Fichtner (2023). **HMCLab: a framework for solving diverse geophysical inverse problems using the Hamiltonian Monte Carlo method**. in prep.

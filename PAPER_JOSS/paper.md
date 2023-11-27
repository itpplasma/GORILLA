---
title: 'GORILLA: Guiding-center ORbit Integration with Local Linearization Approach'
tags:
  - Fortran
  - Theoretical physics
  - Computational physics
  - Plasma physics
  - Plasma dynamics
  - Numerical integration
  - Hamiltonian mechanics
  - Kinetic theory
  - Magnetic confinement fusion
  - Fusion reactors
  - Tokamaks
  - Stellarators
authors:
  - name: Michael Eder^[corresponding author]
    orcid: 0000-0002-1392-1664
    affiliation: 1
  - name: Lukas M. P. Bauer
    orcid: 0000-0003-3341-4085
    affiliation: 1
  - name: Daniel Forstenlechner
    affiliation: 1
  - name: Georg S. Graßler
    affiliation: 1
  - name:  Sergei V. Kasilov
    affiliation: "1, 2, 3"
  - name:  Winfried Kernbichler
    affiliation: 1
  - name:  Markus Meisterhofer
    affiliation: 1
  - name:  Michael Scheidt
    affiliation: 1
  - name: Christopher G. Albert
    orcid: 0000-0003-4773-416X
    affiliation: 1

affiliations:
 - name: Institut für Theoretische Physik - Computational Physics, Technische Universität Graz, Petersgasse 16, 8010 Graz, Austria
   index: 1
 - name: Institute of Plasma Physics, National Science Center, “Kharkov Institute of Physics and Technology,” Akademicheskaya str. 1, 61108 Kharkov, Ukraine
   index: 2
 - name: Department of Applied Physics and Plasma Physics, V. N. Karazin Kharkov National University, Svobody sq. 4, 61022 Kharkov, Ukraine
   index: 3
date: 10 May 2022
bibliography: paper.bib

---
# Introduction
Extremely hot plasmas with a temperature of the order of hundred million degrees Celsius are needed to produce energy from nuclear fusion. Under these conditions, hydrogen isotopes are fused, and energy is released. The energy release from 1 kg of fusion fuel corresponds approximately to that of 10000 tons of coal. This high energy density has made fusion the subject of worldwide research projects. The confinement of such hot plasmas, however, poses major physical and technological problems for researchers. In particular, complex numerical methods are necessary to understand the physics of such plasmas in complicated toroidal magnetic fields.

An important kinetic approach for simulating the collective behavior of a plasma uses direct modeling of particle orbits. A well-known approximation for computing the motion of electrically charged particles in slowly varying electromagnetic fields is to reduce the dynamical equations by separating the relatively fast circular motion around a point called the guiding-center, and primarily treat the relatively slow drift motion of this point. This drift motion is described by the guiding-center equations in non-canonical Hamiltonian form [@boozer_guiding_1980; @littlejohn_variational_1983; @cary_hamiltonian_2009].

Here, we provide `GORILLA`: an efficient code for the purpose of solving the guiding-center equations. This code is a numerical implementation of the novel, quasi-geometric integration method in 3D geometry described by @eder_quasi-geometric_2020. A related method for axisymmetric geometry has been realized in the code K2D [@kasilov_geometric_2016].
Apart from this quasi-geometric approach there exists a variety of codes for tracing guiding-center drift motion including ORBIT [@white_hamiltonian_1984], DCOM [@wakasa2001monte], MOCA [@tribaldos_monte_2001], FORTEC-3D [@satake_non-local_2005], VENUS [@isaev_venusf_2006], NEO-MC [@allmaier_variance_2008], XGC [@ku_full-f_2009], ANTS [@drevlak_thermal_2009], ASCOT [@kurki-suonio_ascot_2009], VENUS-LEVIS [@pfefferle_venus-levis_2014], SIMPLE [@albert_symplectic_2020]; see also @beidler_benchmarking_2011. Most of these codes solve the guiding-center equations primarily in the plasma core of toroidal fusion devices. Due to its formulation in general curvilinear coordinates, `GORILLA` is not limited by the field topology. That means that the computation domain of `GORILLA` covers both the closed field line region (i.e., the plasma core) and the open field line region (i.e., the scrape-off layer).

# Summary

`GORILLA` is a Fortran code that computes guiding-center orbits for charged particles of given mass, charge and energy in toroidal fusion devices with three-dimensional field geometry.
Conventional methods for integrating the guiding-center equations use high order interpolation of the electromagnetic field in space.
In `GORILLA`, a linear interpolation of a specific representation of quantities in the guiding-center system employing a spatial mesh is used for the discretization of the electromagnetic field.
This leads to locally linear equations of motion with piecewise constant coefficients.
As shown by @eder_quasi-geometric_2020, this local linearization approach retains the Hamiltonian structure of the guiding-center equations. For practical purposes this means that the total energy, the magnetic moment and the phase-space volume are conserved.
Furthermore, the approach reduces computational effort and sensitivity to noise in the electromagnetic field. In `GORILLA`, guiding-center orbits are computed without taking into account collisions between particles. Two examples of guiding-center orbits obtained with `GORILLA` can be seen in \autoref{fig:example} where the magnetic field of a real-world fusion device is used, specifically the tokamak “ASDEX Upgrade”.

# Statement of need

`GORILLA` is designed to be used by researchers in scientific plasma physics simulations in the field of magnetic confinement fusion.
In such complex simulations, a simple interface for the efficient integration of the guiding-center equations is needed. Specifically, the initial position in five-dimensional phase-space is provided in terms of guiding-center coordinates (i.e., guiding-center position, parallel and perpendicular velocity) and the main interest is to retrieve the phase-space position after a prescribed time step. Such a pure “orbit time step routine” acting as an interface with a plasma physics simulation is provided.
The integration process itself, however, can also be of great interest. Therefore, a program allowing the detailed analysis of guiding-center orbits, the time evolution of their respective invariants of motion and Poincaré plots is also available. A unique feature of `GORILLA` is the direct way in which fluxes through cell boundaries can be directly obtained while conserving physical invariants. These properties result from the geometric approach of linear approximation in tetrahedal cells. Thus the computation of moments of the distribution function such as current densities from Monte-Carlo sampling can be much more efficient and stable than in usual orbit tracers. Furthermore, the mesh is agnostic to whether an orbit is inside or outside the magnetic separatrix, and does not rely on flux coordinates. Therefore `GORILLA` is a flexible tool that is useful in a variety of applications for non-axisymmetric devices, in particular stellarators and tokamaks with 3D perturbations.

`GORILLA` has already been used by @eder_quasi-geometric_2020 for the application of collisionless guiding-center orbits in an axisymmetric tokamak and a realistic three-dimensional stellarator configuration. There, the code demonstrated stable long-term orbit dynamics conserving invariants.
In the same publication, `GORILLA` was further applied to the Monte Carlo evaluation of transport coefficients. There, the computational efficiency of `GORILLA` was shown to be an order of magnitude higher than with a standard fourth order Runge–Kutta integrator.
Currently, `GORILLA` is part of the “EUROfusion Theory, Simulation, Validation and Verification Task for Impurity Sources, Transport, and Screening”, where it is tested for the kinetic modelling of the impurity ion component. First results of this research project and, in addition, `GORILLA`'s application to the
computation of fusion alpha particle losses in a realistic stellarator configuration have been reported by @eder_integration_2021.
The source code for `GORILLA` has been archived on Zenodo [@eder_gorilla_2021].


# Acknowledgements

The authors would like to thank Michael Drevlak and the ASDEX Upgrade Team for providing the stellarator field configuration and the ASDEX Upgrade MHD equilibrium for shot 26884 at 4300 ms which are used in the code examples.
Further thanks to Martin Heyn, Alex Runov, Philipp Ulbl, Rico Buchholz, Patrick Lainer, and Markus Richter for useful discussions.
This work has been carried out within the framework of the EUROfusion Consortium, funded by the European Union via the Euratom Research and Training Programme (Grant Agreement No 101052200 — EUROfusion). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Commission. Neither the European Union nor the European Commission can be held responsible for them.

![Illustration of (a) passing particle and (b) trapped particle guiding-center orbits of a Deuterium ion with a kinetic energy of 3 keV in the axisymmetric magnetic field configuration of ASDEX Upgrade. The blue transparent area shows the poloidal $\varphi = 0$ plane with blue dots indicating the intersections of the orbit with this plane (Poincaré cut).  Red solid lines represent the guiding-center orbits.\label{fig:example}](figure.png)

# References

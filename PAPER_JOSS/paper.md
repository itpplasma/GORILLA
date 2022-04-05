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
  - name: Christopher G. Albert
    orcid: 0000-0003-4773-416X
    affiliation: "1, 2"
  - name: Lukas M. P. Bauer
    orcid: 0000-0003-3341-4085
    affiliation: 1
  - name: Daniel Forstenlechner
  	affiliation: 1
  - name:  Sergei V. Kasilov
    affiliation: "1, 3, 4"
  - name:  Winfried Kernbichler
    affiliation: 1
  - name:  Markus Meisterhofer
    affiliation: 1
  - name:  Michael Scheidt
    affiliation: 1

affiliations:
 - name: Institut für Theoretische Physik - Computational Physics, Technische Universität Graz, Petersgasse 16, 8010 Graz, Austria
   index: 1
 - name: Max-Planck-Institut für Plasmaphysik, Boltzmannstr. 2, 85748 Garching, Germany
   index: 2
 - name: Institute of Plasma Physics, National Science Center, “Kharkov Institute of Physics and Technology,” Akademicheskaya str. 1, 61108 Kharkov, Ukraine
   index: 3
 - name: Department of Applied Physics and Plasma Physics, V. N. Karazin Kharkov National University, Svobody sq. 4, 61022 Kharkov, Ukraine
   index: 4
date: 16 March 2021
bibliography: paper.bib

---
# Introduction
Extremely hot plasmas with a temperature of the order of hundred million degrees Celsius are needed to produce energy from nuclear fusion. Under these conditions, hydrogen isotopes are fused, and energy is released. The energy release from 1 kg of fusion fuel corresponds approximately to that of 10000 tons of coal. A future use of this energy source is the subject of worldwide research projects. The confinement of such hot plasmas, however, poses major physical and technological problems for researchers. In particular, complex numerical methods are necessary to understand the physics of such plasmas in complicated toroidal magnetic fields.

An important kinetic approach for simulating the collective behavior of a plasma utilizes direct modeling of particle orbits. A well-known approximation for computing the motion of electrically charged particles in slowly varying electromagnetic fields is to reduce the dynamical equations by separating the relatively fast circular motion around a point called the guiding-center, and primarily treat the relatively slow drift motion of this point. This drift motion is described by the guiding-center equations; see, e.g., [@littlejohn_variational_1983], [@boozer_guiding_1980] and [@cary_hamiltonian_2009].

Here, we provide an efficient code for the purpose of solving the guiding-center equations. This code is a numerical implementation of the novel, quasi-geometric integration method described by @eder_quasi-geometric_2020.


# Summary

`GORILLA` is a Fortran code that computes guiding-center orbits for charged particles of given mass, charge and energy in toroidal fusion devices with three-dimensional field geometry. 
Conventional methods for integrating the guiding-center equations utilize high order interpolation of the electromagnetic field in space.
In `GORILLA`, a special linear interpolation employing a spatial mesh is used for the discretization of the electromagnetic field.
This leads to locally linear equations of motion with piecewise constant coefficients. 
As shown by @eder_quasi-geometric_2020, this local linearization approach retains the Hamiltonian structure of the guiding-center equations. For practical purposes this means that the total energy, the magnetic moment and the phase space volume are conserved.
Furthermore, the approach reduces computational effort and sensitivity to noise in the electromagnetic field. In `GORILLA` guiding-center orbits are computed without taking into account collisions in-between particles. Such exemplary guiding-center orbits obtained with `GORILLA` can be seen in \autoref{fig:example} where the magnetic field of a real-world fusion device is used, specifically the tokamak “ASDEX Upgrade”. 

# Statement of need

`GORILLA` is designed to be used by researchers in scientific plasma physics simulations in the field of magnetic confinement fusion. 
In such complex simulations a simple interface for the efficient integration of the guiding-center equations is needed. Specifically, the initial condition in five-dimensional phase space is provided (i.e. guiding-center position, parallel and perpendicular velocity) and the main interest is in the condition after a prescribed time step while the integration process itself is irrelevant. Such a pure “orbit time step routine” acting as an interface with a plasma physics simulation is provided.
The integration process itself, however, can also be of great interest and a program allowing the detailed analysis of guiding-center orbits, the time evolution of their respective invariants of motion and Poincaré plots is thus also provided.

`GORILLA` has already been used by @eder_quasi-geometric_2020 for the application of collisionless guiding-center orbits in an axisymmetric tokamak and a realistic three-dimensional stellarator configuration. There, the code demonstrated stable long-term orbit dynamics conserving invariants.
Further, in the same publication, `GORILLA` was applied to the Monte Carlo evaluation of transport coefficients. There, the computational efficiency of `GORILLA` was shown to be an order of magnitude higher than with a standard fourth order Runge–Kutta integrator.
Currently, `GORILLA` is part of the “EUROfusion Theory, Simulation, Validation and Verification Task for Impurity Sources, Transport, and Screening” where it is tested for the kinetic modelling of the impurity ion component. 
The source code for `GORILLA` has been archived on Zenodo with the linked DOI: [@eder_gorilla_2021]

![Illustration of (a) passing particle and (b) trapped particle guiding-center orbits of a Deuterium ion with a kinetic energy of 3 keV in the axisymmetric magnetic field configuration of ASDEX Upgrade. The blue transparent area shows the poloidal $\varphi = 0$ plane with blue dots indicating the intersections of the orbit with this plane (Poincaré cut).  Red solid lines represent the guiding-center orbits.\label{fig:example}](figure.png)

# Acknowledgements

The authors would like to thank Michael Drevlak and the ASDEX Upgrade Team for providing the stellarator field configuration and the ASDEX Upgrade MHD equilibrium for shot 26884 at 4300 ms which are used in the code examples.
Further thanks to Martin Heyn, Alex Runov, Philipp Ulbl, Rico Buchholz, Patrick Lainer, and Markus Richter for useful discussions.
This work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programmes 2014–2018 and 2019–2020 under Grant Agreement No. 633053.
The views and opinions expressed herein do not necessarily reflect those of the European Commission. The study was supported by the Reduced Complexity Models Grant No. ZT-I-0010 funded by the Helmholtz Association of German Research Centers. Support from NAWI Graz, and from the OeAD under the WTZ Grant Agreement with Ukraine No. UA 04/2017 is gratefully acknowledged.

# References

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
  - name:  Sergei V. Kasilov
    affiliation: "1, 3"
  - name:  Winfried Kernbichler
    affiliation: 1
  - name:  Markus Meisterhofer
    affiliation: 1

affiliations:
 - name: Institut für Theoretische Physik - Computational Physics, Technische Universität Graz, Petersgasse 16, 8010 Graz, Austria
   index: 1
 - name: Max-Planck-Institut für Plasmaphysik, Boltzmannstr. 2, 85748 Garching, Germany
   index: 2
 - name: Institute of Plasma Physics, National Science Center, “Kharkov Institute of Physics and Technology,” Akademicheskaya str. 1, 61108 Kharkov, Ukraine
   index: 3
date: 10 March 2021
bibliography: paper.bib

---
# Introduction
Extremely hot plasmas with a temperature of several million degrees Celsius are needed to produce energy from nuclear fusion. Under these conditions, hydrogen isotopes are fused and energy is released. The energy release from 1 kg of fusion fuel corresponds approximately to that of 10000 tons of coal. A future use of this energy source is the subject of worldwide research projects. However, the confinement of such hot plasmas poses major physical and technological problems for researchers. In particular, complex numerical methods are necessary to understand the physics of such plasmas in complicated toroidal magnetic fields.

The kinetic approach for simulating the collective behavior of a plasma utilizes direct modeling of particle orbits. A well-known approximation for computing the motion of electrically charged particles in slowly varying electromagnetic fields is to reduce the dynamical equations by separating the relatively fast circular motion around a point called the guiding-center, and instead treat the relatively slow drift motion of this point. This drift motion is described by the guiding-center equations, [@littlejohn_variational_1983] and [@boozer_guiding_1980].


# Summary


summary test
[@eder_quasi-geometric_2020]


The source code for `GORILLA` has been archived to Zenodo with the linked DOI: [@eder_gorilla_2021]

# Statement of need

statemend of need

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

citations

# Figures

figure

# Acknowledgements

The authors would like to thank Michael Drevlak and the ASDEX Upgrade Team for providing the stellarator field configuration and the ASDEX Upgrade MHD equilibrium for shot 26884 at 4300 ms.
Further thanks to Martin Heyn, Philipp Ulbl, Rico Buchholz, Patrick Lainer, and Markus Richter for useful discussions.
This work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programmes 2014–2018 and 2019–2020 under Grant Agreement No. 633053.
The views and opinions expressed herein do not necessarily reflect those of the European Commission. The study was supported by the Reduced Complexity Models Grant No. ZT-I-0010 funded by the Helmholtz Association of German Research Centers. Support from NAWI Graz, and from the OeAD under the WTZ Grant Agreement with Ukraine No. UA 04/2017 is gratefully acknowledged.

# References

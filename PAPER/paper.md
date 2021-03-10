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
  - name: Michael Eder^[Author to whom correspondence should be addressed: <eder@tugraz.at>]
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
 - name: Fusion@ÖAW, Institut für Theoretische Physik—Computational Physics, Technische Universität Graz, Petersgasse 16, 8010 Graz, Austria
   index: 1
 - name: Max-Planck-Institut für Plasmaphysik, Boltzmannstr. 2, 85748 Garching, Germany
   index: 2
 - name: Institute of Plasma Physics, National Science Center, “Kharkov Institute of Physics and Technology,” Akademicheskaya str. 1, 61108 Kharkov, Ukraine
   index: 3
date: 10 March 2021
bibliography: paper.bib

---

# Summary

summary

# Statement of need

[@eder_quasi-geometric_2020]

`GORILLA`

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

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures



# Acknowledgements

The authors would like to thank Michael Drevlak and the ASDEX Upgrade Team for providing the stellarator field configuration and the ASDEX Upgrade MHD equilibrium for shot 26884 at 4300 ms, respectively.
Further thanks to Martin Heyn, Philipp Ulbl, Rico Buchholz, Patrick Lainer, and Markus Richter for useful discussions.
This work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programmes 2014–2018 and 2019–2020 under Grant Agreement No. 633053.
The views and opinions expressed herein do not necessarily reflect those of the European Commission. The study was supported by the Reduced Complexity Models Grant No. ZT-I-0010 funded by the Helmholtz Association of German Research Centers. Support from NAWI Graz, and from the OeAD under the WTZ Grant Agreement with Ukraine No. UA 04/2017 is gratefully acknowledged.

# References

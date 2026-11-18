[![DOI](https://zenodo.org/badge/343767192.svg)](https://zenodo.org/badge/latestdoi/343767192)
[![](https://img.shields.io/github/actions/workflow/status/itpplasma/GORILLA/Ubuntu.yml?branch=main)](https://github.com/itpplasma/GORILLA/actions/workflows/Ubuntu.yml)
[![GitHub](https://img.shields.io/github/license/itpplasma/GORILLA?style=flat)](https://github.com/itpplasma/GORILLA/blob/main/LICENSE)
[![status](https://joss.theoj.org/papers/842b83a430e4c2985d202c79ed626587/status.svg)](https://joss.theoj.org/papers/842b83a430e4c2985d202c79ed626587)

# GORILLA
**G**uiding-center **OR**bit **I**ntegration with **L**ocal **L**inearization **A**pproach

![alt text](https://github.com/itpplasma/GORILLA/blob/main/DOCUMENTATION/LATEX/figures/title_image.png "Guiding-center orbit Poincaré plot with GORILLA")

GORILLA computes collisionless guiding-center orbits of charged particles in
toroidal fusion devices with three-dimensional field geometry. It is meant to be
used as the orbit integrator inside larger plasma simulations, for example
global kinetic calculations of quasi-steady plasma parameters or fast-particle
loss estimates.

The electromagnetic fields are linearly interpolated on a tetrahedral grid, so
the guiding-center equations of motion become locally linear with piecewise
constant coefficients and are solved exactly in each cell (Ref. [1]). This
conserves total energy, magnetic moment and phase-space volume. Because the
method works in general curvilinear coordinates, the computational domain can
include both the closed field line region and the scrape-off layer.

## Features

- Orbit time-step interface `orbit_timestep_gorilla`: given a guiding-center
  position and parallel and perpendicular velocity, return the state after a
  prescribed time step.
- Orbit analysis program `gorilla_plot`: Poincaré plots, full orbits and the
  time evolution of invariants of motion.
- Field input: axisymmetric g-file (EFIT) and 3D VMEC NetCDF equilibria,
  extended Boozer chartmap NetCDF files, SOLEDGE3X-EIRENE meshes for the edge,
  and an analytic circular tokamak.
- Orbit pushers within each cell: polynomial series solution of order
  K = 2, 3, 4 (with optional adaptive sub-steps to bound energy fluctuations),
  or numerical integration with Runge-Kutta 4 or adaptive ODE45.
- OpenMP parallelisation over particles on a single node.

## Building

Requirements: GNU Fortran, CMake, Ninja, NetCDF (C and Fortran), BLAS and
LAPACK. On Ubuntu:

```bash
sudo apt install git gfortran make cmake ninja-build liblapack-dev libnetcdff-dev
```

Then, from the repository root:

```bash
make            # configures with CMake and builds BUILD/test_gorilla_main.x
make test       # runs the registered CTest checks
```

See [docs/building.md](docs/building.md) for macOS, plain CMake, coverage
builds with pFUnit, and notes for the CINECA Leonardo cluster.

## Running an example

Eight examples are provided in `EXAMPLES/example_1` to `EXAMPLES/example_8`.
To build, run example 1 and plot the result with Python
(`pip install f90nml numpy matplotlib`):

```bash
make example EXAMPLE=1
```

Each example folder contains the input files and a link to the executable, so
an example can also be run directly with `./test_gorilla_main.x` in its folder.
Set `OMP_NUM_THREADS` to the number of physical cores.

## Documentation

- [docs/usage.md](docs/usage.md): input files, equilibria and grids.
- [docs/examples.md](docs/examples.md): the examples and the MATLAB and Python
  plotting tutorials.
- [docs/building.md](docs/building.md): build options, tests and coverage.
- [docs/contributing.md](docs/contributing.md): reporting bugs and contributing.
- `DOCUMENTATION/GORILLA_DOC.pdf`: detailed description of the method. Theses
  and the preprint of Ref. [1] are in `DOCUMENTATION/SUPPLEMENTAL_MATERIAL`.

## Contact

Questions: Michael Eder, <eder@tugraz.at>, or the
[issue tracker](https://github.com/itpplasma/GORILLA/issues).

## License

MIT License, see [LICENSE](LICENSE). The Runge-Kutta-Fehlberg routines in
`SRC/contrib/rkf45.f90` are from
https://people.sc.fsu.edu/~jburkardt/f_src/rkf45/rkf45.html under the GNU LGPL.
The polynomial root solver in `SRC/contrib/cmplx_roots_sg.f90` is under the
Apache 2.0 licence ([SRC/contrib/LICENSE.cmplx_roots_sg](SRC/contrib/LICENSE.cmplx_roots_sg)).

## Authors

Michael Eder, Christopher G. Albert, Lukas M. P. Bauer, Georg S. Graßler,
Daniel Forstenlechner, Sergei V. Kasilov, Winfried Kernbichler, Markus
Meisterhofer, Michael Scheidt (Fusion@OEAW, Institut für Theoretische Physik -
Computational Physics, Technische Universität Graz). S. V. Kasilov is also with
the Institute of Plasma Physics, National Science Center "Kharkov Institute of
Physics and Technology", and the Department of Applied Physics and Plasma
Physics, V. N. Karazin Kharkov National University.

## Acknowledgments
The authors would like to thank Michael Drevlak, the ASDEX Upgrade Team and the SOLEDGE3X-EIRENE Team for respectively providing the stellarator field configuration, the ASDEX Upgrade MHD equilibrium for shot 26884 at 4300 ms and the WEST MHD equilibrium for shot 54903 at 8 s.
Further thanks to Martin Heyn, Philipp Ulbl, Rico Buchholz, Patrick Lainer, and Markus Richter for useful discussions.
This work has been carried out within the framework of the EUROfusion Consortium, funded by the European Union via the Euratom Research and Training Programme (Grant Agreement No 101052200 — EUROfusion). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Commission. Neither the European Union nor the European Commission can be held responsible for them.


## References for GORILLA
When using this code for scientific publications, please cite both references [1] and [2]:

[1] M. Eder, C.G. Albert, L.M.P. Bauer, S.V. Kasilov and W. Kernbichler
“Quasi-geometric integration of guiding-center orbits in piecewise linear toroidal fields”
Physics of Plasmas 27, 122508 (2020)
<https://doi.org/10.1063/5.0022117>
Preprint: <https://arxiv.org/abs/2007.08151>

[2] M. Eder, C.G. Albert, L.M.P. Bauer, Georg S. Graßler, S.V. Kasilov, W. Kernbichler, M. Meisterhofer and M.Scheidt
“GORILLA: Guiding-center ORbit Integration with Local Linearization Approach”
submitted to Journal of Open Source Software
Preprint: <https://github.com/openjournals/joss-papers/blob/joss.03116/joss.03116/10.21105.joss.03116.pdf>


## References for provided MHD equilibria
[3] M. F. Heyn, I. B. Ivanov, S. V. Kasilov, W. Kernbichler, P. Leitner, and V. Nemov,
“Quasilinear modelling of RMP interaction with a tokamak plasma: application to ASDEX Upgrade ELM mitigation experiments”
Nucl. Fusion 54, 064005 (2014).
<https://doi.org/10.1088/0029-5515/54/6/064005>

[4] M. Drevlak, C. D. Beidler, J. Geiger, P. Helander, and Y. Turkin
“Quasi-Isodynamic Configuration with Improved Confinement”
41st EPS Conference on Plasma Physics ECA (2014), Vol. 38F, p. P1.070.
<http://ocs.ciemat.es/EPS2014PAP/pdf/P1.070.pdf>

[5] H. Bufferand et *al*
“Progress in edge plasma turbulence modelling—hierarchy of models from 2D transport application to 3D fluid simulations in realistic tokamak geometry”
Nucl. Fusion 61, 116052 (2021)
<https://doi.org/10.1088/1741-4326/ac2873>

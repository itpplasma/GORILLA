# GORILLA
**G**uiding-center **OR**bit **I**ntegration with **L**ocal **L**inearization **A**pproach

GORILLA computes guiding-center orbits for charged particles of given mass, charge and energy in toroidal fusion devices with three-dimensional field geometry.
The guiding-center orbits are traced via a quasi-geometric integration method described in Ref. [1].
There, high order interpolation of electromagnetic fields in space is replaced by a special linear interpolation, leading to locally linear Hamiltonian equations of motion with piecewise constant coefficients. The underlying formulation treats the motion in the piecewise linear fields exactly. This further leads to conservation of total energy, magnetic moment and phase space volume. Furthermore, the approach reduces computational effort and noise sensitivity.
Guiding-center orbits are computed without taking collisions into account.

For various simulations in magnetic confinement fusion, direct modeling of guiding-center particle orbits is utilized, e.g. global kinetic computations of quasi-steady plasma parameters or fast alpha particle loss estimation for stellarator optimization. In such complex simulations a simple interface for the guiding-center orbit integration part is needed. Namely, the initial condition of a particle in six-dimensional phase space is provided (e.g. particle position, magnetic moment, parallel and perpendicular velocity) and the main interest is in the condition after a prescribed time step while the integration process itself is irrelevant. Such a pure “orbit time step routine” acting as an interface with a plasma physics simulation is provided (“orbit_timestep_gorilla”).
However, the integration process itself can be of high interest as well, thus, a program allowing the detailed analysis of guiding-center orbits, the time evolution of their respective invariants of motion and Poincaré plots is at disposal as well (“gorilla_plot”).
Both applications are realized for demonstration in the program (“test_gorilla_main”).

The magnetic field can be provided by magnetohydrodynamics (MHD) equilibria with nested magnetic flux surfaces either in 2D (e.g. EFIT) or in 3D (e.g. VMEC). Supported equilibria are in the g-file or NetCDF format, respectively.
VMEC NetCDF equlibrium (`netcdf_file_for_test.nc`) was provided by Michael Drevlak for testing purposes and corresponds to the stellarator field configuration described in Ref. [3], namely, a quasi-isodynamic reactor-scale device with five toroidal field periods and a major radius of 25 m. 

The code is free to use and modify under the MIT License and links to Runge-Kutta-Fehlberg routines in
`SRC/contrib/rkf45.f90` from https://people.sc.fsu.edu/~jburkardt/f_src/rkf45/rkf45.html under the GNU LGPL License.

## Building

GORILLA can be built with make from Gorilla.mk.

Required libraries:
* NetCDF
* LAPACK/BLAS

Supported compilers:
* GNU Fortan 
* Intel Fortran

Include external library:
N. Flocke, “Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver for Physical Applications”
<https://doi.org/10.1145/2699468>

* Download supplemental material `954.zip` from above webpage.
* Copy `954/F90/Src/Polynomial234RootSolvers.f90` to `GORILLA/SRC/` and replace existing file with indentical name.
  (Existing file with identical name is a placeholder which is necessary for compilation.)


Building: 
```bash
cd /path/to/GORILLA
make -f Gorilla.mk
```
This will produce `test_gorilla_main.x` required to run the code.

## Usage

GORILLA currently runs on a single node with OpenMP shared memory parallelization with one particle per thread and background fields residing in main memory. 

The main executable is `test_gorilla_main.x`. 
As an input it takes .... 

... the following input files which can be found in the folder `INPUT/`
* `tetra_grid.inp`
* `gorilla.inp`
* `gorilla_plot.inp`
* `field_divB0.inp`
* `preload_for_SYNCH.inp`

... and the MHD equilibrium files which can be found in the folder `MHD_EQUILIBRIA/`
    * VMEC NetCDF equlibrium (File name can be specified in `tetra_grid.inp`.)
    * g-file equilibrium (File name can be specified in `tetra_grid.inp`.)

Five example input files with explanation of each parameter can be found in `RUN/example_1` - `RUN/example_5`.
After appropriate compilation, the code can be executed in all of these 5 example folders.

Moreover, a detailed explanation of all examples including the generation of the appropriate input files and plotting of the results with MATLAB can be found in the folder `matlab_test_cases`. There, a step-by-step tutorial in form of a MATLAB Live Script with the name `plot_poincare.mlx` is at disposal as well.


## References
When using this code for scientific publications, please cite the references [1] and [2]:

[1] M. Eder, C.G. Albert, L.M.P. Bauer, S.V. Kasilov and W. Kernbichler
“Quasi-geometric integration of guiding-center orbits in piecewise linear toroidal fields”
Physics of Plasmas 27, 122508 (2020)
<https://doi.org/10.1063/5.0022117>
Preprint: <https://arxiv.org/abs/2007.08151>

[2] M. Eder
“Placeholder for JOSS - GORILLA: Guiding-center ORbit Integration with Local Linearization Approach”


## References for MHD equilibria
[3] M. Drevlak, C. D. Beidler, J. Geiger, P. Helander, and Y. Turkin
“Quasi-Isodynamic Configuration with Improved Confinement”
41st EPS Conference on Plasma Physics ECA (2014), Vol. 38F, p. P1.070.
<http://ocs.ciemat.es/EPS2014PAP/pdf/P1.070.pdf>

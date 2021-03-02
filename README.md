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
* Copy `954/F90/Src/Polynomial234RootSolvers.f90` and `954/F90/Src/SetWorkingPrecision.f90` to `GORILLA/SRC/`. (Placeholder files should be re-placed with these two files from external source.)


Building: 
```bash
cd /path/to/GORILLA
make -f Gorilla.mk
```
This will produce `test_gorilla_main.x` required to run the code.

## Usage

GORILLA currently runs on a single node with OpenMP shared memory parallelization with one particle per thread and background fields residing in main memory. 

The main executable is `test_gorilla_main.x`. As an input it takes
* `tetra_grid.inp`
* `gorilla.inp`
* `gorilla_plot.inp`
* MHD equilibrium file
    * VMEC NetCDF equlibrium (wout.nc) file with name specified in `simple.in`
    * g-file equilibrium with name specified in 

Five example input files with explanation of each parameter can be found in `RUN/example_1` - `RUN/example_5`.
After appropriate compilation, the code can be executed in all of these 5 example folders.

Moreover, a detailed explanation of all examples including the generation of the appropriate input files and plotting of the results with MATLAB can be found in the folder `matlab_test_cases`. There, a step-by-step tutorial in form of a MATLAB Live Script with the name `plot_poincare.mlx` is at disposal as well.


## References
When using this code for scientific publications, please cite the according references:

[1] M. Eder, C.G. Albert, L.M.P. Bauer, S.V. Kasilov and W. Kernbichler
“Quasi-geometric integration of guiding-center orbits in piecewise linear toroidal fields”
Physics of Plasmas 27, 122508 (2020)
<https://doi.org/10.1063/5.0022117>

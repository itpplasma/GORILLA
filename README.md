# GORILLA
**G**uiding-center **OR**bit **I**ntegration with **L**ocal **L**inearization **A**pproach

GORILLA computes guiding-center orbits for charged particles of given mass, charge and energy in toroidal fusion devices with three-dimensional field geometry.
The guiding-center orbits are traced via a quasi-geometric integration method described in Ref. [1].
There, high order interpolation of electromagnetic fields in space is replaced by a special linear interpolation, leading to locally linear Hamiltonian equations of motion with piecewise constant coefficients. The underlying formulation treats the motion in the piecewise linear fields exactly. This further leads to conservation of total energy, magnetic moment and phase space volume. Furthermore, the approach reduces computational effort and noise sensitivity.
Guiding-center orbits are computed without taking collisions into account.

For various simulations in magnetic confinement fusion, direct modeling of guiding-center particle orbits is utilized, e.g. global kinetic computations of quasi-steady plasma parameters or fast alpha particle loss estimation for stellarator optimization. In such complex simulations a simple interface for the guiding-center orbit integration part is needed. Namely, the initial condition of a particle in six-dimensional phase space is provided (e.g. particle position, magnetic moment, parallel and perpendicular velocity) and the main interest is in the condition after a prescribed time step while the integration process itself is irrelevant. Such a pure “orbit time step routine” acting as an interface with a plasma physics simulation is provided (“orbit_timestep_gorilla”).
However, the integration process itself can be of high interest as well, thus, a program allowing the detailed analysis of guiding-center orbits, the time evolution of their respective invariants of motion and Poincaré plots is at disposal as well (“gorilla_plot”).
Both applications are realized for demonstration in the program (“test_gorilla_main”).

The code is free to use and modify under the MIT License and links to Runge-Kutta-Fehlberg routines in
`SRC/contrib/rkf45.f90` from https://people.sc.fsu.edu/~jburkardt/f_src/rkf45/rkf45.html under the GNU LGPL License.

### Magnetic field input
The magnetic field can be provided by magnetohydrodynamics (MHD) equilibria with nested magnetic flux surfaces either in 2D (e.g. EFIT) or in 3D (e.g. VMEC). Supported equilibria are in the g-file or NetCDF format, respectively.
For both equilibria formats, test files for the limited purpose of computing guiding-center orbits are provided. 
The g-file test equilibrium (`g_file_for_test`) was provided by the ASDEX Upgrade Team for testing purposes and corresponds to the axisymmetric tokamak field configuration of ASDEX Upgrade (shot 26884 at 4300 ms) described in Ref. [3]. 
The VMEC NetCDF test equlibrium (`netcdf_file_for_test.nc`) was provided by Michael Drevlak for testing purposes and corresponds to the stellarator field configuration described in Ref. [4], namely, a quasi-isodynamic reactor-scale device with five toroidal field periods and a major radius of 25 m. 

## Documentation
A detailed description of the working principle of GORILLA can be found in `DOCUMENTATION/GORILLA_DOC.pdf`.
The following supplemental material is available in `DOCUMENTATION/SUPPLEMENTAL_MATERIAL`:
* arXiv preprint of Ref. [1]
* Master's thesis of M. Eder (preliminary work with some detailed explanations referenced in `GORILLA_DOC.pdf`)
* Master's thesis of L. Bauer (preliminary work with some detailed explanations referenced in `GORILLA_DOC.pdf`)

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
* GORILLA can be run without this external library. The computation of guiding-center orbits is then limited to the numerical Runge-Kutta option of GORILLA.


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
* `tetra_grid.inp`                       (Input file for settings of the tetrahedronal grid used in GORLLA)
* `gorilla.inp`                             (Input file for settings of GORILLA)
* `gorilla_plot.inp`                   (Input file for the program for the analysis of guiding-center orbits)
* `field_divB0.inp`                     (Input file for loading g-file equilibria - Do not change this file.)
* `preload_for_SYNCH.inp`         (Input file for splining magnetic field data of g-file equilibria - Do not change this file.)

... and the MHD equilibrium files which can be found in the folder `MHD_EQUILIBRIA/`
* `netcdf_file_for_test.nc`: VMEC NetCDF equlibrium (File name can be specified in `tetra_grid.inp`.)
* `g_file_for_test`: g-file equilibrium (File name can be specified in `tetra_grid.inp`.)

## Examples

Expamles are realized in MATLAB and Python. For plotting with MATLAB, see below.
Five examples for plotting Poincaré cuts, full guiding-center orbits and the appropriate time evolution of invariants of motion can be found in `EXAMPLES/example_1` - `EXAMPLES/example_5`. There, the necessary soft links are already created and the input files are given. Detailed descriptions of the respective input files can be found in `INPUT`.
After appropriate compilation of GORILLA, the code can be executed in all of these 5 example folders, respectively.
For the visualization of the output of these five examples, appropriate plotting methods are at disposal at. `PYTHON/plot_example_1.py` - `PYTHON/plot_example_5.py`.

### Example 1
* Compute a collisionless guiding-center orbit with GORILLA for a trapped Deuterium particle.
* Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
* Use the GORILLA polynomial option with order K = 4.
* Create a figure with the Poincaré sections ($v_\parallel = 0$) in cylindrical and symmetry flux coordinates.
* Compute the normalized parallel adiabatic invariant as a function of banana bounces.

### Example 2
* Compute a collisionless guiding-center orbit with GORILLA for a passing Deuterium particle.
* Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
* Use the GORILLA Runge-Kutta 4 option.
* Create a figure with the Poincaré plots ($\varphi = 0$) in cylindrical and symmetry flux coordinates.
* Compute the normalized total energy as a function of toroidal mappings. 

### Example 3
* Compute a collisionless guiding-center orbit with GORILLA for a passing Deuterium particle.
* Use a field-aligned grid for an axisymmetric tokamak equilibrium (g-file)
* Use the GORILLA polynomial option with order K = 3.
* Create a figure with the Poincaré plots (\varphi = 0) in cylindrical and symmetry flux coordinates.
* Compute the normalized toroidal angular momentum as a function of toroidal mappings.

### Example 4
* Compute collisionless guiding-center orbits with GORILLA for two passing and one trapped Deuterium particle.
* Use a field-aligned grid for an axisymmetric tokamak equilibrium (g-file)
* Create a figure with the Poincaré plots ($\varphi = 0$) in cylindrical and symmetry flux coordinates.

### Example 5 (only MATLAB)
* Compute collisionless guiding-center orbits with GORILLA for a passing and a trapped Deuterium particle.
* Use a field-aligned grid for an axisymmetric tokamak equilibrium (g-file).
* Plot the plasma boundary, the guiding-center orbits, and the resulting Poincare plot ($\varphi = 0$) for both orbits.

### MATLAB: Generation of input files and plotting
A detailed explanation of all examples including the generation of the appropriate input files (including the example folders in `EXAMPLES/MATLAB_RUN`) and plotting of the results with MATLAB can be found in the folder `MATLAB`.
Here, the results of GORILLA with different polynominal orders K=2,3,4 and Runge-Kutta 4 are compared.

### MATLAB:  Step-by-step plotting tutorial
* MATLAB Live Script with the name `plotting_tutorial.mlx` is at disposal as a step-by-step tutorial for all plotting features of GORILL.


## Acknowledgments
The authors would like to thank Michael Drevlak and the ASDEX Upgrade Team for providing the stellarator field configuration and the ASDEX Upgrade MHD equilibrium for shot 26884 at 4300 ms.
Further thanks to Martin Heyn, Philipp Ulbl, Rico Buchholz, Patrick Lainer, and Markus Richter for useful discussions.
This work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programmes 2014–2018 and 2019–2020 under Grant Agreement No. 633053.
The views and opinions expressed herein do not necessarily reflect those of the European Commission. The study was supported by the Reduced Complexity Models Grant No. ZT-I-0010 funded by the Helmholtz Association of German Research Centers. Support from NAWI Graz, and from the OeAD under the WTZ Grant Agreement with Ukraine No. UA 04/2017 is gratefully acknowledged.


## References for GORILLA
When using this code for scientific publications, please cite both references [1] and [2]:

[1] M. Eder, C.G. Albert, L.M.P. Bauer, S.V. Kasilov and W. Kernbichler
“Quasi-geometric integration of guiding-center orbits in piecewise linear toroidal fields”
Physics of Plasmas 27, 122508 (2020)
<https://doi.org/10.1063/5.0022117>
Preprint: <https://arxiv.org/abs/2007.08151>

[2] M. Eder
“Placeholder for JOSS - GORILLA: Guiding-center ORbit Integration with Local Linearization Approach”


## References for provided MHD equilibria
[3] M. F. Heyn, I. B. Ivanov, S. V. Kasilov, W. Kernbichler, P. Leitner, and V. Nemov,
“Quasilinear modelling of RMP interaction with a tokamak plasma: application to ASDEX Upgrade ELM mitigation experiments”
Nucl. Fusion 54, 064005 (2014).
<https://doi.org/10.1088/0029-5515/54/6/064005>

[4] M. Drevlak, C. D. Beidler, J. Geiger, P. Helander, and Y. Turkin
“Quasi-Isodynamic Configuration with Improved Confinement”
41st EPS Conference on Plasma Physics ECA (2014), Vol. 38F, p. P1.070.
<http://ocs.ciemat.es/EPS2014PAP/pdf/P1.070.pdf>

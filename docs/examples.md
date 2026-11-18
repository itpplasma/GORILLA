# Examples and tutorials

## Running the examples

The folders `EXAMPLES/example_1` to `EXAMPLES/example_8` contain the input
files and a link to `BUILD/test_gorilla_main.x`. After building, run an
example in its folder with

```bash
./test_gorilla_main.x
```

or build, run, check the output files and plot in one step from the repository
root with

```bash
make example EXAMPLE=<number>
```

which calls `EXAMPLES/run_example.sh` and then `PYTHON/plot_example_<number>.py`.
The figure is written to `PYTHON/example_<number>.png`. Limit
`OMP_NUM_THREADS` to the number of physical cores to avoid hyperthreading
slowdowns.

The plotting scripts need Python 3 with [f90nml](https://github.com/marshallward/f90nml),
NumPy and Matplotlib (example 8 also uses SciPy):

```bash
pip install f90nml numpy matplotlib scipy
```

GORILLA itself does not need Python.

## Example descriptions

All examples compute collisionless guiding-center orbits of deuterium unless
stated otherwise.

1. Trapped particle, field-aligned grid for the VMEC stellarator equilibrium,
   polynomial pusher with K = 4. Poincaré sections at v_par = 0 in
   cylindrical and symmetry flux coordinates, and the normalised parallel
   adiabatic invariant over banana bounces.
2. Passing particle, VMEC equilibrium, Runge-Kutta 4 pusher. Poincaré plots
   at phi = 0 and the normalised total energy over toroidal mappings.
3. Passing particle, field-aligned grid for the axisymmetric ASDEX Upgrade
   g-file, polynomial pusher with K = 3. Poincaré plots at phi = 0 and the
   normalised toroidal angular momentum over toroidal mappings.
4. Two passing and one trapped particle, g-file equilibrium. Poincaré plots
   at phi = 0.
5. Trapped particle, g-file equilibrium. Plasma boundary, the orbit and its
   Poincaré plot at phi = 0.
6. Passing particles in the scrape-off layer of WEST, using a toroidal
   extension of the SOLEDGE3X-EIRENE 2D mesh. 2D projection of the orbits on
   the original mesh.
7. Trapped particle, VMEC equilibrium, adaptive scheme. Poincaré plots at
   phi = 0, the parallel adiabatic invariant over bounces, and the
   fluctuation and evolution of energy.
8. Analytic large-aspect-ratio circular tokamak (`grid_kind = 5`, no
   equilibrium file). Poincaré plot at phi = 0 of two deeply trapped, one
   marginally trapped and one passing orbit, compared with banana widths
   predicted from conservation of toroidal canonical momentum.

The MATLAB script `MATLAB/example_8.m` is a different example: it derives an
electric field from a potential proportional to the poloidal flux (ASDEX
Upgrade on the SOLEDGE3X-EIRENE mesh) and compares runs of fully ionised
tungsten with and without the additional terms of the strong electric field
Lagrangian. It plots Poincaré sections, the fluctuations of total energy and
toroidal momentum (each defined consistently with its Lagrangian), the
electric field and the ExB drift.

## Generating inputs and plotting in MATLAB and Python

The scripts `MATLAB/example_<n>.m` and `PYTHON/example_<n>.py` (examples 1-7)
generate the input files, run GORILLA in `EXAMPLES/MATLAB_RUN` or
`EXAMPLES/PYTHON_RUN` and plot the results. For examples 1-3 they compare
polynomial orders K = 2, 3, 4 with Runge-Kutta 4. For examples 5 and 6 they
compute both trapped and passing orbits. For example 7 they compare the
adaptive and non-adaptive schemes in detail.

## Tutorials

Step-by-step tutorials for all plotting features are available in both
languages:

- MATLAB Live Script: `MATLAB/plotting_tutorial.mlx`
- Jupyter Notebook: `PYTHON/plotting_tutorial.ipynb`

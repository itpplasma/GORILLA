# Using GORILLA

GORILLA runs on a single node. Particles are distributed over OpenMP threads,
one particle per thread, and the background fields are held in main memory.

## Programs and interfaces

- `orbit_timestep_gorilla` (module `orbit_timestep_gorilla_mod`,
  `SRC/orbit_timestep_gorilla.f90`): the interface for use inside a plasma
  simulation. It takes the guiding-center position and the parallel and
  perpendicular velocity and returns the state after a prescribed time step.
  The integration itself is hidden from the caller.
- `gorilla_plot` (`SRC/gorilla_plot_mod.f90`): integrates single orbits or
  sets of orbits and writes Poincaré sections, full orbits and the time
  evolution of the invariants of motion.
- `test_gorilla_main.x`: the main executable, which demonstrates both.

## Input files

The executable reads the following files from the working directory.
Templates with comments on every option are in `INPUT/`:

- `tetra_grid.inp`: tetrahedral grid (grid kind, resolution, equilibrium file
  names).
- `gorilla.inp`: GORILLA settings (pusher, polynomial order, time tracing,
  adaptive scheme, electric field options).
- `gorilla_plot.inp`: settings for `gorilla_plot`.
- `field_divB0.inp`: settings for loading g-file equilibria. Do not change.
- `preload_for_SYNCH.inp`: settings for splining g-file field data. Do not
  change.

## Magnetic field input and grids

The grid is selected with `grid_kind` in `tetra_grid.inp`:

1. rectangular grid for an axisymmetric EFIT equilibrium (g-file),
2. field-aligned grid for an axisymmetric EFIT equilibrium (g-file),
3. field-aligned grid for a non-axisymmetric VMEC equilibrium (NetCDF),
4. SOLEDGE3X-EIRENE grid (2D mesh extended toroidally),
5. analytic circular tokamak (no input file),
6. field-aligned grid for an extended Boozer chartmap (NetCDF).

The extended Boozer chartmap format contains `rho`, `s`, `theta`, `zeta`,
Cartesian geometry, `A_phi`, `B_theta`, `B_phi`, `Bmod`, `num_field_periods`
and the global `torflux` attribute, in the chartmap units used by
[SIMPLE](https://github.com/itpplasma/SIMPLE).

## Test equilibria

`MHD_EQUILIBRIA/` contains equilibria for the limited purpose of computing
guiding-center orbits. File names are set in `tetra_grid.inp`.

- `g_file_for_test`: ASDEX Upgrade, shot 26884 at 4300 ms, provided by the
  ASDEX Upgrade Team (Ref. [3] in the README).
- `netcdf_file_for_test.nc`: VMEC equilibrium of a quasi-isodynamic
  reactor-scale stellarator with five field periods and 25 m major radius,
  provided by Michael Drevlak (Ref. [4]).
- `g_file_for_test_WEST`: WEST, shot 54903 at 8 s, provided by the
  SOLEDGE3X-EIRENE Team (Ref. [5]).
- `convex_wall_for_test.dat`, `convex_wall_for_test_WEST.dat`: convex walls
  for the two g-files.

## SOLEDGE3X-EIRENE meshes

For WEST geometry (`grid_kind = 4`) the original 2D mesh is read from
`MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE`:

- `knots_for_test.dat`: coordinates (R, Z) of the mesh vertices,
- `triangles_for_test.dat`: vertex triples forming the triangles.

Both file names are set in `tetra_grid.inp`. MATLAB scripts to produce these
files and the g-file from SOLEDGE3X-EIRENE output are in
`PREPROCESSING/SOLEDGE3X_EIRENE/MESH` and
`PREPROCESSING/SOLEDGE3X_EIRENE/MHD_EQUILIBRIUM`.

## Further reading

`DOCUMENTATION/GORILLA_DOC.pdf` describes the method and the code in detail.
`DOCUMENTATION/SUPPLEMENTAL_MATERIAL` contains the arXiv preprint of Ref. [1],
the Master's theses of M. Eder and L. Bauer (referenced in `GORILLA_DOC.pdf`)
and the Bachelor's thesis of D. Forstenlechner on unit testing with pFUnit.

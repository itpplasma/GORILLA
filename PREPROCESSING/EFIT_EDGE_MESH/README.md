# EFIT mesh extended to the limiter (grid_kind 4 + EFIT)

Extends the flux-aligned `grid_kind = 2` mesh of an EFIT equilibrium (closed surfaces up
to psi_N ~ 0.9995) out to the limiter contour stored in the g-file.

* Core: the `grid_kind = 2` vertices and triangles are kept unchanged.
* Edge: vertices lie on psi = const contours (SOL, separatrix legs, private-flux region)
  and on the limiter. They are connected by a constrained Delaunay triangulation
  (Shewchuk's Triangle), using the contours as constraint segments, so triangle edges
  follow flux surfaces wherever a contour exists.
* The result is written in the knots/triangles format of the `grid_kind = 4` reader and
  is used with `boole_ext_mesh_efit = .true.` (EFIT g-file instead of the WEST format).
  Only `coord_system = 1` (cylindrical) is possible: there are no flux coordinates
  outside the separatrix.

## Important: use `theta_geom_flux = 2` for the core

With `theta_geom_flux = 1`, the `grid_kind = 2` points crowd toward the X-point near
the edge. The ring-to-ring connection then FOLDS in (R,Z): there are 97 folded edges
for 30x60 and 232 for the example_4 grid (100x30). GORILLA's own check reports
"Some normal vectors are invalid" when such a grid is used in cylindrical coordinates.
With `theta_geom_flux = 2` the core has no folds. The generator refuses a folded core;
`--core cdt` re-triangulates the core instead and keeps its vertices.

## Recipe

1. Build GORILLA, then compile the dump driver against it:
   `gfortran -O2 $GORILLA_INC -I BUILD/SRC dump_mesh_2d.f90 -o dump_mesh_2d.x BUILD/SRC/libGORILLA.a $GORILLA_LINK -lnetcdff -lnetcdf -llapack -fopenmp`
2. In a run directory with `field_divB0.inp` and `preload_for_SYNCH.inp` (as in
   EXAMPLES/example_4) and a `tetra_grid.inp` with `grid_kind = 2`, `theta_geom_flux = 2`,
   and the n1/n3 wanted for the core (n2 does not matter), run `dump_mesh_2d.x`. It writes
   `verts_2d.dat` and `tets_slice1.dat`.
3. `python aug_edge_mesh.py DUMPDIR GFILE OUTDIR [--plot]`. Needs numpy, scipy,
   matplotlib, contourpy and triangle (pip). It writes `knots_aug_edge.dat`,
   `triangles_aug_edge.dat`, `triangle_region.dat` (0 = core, 1 = edge) and
   `mesh_report.txt`. It checks for folds, non-manifold edges, and that the total area
   equals the limiter area minus the core hole.
4. tetra_grid.inp for GORILLA (file names are limited to 64 characters):
   ```
   grid_kind = 4 , boole_ext_mesh_efit = .true. , n2 = <toroidal slices> ,
   g_file_filename = 'MHD_EQUILIBRIA/g_file_for_test' ,
   knots_SOLEDGE3X_EIRENE_filename = 'mesh/knots_aug_edge.dat' ,
   triangles_SOLEDGE3X_EIRENE_filename = 'mesh/triangles_aug_edge.dat' ,
   ```
   and gorilla.inp with `coord_system = 1`.

Options: `--d1/--growth/--dmax` set the contour spacing. The first SOL contour lies
d1 from the last ring at the outboard midplane; PFR contours are spaced downward from
the X-point. `--hmin/--hmax/--hgrad` set the vertex spacing along contours, finest at
the X-point. `--minangle` enables Triangle's quality refinement. Its Steiner points
are not exactly on contours, so it is off by default.

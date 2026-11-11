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

## Use from a grid_kind 2 run (recommended)

Set in `tetra_grid.inp` (with `coord_system = 1` in `gorilla.inp`):

    grid_kind = 2 , n1 = 100 , n2 = 60 , n3 = 120 , sfc_s_min = 1.d-3 ,
    theta_geom_flux = 2 , theta0_at_xpoint = .true. ,
    theta_xpoint_clustering = 0.9 , theta_xpoint_clustering_width = 0.1 ,   ! optional
    boole_extend_to_wall = .true. ,
    edge_mesh_dir = 'mesh' ,
    edge_mesh_options = '--fill-axis --xref-radius 12' ,

The first run stops and prints the `make_edge_mesh.sh` command that builds a matching mesh into
`edge_mesh_dir` (about 10 s). Afterwards every run checks that `edge_mesh_dir/edge_mesh.nml` matches
the grid_kind 2 settings, the g-file (first line and size) and `edge_mesh_options`, and uses the mesh
through the external-mesh path (internally grid_kind 4 with the EFIT g-file). Several runs can share
one mesh directory. In the run, `tetra_grid_mod` provides `verts_s_theta_core` (s, theta of the core
vertices, -1 outside) and `tetra_in_core` (tetrahedra of the flux-aligned core).

Files in the mesh directory: `knots_aug_edge.dat`, `triangles_aug_edge.dat` (mesh),
`edge_mesh_core.dat` (core s, theta and triangle regions), `edge_mesh.nml` (settings record),
`mesh_report.txt`, `generator.log`, `core/` (the grid_kind 2 core run).

## One step

`make_edge_mesh.sh OUTDIR N1 N3 [SFC_S_MIN] [-- generator options]` builds the core with
GORILLA (grid_kind 2, theta_geom_flux 2) and runs the generator. It needs `dump_mesh_2d.x`
(built with GORILLA as `BUILD/dump_mesh_2d.x`; other path via `DUMP_X`) and a python with the packages below (`PYTHON`). Example:

    make_edge_mesh.sh mesh_100x120 100 120 1e-3 -- --fill-axis --xref-radius 12 --plot

Timing for AUG g26884 on one core: 30x60 3 s + 3 s, 100x120 5 s + 5 s, 200x240 4 s + 10 s
(core + edge generator).

## Recipe (individual steps)

1. Build GORILLA; this also builds the dump driver `BUILD/dump_mesh_2d.x`.
2. In a run directory with `field_divB0.inp` and `preload_for_SYNCH.inp` (as in
   EXAMPLES/example_4) and a `tetra_grid.inp` with `grid_kind = 2`, `theta_geom_flux = 2`,
   and the n1/n3 wanted for the core (n2 does not matter), run `dump_mesh_2d.x`. It writes
   `verts_2d.dat`, `tets_slice1.dat` and `core_settings.nml`.
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

Options:
* `--d1/--growth/--dmax`: contour spacing. The first SOL contour lies d1 from the last ring
  at the outboard midplane (default: the last core ring spacing); the spacing then grows
  geometrically. PFR levels are spaced from the X-point down a vertical line.
* `--leg-aspect` (0.5), `--leg-cuts`: on cuts normal to each separatrix leg, neighbouring
  levels are at least `max(d_k, leg_aspect*h)` apart. This stops SOL/PFR contours from
  bunching along the legs.
* `--gap-min` (0.07): each contour keeps a gap of at least `gap_min*h` to its neighbour level
  everywhere along it (about 8 deg minimum angle between contours).
* `--hmin/--hmax/--hgrad`: vertex spacing along contours, finest at the X-point.
* `--snap` (0.3): contour end points closer than `snap*h` to a limiter corner are moved onto it.
* `--xref-radius R`, `--xref-h`: X-point refinement. Free points (not on contours) on circles
  around the X-point with spacing `xref_h + hgrad*r` for r < R, at least half a spacing away
  from every contour and vertex. The core ring is never split, so the refinement stops at the
  core boundary.
* `--fill-axis`: closes the hole inside the innermost grid_kind 2 ring (axis point plus
  triangles, region 2). Only meaningful in cylindrical coordinates. With a small `sfc_s_min`
  (e.g. 1e-3) GORILLA's `n_extra_rings` keep the near-axis rings fine and the filled hole is
  tiny.
* `--minangle`: Triangle's global quality refinement. Its Steiner points are not on contours,
  so it is off by default.

## Tests

`ctest` in the build directory runs `efit_grid_*` (driver `SRC/TESTS/test_efit_grid_driver.f90`, AUG test
equilibrium): the stop for folded grid_kind 2 grids (both coordinate systems), the X-point clustering, and the
missing-mesh stop. If CMake finds a python with numpy, scipy, matplotlib, contourpy and triangle
(choose it with `-DPython3_EXECUTABLE=...`), it also builds a small extended mesh with `make_edge_mesh.sh`
and checks the run with `boole_extend_to_wall` (core flags, valid faces, magnetic axis from the g-file) and
the stop for a mesh built with other settings. About 20 s in total.

## License note

The generator uses the `triangle` python package, a wrapper of J. R. Shewchuk's Triangle. Triangle is free
for research and other non-commercial use; commercial use needs a license from its author. It is not part
of GORILLA and is installed separately (`pip install triangle`).


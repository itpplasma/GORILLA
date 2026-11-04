# JOREK preprocessing

`jorek2gorilla.py` converts an axisymmetric JOREK **model-711** restart
(`jorek??????.h5`) into a netCDF file with a triangular poloidal mesh and all
quantities GORILLA needs at the mesh vertices.

The JOREK finite elements are evaluated with Haowei Zhang's reader
[`plot_jorek_model711_h5.py`](https://github.com/changhw/plot_jorek_model711_h5),
which is imported from its own directory, not copied here.

## Usage

```bash
python jorek2gorilla.py jorek050000.h5 jorek2D.nc --n-div 4 \
    --reader-dir /path/to/plot_jorek_model711_h5
```

Requirements: `numpy`, `scipy`, `h5py`, `matplotlib` (plus `tqdm`, which the reader imports).

| Option | Meaning |
|---|---|
| `--n-div N` | subdivide every JOREK element N x N times in its local coordinates (default 2) |
| `--reader-dir DIR` | directory of the reader (default: `$JOREK711_READER`) |
| `--n0-only` | take only the n = 0 part of a 3D restart |
| `--fine-div`, `--grid` | resolution of the integration that gives A_Z (defaults are fine) |

## Use in GORILLA (`grid_kind = 7`)

```fortran
&TETRA_GRID_NML
grid_kind = 7 ,
vertex_data_filename = 'jorek2D.nc' ,
n2 = 60 ,               ! number of toroidal slices
boole_n_field_periods = .true. ,
...
/
```

GORILLA reads the triangles and the vertex data (`SRC/vertex_data_mod.f90`),
extrudes the triangles toroidally into prisms and splits them into tetrahedra
(the same algorithm as for the SOLEDGE3X-EIRENE grid). Only cylindrical
coordinates (`coord_system = 1`) are supported. The electrostatic potential is
`Phi` from the file plus the usual `eps_Phi * A_phi`.

## Mesh

Each JOREK element is subdivided in (s,t) and every sub-quadrilateral is split
into two triangles. Coinciding points on shared element edges and the collapsed
points at the magnetic axis are merged, and degenerate triangles at the axis are
dropped. All values are evaluated exactly at the vertices with the JOREK
finite-element representation. Where the geometry Jacobian is singular (at the
axis), B is evaluated at a local offset of 1e-5.

## Output file

CGS units, temperatures in eV, **right-handed** (R, phi, Z). JOREK uses
x = R cos(phi), y = -R sin(phi), so every phi component is negated.

| Variable | Units | Content |
|---|---|---|
| `triangles` | | one-based vertex indices, counterclockwise in (R,Z) |
| `R`, `Z` | cm | vertex positions |
| `A_R`, `A_phi`, `A_Z` | G cm, G cm^2, G cm | covariant vector potential, B = curl A |
| `B_R`, `B_phi`, `B_Z`, `Bmod` | G | physical magnetic field |
| `Phi` | statV | electrostatic potential (not stored by model 711; currently zero) |
| `n` | cm^-3 | ion number density (rho n0) |
| `T_i`, `T_e` | eV | temperatures |
| `u_R`, `u_phi`, `u_Z` | cm/s | fluid velocity |
| `mesh_label` | 1 | radial ring index of the vertex; orients the prisms in GORILLA |
| `psi_eq` | G cm^2 | initial equilibrium flux, same sign convention as `A_phi` |

Global attributes record the source file and its SHA-256, the time, the
normalisation, the options and the command.

### Gauge of the vector potential

`A_phi = -A3` (poloidal flux). In the axisymmetric case JOREK's `A_R`, `A_Z`
and the equilibrium `F/R` contribute only to B_phi. JOREK's `A_R` and `A_Z` have
large, oscillating gauge parts near the edge that a piecewise-linear
interpolation cannot represent, so they are replaced by the equivalent gauge
`A_R = 0`, `dA_Z/dR = -B_phi`, integrated outward from the magnetic axis.

## Checks printed by the converter

* number of vertices against the structured-grid expectation, spread of the
  values at merged vertices (round-off for all stored fields; ~1e-3 for B,
  which comes from first derivatives);
* curl of the piecewise-linear A per triangle against the vertex-averaged B.
  The median error converges at first order in 1/n_div. The largest errors
  appear where the JOREK solution varies on the scale of the mesh, typically
  near the boundary.

### Why `mesh_label`

GORILLA orients each prism by a label at the triangle vertices: the two vertices
with the closest label form the base. With the poloidal flux as label this fails
where the flux surfaces have moved away from the JOREK grid during the JOREK
run (some prisms then cannot be connected). The exact ring index of the
polar JOREK grid makes the classification unambiguous. For a non-polar grid the
converter falls back to `psi_eq`.

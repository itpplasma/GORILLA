#!/usr/bin/env python3
"""Convert an axisymmetric JOREK model-711 restart into a GORILLA vertex-data file.

The poloidal plane is triangulated by subdividing every JOREK Bezier element
n_div x n_div times in its local coordinates (s,t) and splitting each
sub-quadrilateral into two triangles.  All quantities are evaluated exactly at
the triangle vertices with the JOREK finite-element representation, using
Haowei Zhang's reader plot_jorek_model711_h5.py (imported, not copied).

Conventions of the output file (GORILLA conventions):
  * CGS units (cm, G, G cm, G cm^2, cm^-3, cm/s); temperatures in eV.
  * Right-handed cylindrical coordinates (R, phi, Z).  JOREK uses
    x = R cos(phi), y = -R sin(phi), so every phi component changes sign.
  * A_R, A_Z, A_phi are covariant components of the vector potential with
    B = curl A.  A_phi = -JOREK A3 (poloidal flux).  In the axisymmetric case
    JOREK's A_R, A_Z and the equilibrium F/R only provide B_phi; they are
    replaced by the equivalent gauge A_R = 0, dA_Z/dR = -B_phi, which is
    smooth (JOREK's A_R, A_Z have large, oscillating gauge parts at the edge
    that a piecewise-linear interpolation cannot represent).
  * B_R, B_phi, B_Z are physical components.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import importlib
import os
import sys
from pathlib import Path

import h5py
import matplotlib.tri as mtri
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf_file
from scipy.spatial import cKDTree

CM_PER_M = 1.0e2
GAUSS_PER_TESLA = 1.0e4

# Offset from the collapsed axis side of the axis elements, used where the
# geometry Jacobian is singular and derivatives cannot be evaluated exactly.
AXIS_OFFSET_S = 1.0e-5

# Vertices closer than this (in m) are treated as the same point.
MERGE_TOLERANCE = 1.0e-9


def import_reader(reader_dir: Path):
    sys.path.insert(0, str(reader_dir))
    return importlib.import_module("plot_jorek_model711_h5")


def local_samples(n_div: int) -> tuple[np.ndarray, np.ndarray]:
    """Local (s,t) of the subdivision vertices, s-major ordering."""
    sample = np.linspace(0.0, 1.0, n_div + 1)
    s, t = np.meshgrid(sample, sample, indexing="ij")
    return s.ravel(), t.ravel()


def evaluate_elements(reader, restart, fourier_basis, s, t) -> dict[str, np.ndarray]:
    """Evaluate geometry and all reader fields at local (s,t) of every element.

    Returns arrays of shape (element, point).
    """
    poloidal_basis = reader.basis_functions(restart["order"], s, t)
    geometry = reader.calculate_geometry(restart, poloidal_basis, fourier_basis)
    fields = reader.calculate_fields(restart, poloidal_basis, fourier_basis, geometry)
    shape = geometry["R"].shape
    result = {name: value.reshape(shape) for name, value in fields.items()}
    result["R"] = geometry["R"]
    result["Z"] = geometry["Z"]
    return result


def evaluate_with_axis_fix(reader, restart, fourier_basis, n_div):
    """Evaluate on the subdivision vertices; repair B where the Jacobian is singular.

    At the magnetic axis the axis elements collapse one side, so the
    derivatives (and thus B) are undefined there.  B is then taken at a small
    offset AXIS_OFFSET_S into the element instead.
    """
    s, t = local_samples(n_div)
    values = evaluate_elements(reader, restart, fourier_basis, s, t)

    singular = ~np.isfinite(values["BR"])
    if np.any(singular):
        s_shifted = np.where(s == 0.0, AXIS_OFFSET_S, s)
        s_shifted = np.where(s == 1.0, 1.0 - AXIS_OFFSET_S, s_shifted)
        shifted = evaluate_elements(reader, restart, fourier_basis, s_shifted, t)
        for name in ("BR", "BZ", "B_phi"):
            values[name] = np.where(singular, shifted[name], values[name])
    n_singular = int(np.count_nonzero(singular))
    n_unrepaired = int(np.count_nonzero(~np.isfinite(values["BR"])))
    return values, n_singular, n_unrepaired


def merge_vertices(R: np.ndarray, Z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Map every sample to a unique vertex.  Returns (owner index, representative sample)."""
    points = np.column_stack((R, Z))
    pairs = cKDTree(points).query_pairs(MERGE_TOLERANCE, output_type="ndarray")

    parent = np.arange(len(points))

    def find(i):
        root = i
        while parent[root] != root:
            root = parent[root]
        while parent[i] != root:
            parent[i], i = root, parent[i]
        return root

    for i, j in pairs:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    roots = np.array([find(i) for i in range(len(points))])
    representative, vertex_of_sample = np.unique(roots, return_inverse=True)
    return vertex_of_sample, representative


def build_triangles(n_elements: int, n_div: int, vertex_of_sample: np.ndarray) -> np.ndarray:
    """Two triangles per sub-quadrilateral; degenerate triangles at the axis are dropped."""
    n_side = n_div + 1
    i, j = np.meshgrid(np.arange(n_div), np.arange(n_div), indexing="ij")
    p00 = (i * n_side + j).ravel()
    p10 = ((i + 1) * n_side + j).ravel()
    p11 = ((i + 1) * n_side + j + 1).ravel()
    p01 = (i * n_side + j + 1).ravel()
    local = np.concatenate(
        (np.column_stack((p00, p10, p11)), np.column_stack((p00, p11, p01)))
    )
    offsets = np.arange(n_elements)[:, None, None] * n_side**2
    samples = (local[None, :, :] + offsets).reshape(-1, 3)
    triangles = vertex_of_sample[samples]
    distinct = (
        (triangles[:, 0] != triangles[:, 1])
        & (triangles[:, 1] != triangles[:, 2])
        & (triangles[:, 0] != triangles[:, 2])
    )
    return triangles[distinct]


def ring_labels(n_elements: int, n_div: int, vertex_of_sample: np.ndarray, n_vertices: int):
    """Radial ring index of every vertex for a polar JOREK grid, or None.

    Elements whose s=0 side collapses onto one vertex touch the axis (ring
    index 0). The radial neighbour of an element is the one whose s=0 side
    consists of the same vertices as its s=1 side. GORILLA uses the label to
    orient the prisms: the two vertices of a triangle on the same ring form
    its base.
    """
    sides = vertex_of_sample.reshape(n_elements, n_div + 1, n_div + 1)
    inner_side = {frozenset(sides[e, 0, :]): e for e in range(n_elements)}
    element_ring = np.full(n_elements, -1)
    frontier = [e for e in range(n_elements) if len(set(sides[e, 0, :])) == 1]
    element_ring[frontier] = 0
    while frontier:
        following = []
        for e in frontier:
            outer = inner_side.get(frozenset(sides[e, n_div, :]))
            if outer is not None and element_ring[outer] < 0:
                element_ring[outer] = element_ring[e] + 1
                following.append(outer)
        frontier = following
    if np.any(element_ring < 0):
        return None
    sample_ring = (element_ring[:, None, None] * n_div + np.arange(n_div + 1)[None, :, None]
                   + np.zeros((1, 1, n_div + 1), dtype=int)).ravel()
    label = np.full(n_vertices, -1)
    label[vertex_of_sample] = sample_ring
    if np.any(label[vertex_of_sample] != sample_ring):
        return None
    return label.astype(float)


def signed_areas(R, Z, triangles):
    r = R[triangles]
    z = Z[triangles]
    return 0.5 * ((r[:, 1] - r[:, 0]) * (z[:, 2] - z[:, 0]) - (r[:, 2] - r[:, 0]) * (z[:, 1] - z[:, 0]))


def toroidal_field_potential(R_vert, Z_vert, F_fine, R_fine, Z_fine, tri_fine, F_ref, R_ref, n_grid):
    """Covariant A_Z (T m) with A_R = 0 whose curl is the toroidal field F/R.

    F = R B_phi is given with JOREK's left-handed phi.  In right-handed
    coordinates B_phi = dA_R/dZ - dA_Z/dR = -F/R, hence with A_R = 0
        A_Z = F_ref ln R + int_{R_ref}^{R} (F - F_ref)/R' dR'.
    The deviation F - F_ref is interpolated linearly from a fine sample mesh
    onto a rectangular grid (zero outside the domain) and integrated along R.
    Integrating outward from R_ref (the magnetic axis) keeps the noisy edge
    values from entering A_Z in the interior: a Z-dependent offset collected
    at the edge would not change dA_Z/dR, but piecewise-linear triangles
    alias it into a spurious B_phi.
    """
    interpolator = mtri.LinearTriInterpolator(
        mtri.Triangulation(R_fine, Z_fine, tri_fine), F_fine - F_ref
    )
    margin = 0.01
    R_axis = np.linspace(R_fine.min() - margin, R_fine.max() + margin, n_grid[0])
    Z_axis = np.linspace(Z_fine.min() - margin, Z_fine.max() + margin, n_grid[1])
    RR, ZZ = np.meshgrid(R_axis, Z_axis, indexing="ij")
    deviation = np.ma.filled(interpolator(RR, ZZ), 0.0)
    integrand = deviation / RR
    integral = np.concatenate(
        (
            np.zeros((1, len(Z_axis))),
            np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(R_axis)[:, None], axis=0),
        )
    )
    for k in range(len(Z_axis)):
        integral[:, k] -= np.interp(R_ref, R_axis, integral[:, k])
    integral_at_vertices = RegularGridInterpolator((R_axis, Z_axis), integral)(
        np.column_stack((R_vert, Z_vert))
    )
    return F_ref * np.log(R_vert) + integral_at_vertices


def curl_check(R, Z, triangles, A_R, A_phi, A_Z, B_R, B_phi, B_Z):
    """Compare curl of the piecewise-linear A with the vertex-averaged B per triangle."""
    r = R[triangles]
    z = Z[triangles]
    det = (r[:, 1] - r[:, 0]) * (z[:, 2] - z[:, 0]) - (r[:, 2] - r[:, 0]) * (z[:, 1] - z[:, 0])

    def gradient(f):
        fv = f[triangles]
        d_R = ((fv[:, 1] - fv[:, 0]) * (z[:, 2] - z[:, 0]) - (fv[:, 2] - fv[:, 0]) * (z[:, 1] - z[:, 0])) / det
        d_Z = ((r[:, 1] - r[:, 0]) * (fv[:, 2] - fv[:, 0]) - (r[:, 2] - r[:, 0]) * (fv[:, 1] - fv[:, 0])) / det
        return d_R, d_Z

    R_c = r.mean(axis=1)
    dAphi_dR, dAphi_dZ = gradient(A_phi)
    dAR_dR, dAR_dZ = gradient(A_R)
    dAZ_dR, dAZ_dZ = gradient(A_Z)
    curl = {
        "B_R": -dAphi_dZ / R_c,
        "B_phi": dAR_dZ - dAZ_dR,
        "B_Z": dAphi_dR / R_c,
    }
    reference = {"B_R": B_R, "B_phi": B_phi, "B_Z": B_Z}
    B_pol = np.sqrt(B_R**2 + B_Z**2)[triangles].mean(axis=1)
    B_mod = np.sqrt(B_R**2 + B_phi**2 + B_Z**2)[triangles].mean(axis=1)
    report = {}
    for name in curl:
        difference = np.abs(curl[name] - reference[name][triangles].mean(axis=1))
        scale = B_mod if name == "B_phi" else B_pol
        relative = difference / scale
        report[name] = (np.median(relative), np.percentile(relative, 99), relative.max())
    return report


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_netcdf(path, vertices, triangles, attributes):
    with netcdf_file(path, "w", version=2) as nc:
        nc.createDimension("n_vertices", len(vertices["R"][0]))
        nc.createDimension("n_triangles", len(triangles))
        nc.createDimension("three", 3)
        variable = nc.createVariable("triangles", "i4", ("n_triangles", "three"))
        variable[:] = triangles + 1
        variable.description = "one-based vertex indices, counterclockwise in (R,Z)"
        for name, (data, units, description) in vertices.items():
            variable = nc.createVariable(name, "f8", ("n_vertices",))
            variable[:] = data
            variable.units = units
            variable.description = description
        for name, value in attributes.items():
            setattr(nc, name, value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("restart", type=Path, help="JOREK model-711 restart (jorek??????.h5)")
    parser.add_argument("output", type=Path, help="output netCDF file")
    parser.add_argument("--n-div", type=int, default=2, help="subdivisions per element direction (default 2)")
    parser.add_argument(
        "--reader-dir",
        type=Path,
        default=os.environ.get("JOREK711_READER"),
        help="directory containing plot_jorek_model711_h5.py (default: $JOREK711_READER)",
    )
    parser.add_argument(
        "--n0-only",
        action="store_true",
        help="use only the n=0 part of a 3D restart (otherwise n_tor must be 1)",
    )
    parser.add_argument("--fine-div", type=int, default=8, help="subdivisions of the mesh used for R B_phi (default 8)")
    parser.add_argument("--grid", type=int, nargs=2, default=(1000, 1400), metavar=("NR", "NZ"),
                        help="rectangular grid for integrating B_phi (default 1000 1400)")
    args = parser.parse_args()

    if args.reader_dir is None:
        parser.error("give --reader-dir or set JOREK711_READER")
    if args.n_div < 1:
        parser.error("--n-div must be at least 1")
    reader = import_reader(args.reader_dir)

    restart = reader.load_restart(args.restart)
    if restart["n_tor"] > 1 and not args.n0_only:
        sys.exit(f"{args.restart} has n_tor={restart['n_tor']}; use --n0-only for its axisymmetric part")
    fourier_basis = (np.eye(restart["n_tor"])[0], np.zeros(restart["n_tor"]))
    with h5py.File(args.restart, "r") as h5:
        F_ref = float(h5["F0"][0]) if "F0" in h5 else float(restart["fprof_eq"][0, 0])

    si_scale = reader.enable_si_units(restart)
    n_elements = restart["vertex"].shape[1]

    values, n_singular, n_unrepaired = evaluate_with_axis_fix(reader, restart, fourier_basis, args.n_div)
    if n_unrepaired:
        sys.exit(f"B is undefined at {n_unrepaired} samples after the axis repair")

    R_samples = values["R"].ravel()
    Z_samples = values["Z"].ravel()
    vertex_of_sample, representative = merge_vertices(R_samples, Z_samples)
    n_vertices = len(representative)

    # Average every field over the samples that share a vertex and record the spread.
    counts = np.bincount(vertex_of_sample, minlength=n_vertices)
    averaged = {}
    spread = {}
    for name, data in values.items():
        flat = data.ravel()
        mean = np.bincount(vertex_of_sample, weights=flat, minlength=n_vertices) / counts
        averaged[name] = mean
        deviation = np.abs(flat - mean[vertex_of_sample])
        spread[name] = deviation.max() / max(np.abs(flat).max(), np.finfo(float).tiny)

    triangles = build_triangles(n_elements, args.n_div, vertex_of_sample)
    mesh_label = ring_labels(n_elements, args.n_div, vertex_of_sample, n_vertices)
    if mesh_label is None:
        print("grid is not polar: using psi_eq as mesh label")
        mesh_label = -averaged["psi_eq"]
    R = averaged["R"]
    Z = averaged["Z"]
    areas = signed_areas(R, Z, triangles)
    triangles[areas < 0] = triangles[areas < 0][:, ::-1]
    areas = np.abs(areas)

    # Fine sample mesh of F = R B_phi for the toroidal-field potential.
    fine, _, fine_unrepaired = evaluate_with_axis_fix(reader, restart, fourier_basis, args.fine_div)
    if fine_unrepaired:
        sys.exit(f"B is undefined at {fine_unrepaired} fine samples after the axis repair")
    R_fine = fine["R"].ravel()
    Z_fine = fine["Z"].ravel()
    vertex_fine, representative_fine = merge_vertices(R_fine, Z_fine)
    tri_fine = build_triangles(n_elements, args.fine_div, vertex_fine)
    A_Z_toroidal = toroidal_field_potential(
        R, Z,
        (R_fine * fine["B_phi"].ravel())[representative_fine],
        R_fine[representative_fine],
        Z_fine[representative_fine],
        tri_fine, F_ref, R[np.argmax(np.abs(averaged["A3/psi"]))], args.grid,
    )

    # GORILLA conventions: right-handed phi, covariant A, physical B (SI first).
    A_R = np.zeros(n_vertices)
    A_Z = A_Z_toroidal
    A_phi = -averaged["A3/psi"]
    B_R = averaged["BR"]
    B_Z = averaged["BZ"]
    B_phi = -averaged["B_phi"]
    report = curl_check(R, Z, triangles, A_R, A_phi, A_Z, B_R, B_phi, B_Z)

    density = averaged["rho"] * si_scale["rho"] / CM_PER_M**3
    vertices = {
        "R": (R * CM_PER_M, "cm", "major radius"),
        "Z": (Z * CM_PER_M, "cm", "vertical position"),
        "A_R": (A_R * GAUSS_PER_TESLA * CM_PER_M, "G cm", "covariant vector potential, R (zero in this gauge)"),
        "A_phi": (A_phi * GAUSS_PER_TESLA * CM_PER_M**2, "G cm^2",
                  "covariant vector potential, phi (right-handed) = -JOREK A3"),
        "A_Z": (A_Z * GAUSS_PER_TESLA * CM_PER_M, "G cm",
                "covariant vector potential, Z (gauge A_R = 0, gives B_phi)"),
        "B_R": (B_R * GAUSS_PER_TESLA, "G", "physical magnetic field, R"),
        "B_phi": (B_phi * GAUSS_PER_TESLA, "G", "physical magnetic field, phi (right-handed)"),
        "B_Z": (B_Z * GAUSS_PER_TESLA, "G", "physical magnetic field, Z"),
        "Bmod": (np.sqrt(B_R**2 + B_phi**2 + B_Z**2) * GAUSS_PER_TESLA, "G", "magnetic field modulus"),
        "Phi": (np.zeros(n_vertices), "statV", "electrostatic potential (not stored by model 711; zero)"),
        "n": (density, "cm^-3", "ion number density rho*n0 (= electron density for Z=1)"),
        "T_i": (averaged["T_i"] * si_scale["T_i"], "eV", "ion temperature"),
        "T_e": (averaged["T_e"] * si_scale["T_e"], "eV", "electron temperature"),
        "u_R": (averaged["uR"] * si_scale["uR"] * CM_PER_M, "cm/s", "fluid velocity, R"),
        "u_phi": (-averaged["u_phi"] * si_scale["u_phi"] * CM_PER_M, "cm/s",
                  "fluid velocity, phi (right-handed)"),
        "u_Z": (averaged["uZ"] * si_scale["uZ"] * CM_PER_M, "cm/s", "fluid velocity, Z"),
        "mesh_label": (mesh_label, "1",
                       "label constant along the poloidal mesh lines (radial ring index), orients the prisms"),
        "psi_eq": (-averaged["psi_eq"] * GAUSS_PER_TESLA * CM_PER_M**2, "G cm^2",
                   "initial equilibrium flux, same sign convention as A_phi"),
    }
    attributes = {
        "title": "JOREK model-711 axisymmetric field on a triangular mesh for GORILLA",
        "source_file": str(args.restart.resolve()),
        "source_sha256": sha256(args.restart),
        "source_step": restart["step"],
        "source_time_ms": restart["time"] * restart["time_scale"] * 1.0e3,
        "central_density_1e20_m3": restart["central_density"],
        "central_mass_mp": restart["central_mass"],
        "F_ref_T_m": F_ref,
        "n_div": args.n_div,
        "phi_convention": "right-handed (R,phi,Z); JOREK phi components were negated",
        "reader": str((Path(args.reader_dir) / "plot_jorek_model711_h5.py").resolve()),
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
    }
    write_netcdf(args.output, vertices, triangles, attributes)

    # Polar grid: n_pol elements touch the axis (their s=0 side collapses), n_rad rings.
    s_zero = local_samples(args.n_div)[0] == 0.0
    n_pol = int(np.count_nonzero(np.ptp(values["R"][:, s_zero], axis=1) < MERGE_TOLERANCE))
    expected = 1 + (n_elements // n_pol) * args.n_div * n_pol * args.n_div if n_pol else "n/a"
    print(f"{args.restart.name}: t = {attributes['source_time_ms']:.3f} ms, {n_elements} elements")
    print(f"vertices {n_vertices} (structured-grid expectation {expected}), triangles {len(triangles)}")
    print(f"samples with singular Jacobian repaired at s offset {AXIS_OFFSET_S}: {n_singular}")
    print(f"mesh area {areas.sum():.6f} m^2, smallest triangle {areas.min():.3e} m^2")
    print("max relative spread of shared-vertex values:",
          ", ".join(f"{name} {value:.1e}" for name, value in spread.items() if name not in ("R", "Z")))
    print("curl(A) vs B per triangle, relative error (median, 99%, max):")
    for name, (median, p99, worst) in report.items():
        print(f"  {name:6s} {median:.2e} {p99:.2e} {worst:.2e}")
    print(f"written {args.output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Add flux-function gradients of n, T_i, T_e to a GORILLA vertex-data file.

The delta-f source of GORILLA_APPLETS needs d(ln f_M)/dpsi with psi = A_phi.
In an axisymmetric state n, T_i and T_e are (close to) flux functions of A_phi.
Their derivatives are therefore taken from smooth 1D fits n(psi_N), T(psi_N)
(least-squares cubic splines through binned medians), not from per-triangle
gradients, which would amplify any small-scale noise in the input data.

The derivatives are multiplied by a taper that is 1 for psi_N < psin_cut - taper
and falls smoothly (cos^2) to 0 at psi_N = psin_cut, so that the source can be
switched off in an edge region where the input is not trusted.

psi_N = (A_phi - A_phi_axis) / (A_phi_boundary - A_phi_axis), with the axis at
the extremum of A_phi and the boundary value taken from the outermost mesh ring
(an exact flux surface of the fixed-boundary JOREK run).

Added variables (per vertex): psi_N [1], dn_dA_phi [cm^-3 / (G cm^2)],
dT_i_dA_phi and dT_e_dA_phi [eV / (G cm^2)], source_taper [1].
"""

from __future__ import annotations

import argparse
import datetime
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import LSQUnivariateSpline
from scipy.io import netcdf_file


def read_all(path):
    with netcdf_file(path, "r", mmap=False) as nc:
        dims = {name: size for name, size in nc.dimensions.items()}
        variables = {}
        for name, var in nc.variables.items():
            variables[name] = (var.dimensions, var.data.copy(), dict(var._attributes))
        attributes = dict(nc._attributes)
    return dims, variables, attributes


def psi_normalised(A_phi, mesh_label):
    outer = mesh_label == mesh_label.max()
    A_boundary = A_phi[outer].mean()
    # axis = extremum of A_phi on the opposite side from the boundary value
    A_axis = A_phi.min() if A_boundary > A_phi.mean() else A_phi.max()
    spread = np.ptp(A_phi[outer]) / abs(A_boundary - A_axis)
    return (A_phi - A_axis) / (A_boundary - A_axis), A_axis, A_boundary, spread


def fit_profile(psin, values, psin_max, knot_spacing, bin_width):
    """Least-squares cubic spline through binned medians on [0, psin_max]."""
    edges = np.arange(0.0, psin_max + bin_width, bin_width)
    centres, medians = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (psin >= lo) & (psin < hi)
        if np.count_nonzero(m) >= 3:
            centres.append(np.median(psin[m]))
            medians.append(np.median(values[m]))
    centres, medians = np.array(centres), np.array(medians)
    knots = np.arange(knot_spacing, centres[-1] - 0.5 * knot_spacing, knot_spacing)
    spline = LSQUnivariateSpline(centres, medians, knots, k=3)
    residual = np.abs(spline(centres) - medians) / np.abs(medians)
    return spline, residual.max()


def taper(psin, psin_cut, width):
    x = np.clip((psin - (psin_cut - width)) / width, 0.0, 1.0)
    return np.cos(0.5 * np.pi * x) ** 2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=Path, help="vertex-data netCDF written by jorek2gorilla.py")
    parser.add_argument("output", type=Path, help="output netCDF (input plus gradient variables)")
    parser.add_argument("--psin-cut", type=float, default=0.9, help="source is zero beyond this psi_N (default 0.9)")
    parser.add_argument("--taper", type=float, default=0.05, help="width of the cos^2 taper in psi_N (default 0.05)")
    parser.add_argument("--knot-spacing", type=float, default=0.05, help="spline knot spacing in psi_N (default 0.05)")
    parser.add_argument("--bin-width", type=float, default=0.005, help="psi_N bin width for medians (default 0.005)")
    parser.add_argument("--plot", type=Path, help="optional check plot (png)")
    args = parser.parse_args()

    dims, variables, attributes = read_all(args.input)
    A_phi = variables["A_phi"][1]
    mesh_label = variables["mesh_label"][1]
    psin, A_axis, A_boundary, spread = psi_normalised(A_phi, mesh_label)
    dpsin_dA = 1.0 / (A_boundary - A_axis)
    fit_range = min(1.0, args.psin_cut + 0.03)
    w = taper(psin, args.psin_cut, args.taper)

    added = {
        "psi_N": (psin, "1", "normalised poloidal flux (A_phi - A_phi_axis)/(A_phi_boundary - A_phi_axis)"),
        "source_taper": (w, "1", f"cos^2 taper, 1 below psi_N = {args.psin_cut - args.taper}, 0 above {args.psin_cut}"),
    }
    fits = {}
    for name, unit in (("n", "cm^-3"), ("T_i", "eV"), ("T_e", "eV")):
        spline, residual = fit_profile(psin, variables[name][1], fit_range, args.knot_spacing, args.bin_width)
        fits[name] = spline
        derivative = spline.derivative()(np.minimum(psin, fit_range)) * dpsin_dA * w
        added[f"d{name}_dA_phi"] = (derivative, f"{unit} / (G cm^2)",
                                    f"tapered d{name}/dA_phi from a smooth spline fit {name}(psi_N)")
        print(f"{name}: max relative fit residual of the bin medians {residual:.1e}")
    print(f"A_phi axis {A_axis:.6e}, boundary {A_boundary:.6e} G cm^2 (relative spread on the boundary ring {spread:.1e})")

    with netcdf_file(args.output, "w", version=2) as nc:
        for name, size in dims.items():
            nc.createDimension(name, size)
        for name, (dimensions, data, attrs) in variables.items():
            var = nc.createVariable(name, data.dtype, dimensions)
            var[:] = data
            for key, value in attrs.items():
                setattr(var, key, value)
        for name, (data, units, description) in added.items():
            var = nc.createVariable(name, "f8", ("n_vertices",))
            var[:] = data
            var.units = units
            var.description = description
        for key, value in attributes.items():
            setattr(nc, key, value)
        nc.profile_gradients_command = " ".join(sys.argv)
        nc.profile_gradients_created = datetime.datetime.now().isoformat(timespec="seconds")
    print(f"written {args.output}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        grid = np.linspace(0.0, fit_range, 400)
        for k, name in enumerate(("n", "T_i", "T_e")):
            axes[0, k].plot(psin, variables[name][1], ",", alpha=0.3, label="vertex data")
            axes[0, k].plot(grid, fits[name](grid), "k-", label="spline fit")
            axes[0, k].set_title(name); axes[0, k].set_xlabel("psi_N"); axes[0, k].legend(markerscale=20)
            axes[1, k].plot(psin, added[f"d{name}_dA_phi"][0], ",")
            axes[1, k].set_title(f"tapered d{name}/dA_phi"); axes[1, k].set_xlabel("psi_N")
        fig.tight_layout(); fig.savefig(args.plot, dpi=100)
        print(f"plot {args.plot}")


if __name__ == "__main__":
    main()

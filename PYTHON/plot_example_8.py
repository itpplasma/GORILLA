#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example 8:
 * Collisionless guiding-center orbits in an analytic large-aspect-ratio circular tokamak.
 * grid_kind=5: no external equilibrium file required.
 * Poincare plot (phi=0) with predicted banana widths from p_phi conservation:
     Delta R_b = 2 * q(r) * rho_L * lambda / epsilon
 * Four orbits: two deeply trapped (lambda=0.10), one passing (lambda=0.70),
   one marginally trapped (lambda=0.48 ~ sqrt(2eps/(1+eps)) - 0.03).
"""

import f90nml
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

mpl.rcdefaults()

# ── Physical constants (Gaussian CGS) ─────────────────────────────────────────
m_D = 2.0 * 1.67262e-24  # deuterium mass [g]
e_cgs = 4.80326e-10  # elementary charge [esu]
c_cgs = 2.99792e10  # speed of light [cm/s]
eV_erg = 1.60218e-12  # 1 eV in erg

# ── Tokamak / run parameters read from input files ────────────────────────────
_grid = f90nml.read("../EXAMPLES/example_8/tetra_grid.inp")["tetra_grid_nml"]
_plot = f90nml.read("../EXAMPLES/example_8/gorilla_plot.inp")["gorilla_plot_nml"]
R0 = _grid["r0_analytic_circ"]  # major radius [cm]
a = _grid["a_analytic_circ"]  # minor radius [cm]
B0 = _grid["b0_analytic_circ"]  # toroidal field at R=R0 [G]
q0 = _grid["q0_analytic_circ"]  # safety factor on axis
q1 = _grid["q1_analytic_circ"]  # safety factor gradient
E_eV = _plot["energy_ev_start"]  # kinetic energy [eV]


def q_of_r(r):
    return q0 + q1 * (r / a) ** 2


v = np.sqrt(2.0 * E_eV * eV_erg / m_D)
rho_L = m_D * v * c_cgs / (e_cgs * B0)

# ── Load Poincare data ─────────────────────────────────────────────────────────
poincare = np.genfromtxt("../EXAMPLES/example_8/poincare_plot_phi_0_rphiz.dat")
R_all = poincare[:, 0]
Z_all = poincare[:, 2]

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, axs = plt.subplot_mosaic([["poinc", "hist"], ["poinc", "bar"]], figsize=(14, 6))
ax_poinc = axs["poinc"]
ax_hist = axs["hist"]
ax_bar = axs["bar"]
fig.suptitle(
    "Analytic circular tokamak — Poincaré plot ($\\phi=0$)\n"
    "3 keV deuterium: two deeply trapped ($\\lambda=0.10$), "
    "one passing ($\\lambda=0.70$), one marginally trapped ($\\lambda=0.48$)",
    fontsize=10,
)

# ── Left: Poincare scatter ────────────────────────────────────────────────────
ax_poinc.plot(
    R_all,
    Z_all,
    ".",
    markersize=1.0,
    color="steelblue",
    alpha=0.35,
    rasterized=True,
    label="Poincaré points",
)

theta = np.linspace(0, 2 * np.pi, 500)
ax_poinc.plot(
    R0 + a * np.cos(theta), a * np.sin(theta), "k--", linewidth=1, label="limiter"
)


# ── Orbit definitions ──────────────────────────────────────────────────────────
orbit_start_data = np.genfromtxt("../EXAMPLES/example_8/orbit_start_rphizlambda.dat")
if orbit_start_data.ndim == 1:
    orbit_start_data = orbit_start_data[np.newaxis, :]
orbits = [
    dict(
        R_start=row[0],
        Z_start=row[2],
        lam=row[3],
        color=f"C{i + 1}",
        label=f"Orbit {i + 1} (λ={row[3]:.2f})",
    )
    for i, row in enumerate(orbit_start_data)
]

for orb in orbits:
    r = np.sqrt((orb["R_start"] - R0) ** 2 + orb["Z_start"] ** 2)
    q = q_of_r(r)
    eps = r / R0
    lam = orb["lam"]
    orb["r"] = r
    orb["q"] = q
    orb["eps"] = eps
    # exact trapping boundary at this flux surface
    B_min = B0 * (1 - eps)
    B_max = B0 * (1 + eps)
    orb["lam_trap"] = np.sqrt(1 - B_min / B_max)  # exact: sqrt(2eps/(1+eps))
    orb["passing"] = orb["lam"] >= orb["lam_trap"]
    if not orb["passing"]:
        orb["dR_b"] = 2.0 * q * rho_L * lam / eps

# ── Console summary ────────────────────────────────────────────────────────────
print(f"Toroidal Larmor radius  rho_L = {rho_L:.4f} cm")
print(f"Speed                   v     = {v:.3e} cm/s")
print()
print(
    f"{'Orbit':<30} {'r':>6} {'eps':>6} {'q':>6} {'lam_trap(exact)':>16} {'lam':>6} {'margin':>8} {'dR_b':>8}"
)
print("-" * 96)
for orb in orbits:
    trapped_str = "passing" if orb.get("passing") else "trapped"
    margin = "" if orb.get("passing") else f"{orb['lam_trap'] - orb['lam']:+.4f}"
    dR_str = "" if orb.get("passing") else f"{orb['dR_b']:.3f}"
    print(
        f"{orb['label']:<30} {orb['r']:>6.2f} {orb['eps']:>6.4f} {orb['q']:>6.4f} "
        f"{orb['lam_trap']:>16.4f} {orb['lam']:>6.2f} {margin:>8} {dR_str:>8}"
    )


for orb in orbits:
    if orb.get("passing"):
        continue
    R1 = orb["R_start"]
    R2 = R1 + orb["dR_b"]
    c = orb["color"]
    ax_poinc.axvline(R1, color=c, ls="-", lw=1.4)
    ax_poinc.axvline(
        R2,
        color=c,
        ls="--",
        lw=1.4,
        label=f"{orb['label']}: $\\Delta R_b={orb['dR_b']:.2f}$ cm",
    )


def get_middle_bin(bins: np.ndarray, loc: int) -> float:
    return (bins[loc] + bins[loc + 1]) * 0.5


trapped_orbits = [orb for orb in orbits if not orb.get("passing")]
passing_R_max = max(
    (orb["R_start"] for orb in orbits if orb.get("passing")), default=0.0
)
R_near_0 = R_all[(np.abs(Z_all) < 0.2) & (R_all > passing_R_max + 5.0)]
R_hist_lo = min(orb["R_start"] for orb in trapped_orbits) - 5.0
R_hist_hi = max(orb["R_start"] + orb["dR_b"] for orb in trapped_orbits) + 5.0
counts, bins = np.histogram(R_near_0, bins=200, range=(R_hist_lo, R_hist_hi))
banana_locs = find_peaks(counts)
peak_Rs = np.array([get_middle_bin(bins, loc) for loc in banana_locs[0]])
prediced = [orb["dR_b"] for orb in trapped_orbits]
measured = [
    peak_Rs[np.argmin(np.abs(peak_Rs - (orb["R_start"] + orb["dR_b"])))]
    - peak_Rs[np.argmin(np.abs(peak_Rs - orb["R_start"]))]
    for orb in trapped_orbits
]

ax_poinc.vlines(
    peak_Rs,
    -5,
    5,
    color="gray",
    ls=":",
    lw=1.0,
    label="Estimated banana edges (histogram peaks)",
)

x = np.arange(len(prediced))
width = 0.35
ax_bar.bar(
    x - width / 2,
    prediced,
    width,
    label="Predicted $\\Delta R_b$",
    color=[o["color"] for o in trapped_orbits],
    alpha=0.8,
)
ax_bar.bar(
    x + width / 2,
    measured,
    width,
    label="Measured $\\Delta R_b$ (Poincaré data)",
    color=[o["color"] for o in trapped_orbits],
    alpha=0.4,
    hatch="//",
)
ax_bar.set_xticks(x)
ax_bar.set_xticklabels([o["label"] for o in trapped_orbits], fontsize=9)
ax_bar.set_ylabel("Banana width $\\Delta R_b$ [cm]")
ax_bar.set_title("Predicted vs measured banana width")
ax_bar.legend(fontsize=8)
for i, (p, m) in enumerate(zip(prediced, measured)):
    ax_bar.text(
        i - width / 2, p + 0.05, f"{p:.2f}", ha="center", va="bottom", fontsize=8
    )
    ax_bar.text(
        i + width / 2, m + 0.05, f"{m:.2f}", ha="center", va="bottom", fontsize=8
    )

ax_hist.stairs(
    counts,
    bins,
    label="Poincaré data near Z=0",
    color="steelblue",
    alpha=0.7,
    fill=True,
)  # histogram of R near Z=0
ax_hist.set_xlabel("R [cm]")
ax_hist.set_ylabel("Count")
ax_hist.set_title("Histogram of R near Z=0 (banana width estimate)")
ax_hist.legend(fontsize=8)

ax_poinc.set_xlabel("R [cm]")
ax_poinc.set_ylabel("Z [cm]")
ax_poinc.set_aspect("equal")
ax_poinc.legend(fontsize=7, loc="upper left", framealpha=0.85)
ax_poinc.set_xlim(118, 222)
ax_poinc.set_ylim(-55, 55)

plt.tight_layout()
plt.savefig("example_8.png", dpi=150)
print("\nPlot saved to example_8.png")
plt.show()

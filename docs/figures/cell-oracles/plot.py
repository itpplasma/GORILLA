#!/usr/bin/env python3
"""Plot exported native analytic-oracle values without filtering failures."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
with (root / "cpu_trace.csv").open() as source:
    cpu = list(csv.DictReader(source))
with (root / "cuda_trace.csv").open() as source:
    cuda = list(csv.DictReader(source))
rows = [r for r in cpu if int(r["component"]) == 1]
tau = np.array([float(r["tau"]) for r in rows])
baseline = np.array([float(r["baseline"]) for r in rows])
generated = np.array([float(r["generated"]) for r in rows])
reference = np.array([float(r["reference"]) for r in rows])
if not np.allclose(reference, 1 + .7*tau + .75*tau**2, rtol=0, atol=1e-14):
    raise ValueError("CSV does not describe the native nilpotent oracle")

names = [("partial position", "Partial"), ("backward position", "Backward"),
         ("constant acceleration position", "Acceleration"),
         ("first-cell exit", "First exit"),
         ("second-cell final position", "Second cell"),
         ("periodic final position", "Periodic angle")]
selected = [next(r for r in cuda if r["case"] == name) for name, _ in names]
epsilon = np.finfo(float).eps
errors = [abs(float(r["actual"]) - float(r["reference"])) /
          (epsilon * max(1., abs(float(r["reference"])))) for r in selected]

plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none",
                     "svg.hashsalt": "gorilla-cell-oracles"})
fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), constrained_layout=True)
curve = np.linspace(tau.min(), tau.max(), 300)
axes[0].plot(curve, 1 + .7*curve + .75*curve**2, "--", color="#333333",
             label="Analytic reference", linewidth=1.5)
axes[0].plot(tau, baseline, "o", color="#0072B2", markerfacecolor="none",
             markersize=9, markeredgewidth=1.5, label="Legacy Taylor expression")
axes[0].plot(tau, generated, "x", color="#D55E00", markersize=7,
             markeredgewidth=1.5, label="Generated CPU kernel")
axes[0].set(xlabel=r"Evolution parameter $\tau$ (dimensionless)",
            ylabel=r"State component $z_1$ (dimensionless)",
            title="CPU: nilpotent affine system")
axes[0].legend(loc="upper left", fontsize=9.5)
axes[1].axhline(0, linestyle="--", color="#333333", linewidth=1.3,
                label="Analytic reference")
axes[1].plot(range(len(errors)), errors, "D", color="#D55E00",
             markersize=7, label="CUDA kernel")
axes[1].set_xticks(range(len(errors)), [label for _, label in names], rotation=25, ha="right")
axes[1].set(xlabel="Synthetic cell trajectory",
            ylabel=r"Normalized coordinate error / $\epsilon$ (dimensionless)",
            title="CUDA: independent analytic coordinates",
            ylim=(-.25, max(1., 1.25*max(errors))))
axes[1].legend(loc="upper left", fontsize=9.5)
for axis in axes:
    axis.grid(axis="y", color=".88", linewidth=.7)
fig.savefig(root / "cell-oracles.svg", metadata={"Date": None})
fig.savefig(root / "cell-oracles.png", dpi=180)

svg = root / "cell-oracles.svg"
svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")

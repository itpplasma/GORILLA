#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example 8:
 * Compute collisionless guiding-center orbits with GORILLA for two passing and one trapped
   Deuterium particle in an analytic large-aspect-ratio circular tokamak.
 * Use a rectangular grid driven by analytic field formulas (grid_kind=5) — no external
   equilibrium file required.
 * Create a figure with the Poincare plot (phi=0) in cylindrical coordinates (R, Z).
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()

# Analytic tokamak parameters (must match tetra_grid.inp)
R0 = 170.0  # major radius [cm]
a  =  50.0  # minor radius [cm]

# Load Poincare cut data from example_8 in cylindrical coordinates (R, phi, Z)
poincare_rphiz = np.genfromtxt("../EXAMPLES/example_8/poincare_plot_phi_0_rphiz.dat")

fig, ax = plt.subplots(figsize=(7, 7))
fig.suptitle('Analytic circular tokamak: Poincare plot (phi = 0)\n'
             'Two trapped (lambda=0.10) and one passing (lambda=0.70) deuterium particle')

ax.plot(poincare_rphiz[:, 0], poincare_rphiz[:, 2], '.', markersize=2)

# Draw the limiter circle
theta = np.linspace(0, 2 * np.pi, 400)
ax.plot(R0 + a * np.cos(theta), a * np.sin(theta), 'k--', linewidth=1, label='limiter')

ax.set_xlabel('R [cm]')
ax.set_ylabel('Z [cm]')
ax.set_aspect('equal')
ax.legend()
plt.tight_layout()

plt.savefig('example_8.png', dpi=150)
print("Plot saved to example_8.png")

plt.show()

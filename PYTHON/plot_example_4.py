#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:48:16 2021

Example 4:
 * Compute collisionless guiding-center orbits with GORILLA for two passing and one trapped Deuterium particle.
 * Use a field-aligned grid for an axisymmetric tokamak equilibrium (g-file)
 * Create a figure with the Poincaré plots (\varphi = 0) in cylindrical and symmetry flux coordinates.

@author: Michael Eder
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# Load Poincaré cut data from example_4 in symmetry flux coordinates (s,\vartheta,\varphi)
poincare_sthetaphi = np.genfromtxt("../EXAMPLES/example_4/poincare_plot_phi_0_sthetaphi.dat")

# Load Poincaré cut data from example_4 in cylindrical coordinates (R,\varphi,Z)
poincare_rphiz = np.genfromtxt("../EXAMPLES/example_4/poincare_plot_phi_0_rphiz.dat")



# Plot Poincaré cuts and for two passing and one trapped Deuterium particle.
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Tokamak: Poincaré plots for two passing and one trapped Deuterium particle.')


ax1.plot(poincare_rphiz[:,0],poincare_rphiz[:,2],'.')
ax1.set_xlabel('R [cm]')
ax1.set_ylabel('Z [cm]')
ax1.set_title('Poincaré phi = 0')


ax2.plot(poincare_sthetaphi[:,1],poincare_sthetaphi[:,0],'.')
ax2.set_xlabel('theta')
ax2.set_ylabel('s')
ax2.set_title('Poincaré phi = 0')


plt.tight_layout()
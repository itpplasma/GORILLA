#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:32:20 2021

Example 3:
 * Compute a collisionless guiding-center orbit with GORILLA for a passing Deuterium particle. (Manually execute ''test_gorilla_main.x'')
 * Use a field-aligned grid for an axisymmetric tokamak equilibrium (g-file)
 * Load the results of GORILLA with polynominal order 3.
 * Create a figure with the Poincaré plots (\varphi = 0) in cylindrical and symmetry flux coordinates.
 * Compute the normalized toroidal angular momentum as a function of toroidal mappings.

@author: Michael Eder
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# Load Poincaré cut data from example_3 in symmetry flux coordinates (s,\vartheta,\varphi)
poincare_sthetaphi = np.genfromtxt("../EXAMPLES/example_3/poincare_plot_phi_0_sthetaphi.dat")

# Load Poincaré cut data from example_3 in cylindrical coordinates (R,\varphi,Z)
poincare_rphiz = np.genfromtxt("../EXAMPLES/example_3/poincare_plot_phi_0_rphiz.dat")

# Load toroidal angular momentum data from example_3
p_phi = np.genfromtxt("../EXAMPLES/example_3/p_phi.dat")

# Normalize the toroidal angular momentum to the value at t = 0
p_phi[:,1] = p_phi[:,1] / p_phi[0,1] 

# Plot Poincaré cuts and evolution of the normalized toroidal angular momentum as a function of toroidal mappings
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.suptitle('Tokamak: Poincaré plots and toroidal angular momentum of passing particle')


ax1.plot(poincare_rphiz[:,0],poincare_rphiz[:,2],'.')
ax1.set_xlabel('R [cm]')
ax1.set_ylabel('Z [cm]')
ax1.set_title('Poincaré phi = 0')


ax2.plot(poincare_sthetaphi[:,1],poincare_sthetaphi[:,0],'.')
ax2.set_xlabel('theta')
ax2.set_ylabel('s')
ax2.set_title('Poincaré phi = 0')


ax3.plot(abs(p_phi[:,0]),p_phi[:,1],'.')
ax3.set_xlabel('N_mappings')
ax3.set_ylabel('p_phi / p_phi(t=0)')
ax3.set_title('Toroidal angular momentum')

plt.tight_layout()
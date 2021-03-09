#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:09:29 2021

Example 2:
 * Compute a collisionless guiding-center orbit with GORILLA for a passing Deuterium particle. (Manually execute ''test_gorilla_main.x'')
 * Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
 * Load the results of GORILLA with Runge-Kutta 4.
 * Create a figure with the Poincaré plots (\varphi = 0) in cylindrical and symmetry flux coordinates.
 * Compute the normalized total energy as a function of toroidal mappings. 

@author: Michael Eder
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# Load Poincaré cut data from example_2 in symmetry flux coordinates (s,\vartheta,\varphi)
poincare_sthetaphi = np.genfromtxt("../EXAMPLES/example_2/poincare_plot_phi_0_sthetaphi.dat")

# Load Poincaré cut data from example_2 in cylindrical coordinates (R,\varphi,Z)
poincare_rphiz = np.genfromtxt("../EXAMPLES/example_2/poincare_plot_phi_0_rphiz.dat")

# Load total energy data from example_2
E_tot = np.genfromtxt("../EXAMPLES/example_2/e_tot.dat")

# Normalize total energy to the value at t = 0
E_tot[:,1] = E_tot[:,1] / E_tot[0,1] 

# Plot Poincaré cuts and evolution of total energy as a function of toroidal mappings
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.suptitle('Stellarator: Poincaré plots and total energy of passing particle')


ax1.plot(poincare_rphiz[:,0],poincare_rphiz[:,2],'.')
ax1.set_xlabel('R [cm]')
ax1.set_ylabel('Z [cm]')
ax1.set_title('Poincaré phi = 0')


ax2.plot(poincare_sthetaphi[:,1],poincare_sthetaphi[:,0],'.')
ax2.set_xlabel('theta')
ax2.set_ylabel('s')
ax2.set_title('Poincaré phi = 0')


ax3.plot(E_tot[:,0],E_tot[:,1],'.')
ax3.set_xlabel('N_mappings')
ax3.set_ylabel('E_tot / E_tot(t=0)')
ax3.set_title('Total energy')

plt.tight_layout()
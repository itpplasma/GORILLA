#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 09:18:31 2021
Example 1:
 * Compute a collisionless guiding-center orbit with GORILLA for a trapped Deuterium particle. (Manually execute ''test_gorilla_main.x'')
 * Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
 * Load the results of GORILLA with polynominal order 4.
 * Create a figure with the Poincaré sections (v_\parallel = 0) in cylindrical and symmetry flux coordinates.
 * Compute the normalized parallel adiabatic invariant as a function of banana bounces.
@author: Michael Eder
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# Load Poincaré cut data from example_1 in symmetry flux coordinates (s,\vartheta,\varphi)
poincare_sthetaphi = np.genfromtxt("../EXAMPLES/example_1/poincare_plot_vpar_0_sthetaphi.dat")

# Load Poincaré cut data from example_1 in cylindrical coordinates (R,\varphi,Z)
poincare_rphiz = np.genfromtxt("../EXAMPLES/example_1/poincare_plot_vpar_0_rphiz.dat")

# Load parallel adiabatic invariant data from example_1
J_par = np.genfromtxt("../EXAMPLES/example_1/J_par.dat")

# Normalize parallel adiabatic invariant to the value at t = 0
J_par[:,1] = J_par[:,1] / J_par[0,1] 

# Plot Poincaré sections and evolution of the normalized parallel adiabatic invariant as a function of banana bounces
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.suptitle('Stellarator: Poincaré sections and parallel adiabatic invariant of trapped particle')


ax1.plot(poincare_rphiz[:,0],poincare_rphiz[:,2],'.')
ax1.set_xlabel('R [cm]')
ax1.set_ylabel('Z [cm]')
ax1.set_title('Poincaré v_par = 0')


ax2.plot(poincare_sthetaphi[:,1],poincare_sthetaphi[:,0],'.')
ax2.set_xlabel('theta')
ax2.set_ylabel('s')
ax2.set_title('Poincaré v_par = 0')


ax3.plot(J_par[:,0],J_par[:,1],'.')
ax3.set_xlabel('N_mappings')
ax3.set_ylabel('J_par / J_par(t=0)')
ax3.set_title('Parallel adiabatic invariant')

plt.tight_layout()

plt.show()
#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 11:17:37 2022

Example 7:
 * Compute collisionless guiding-center orbit with GORILLA for a trapped Deuterium particle with adaptive scheme.
 * Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
 * Create a figure with the Poincar� plots (\varphi = 0) in cylindrical and symmetry flux coordinates.
 * Compute the normalized parallel adiabatic invariant as a function of banana bounces.
 * Plot fluctuation and evolution of energy over the bounces.

@author: Georg Graßler
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# Load vpar = 0 poincare data in both coordinate systems from example_7
poincare_vpar_0_rphiz_adaptive = np.genfromtxt('../EXAMPLES/example_7/poincare_plot_vpar_0_rphiz_adaptive.dat')
poincare_vpar_0_sthetaphi_adaptive = np.genfromtxt('../EXAMPLES/example_7/poincare_plot_vpar_0_sthetaphi_adaptive.dat')

# Load parallel adiabatic invariant data from example_7
J_par_adaptive = np.genfromtxt('../EXAMPLES/example_7/J_par_adaptive.dat') 

# Load total particle energy data from example_7
e_tot_adaptive = np.genfromtxt('../EXAMPLES/example_7/e_tot_adaptive.dat')

# Plotting Options
adaptiveColor = np.array([217,95,2])/256

MarkerPoincare = '.'
MarkerJpar = 'o'
MarkerEnergy = '-'
MarkerEnergyFluctuation = '.'

MarkerSizePoincare = 15
MarkerSizeJpar = 5

n_skip_poincare_sthetaphi = 1
n_skip_poincare_rphiz = 1
n_skip_J_par = 1
n_skip_e_tot_fluctuation = 1
n_skip_e_tot = 1

CentralLineWidth = 2
figsize = (12.5,6)

# Plot vpar = 0 poincare plots in both coordinate systems
poincare_fig = plt.figure(figsize=figsize)
ax1 = poincare_fig.add_subplot(121)
poincare_rphiz_plot = ax1.plot(poincare_vpar_0_rphiz_adaptive[0:-1:n_skip_poincare_rphiz,0],poincare_vpar_0_rphiz_adaptive[0:-1:n_skip_poincare_rphiz,2],
                    MarkerPoincare, markersize = MarkerSizePoincare, color = adaptiveColor, label = 'order 2 (adaptive)')
ax1.set_xlabel(r'$R$ [cm]')
ax1.set_ylabel(r'$Z$ [cm]')
ax1.set_title(r'Poincaré cut ($v_{\parallel}=0$)')
ax1.legend()
ax2 = poincare_fig.add_subplot(122)
poincare_sthetaphi_plot = ax2.plot(poincare_vpar_0_sthetaphi_adaptive[0:-1:n_skip_poincare_sthetaphi,1],poincare_vpar_0_sthetaphi_adaptive[0:-1:n_skip_poincare_sthetaphi,0],
                    MarkerPoincare, markersize = MarkerSizePoincare, color = adaptiveColor, label = 'order 2 (adaptive)')
ax2.set_xlabel(r'$\vartheta$')
ax2.set_ylabel(r'$s$')
ax2.set_title(r'Poincaré cut ($v_{\parallel}=0$)')
ax2.legend()

# Plot J_par over banana bounces
J_par_fig = plt.figure(figsize=figsize)
ax1 = J_par_fig.add_subplot(111)
J_par_plot = ax1.plot(J_par_adaptive[0:-1:n_skip_J_par,0],J_par_adaptive[0:-1:n_skip_J_par,1]/J_par_adaptive[0,1],
                    MarkerJpar, markersize = MarkerSizeJpar, color = adaptiveColor, label = 'order 2 (adaptive)')
central_line = ax1.plot(J_par_adaptive[0:-1:n_skip_J_par,0],np.ones(np.shape(J_par_adaptive[0:-1:n_skip_J_par,0])),
                    '-k', linewidth = CentralLineWidth)
plt.xlim(np.min(J_par_adaptive[0:-1:n_skip_J_par,0]),np.max(J_par_adaptive[0:-1:n_skip_J_par,0]))
ax1.set_xlabel(r'$N_{Bounces}$')
ax1.set_ylabel(r'$J_{\parallel}/J_{\parallel}(\tau = 0)$')
ax1.set_title(r'$J_{\parallel}$ over number of bounces')
ax1.legend()

# Plot e_tot over banana bounces
e_tot_fig = plt.figure(figsize=figsize)
ax1 = e_tot_fig.add_subplot(111)
e_tot_plot = ax1.plot(e_tot_adaptive[0:-1:n_skip_e_tot,0],e_tot_adaptive[0:-1:n_skip_e_tot,1]/e_tot_adaptive[0,1] - 1,
                    MarkerEnergy, color = adaptiveColor, label = 'order 2 (adaptive)')
ax1.set_xlabel(r'$N_{Bounces}$')
ax1.set_ylabel(r'$E/E(\tau = 0) - 1$')
ax1.set_title(r'Energy evolution (adaptive/order 2)')

# Plot flucutuations of e_tot over banana bounces (the difference of energy between successive bounces)
e_tot_fluctuation_adaptive = np.zeros(np.shape(e_tot_adaptive))
e_tot_fluctuation_adaptive = np.delete(e_tot_fluctuation_adaptive,-1,0)
e_tot_fluctuation_adaptive[:,0] = e_tot_adaptive[1:,0]
e_tot_fluctuation_adaptive[:,1] = abs(e_tot_adaptive[1:,1]/e_tot_adaptive[0:-1,1] -1)
e_tot_fluctuation_fig = plt.figure(figsize=figsize)
ax1 = e_tot_fluctuation_fig.add_subplot(111)
e_tot_fluctuation_plot = ax1.plot(e_tot_fluctuation_adaptive[0:-1:n_skip_e_tot_fluctuation,0],e_tot_fluctuation_adaptive[0:-1:n_skip_e_tot_fluctuation,1],
                    MarkerEnergy, color = adaptiveColor, label = 'order 2 (adaptive)')
ax1.set_xlabel(r'$N_{Bounces}$')
ax1.set_ylabel(r'$\delta$')
ax1.set_title(r'Energy fluctuation between bounces (adaptive/order 2)')

plt.tight_layout()
plt.show()
#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 11:17:37 2022

Example 7:
 * Compute collisionless guiding-center orbit with GORILLA for a trapped Deuterium particle with/without adaptive scheme.
 * Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
 * Create a figure with the Poincar� plots (\varphi = 0) in cylindrical and symmetry flux coordinates.
 * Compute the normalized parallel adiabatic invariant as a function of banana bounces.
 * Compare energy fluctuation as function of dwell time and evolution of energy between the two schemes.

@author: Georg Graßler
"""

import f90nml
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# For Construction of Demonstration
# -----------------------------------------------------------------------------------------------------------------------

# Predefinig options for calculation
total_orbit_time = -1E1
desired_delta_energy = 1E-10

# Initialize used paths and files
# -----------------------------------------------------------------------------------------------------------------------

# Name of the current calculation to create folder
name_test_case = 'example_7'

# path of PYTHON script
path_script = os.getcwd()

# main path of GORILLA
path_main= path_script + '/..'

# path to RUN folder
path_RUN = path_main + '/EXAMPLES/PYTHON_RUN/' + name_test_case
if not os.path.exists(path_RUN):
  os.makedirs(path_RUN)

# path of input files (blueprints)
path_inp_files = path_main + '/INPUT'

# define path for plots and create a new folder
name_data_folder ='/data_plots/' + name_test_case
path_data_plots = path_script + name_data_folder
if not os.path.exists(path_data_plots):
  os.makedirs(path_data_plots)

# Change to RUN folder and clean it
os. chdir(path_RUN)
if os.listdir('../' + name_test_case):
  os.system('rm -r ../' + name_test_case + '/*')

# load the namelists from the inputfiles with proper comma setting
gorilla = f90nml.read(path_inp_files + '/gorilla.inp')
gorilla.end_comma = True
tetra_grid = f90nml.read(path_inp_files + '/tetra_grid.inp')
tetra_grid.end_comma = True
gorilla_plot = f90nml.read(path_inp_files + '/gorilla_plot.inp')
gorilla_plot.end_comma = True


# Set input variables for GORILLA
# -----------------------------------------------------------------------------------------------------------------------

# Input file gorilla.inp
# All necessary variables for current calculation

# Change in the electrostatic potential within the plasma volume in Gaussian units
gorilla['gorillanml']['eps_Phi'] = 0

# Coordinate system
# 1 ... (R,phi,Z) cylindrical coordinate system
# 2 ... (s,theta,phi) symmetry flux coordinate system
gorilla['gorillanml']['coord_system'] = 2

# particle species
# 1 ... electron, 2 ... deuterium ion, 3 ... alpha particle
gorilla['gorillanml']['ispecies'] = 2

# Switch for initial periodic coordinate particle re-location (modulo operation)
# true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
# false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
gorilla['gorillanml']['boole_periodic_relocation'] = False

# 1 ... numerical RK pusher, 2 ... polynomial pusher
gorilla['gorillanml']['ipusher'] = 2

# true ... Face guessing algorithm is used, false ... NOT used
gorilla['gorillanml']['boole_guess'] = True

# Polynomial order for orbit pusher (from 2 to 4)
gorilla['gorillanml']['poly_order'] = 2

# Switches to calculate optional quantities
gorilla['gorillanml']['boole_time_Hamiltonian'] = False
gorilla['gorillanml']['boole_gyrophase'] = False

# Switch for adaptive time step scheme
# false ... no adaptive scheme
# true ... adaptive scheme to ensure energy conservation up to specified fluctuation
gorilla['gorillanml']['boole_adaptive_time_steps'] = False

# Allowed relative fluctuation of energy between entry and exit of a tetrahedron
gorilla['gorillanml']['desired_delta_energy']  = desired_delta_energy

# Maximum number of intermediate steps allowed for splitting up an orbit section
gorilla['gorillanml']['max_n_intermediate_steps'] = 10000


# Input file tetra_grid.inp
# All necessary variables for current calculation

# Grid Size
# Rectangular: nR, Field-aligned: ns
tetra_grid['tetra_grid_nml']['n1'] = 100
# Rectangular: nphi, Field-aligned: nphi
tetra_grid['tetra_grid_nml']['n2'] = 30
# Rectangular: nZ, Field-aligned: ntheta
tetra_grid['tetra_grid_nml']['n3'] = 30

# Grid kind
# 1 ... rectangular grid for axisymmetric EFIT data
# 2 ... field-aligned grid for axisymmetric EFIT data
# 3 ... field-aligned grid for non-axisymmetric VMEC
# 4 ... SOLEDGE3X_EIRENE grid
tetra_grid['tetra_grid_nml']['grid_kind'] = 3

# Switch for selecting number of field periods automatically or manually
# .true. ... number of field periods is selected automatically (Tokamak = 1, Stellarator depending on VMEC equilibrium)
# .false. ... number of field periods is selected manually (see below)
tetra_grid['tetra_grid_nml']['boole_n_field_periods'] = True


# Input file gorilla_plot.inp
# All necessary variables for current calculation

# Switch for options
# 1 ... Single orbit - Starting positions and pitch for the orbit are taken from file (see below) [First Line]
# 2 ... Single orbit - Starting positions and pitch for the orbit are taken from starting drift surfaces (see below)
# 3 ... Multiple orbits - Starting positions and pitch for orbits are taken from file (see below) [Every Line New Starting position]
# 4 ... Multiple orbits - Starting positions and pitch for orbits are taken from drift surfaces with regular spacing (see below)
gorilla_plot['gorilla_plot_nml']['i_orbit_options'] = 2

# Total individual orbit flight time for plotting
gorilla_plot['gorilla_plot_nml']['total_orbit_time'] = total_orbit_time

# Total Energy of particle in eV
gorilla_plot['gorilla_plot_nml']['energy_eV_start'] = 3000

# Switch for plotting Poincaré cuts at toroidal variable $\varphi$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_phi_0'] = False

# Switch for plotting Poincaré cuts at parallel velocity $v_\parallel$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_vpar_0'] = True

# Number of skipped (non-printed) Poincaré cuts at parallel velocity $v_\parallel$ = 0
gorilla_plot['gorilla_plot_nml']['n_skip_vpar_0'] = 1

# Filename for Poincaré cuts at parallel velocity $v_\parallel$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['filename_poincare_vpar_0_rphiz'] = 'poincare_plot_vpar_0_rphiz.dat'

# Filename for Poincaré cuts at parallel velocity $v_\parallel$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
gorilla_plot['gorilla_plot_nml']['filename_poincare_vpar_0_sthetaphi'] = 'poincare_plot_vpar_0_sthetaphi.dat'

# Switch for plotting full orbit
gorilla_plot['gorilla_plot_nml']['boole_full_orbit'] = False

# Number of skipped (non-printed tetrahedra passings) full orbit
gorilla_plot['gorilla_plot_nml']['n_skip_full_orbit'] = 1

# Plot invariances of motion (ONLY for single orbits)

# Switch for plotting total particle energy
gorilla_plot['gorilla_plot_nml']['boole_e_tot'] = True

# Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
gorilla_plot['gorilla_plot_nml']['boole_p_phi'] = False

# Switch for parallel adiabatic invariant $J_\parallel$
gorilla_plot['gorilla_plot_nml']['boole_J_par'] = True

# Filename for parallel adiabatic invariant $J_\parallel$
gorilla_plot['gorilla_plot_nml']['filename_J_par'] = 'J_par.dat'

# Filename for total particle energy
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot.dat'

# Single orbit from starting drift surface (i_orbit_options = 2)

# Starting drift surface
# = s for (s,$\vartheta$,$\varphi$)
# = R for (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['start_pos_x1_beg'] = 0.5

# Starting value for toroidal variable
# = $\vartheta$ for (s,$\vartheta$,$\varphi$)
# = $\varphi$ for (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['start_pos_x2']  = 0

# Starting value for poloidal variable $\vartheta$
# = $\varphi$ for (s,$\vartheta$,$\varphi$)
# = Z for (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['start_pos_x3']  = 0.63

# Pitch parameter $\lambda$ = $v_\parallel$ / vmod
gorilla_plot['gorilla_plot_nml']['start_pitch_parameter'] = 0.3


# Run GORILLA
# -----------------------------------------------------------------------------------------------------------------------

# write Input files for GORILLA
gorilla.write(path_RUN + '/gorilla.inp', force = True)
tetra_grid.write(path_RUN + '/tetra_grid.inp', force = True)
gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)

# Create softlinks for used files
if os.path.exists(path_main + '/test_gorilla_main.x'):
  os.system('ln -s ../../../test_gorilla_main.x .')
elif os.path.exists(path_main + '/BUILD/SRC/test_gorilla_main.x'):
  os.system('ln -s ../../../BUILD/SRC/test_gorilla_main.x  .') 
else:
  print('GORILLA not built, exiting the PYTHON script')
  raise SystemExit(0)

os.system('ln -s ../../../MHD_EQUILIBRIA .')
os.system('ln -s ../../../INPUT/field_divB0.inp .')
os.system('ln -s ../../../INPUT/preload_for_SYNCH.inp .')

# Run GORILLA code
os.system('./test_gorilla_main.x')

# Setup for comparison run of polynominal order = 4
gorilla['gorillanml']['poly_order'] = 4
gorilla_plot['gorilla_plot_nml']['filename_J_par'] = 'J_par_control_order4.dat'
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot_control_order4.dat'
gorilla_plot['gorilla_plot_nml']['filename_full_orbit_rphiz'] = 'full_orbit_plot_rphiz_control_order4.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_vpar_0_rphiz'] = 'poincare_plot_vpar_0_rphiz_control_order4.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_vpar_0_sthetaphi'] = 'poincare_plot_vpar_0_sthetaphi_control_order4.dat'
gorilla.write(path_RUN + '/gorilla.inp', force = True)
gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)

# Commparison run of 4th order
os.system('./test_gorilla_main.x')

# Changes for adaptive RUN
gorilla['gorillanml']['poly_order'] = 2
gorilla['gorillanml']['boole_adaptive_time_steps'] = True
gorilla_plot['gorilla_plot_nml']['filename_poincare_vpar_0_rphiz'] = 'poincare_plot_vpar_0_rphiz_adaptive.dat'
gorilla_plot['gorilla_plot_nml']['filename_full_orbit_rphiz'] = 'full_orbit_plot_rphiz_adaptive.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_vpar_0_sthetaphi'] = 'poincare_plot_vpar_0_sthetaphi_adaptive.dat'
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot_adaptive.dat'
gorilla_plot['gorilla_plot_nml']['filename_J_par'] = 'J_par_adaptive.dat'
gorilla.write(path_RUN + '/gorilla.inp', force = True)
gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)

# Run with energy control
os.system('./test_gorilla_main.x')

# Setup to observe energy fluctuation in detail
gorilla['gorillanml']['poly_order'] = 2
gorilla['gorillanml']['boole_adaptive_time_steps'] = False
gorilla_plot['gorilla_plot_nml']['boole_poincare_vpar_0'] = False
gorilla_plot['gorilla_plot_nml']['boole_full_orbit'] = True
gorilla_plot['gorilla_plot_nml']['boole_J_par'] = False
gorilla_plot['gorilla_plot_nml']['total_orbit_time'] = total_orbit_time/100
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot_full_orbit.dat'
gorilla.write(path_RUN + '/gorilla.inp', force = True)
gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)

# Run for full orbit analysis (shorter orbit)
os.system('./test_gorilla_main.x')

# Setup to observe energy fluctuation in detail (adaptive)
gorilla['gorillanml']['poly_order'] = 2
gorilla['gorillanml']['boole_adaptive_time_steps'] = True
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot_full_orbit_adaptive.dat'
gorilla.write(path_RUN + '/gorilla.inp', force = True)
gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)

# Run for full orbit analysis (shorter orbit) witch energy control
os.system('./test_gorilla_main.x')


# Create plots of generated data
# -----------------------------------------------------------------------------------------------------------------------

# Here the absolut path to the RUN folder is needed
path_script2RUN = path_RUN
extension = '.png'

# Load Poincaré cut data in cylindrical coordinates (R,\varphi,Z)
poincare_vpar_0_rphiz = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_vpar_0_rphiz.dat')
poincare_vpar_0_rphiz_adaptive = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_vpar_0_rphiz_adaptive.dat')
poincare_vpar_0_rphiz_control_order4 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_vpar_0_rphiz_control_order4.dat')

# Load Poincaré cut data in symmetry flux coordinates (s,\vartheta,\varphi)
poincare_vpar_0_sthetaphi = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_vpar_0_sthetaphi.dat')
poincare_vpar_0_sthetaphi_adaptive = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_vpar_0_sthetaphi_adaptive.dat')
poincare_vpar_0_sthetaphi_control_order4 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_vpar_0_sthetaphi_control_order4.dat')

# Load parallel adiabatic invariant data from example_7
J_par = np.genfromtxt(path_script2RUN + '/' + 'J_par.dat')
J_par_adaptive = np.genfromtxt(path_script2RUN + '/' + 'J_par_adaptive.dat')
J_par_control_order4 = np.genfromtxt(path_script2RUN + '/' + 'J_par_control_order4.dat')

# Load total particle energy data from example_7
e_tot = np.genfromtxt(path_script2RUN + '/' + 'e_tot.dat')
e_tot_adaptive = np.genfromtxt(path_script2RUN + '/' + 'e_tot_adaptive.dat')

#Load total particle energy data after each push for fluctuation visualization
e_tot_full_orbit = np.genfromtxt(path_script2RUN + '/' + 'e_tot_full_orbit.dat')
e_tot_full_orbit_adaptive = np.genfromtxt(path_script2RUN + '/' + 'e_tot_full_orbit_adaptive.dat')

# Plotting Options
order2Color = np.array([27,158,119])/256
adaptiveColor = np.array([217,95,2])/256
order4Color = np.array([117,112,179])/256

MarkerPoincare = '.'
MarkerJpar = 'o'
MarkerEnergy = '-'
MarkerEnergyFluctuation = '.'

MarkerSizePoincare = 6
MarkerSizeJpar = 3
MarkerSizeEnergyFluctuation = 10

n_skip_poincare_sthetaphi = 5
n_skip_poincare_rphiz = 5
n_skip_J_par = 20
n_skip_e_tot_fluctuation = 1
n_skip_e_tot = 1

CentralLineWidth = 2
figsize = (12.5,6)

# Plot vpar = 0 poincare plots in both coordinate systems
poincare_fig = plt.figure(figsize=figsize)
ax1 = poincare_fig.add_subplot(121)
poincare_rphiz_plot = ax1.plot(poincare_vpar_0_rphiz[0:-1:n_skip_poincare_rphiz,0],poincare_vpar_0_rphiz[0:-1:n_skip_poincare_rphiz,2],
                    MarkerPoincare, markersize = MarkerSizePoincare, color = order2Color, label = 'order 2')
poincare_rphiz_plot_control_order4 = ax1.plot(poincare_vpar_0_rphiz_control_order4[0:-1:n_skip_poincare_rphiz,0],poincare_vpar_0_rphiz_control_order4[0:-1:n_skip_poincare_rphiz,2],
                    MarkerPoincare, markersize = MarkerSizePoincare, color = order4Color, label = 'order 4')
poincare_rphiz_plot_adaptive = ax1.plot(poincare_vpar_0_rphiz_adaptive[0:-1:n_skip_poincare_rphiz,0],poincare_vpar_0_rphiz_adaptive[0:-1:n_skip_poincare_rphiz,2],
                    MarkerPoincare, markersize = MarkerSizePoincare, color = adaptiveColor, label = 'order 2 (adaptive)')
ax1.set_xlabel(r'$R$ [cm]')
ax1.set_ylabel(r'$Z$ [cm]')
ax1.set_title(r'Poincaré cut ($v_{\parallel}=0$)')
ax1.legend()
ax2 = poincare_fig.add_subplot(122)
poincare_sthetaphi_plot = ax2.plot(poincare_vpar_0_sthetaphi[0:-1:n_skip_poincare_sthetaphi,1],poincare_vpar_0_sthetaphi[0:-1:n_skip_poincare_sthetaphi,0],
                    MarkerPoincare, markersize = MarkerSizePoincare, color = order2Color, label = 'order 2')
poincare_sthetaphi_plot_control_order4 = ax2.plot(poincare_vpar_0_sthetaphi_control_order4[0:-1:n_skip_poincare_sthetaphi,1],poincare_vpar_0_sthetaphi_control_order4[0:-1:n_skip_poincare_sthetaphi,0],
                    MarkerPoincare, markersize = MarkerSizePoincare, color = order4Color, label = 'order 4')
poincare_sthetaphi_plot_adaptive = ax2.plot(poincare_vpar_0_sthetaphi_adaptive[0:-1:n_skip_poincare_sthetaphi,1],poincare_vpar_0_sthetaphi_adaptive[0:-1:n_skip_poincare_sthetaphi,0],
                    MarkerPoincare, markersize = MarkerSizePoincare, color = adaptiveColor, label = 'order 2 (adaptive)')
ax2.set_xlabel(r'$\vartheta$')
ax2.set_ylabel(r'$s$')
ax2.set_title(r'Poincaré cut ($v_{\parallel}=0$)')
ax2.legend()

# Plot J_par over banana bounces
J_par_fig = plt.figure(figsize=figsize)
ax1 = J_par_fig.add_subplot(111)
J_par_plot = ax1.plot(J_par[0:-1:n_skip_J_par,0],J_par[0:-1:n_skip_J_par,1]/J_par[0,1],
                    MarkerJpar, markersize = MarkerSizeJpar, color = order2Color, label = 'order 2')
J_par_plot_control_order4 = ax1.plot(J_par_control_order4[0:-1:n_skip_J_par,0],J_par_control_order4[0:-1:n_skip_J_par,1]/J_par_control_order4[0,1],
                    MarkerJpar, markersize = MarkerSizeJpar, color = order4Color, label = 'order 4')
J_par_plot_adaptive = ax1.plot(J_par_adaptive[0:-1:n_skip_J_par,0],J_par_adaptive[0:-1:n_skip_J_par,1]/J_par_adaptive[0,1],
                    MarkerJpar, markersize = MarkerSizeJpar, color = adaptiveColor, label = 'order 2 (adaptive)')
central_line = ax1.plot(J_par[0:-1:n_skip_J_par,0],np.ones(np.shape(J_par[0:-1:n_skip_J_par,0])),
                    '-k', linewidth = CentralLineWidth)
plt.xlim(np.min(J_par_adaptive[0:-1:n_skip_J_par,0]),np.max(J_par_adaptive[0:-1:n_skip_J_par,0]))
ax1.set_xlabel(r'$N_{Bounces}$')
ax1.set_ylabel(r'$J_{\parallel}/J_{\parallel}(\tau = 0)$')
ax1.set_title(r'$J_{\parallel}$ over number of bounces')
ax1.legend()

# Plot e_tot over banana bounces
e_tot_fig = plt.figure(figsize=figsize)
ax1 = e_tot_fig.add_subplot(121)
e_tot_plot = ax1.plot(e_tot[0:-1:n_skip_e_tot,0],e_tot[0:-1:n_skip_e_tot,1]/e_tot[0,1] - 1,
                    MarkerEnergy, color = order2Color, label = 'order 2')
ax1.set_xlabel(r'$N_{Bounces}$')
ax1.set_ylabel(r'$E/E(\tau = 0) - 1$')
ax1.set_title(r'Energy evolution (order 2)')

ax2 = e_tot_fig.add_subplot(122)
e_tot_plot_adaptive = ax2.plot(e_tot_adaptive[0:-1:n_skip_e_tot,0],e_tot_adaptive[0:-1:n_skip_e_tot,1]/e_tot_adaptive[0,1] - 1,
                    MarkerEnergy, color = adaptiveColor, label = 'order 2 (adaptive)')
ax2.set_xlabel(r'$N_{Bounces}$')
ax2.set_ylabel(r'$E/E(\tau = 0) - 1$')
ax2.set_title(r'Energy evolution (adaptive/order 2)')

# Plot flucutuations of e_tot as function of dwell time (the difference of energy between successive pushings)
e_tot_fluctuation_fig = plt.figure(figsize=figsize)
ax1 = e_tot_fluctuation_fig.add_subplot(121)
e_tot_fluctuation = np.zeros(np.shape(e_tot_full_orbit))
e_tot_fluctuation = np.delete(e_tot_fluctuation,-1,0)
e_tot_fluctuation[:,0] = abs(e_tot_full_orbit[1:,0] - e_tot_full_orbit[0:-1,0])
e_tot_fluctuation[:,1] = abs(e_tot_full_orbit[1:,1]/e_tot_full_orbit[0:-1,1] -1)
e_tot_fluctuation_plot = ax1.plot(e_tot_fluctuation[0:-1:n_skip_e_tot_fluctuation,0],e_tot_fluctuation[0:-1:n_skip_e_tot_fluctuation,1],
                    MarkerEnergyFluctuation, markersize = MarkerSizeEnergyFluctuation, color = order2Color, label = 'order 2')
ax1.set_xlabel(r'$\Delta{}t$')
ax1.set_ylabel(r'$\delta$')
ax1.set_title(r'Energy fluctuation (order 2)')

ax2 = e_tot_fluctuation_fig.add_subplot(122)
e_tot_fluctuation_adaptive = np.zeros(np.shape(e_tot_full_orbit_adaptive))
e_tot_fluctuation_adaptive = np.delete(e_tot_fluctuation_adaptive,-1,0)
e_tot_fluctuation_adaptive[:,0] = abs(e_tot_full_orbit_adaptive[1:,0] - e_tot_full_orbit_adaptive[0:-1,0])
e_tot_fluctuation_adaptive[:,1] = abs(e_tot_full_orbit_adaptive[1:,1]/e_tot_full_orbit_adaptive[0:-1,1] -1)
e_tot_fluctuation_plot_adaptive = ax2.plot(e_tot_fluctuation_adaptive[0:-1:n_skip_e_tot_fluctuation,0],e_tot_fluctuation_adaptive[0:-1:n_skip_e_tot_fluctuation,1],
                    MarkerEnergyFluctuation, markersize = MarkerSizeEnergyFluctuation, color = adaptiveColor, label = 'order 2 (adaptive)')
ax2.set_xlabel(r'$\Delta{}t$')
ax2.set_ylabel(r'$\delta$')
ax2.set_title(r'Energy fluctuation (adaptive/order 2)')


plt.tight_layout()

# Save figure to file and show
poincare_fig.savefig(path_data_plots + '/results' + name_test_case + 'poincare' + extension)
J_par_fig.savefig(path_data_plots + '/results' + name_test_case + 'Jpar' + extension)
e_tot_fig.savefig(path_data_plots + '/results' + name_test_case + 'etot' + extension)
e_tot_fluctuation_fig.savefig(path_data_plots + '/results' + name_test_case + 'etotfluctuation' + extension)
plt.show()

# Go back to path of PYTHON script
os.chdir(path_script)

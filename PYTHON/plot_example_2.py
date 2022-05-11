#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:09:29 2021

Example 2:
 * Compute a collisionless guiding-center orbit with GORILLA for a passing Deuterium particle.
 * Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
 * Compare the results of GORILLA with different polynominal orders and Runge-Kutta 4.
 * Create a figure with the Poincaré plots (\varphi = 0) in cylindrical and symmetry flux coordinates.
 * Compute the normalized total energy as a function of toroidal mappings. 

@author: Michael Eder
"""

from tracemalloc import stop
import f90nml
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# Initialize used paths and files
# -----------------------------------------------------------------------------------------------------------------------

# Name of the current calculation to create folder
name_test_case = 'example_2'

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

# Polynomial order for orbit pusher (from 2 to 4)
gorilla['gorillanml']['poly_order'] = 2


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
gorilla_plot['gorilla_plot_nml']['total_orbit_time'] = 2

# Total Energy of particle in eV
gorilla_plot['gorilla_plot_nml']['energy_eV_start'] = 3000

# Switch for plotting Poincaré cuts at toroidal variable $\varphi$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_phi_0'] = True

# Number of skipped (non-printed) Poincaré cuts at parallel velocity $v_\parallel$ = 0
gorilla_plot['gorilla_plot_nml']['n_skip_phi_0'] = 10

# Filename for Poincaré cuts at parallel velocity $v_\parallel$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_rphiz'] = 'poincare_plot_phi_0_rphiz_order2.dat'

# Filename for Poincaré cuts at parallel velocity $v_\parallel$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_sthetaphi'] = 'poincare_plot_phi_0_sthetaphi_order2.dat'

# Switch for plotting Poincaré cuts at parallel velocity $v_\parallel$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_vpar_0'] = False

# Switch for plotting full orbit
gorilla_plot['gorilla_plot_nml']['boole_full_orbit'] = False

# Plot invariances of motion (ONLY for single orbits)

# Switch for plotting total particle energy
gorilla_plot['gorilla_plot_nml']['boole_e_tot'] = True

# Filename for total particle energy
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot_order2.dat'

# Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
gorilla_plot['gorilla_plot_nml']['boole_p_phi'] = False

# Switch for parallel adiabatic invariant $J_\parallel$
gorilla_plot['gorilla_plot_nml']['boole_J_par'] = False

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
gorilla_plot['gorilla_plot_nml']['start_pitch_parameter'] = 0.7


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
  stop

os.system('ln -s ../../../MHD_EQUILIBRIA .')
os.system('ln -s ../../../INPUT/field_divB0.inp .')
os.system('ln -s ../../../INPUT/preload_for_SYNCH.inp .')

# Run GORILLA code
os.system('./test_gorilla_main.x')

# Repeat orbit calculation for different polynominal orders and Runge-Kutta 4
# polynominal order = 3
gorilla['gorillanml']['poly_order'] = 3
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot_order3.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_rphiz'] = 'poincare_plot_phi_0_rphiz_order3.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_sthetaphi'] = 'poincare_plot_phi_0_sthetaphi_order3.dat'
gorilla.write(path_RUN + '/gorilla.inp', force = True)
gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)
os.system('./test_gorilla_main.x')

# polynominal order = 4
gorilla['gorillanml']['poly_order'] = 4
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot_order4.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_rphiz'] = 'poincare_plot_phi_0_rphiz_order4.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_sthetaphi'] = 'poincare_plot_phi_0_sthetaphi_order4.dat'
gorilla.write(path_RUN + '/gorilla.inp', force = True)
gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)
os.system('./test_gorilla_main.x')

# numerical RK4
gorilla['gorillanml']['ipusher'] = 1
gorilla['gorillanml']['boole_pusher_ode45'] = False
gorilla_plot['gorilla_plot_nml']['filename_e_tot'] = 'e_tot_rk4.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_rphiz'] = 'poincare_plot_phi_0_rphiz_rk4.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_sthetaphi'] = 'poincare_plot_phi_0_sthetaphi_rk4.dat'
gorilla.write(path_RUN + '/gorilla.inp', force = True)
gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)
os.system('./test_gorilla_main.x')


# Create plots of generated data
# -----------------------------------------------------------------------------------------------------------------------

# Here the absolut path to the RUN folder is needed
path_script2RUN = path_RUN
extension = '.png'

# Load Poincaré cut data in symmetry flux coordinates (s,\vartheta,\varphi)
poincare_sthetaphi_order2 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_sthetaphi_order2.dat')
poincare_sthetaphi_order3 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_sthetaphi_order3.dat')
poincare_sthetaphi_order4 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_sthetaphi_order4.dat')
poincare_sthetaphi_rk4 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_sthetaphi_rk4.dat')

# Load Poincaré cut data in cylindrical coordinates (R,\varphi,Z)
poincare_rphiz_order2 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_rphiz_order2.dat')
poincare_rphiz_order3 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_rphiz_order3.dat')
poincare_rphiz_order4 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_rphiz_order4.dat')
poincare_rphiz_rk4 = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_rphiz_rk4.dat')

# Load total particle energy data 
e_tot_order2 = np.genfromtxt(path_script2RUN + '/' + 'e_tot_order2.dat')
e_tot_order3 = np.genfromtxt(path_script2RUN + '/' + 'e_tot_order3.dat')
e_tot_order4 = np.genfromtxt(path_script2RUN + '/' + 'e_tot_order4.dat')
e_tot_rk4 = np.genfromtxt(path_script2RUN + '/' + 'e_tot_rk4.dat')

# Normalize total particle energy to the value at t = 0
e_tot_order2[:,1] = e_tot_order2[:,1] / e_tot_order2[0,1]
e_tot_order3[:,1] = e_tot_order3[:,1] / e_tot_order3[0,1] 
e_tot_order4[:,1] = e_tot_order4[:,1] / e_tot_order4[0,1] 
e_tot_rk4[:,1] = e_tot_rk4[:,1] / e_tot_rk4[0,1]  

# Plot Poincaré sections and evolution of the normalized total particle energy as a function of toroidal mappings
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.set_size_inches(18.5,10.5)
fig.suptitle('Stellarator: Poincaré sections and total energy of passing particle')


ax1.plot(poincare_rphiz_order2[:,0],poincare_rphiz_order2[:,2],'s',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax1.plot(poincare_rphiz_order3[:,0],poincare_rphiz_order3[:,2],'d',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax1.plot(poincare_rphiz_order4[:,0],poincare_rphiz_order4[:,2],'v',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax1.plot(poincare_rphiz_rk4[:,0],poincare_rphiz_rk4[:,2],'^',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax1.set_xlabel('$R$ [cm]')
ax1.set_ylabel('$Z$ [cm]')
ax1.legend(['GORILLA Poly2','GORILLA Poly3','GORILLA Poly4','GORILLA RK4'])
ax1.set_title(r'Poincaré $\varphi = 0$')


ax2.plot(poincare_sthetaphi_order2[:,1],poincare_sthetaphi_order2[:,0],'s',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax2.plot(poincare_sthetaphi_order3[:,1],poincare_sthetaphi_order3[:,0],'d',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax2.plot(poincare_sthetaphi_order4[:,1],poincare_sthetaphi_order4[:,0],'v',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax2.plot(poincare_sthetaphi_rk4[:,1],poincare_sthetaphi_rk4[:,0],'^',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax2.set_xlabel(r'$\vartheta$')
ax2.set_ylabel('$s$')
ax2.legend(['GORILLA Poly2','GORILLA Poly3','GORILLA Poly4','GORILLA RK4'])
ax2.set_title(r'Poincaré $\varphi = 0$')


ax3.plot(np.abs(e_tot_order2[:,0]),e_tot_order2[:,1],'s',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax3.plot(np.abs(e_tot_order3[:,0]),e_tot_order3[:,1],'d',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax3.plot(np.abs(e_tot_order4[:,0]),e_tot_order4[:,1],'v',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax3.plot(np.abs(e_tot_rk4[:,0]),e_tot_rk4[:,1],'^',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax3.set_xlabel('$N_{mappings}$')
ax3.set_ylabel('$E_{tot} / E_{tot}(t=0)$')
ax3.legend(['GORILLA Poly2','GORILLA Poly3','GORILLA Poly4','GORILLA RK4'])
ax3.set_title('Normalized total energy')

plt.tight_layout()

# Save figure to file and show
fig.savefig(path_data_plots + '/results' + name_test_case + extension)
plt.show()

# Go back to path of PYTHON script
os.chdir(path_script)

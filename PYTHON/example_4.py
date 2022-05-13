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
name_test_case = 'example_4'

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
gorilla['gorillanml']['poly_order'] = 4


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
tetra_grid['tetra_grid_nml']['grid_kind'] = 2

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
gorilla_plot['gorilla_plot_nml']['i_orbit_options'] = 3

# Total individual orbit flight time for plotting
gorilla_plot['gorilla_plot_nml']['total_orbit_time'] = 2

# Total Energy of particle in eV
gorilla_plot['gorilla_plot_nml']['energy_eV_start'] = 3000

# Switch for plotting Poincaré cuts at toroidal variable $\varphi$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_phi_0'] = True

# Number of skipped (non-printed) Poincaré cuts at parallel velocity $v_\parallel$ = 0
gorilla_plot['gorilla_plot_nml']['n_skip_phi_0'] = 100

# Filename for Poincaré cuts at parallel velocity $v_\parallel$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_rphiz'] = 'poincare_plot_phi_0_rphiz.dat'

# Filename for Poincaré cuts at parallel velocity $v_\parallel$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_sthetaphi'] = 'poincare_plot_phi_0_sthetaphi.dat'

# Switch for plotting Poincaré cuts at parallel velocity $v_\parallel$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_vpar_0'] = False

# Switch for plotting full orbit
gorilla_plot['gorilla_plot_nml']['boole_full_orbit'] = False

# Plot invariances of motion (ONLY for single orbits)

# Switch for plotting total particle energy
gorilla_plot['gorilla_plot_nml']['boole_e_tot'] = False

# Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
gorilla_plot['gorilla_plot_nml']['boole_p_phi'] = False

# Switch for parallel adiabatic invariant $J_\parallel$
gorilla_plot['gorilla_plot_nml']['boole_J_par'] = False

# Starting positions of particles (BOTH single and multiple orbits) and starting pitch parameter
# File (either rphiz or sthetaphi) is automatically chosen dependent on coord_system in 'gorilla.inp')

# Filename for list of starting position(s) of particle(s) in cylindrical coordinates (R,$\varphi$,Z) and pitch ($\lambda$)
gorilla_plot['gorilla_plot_nml']['filename_orbit_start_pos_rphiz'] = 'orbit_start_rphizlambda.dat'

# Filename for list of starting position(s) of particle(s) in symmetry flux coordinates (s,$\vartheta$,$\varphi$) and pitch ($\lambda$)
gorilla_plot['gorilla_plot_nml']['filename_orbit_start_pos_sthetaphi'] = 'orbit_start_sthetaphilambda.dat'

# Create inputfile for starting positions
# s,theta,phi,lambda
start_pos_1=[0.8,0.0,0.63,0.7]
start_pos_2=[0.6,0.0,0.63,0.7]
start_pos_3=[0.4,0.0,0.63,0.2]
start_pos=[start_pos_1,start_pos_2,start_pos_3]
fileID = open(path_RUN + '/' + gorilla_plot['gorilla_plot_nml']['filename_orbit_start_pos_sthetaphi'],'w')
for row in start_pos:
    fileID.write(" ".join(map(lambda n: '%.8f'%n, row)))
    fileID.write("\n")
fileID.close()


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


# Create plots of generated data
# -----------------------------------------------------------------------------------------------------------------------

# Here the absolut path to the RUN folder is needed
path_script2RUN = path_RUN
extension = '.png'

# Load Poincaré cut data in symmetry flux coordinates (s,\vartheta,\varphi)
poincare_sthetaphi = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_sthetaphi.dat')

# Load Poincaré cut data in cylindrical coordinates (R,\varphi,Z)
poincare_rphiz = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_rphiz.dat') 

# Plot Poincaré cuts and for two passing and one trapped Deuterium particle.
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.set_size_inches(18.5,10.5)
fig.suptitle('Tokamak: Poincaré plots for two passing and one trapped Deuterium particle.')


ax1.plot(poincare_rphiz[:,0],poincare_rphiz[:,2],'s',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax1.set_xlabel('$R$ [cm]')
ax1.set_ylabel('$Z$ [cm]')
ax1.set_title(r'Poincaré $\varphi = 0$')


ax2.plot(poincare_sthetaphi[:,1],poincare_sthetaphi[:,0],'s',markersize=3,markeredgewidth=0.5,markerfacecolor="None")
ax2.set_xlabel(r'$\vartheta$')
ax2.set_ylabel('$s$')
ax2.set_title(r'Poincaré $\varphi = 0$')
ax2.set_xlim(0,2*np.pi)


plt.tight_layout()

# Save figure to file and show
fig.savefig(path_data_plots + '/results' + name_test_case + extension)
plt.show()

# Go back to path of PYTHON script
os.chdir(path_script)

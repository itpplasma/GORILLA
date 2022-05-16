#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 15:34:30 2022

Example 6:
 * Compute collisionless guiding-center orbits with GORILLA for passing and trapped Deuterium particles (WEST geometry).
 * Construct a 3D extension of the Soledge3x-EIRENE 2D-mesh for an axisymmetric tokamak equilibrium (g-file)
 * Plot the 2D projection of the guiding-center orbits on the original Soledge3x-EIRENE grid for both particle types.

@author: Georg Graßler
"""

from tracemalloc import stop
from auxillary_functions import boundary_search
import f90nml
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# Initialize used paths and files
# -----------------------------------------------------------------------------------------------------------------------

# Name of the current calculation to create folder
name_test_case = 'example_6'

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
gorilla['gorillanml']['eps_Phi'] = -1e-7

# Coordinate system
# 1 ... (R,phi,Z) cylindrical coordinate system
# 2 ... (s,theta,phi) symmetry flux coordinate system
gorilla['gorillanml']['coord_system'] = 1

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
tetra_grid['tetra_grid_nml']['n2'] = 60
# Rectangular: nZ, Field-aligned: ntheta
tetra_grid['tetra_grid_nml']['n3'] = 60

# Grid kind
# 1 ... rectangular grid for axisymmetric EFIT data
# 2 ... field-aligned grid for axisymmetric EFIT data
# 3 ... field-aligned grid for non-axisymmetric VMEC
# 4 ... SOLEDGE3X_EIRENE grid
tetra_grid['tetra_grid_nml']['grid_kind'] = 4

# MHD equilibrium filename
tetra_grid['tetra_grid_nml']['g_file_filename'] = 'MHD_EQUILIBRIA/g_file_for_test_WEST'
tetra_grid['tetra_grid_nml']['convex_wall_filename'] = 'MHD_EQUILIBRIA/convex_wall_for_test_WEST.dat'


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
gorilla_plot['gorilla_plot_nml']['total_orbit_time'] = 0.004

# Total Energy of particle in eV
gorilla_plot['gorilla_plot_nml']['energy_eV_start'] = 3000

# Switch for plotting Poincaré cuts at toroidal variable $\varphi$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_phi_0'] = False

# Switch for plotting Poincaré cuts at parallel velocity $v_\parallel$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_vpar_0'] = False

# Switch for plotting full orbit
gorilla_plot['gorilla_plot_nml']['boole_full_orbit'] = True

# Number of skipped (non-printed tetrahedra passings) full orbit
gorilla_plot['gorilla_plot_nml']['n_skip_full_orbit'] = 1

# Filename for full orbit in cylindrical coordinates (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['filename_full_orbit_rphiz'] = 'full_orbit_plot_rphiz_passing.dat'

# Filename for full orbit in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
gorilla_plot['gorilla_plot_nml']['filename_full_orbit_sthetaphi'] = 'full_orbit_plot_sthetaphi_passing.dat'


# Plot invariances of motion (ONLY for single orbits)

# Switch for plotting total particle energy
gorilla_plot['gorilla_plot_nml']['boole_e_tot'] = False

# Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
gorilla_plot['gorilla_plot_nml']['boole_p_phi'] = False

# Switch for parallel adiabatic invariant $J_\parallel$
gorilla_plot['gorilla_plot_nml']['boole_J_par'] = False

# Multiple orbits: starting positions taken from file (i_orbit_options = 3)

# Filename for list of starting position(s) of particle(s) in cylindrical coordinates (R,$\varphi$,Z) and pitch ($\lambda$)
gorilla_plot['gorilla_plot_nml']['filename_orbit_start_pos_rphiz'] = 'orbit_start_rphizlambda_passing.dat'

# Create inputfile for starting positions
# R,phi,Z,lambda
R_left_values = np.linspace(192.245,205.714,7).tolist()
Z_left_values = np.linspace(57.143,40.000,7).tolist()
R_inner_left_values = [207.959]
Z_inner_left_values = [37.143]
R_inner_right_values = [266.327,268.571]
Z_inner_right_values = [-37.143,-40.000]
R_right_values = np.linspace(270.816,291.020,10).tolist()
Z_right_values = np.linspace(-42.857,-68.571,10).tolist()
R_values = R_left_values + R_inner_left_values + R_inner_right_values + R_right_values
Z_values = Z_left_values + Z_inner_left_values + Z_inner_right_values + Z_right_values

lambda_values = np.empty(2)
lambda_values[0] = 0.9

fileID = open(path_RUN + '/' + gorilla_plot['gorilla_plot_nml']['filename_orbit_start_pos_rphiz'],'w')
for k in range(len(R_values)):
  entry = [R_values[k],0,Z_values[k],lambda_values[0]]
  fileID.write(" ".join(map(lambda n: '%.8f'%n, entry)))
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


# Second run for trapped particle
# -----------------------------------------------------------------------------------------------------------------------

# Changes in Input files
gorilla_plot['gorilla_plot_nml']['filename_full_orbit_rphiz'] = 'full_orbit_plot_rphiz_trapped.dat'
gorilla_plot['gorilla_plot_nml']['filename_full_orbit_sthetaphi'] = 'full_orbit_plot_sthetaphi_trapped.dat'
gorilla_plot['gorilla_plot_nml']['filename_orbit_start_pos_rphiz'] = 'orbit_start_rphizlambda_trapped.dat'
lambda_values[1] = 0.3
R_values[9] +=  0.5
del R_values[10]
del R_values[7]
Z_values[9] += -0.5
del Z_values[10]
del Z_values[7]

fileID = open(path_RUN + '/' + gorilla_plot['gorilla_plot_nml']['filename_orbit_start_pos_rphiz'],'w')
for k in range(len(R_values)):
  entry = [R_values[k],0,Z_values[k],lambda_values[1]]
  fileID.write(" ".join(map(lambda n: '%.8f'%n, entry)))
  fileID.write("\n")
fileID.close()

gorilla_plot.write(path_RUN + '/gorilla_plot.inp', force = True)

# Run GORILLA code
os.system('./test_gorilla_main.x')


# Create plots of generated data
# -----------------------------------------------------------------------------------------------------------------------

# Here the absolut path to the RUN folder is needed
path_script2RUN = path_RUN
extension = '.png'

# Names for particle calculation in plot
particle_type = ['Passing','Trapped']

# Read in SOLEDGE3X-EIRENE mesh data
coordinates = []
fileID = open('MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/knots_for_test.dat','r')
desc = fileID.readline()
s = desc.split()
n_vertex = int(s[0])
for i in range(n_vertex):
  desc = fileID.readline()
  s = desc.split()
  for k in s:
    coordinates.append(float(k))
fileID.close()
coordinates = np.array(coordinates)
coordinates = np.reshape(coordinates,(n_vertex,2))

triangles = []
fileID = open('MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/triangles_for_test.dat','r')
desc = fileID.readline()
s = desc.split()
n_triangles = int(s[0])
for i in range(n_triangles):
  desc = fileID.readline()
  s = desc.split()
  for k in s:
    triangles.append(int(k))
fileID.close()
triangles = np.array(triangles)
triangles = np.reshape(triangles,(n_triangles,3))

# Prepare 2D grid for plotting (NAN seperating the individual triangles -> one single grafic object)
grid = float("NAN")*np.ones((n_triangles*5,2))
for t in range(n_triangles):
    for k in range(3):
        grid[t*5 + k,:] = coordinates[triangles[t,k]-1,:]
    grid[t*5 + 3,:] = coordinates[triangles[t,0]-1,:]
    grid[t*5 + 4,:] = float("NAN")

# Colorscheme
grid_color = np.array([204,204,204])/256
grid_thickness = 0.3
orbit_color = np.array([4,90,141])/256
orbit_thickness = 2

# Loop over GORILLA calculations
for i in range(2):
  fig = plt.figure(figsize=(15, 12))
  ax = fig.add_subplot(111)

  # Data for passing or trapped paricle
  if (i == 0):
      full_orbit_rphiz = np.genfromtxt(path_script2RUN + '/' + 'full_orbit_plot_rphiz_passing.dat') 
      start_pos = np.genfromtxt(path_script2RUN + '/' + 'orbit_start_rphizlambda_passing.dat') 
  elif (i==1):
      full_orbit_rphiz = np.genfromtxt(path_script2RUN + '/' + 'full_orbit_plot_rphiz_trapped.dat') 
      start_pos = np.genfromtxt(path_script2RUN + '/' + 'orbit_start_rphizlambda_trapped.dat')  

  # Plot the grid
  grid_plot = ax.plot(grid[:,0],grid[:,1],color = grid_color,linewidth = grid_thickness, 
                      label = 'SOLEDGE3X-EIRENE mesh')

  # Plot the guiding-center orbit projection
  orbit_plot = ax.plot(full_orbit_rphiz[:,0],full_orbit_rphiz[:,2],'.',color = orbit_color,
                        markersize = orbit_thickness, label = 'Guiding-center orbit-projection')

  # Plot the starting positions
  start_pos_plot = ax.plot(start_pos[:,0],start_pos[:,2],'o',markerfacecolor = orbit_color,
                            markeredgecolor = orbit_color, label = 'Guiding-center starting postions')

  # Legend, labels and title
  lh = ax.legend(loc = 'lower left')
  ax.set_xlabel('R [cm]')
  ax.set_ylabel('Z [cm]')
  ax.set_title(r'SOLEDGE3X-EIRENE: Toroidal projection ($\varphi=0$) of guiding-center orbits ($\lambda =$ ' + str(lambda_values[i]) + ')')

  # Save figure to file and show
  fig.savefig(path_data_plots + '/results' + name_test_case + particle_type[i] + extension)

plt.show()

# Go back to path of PYTHON script
os.chdir(path_script)
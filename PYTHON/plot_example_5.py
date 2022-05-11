#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 17:32:30 2022

Example 5:
 * Compute collisionless guiding-center orbits with GORILLA for a passing and a trapped Deuterium particle.
 * Use a field-aligned grid for an axisymmetric tokamak equilibrium (g-file).
 * Plot the plasma boundary, the guiding-center orbits, and the resulting Poincare plot (\varphi = 0) for both orbits.

@author: Georg Graßler
"""

from cmath import pi
from tracemalloc import stop
from auxillary_functions import boundary_search
import f90nml
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
mpl.rcdefaults()


# Initialize used paths and files
# -----------------------------------------------------------------------------------------------------------------------

# Name of the current calculation to create folder
name_test_case = 'example_5'

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
tetra_grid['tetra_grid_nml']['n2'] = 40
# Rectangular: nZ, Field-aligned: ntheta
tetra_grid['tetra_grid_nml']['n3'] = 40

# Grid kind
# 1 ... rectangular grid for axisymmetric EFIT data
# 2 ... field-aligned grid for axisymmetric EFIT data
# 3 ... field-aligned grid for non-axisymmetric VMEC
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
gorilla_plot['gorilla_plot_nml']['i_orbit_options'] = 2

# Total individual orbit flight time for plotting
gorilla_plot['gorilla_plot_nml']['total_orbit_time'] = 0.004

# Total Energy of particle in eV
gorilla_plot['gorilla_plot_nml']['energy_eV_start'] = 3000

# Switch for plotting Poincaré cuts at toroidal variable $\varphi$ = 0
gorilla_plot['gorilla_plot_nml']['boole_poincare_phi_0'] = True

# Number of skipped (non-printed) Poincaré cuts at parallel velocity $v_\parallel$ = 0
gorilla_plot['gorilla_plot_nml']['n_skip_phi_0'] = 1

# Filename for Poincaré cuts at parallel velocity $v_\parallel$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_rphiz'] = 'poincare_plot_phi_0_rphiz_passing.dat'

# Filename for Poincaré cuts at parallel velocity $v_\parallel$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_sthetaphi'] = 'poincare_plot_phi_0_sthetaphi_passing.dat'

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

# Single orbit from starting drift surface (i_orbit_options = 2)

# Starting drift surface
# = s for (s,$\vartheta$,$\varphi$)
# = R for (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['start_pos_x1_beg'] = 0.5

# End drift surface
# = s for (s,$\vartheta$,$\varphi$)
# = R for (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['start_pos_x1_end'] = 0.9

# Number of drift surfaces in between start and end
gorilla_plot['gorilla_plot_nml']['n_surfaces'] = 30

# Starting value for toroidal variable
# = $\vartheta$ for (s,$\vartheta$,$\varphi$)
# = $\varphi$ for (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['start_pos_x2'] = 0.1

# Starting value for poloidal variable $\vartheta$
# = $\varphi$ for (s,$\vartheta$,$\varphi$)
# = Z for (R,$\varphi$,Z)
gorilla_plot['gorilla_plot_nml']['start_pos_x3']  = 3.63

# Pitch parameter $\lambda$ = $v_\parallel$ / vmod
gorilla_plot['gorilla_plot_nml']['start_pitch_parameter'] = 0.9


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
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_rphiz'] = 'poincare_plot_phi_0_rphiz_trapped.dat'
gorilla_plot['gorilla_plot_nml']['filename_poincare_phi_0_sthetaphi'] = 'poincare_plot_phi_0_sthetaphi_trapped.dat'
gorilla_plot['gorilla_plot_nml']['filename_full_orbit_rphiz'] = 'full_orbit_plot_rphiz_trapped.dat'
gorilla_plot['gorilla_plot_nml']['filename_full_orbit_sthetaphi'] = 'full_orbit_plot_sthetaphi_trapped.dat'
gorilla_plot['gorilla_plot_nml']['start_pitch_parameter'] = 0.4

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

# Number of points shown from full_orbit_files [passing, trapped]
n_show_points = [290,300]

# Loop over GORILLA calculations
for i in range(2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Data for passing or trapped paricle
    if (i == 0):
        full_orbit_rphiz = np.genfromtxt(path_script2RUN + '/' + 'full_orbit_plot_rphiz_passing.dat') 
        poincare_rphiz = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_rphiz_passing.dat') 
    elif (i==1):
        full_orbit_rphiz = np.genfromtxt(path_script2RUN + '/' + 'full_orbit_plot_rphiz_trapped.dat') 
        poincare_rphiz = np.genfromtxt(path_script2RUN + '/' + 'poincare_plot_phi_0_rphiz_trapped.dat') 

    # XYZ Data of orbit for plot
    full_orbit_xyz = full_orbit_rphiz.copy()
    full_orbit_xyz[:,0] = full_orbit_rphiz[:,0]*np.cos(full_orbit_rphiz[:,1])
    full_orbit_xyz[:,1] = full_orbit_rphiz[:,0]*np.sin(full_orbit_rphiz[:,1])

    poincare_xyz = poincare_rphiz.copy()
    poincare_xyz[:,0] = poincare_rphiz[:,0]*np.cos(poincare_rphiz[:,1])
    poincare_xyz[:,1] = poincare_rphiz[:,0]*np.sin(poincare_rphiz[:,1])

    # Used points of full orbit
    full_orbit_xyz=full_orbit_xyz[-1-n_show_points[i]:-1,:]

    # Plasma boundaries
    filename = path_RUN + '/MHD_EQUILIBRIA/g_file_for_test'
    conversion = 100
    [rbbbs,zbbbs] = boundary_search(filename,'e')
    rbbbs = rbbbs[zbbbs>(np.min(full_orbit_rphiz[:,2])-48)/conversion]
    zbbbs = zbbbs[zbbbs>(np.min(full_orbit_rphiz[:,2])-48)/conversion]
    phibbbs = np.linspace(0,1.5*np.pi,100)
    [R,PHI] = np.meshgrid(rbbbs*conversion, phibbbs)
    Z = np.tile(zbbbs*conversion, (len(phibbbs),1))
    sb = ax.plot_surface(R*np.cos(PHI),R*np.sin(PHI),Z, label = 'Plasma Boundary')
    sb.set_color((0, 0.3, 0.3, 0.3))
    sb.set_edgecolor('none')
    sb._edgecolors2d = sb._edgecolor3d
    sb._facecolors2d = sb._facecolor3d

    # Plane of Poincare plot
    xp = [np.min(poincare_rphiz[:,0])-30,np.max(poincare_rphiz[:,0])+30]
    zp = [np.min(full_orbit_rphiz[:,2])-48,np.max(full_orbit_rphiz[:,2])+30]
    x1 = np.array([ [xp[0], xp[1]], [xp[0], xp[1]]])
    y1 = np.zeros(np.shape(x1))
    z1 = np.array([ [zp[0], zp[0]], [zp[1], zp[1]]])

    v = ax.plot_surface(x1,y1,z1)
    v.set_color((0.5, 0.5, 0.5, 0.2))
    v.set_edgecolor((0.1, 0.1, 0.1, 1))

    t = ax.text(xp[1]-70,0,zp[0]+15, '%s' % (r'$\varphi=0$'), size=15, zorder=1, color='k')

    # Full orbit positons
    p0 = ax.plot(full_orbit_xyz[:,0],full_orbit_xyz[:,1],full_orbit_xyz[:,2],color ='r',linewidth=2)

    # Current particle position
    p1 = ax.scatter(full_orbit_xyz[-1,0],full_orbit_xyz[-1,1],full_orbit_xyz[-1,2],color = 'r',label= particle_type[i] + ' particle')

    # Poincare plot
    p2 = ax.scatter(poincare_xyz[:,0],poincare_xyz[:,1],poincare_xyz[:,2],'b',s = 6,label = 'Poincaré plot')

    # Legend
    lh = ax.legend(handles = [sb,p1,p2])

    # Limits of plot and view
    ax.set_xlim(-xp[1],xp[1])
    ax.set_ylim(-xp[1],xp[1])
    ax.set_zlim(-xp[1],xp[1])
    ax.view_init(elev=24.2848,azim = -57)
    ax.axis('off')
    ax.dist = 6

    # Save figure to file and show
    fig.savefig(path_data_plots + '/results' + name_test_case + particle_type[i] + extension)

plt.show()

# Go back to path of PYTHON script
os.chdir(path_script)
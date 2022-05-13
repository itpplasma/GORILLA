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

from cProfile import label
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


# Names for particle calculation in plot
particle_type = 'Passing'
start_lambda = 0.9

# Load full orbit data from example_6 in cylindrical coordinates (R,\varphi,Z)
full_orbit_rphiz = np.genfromtxt('../EXAMPLES/example_6/full_orbit_plot_rphiz_passing.dat') 

# Load Poincaré cut data from example_6 in cylindrical coordinates (R,\varphi,Z)
start_pos = np.genfromtxt('../EXAMPLES/example_6/orbit_start_rphizlambda_passing.dat')

# Read in SOLEDGE3X-EIRENE mesh data
coordinates = []
fileID = open('../MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/knots_for_test.dat','r')
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
fileID = open('../MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/triangles_for_test.dat','r')
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

# Plot the grid
fig = plt.figure(figsize=(15, 12))
ax = fig.add_subplot(111)
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
ax.set_title(r'SOLEDGE3X-EIRENE: Toroidal projection ($\varphi=0$) of guiding-center orbits ($\lambda =$ ' + str(start_lambda) + ')')

plt.show()
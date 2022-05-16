#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 18:40:14 2022

Example 5:
 * Compute collisionless guiding-center orbits with GORILLA for a trapped Deuterium particle. (Manually execute 'test_gorilla_main.x' in the corresponding example folder)
 * Use a field-aligned grid for an axisymmetric tokamak equilibrium (g-file).
 * Plot the plasma boundary, the guiding-center orbits, and the resulting Poincare plot (\varphi = 0).

@author: Georg Graßler
"""

from auxillary_functions import boundary_search
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
mpl.rcdefaults()


# Names for particle calculation in plot
particle_type = 'Trapped'

# Number of points shown from full_orbit_files [passing, trapped]
n_show_points = 300

# Load full orbit data from example_5 in cylindrical coordinates (R,\varphi,Z)
full_orbit_rphiz = np.genfromtxt('../EXAMPLES/example_5/full_orbit_plot_rphiz_trapped.dat') 

# Load Poincaré cut data from example_5 in cylindrical coordinates (R,\varphi,Z)
poincare_rphiz = np.genfromtxt('../EXAMPLES/example_5/poincare_plot_phi_0_rphiz_trapped.dat') 

# XYZ Data of orbit for plot
full_orbit_xyz = full_orbit_rphiz.copy()
full_orbit_xyz[:,0] = full_orbit_rphiz[:,0]*np.cos(full_orbit_rphiz[:,1])
full_orbit_xyz[:,1] = full_orbit_rphiz[:,0]*np.sin(full_orbit_rphiz[:,1])

poincare_xyz = poincare_rphiz.copy()
poincare_xyz[:,0] = poincare_rphiz[:,0]*np.cos(poincare_rphiz[:,1])
poincare_xyz[:,1] = poincare_rphiz[:,0]*np.sin(poincare_rphiz[:,1])

# Used points of full orbit
full_orbit_xyz=full_orbit_xyz[-1-n_show_points:-1,:]

# Plot plasma boundaries
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')


filename = '../MHD_EQUILIBRIA/g_file_for_test'
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

# Plot plane of Poincare plot
xp = [np.min(poincare_rphiz[:,0])-30,np.max(poincare_rphiz[:,0])+30]
zp = [np.min(full_orbit_rphiz[:,2])-48,np.max(full_orbit_rphiz[:,2])+30]
x1 = np.array([ [xp[0], xp[1]], [xp[0], xp[1]]])
y1 = np.zeros(np.shape(x1))
z1 = np.array([ [zp[0], zp[0]], [zp[1], zp[1]]])


v = ax.plot_surface(x1,y1,z1)
v.set_color((0.5, 0.5, 0.5, 0.2))
v.set_edgecolor((0.1, 0.1, 0.1, 1))


t = ax.text(xp[1]-70,0,zp[0]+15, '%s' % (r'$\varphi=0$'), size=15, zorder=1, color='k')

# Plot full orbit
p0 = ax.plot(full_orbit_xyz[:,0],full_orbit_xyz[:,1],full_orbit_xyz[:,2],color ='r',linewidth=2)

# Plot current particle position
p1 = ax.scatter(full_orbit_xyz[-1,0],full_orbit_xyz[-1,1],full_orbit_xyz[-1,2],color = 'r',label= particle_type + ' particle')

# Plot Poincare plot
p2 = ax.scatter(poincare_xyz[:,0],poincare_xyz[:,1],poincare_xyz[:,2],'b',s = 6,label = 'Poincaré plot')

# Legend, limits and view
lh = ax.legend(handles = [sb,p1,p2])
ax.set_xlim(-xp[1],xp[1])
ax.set_ylim(-xp[1],xp[1])
ax.set_zlim(-xp[1],xp[1])
ax.view_init(elev=24.2848,azim = -57)
ax.axis('off')
ax.dist = 6.5

plt.show()
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,'tools'); from gfile import Equilibrium
d=sys.argv[1]; box=[float(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 and sys.argv[2] else None
P=np.loadtxt(d+'/knots_aug_edge.dat',skiprows=1); T=np.loadtxt(d+'/triangles_aug_edge.dat',skiprows=1,dtype=int)-1
reg=np.loadtxt(d+'/triangle_region.dat',dtype=int)
eq=Equilibrium('GORILLA/MHD_EQUILIBRIA/g_file_for_test')
fig,ax=plt.subplots(figsize=(7,10))
for r,c in ((0,'tab:blue'),(1,'tab:red')):
    ax.add_collection(PolyCollection(P[T[reg==r]],facecolors='none',edgecolors=c,linewidths=0.3))
ax.plot(*eq.lim.T,'k-',lw=0.8); ax.plot(*eq.xpt,'kx'); ax.set_aspect('equal')
if box: ax.set_xlim(box[0],box[1]); ax.set_ylim(box[2],box[3])
else: ax.autoscale_view()
out=sys.argv[3] if len(sys.argv)>3 else d+'/mesh.png'
fig.tight_layout(); fig.savefig(out,dpi=130)

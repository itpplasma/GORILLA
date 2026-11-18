import numpy as np
from collections import defaultdict
def folds(P, T):
    """Interior edges whose two triangles lie on the same side (mesh folds)."""
    E = defaultdict(list)
    for it, t in enumerate(T):
        for k in range(3):
            a, b, c = t[k], t[(k+1)%3], t[(k+2)%3]
            E[tuple(sorted((a, b)))].append(c)
    bad = []
    for (a, b), opp in E.items():
        if len(opp) != 2: continue
        d = P[b]-P[a]
        s = [d[0]*(P[o][1]-P[a][1]) - d[1]*(P[o][0]-P[a][0]) for o in opp]
        if s[0]*s[1] >= 0: bad.append((a, b, opp))
    return bad
if __name__ == '__main__':
    import sys
    d = sys.argv[1]
    P = np.loadtxt(d+'/knots_aug_edge.dat', skiprows=1); T = np.loadtxt(d+'/triangles_aug_edge.dat', skiprows=1, dtype=int)-1
    b = folds(P, T); print('folded edges:', len(b))
    v = np.array([e[0] for e in b]); 
    if len(b): print('first:', b[:5]); print('R,Z:', P[v[:5]])

"""Minimal EFIT g-file reader plus psi spline (lengths in cm, like GORILLA)."""
import re
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import fsolve


def _numbers(lines):
    txt = " ".join(lines)
    # EFIT writes 5e16.9 without separators; split by regex for floats
    return [float(x) for x in re.findall(r"[-+]?\d*\.\d+(?:[eEdD][-+]?\d+)?", txt.replace("D", "E"))]


def read_gfile(path):
    with open(path) as f:
        lines = f.read().splitlines()
    head = lines[0].split()
    nw, nh = int(head[-2]), int(head[-1])
    # Numbers in fixed width 16 chars per entry
    body = "".join(l.ljust(80)[:80] for l in lines[1:])
    vals = []
    pos = 0
    # first parse floats in fixed-width chunks until the limiter counts line
    flt = []
    for l in lines[1:]:
        flt.append(l)
    def fixed(line_iter, n):
        out = []
        while len(out) < n:
            l = next(line_iter)
            for k in range(0, len(l.rstrip()), 16):
                s = l[k:k + 16].strip()
                if s:
                    out.append(float(s.replace("D", "E")))
        return out
    it = iter(lines[1:])
    rdim, zdim, rcentr, rleft, zmid = fixed(it, 5)
    rmaxis, zmaxis, simag, sibry, bcentr = fixed(it, 5)
    current, simag2, xdum, rmaxis2, xdum = fixed(it, 5)
    zmaxis2, xdum, sibry2, xdum, xdum = fixed(it, 5)
    fpol = fixed(it, nw)
    pres = fixed(it, nw)
    ffprim = fixed(it, nw)
    pprime = fixed(it, nw)
    psirz = np.array(fixed(it, nw * nh)).reshape(nh, nw)
    qpsi = fixed(it, nw)
    nbbbs, limitr = [int(x) for x in next(it).split()]
    bbbs = np.array(fixed(it, 2 * nbbbs)).reshape(nbbbs, 2)
    lim = np.array(fixed(it, 2 * limitr)).reshape(limitr, 2)
    R = rleft + np.linspace(0, rdim, nw)
    Z = zmid - zdim / 2 + np.linspace(0, zdim, nh)
    return dict(nw=nw, nh=nh, R=R * 100, Z=Z * 100, psi=psirz, simag=simag, sibry=sibry,
                rmaxis=rmaxis * 100, zmaxis=zmaxis * 100, bbbs=bbbs * 100, lim=lim * 100,
                fpol=np.array(fpol), qpsi=np.array(qpsi), bcentr=bcentr, rcentr=rcentr * 100)


class Equilibrium:
    def __init__(self, path):
        g = read_gfile(path)
        self.g = g
        self.spl = RectBivariateSpline(g["R"], g["Z"], g["psi"].T, kx=3, ky=3)
        self.R, self.Z = g["R"], g["Z"]
        self.lim = g["lim"]
        self.psi_ax = float(self.spl(g["rmaxis"], g["zmaxis"])[0, 0])
        self.raxis, self.zaxis = g["rmaxis"], g["zmaxis"]
        self.xpt = self.find_xpoint()
        self.psi_x = float(self.spl(*self.xpt)[0, 0])

    def psi(self, r, z):
        return self.spl.ev(r, z)

    def psin(self, r, z):
        return (self.spl.ev(r, z) - self.psi_ax) / (self.psi_x - self.psi_ax)

    def grad(self, r, z):
        return self.spl.ev(r, z, dx=1), self.spl.ev(r, z, dy=1)

    def find_xpoint(self):
        # start guess: lowest point of the EFIT boundary
        b = self.g["bbbs"]
        i = np.argmin(b[:, 1])
        x0 = b[i]
        f = lambda x: [self.spl.ev(x[0], x[1], dx=1), self.spl.ev(x[0], x[1], dy=1)]
        x = fsolve(f, x0, xtol=1e-12)
        return np.array(x)

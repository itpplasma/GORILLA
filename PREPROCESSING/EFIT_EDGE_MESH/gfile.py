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
        self.xpt = check_lower_single_null(self)
        self.psi_x = float(self.spl(*self.xpt)[0, 0])

    def psi(self, r, z):
        return self.spl.ev(r, z)

    def psin(self, r, z):
        return (self.spl.ev(r, z) - self.psi_ax) / (self.psi_x - self.psi_ax)

    def grad(self, r, z):
        return self.spl.ev(r, z, dx=1), self.spl.ev(r, z, dy=1)

    def critical_points(self, step=0.5):
        """Points inside the wall contour where grad psi = 0, as (R, Z, kind) with kind 'O' (extremum)
        or 'X' (saddle). Candidates are local minima of |grad psi|^2 on a grid of the given step [cm],
        refined with fsolve."""
        from matplotlib.path import Path
        lim_path = Path(self.lim)
        R = np.arange(self.lim[:, 0].min(), self.lim[:, 0].max(), step)
        Z = np.arange(self.lim[:, 1].min(), self.lim[:, 1].max(), step)
        RR, ZZ = np.meshgrid(R, Z)
        g2 = self.spl.ev(RR, ZZ, dx=1) ** 2 + self.spl.ev(RR, ZZ, dy=1) ** 2
        gscale = np.sqrt(g2.max())
        inner = g2[1:-1, 1:-1]
        is_min = np.ones_like(inner, bool)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di or dj:
                    is_min &= inner <= g2[1 + di:g2.shape[0] - 1 + di, 1 + dj:g2.shape[1] - 1 + dj]
        cand = np.column_stack([RR[1:-1, 1:-1][is_min], ZZ[1:-1, 1:-1][is_min]])
        grad = lambda x: [self.spl.ev(x[0], x[1], dx=1), self.spl.ev(x[0], x[1], dy=1)]
        found = []
        for c in cand:
            x, _, ier, _ = fsolve(grad, c, full_output=True, xtol=1e-12)
            if ier != 1 or not lim_path.contains_point(x) or np.hypot(*grad(x)) > 1e-6 * gscale:
                continue
            if any(np.hypot(x[0] - f[0], x[1] - f[1]) < 2 * step for f in found):
                continue
            det = (self.spl.ev(x[0], x[1], dx=2) * self.spl.ev(x[0], x[1], dy=2)
                   - self.spl.ev(x[0], x[1], dx=1, dy=1) ** 2)
            found.append((float(x[0]), float(x[1]), "X" if det < 0 else "O"))
        return found


class TopologyError(SystemExit):
    pass


def check_lower_single_null(eq, tol_psin=0.01, tol_axis=2.0):
    """The edge generator handles lower single null equilibria only: inside the wall contour exactly one
    O-point (the magnetic axis of the g-file) and exactly one X-point, which lies on the plasma boundary
    (psi_N = 1 with respect to the g-file's boundary flux sibry) below the axis.
    Returns the X-point (R, Z) [cm]; stops with a description of the topology otherwise."""
    pts = eq.critical_points()
    opts = [p for p in pts if p[2] == "O"]
    xpts = [p for p in pts if p[2] == "X"]
    psin_b = lambda p: (float(eq.spl.ev(p[0], p[1])) - eq.psi_ax) / (eq.g["sibry"] - eq.psi_ax)
    desc = "; ".join("%s-point at R=%.1f Z=%.1f cm (psi_N=%.4f)" % (p[2], p[0], p[1], psin_b(p)) for p in pts)
    head = ("Edge mesh: unsupported magnetic topology (only lower single null is implemented).\n"
            "  Found inside the wall contour: %s.\n  " % (desc or "no critical points"))
    near_axis = [p for p in opts if np.hypot(p[0] - eq.raxis, p[1] - eq.zaxis) < tol_axis]
    if len(near_axis) != 1 or len(opts) != 1:
        raise TopologyError(head + "Expected exactly one O-point, at the magnetic axis of the g-file "
                            "(R=%.1f Z=%.1f cm): additional O-points mean islands or a second axis."
                            % (eq.raxis, eq.zaxis))
    on_boundary = [p for p in xpts if abs(psin_b(p) - 1.0) < tol_psin]
    if not on_boundary:
        raise TopologyError(head + "No X-point on the plasma boundary (psi_N = 1) inside the wall: a limited "
                            "plasma, or the X-point lies outside the wall contour.")
    x = on_boundary[0]
    if x[1] > eq.zaxis:
        raise TopologyError(head + "The X-point on the plasma boundary lies above the magnetic axis: upper "
                            "single null.")
    if len(xpts) > 1:
        raise TopologyError(head + "More than one X-point inside the wall: double null, or a second "
                            "separatrix crossing the scrape-off layer.")
    return np.array(x[:2])

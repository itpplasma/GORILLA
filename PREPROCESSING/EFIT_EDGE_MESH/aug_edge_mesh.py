#!/usr/bin/env python
"""Extend a GORILLA grid_kind=2 (flux-aligned, closed surfaces) EFIT mesh out to the
g-file limiter.

Inside the last closed ring the grid_kind=2 vertices and triangles are kept unchanged.
Between that ring and the limiter, vertices are placed on psi = const contours (SOL,
divertor legs and private-flux region) and connected by a constrained Delaunay
triangulation (Shewchuk's Triangle) with the contours as constraint segments, so
triangle edges follow flux surfaces wherever a contour exists.

Output: knots/triangles files in the SOLEDGE3X_EIRENE format read by GORILLA
(grid_kind=4): knots "N 2" + R Z [cm]; triangles "M 3" + 1-based vertex indices.

Usage: aug_edge_mesh.py DUMPDIR GFILE OUTDIR [options]
DUMPDIR holds verts_2d.dat / tets_slice1.dat from tools/dump_mesh_2d.x.
"""
import argparse
import os
import sys

import numpy as np
import triangle
import contourpy
from matplotlib.path import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gfile import Equilibrium  # noqa: E402


# ----------------------------------------------------------------------------- inner mesh
def load_inner(dumpdir):
    with open(os.path.join(dumpdir, "verts_2d.dat")) as f:
        nvs = int(f.readline().split()[0])
    v = np.loadtxt(os.path.join(dumpdir, "verts_2d.dat"), skiprows=1)
    assert len(v) == nvs
    tets = np.loadtxt(os.path.join(dumpdir, "tets_slice1.dat"), dtype=int)
    # The prisms of grid_kind 2 are sheared toroidally, so projecting a tetrahedron onto the
    # poloidal plane is not a triangle of the mesh. Use instead the tetrahedron faces whose
    # three knots all lie in the first slice (phi = 0): these are the prism bottoms.
    tris = set()
    for row in tets:
        for drop in range(4):
            f = [row[k] for k in range(4) if k != drop]
            if max(f) <= nvs:
                tris.add(tuple(sorted(x - 1 for x in f)))
    tris = np.array(sorted(tris))
    RZ = v[:, 1:3]
    s = v[:, 3]
    load_inner.theta = v[:, 4]
    return RZ, s, tris


def boundary_loops(tris):
    """Return the boundary edges of a triangle set grouped into closed loops."""
    from collections import Counter, defaultdict
    cnt = Counter()
    for t in tris:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            cnt[tuple(sorted((a, b)))] += 1
    bad = [e for e, c in cnt.items() if c > 2]
    if bad:
        raise RuntimeError("non-manifold edges: %s" % bad[:5])
    bedges = [e for e, c in cnt.items() if c == 1]
    nb = defaultdict(list)
    for a, b in bedges:
        nb[a].append(b)
        nb[b].append(a)
    loops, seen = [], set()
    for start in nb:
        if start in seen:
            continue
        loop = [start]
        seen.add(start)
        prev, cur = None, start
        while True:
            nxt = [x for x in nb[cur] if x != prev]
            nxt = nxt[0] if prev is not None else nb[cur][0]
            if nxt == start:
                break
            loop.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        loops.append(loop)
    return loops


# ----------------------------------------------------------------------------- geometry
def seg_intersect(p1, p2, q1, q2):
    """Intersection parameter (t along p, u along q) or None."""
    d1, d2 = p2 - p1, q2 - q1
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-14:
        return None
    w = q1 - p1
    t = (w[0] * d2[1] - w[1] * d2[0]) / den
    u = (w[0] * d1[1] - w[1] * d1[0]) / den
    if -1e-12 <= t <= 1 + 1e-12 and -1e-12 <= u <= 1 + 1e-12:
        return t, u
    return None


def clip_polyline(line, poly, path):
    """Split a polyline into the runs that lie inside the closed polygon poly.
    Returns a list of (points, start_hit, end_hit); hit = (edge index, u) on poly or None."""
    inside = path.contains_points(line)
    runs, cur, start_hit = [], None, None
    npoly = len(poly) - 1  # poly is closed (last == first)

    def crossing(a, b):
        best = None
        for k in range(npoly):
            r = seg_intersect(a, b, poly[k], poly[k + 1])
            if r is not None and (best is None or r[0] < best[0]):
                best = (r[0], k, r[1])
        return best

    for i in range(len(line)):
        if inside[i]:
            if cur is None:
                cur = []
                start_hit = None
                if i > 0:
                    c = crossing(line[i - 1], line[i])
                    if c is not None:
                        t, k, u = c
                        cur.append(line[i - 1] + t * (line[i] - line[i - 1]))
                        start_hit = (k, u)
            cur.append(line[i])
        else:
            if cur is not None:
                c = crossing(line[i - 1], line[i])
                end_hit = None
                if c is not None:
                    t, k, u = c
                    cur.append(line[i - 1] + t * (line[i] - line[i - 1]))
                    end_hit = (k, u)
                runs.append((np.array(cur), start_hit, end_hit))
                cur = None
    if cur is not None:
        runs.append((np.array(cur), start_hit, None))
    return runs


def resample(line, hfun):
    """Resample a polyline with local spacing hfun(point); endpoints kept."""
    seg = np.linalg.norm(np.diff(line, axis=0), axis=1)
    sarc = np.concatenate([[0], np.cumsum(seg)])
    L = sarc[-1]
    if L == 0:
        return line[:1]
    # integrate 1/h along the line to place points with local spacing h
    fine = np.linspace(0, L, max(int(L / 0.05), 50))
    pts = np.column_stack([np.interp(fine, sarc, line[:, 0]), np.interp(fine, sarc, line[:, 1])])
    dens = 1.0 / hfun(pts)
    N = np.concatenate([[0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(fine))])
    n = max(int(round(N[-1])), 1)
    targets = np.linspace(0, N[-1], n + 1)
    st = np.interp(targets, N, fine)
    return np.column_stack([np.interp(st, sarc, line[:, 0]), np.interp(st, sarc, line[:, 1])])


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dumpdir")
    ap.add_argument("gfile")
    ap.add_argument("outdir")
    ap.add_argument("--d1", type=float, default=None,
                    help="first SOL contour distance from the last ring at the outboard midplane [cm] "
                         "(default: radial spacing of the last two inner rings there)")
    ap.add_argument("--growth", type=float, default=1.25, help="geometric growth of contour spacing")
    ap.add_argument("--dmax", type=float, default=4.0, help="max contour spacing [cm]")
    ap.add_argument("--hmin", type=float, default=0.8, help="vertex spacing along contours at the X-point [cm]")
    ap.add_argument("--hmax", type=float, default=4.0, help="max vertex spacing along contours [cm]")
    ap.add_argument("--hgrad", type=float, default=0.15, help="growth of spacing with distance from X-point")
    ap.add_argument("--minangle", type=float, default=0.0, help="Triangle quality: min angle [deg] (0 = off)")
    ap.add_argument("--leg-aspect", type=float, default=0.5,
                    help="min contour spacing on the leg cuts in units of the local vertex spacing h")
    ap.add_argument("--leg-cuts", type=float, nargs="+", default=[0.25, 0.5, 0.75, 0.95],
                    help="positions of the cuts along each leg (fraction of the leg length from the X-point)")
    ap.add_argument("--no-legs", action="store_true", help="do not add the separatrix legs as constraints")
    ap.add_argument("--core", choices=["native", "cdt"], default="native",
                    help="native: keep grid_kind 2 triangles inside the last ring; cdt: re-triangulate")
    ap.add_argument("--xref-radius", type=float, default=0.0,
                    help="refine the edge triangles within this distance of the X-point [cm] (0 = off)")
    ap.add_argument("--xref-h", type=float, default=0.5, help="target triangle size at the X-point [cm]")
    ap.add_argument("--fill-axis", action="store_true",
                    help="close the hole inside the innermost grid_kind 2 ring (axis point + triangles)")
    ap.add_argument("--gap-min", type=float, default=0.07,
                    help="every contour keeps a gap >= gap_min*h to its neighbour level everywhere "
                         "(0.07 ~ 8 deg minimum angle between contours)")
    ap.add_argument("--snap", type=float, default=0.3,
                    help="contour end points closer than snap*h to a limiter corner are moved onto it")
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    eq = Equilibrium(a.gfile)
    RZ_in, s_in, tri_in = load_inner(a.dumpdir)
    theta_in = load_inner.theta
    nin = len(RZ_in)

    # ---- outer ring of the inner mesh (= the one with the largest s)
    loops = [None, None]
    ring = np.where(np.isclose(s_in, s_in.max()))[0]
    ringRZ = RZ_in[ring]
    # orient counter-clockwise
    area = 0.5 * np.sum(ringRZ[:, 0] * np.roll(ringRZ[:, 1], -1) - np.roll(ringRZ[:, 0], -1) * ringRZ[:, 1])
    if area < 0:
        ring, ringRZ = ring[::-1], ringRZ[::-1]
    psin_ring = eq.psin(ringRZ[:, 0], ringRZ[:, 1])
    print("inner mesh: %d verts, %d triangles, %d boundary loops; last ring %d pts, psi_N %.5f..%.5f"
          % (nin, len(tri_in), len(loops), len(ring), psin_ring.min(), psin_ring.max()))

    ring_path = Path(np.vstack([ringRZ, ringRZ[:1]]))
    lim = eq.lim.copy()
    if np.linalg.norm(lim[0] - lim[-1]) > 1e-8:
        lim = np.vstack([lim, lim[:1]])
    # remove duplicate consecutive limiter points
    keep = np.concatenate([[True], np.linalg.norm(np.diff(lim, axis=0), axis=1) > 1e-6])
    lim = lim[keep]
    lim_path = Path(lim)
    assert lim_path.contains_points(ringRZ).all(), "last ring not inside limiter"

    xpt = eq.xpt
    print("X-point R,Z = %.3f %.3f cm" % tuple(xpt))

    def hfun(p):
        d = np.linalg.norm(np.atleast_2d(p) - xpt, axis=1)
        return np.clip(a.hmin + a.hgrad * d, a.hmin, a.hmax)

    # ---- contour lines on a fine grid (needed first: the separatrix legs define the cuts)
    Rg = np.linspace(lim[:, 0].min() - 2, lim[:, 0].max() + 2, 900)
    Zg = np.linspace(lim[:, 1].min() - 2, lim[:, 1].max() + 2, 1600)
    RR, ZZ = np.meshgrid(Rg, Zg)
    PN = eq.psin(RR, ZZ)
    gen = contourpy.contour_generator(Rg, Zg, PN, line_type="Separate")

    # separatrix legs: psi_N = 1 below the X-point, away from the X-point, reaching the limiter
    legs = []
    for line in gen.lines(1.0):
        if len(line) < 3:
            continue
        dx = np.linalg.norm(line - xpt, axis=1)
        below = (line[:, 1] < xpt[1] - 0.5 * a.hmin) & (dx > 3 * a.hmin)
        idx = np.where(below)[0]
        if len(idx) == 0:
            continue
        breaks = np.where(np.diff(idx) > 1)[0]
        for part in np.split(idx, breaks + 1):
            sub = line[part]
            if len(sub) < 3:
                continue
            for pts, sh, eh in clip_polyline(sub, lim, lim_path):
                if sh is None and eh is None:
                    continue
                legs.append((pts, sh, eh))

    # ---- contour levels.
    # Each level is placed so that its distance to the previous level is at least
    #   * d_k (d1 * growth^k, capped at dmax) on the outboard midplane (SOL; measured from the
    #     last core ring) and on a vertical line below the X-point (PFR), and
    #   * max(d_k, leg_aspect * h) on cuts normal to each separatrix leg (SOL and PFR side),
    #     h = local vertex spacing along the contours.
    # The leg cuts stop the contours from bunching along the legs (thin triangles).
    zmid = eq.zaxis
    rr = np.linspace(eq.raxis, eq.lim[:, 0].max() + 30, 4000)
    pn_mid = eq.psin(rr, np.full_like(rr, zmid))
    r_sep = np.interp(1.0, pn_mid, rr)
    r_ring = np.interp(psin_ring.mean(), pn_mid, rr)
    # radial spacing of the last two inner rings at the outboard midplane
    sv = np.unique(np.round(s_in, 12))
    s_prev = sv[-2]
    prev_pts = RZ_in[np.isclose(s_in, s_prev)]
    psin_prev = eq.psin(prev_pts[:, 0], prev_pts[:, 1]).mean()
    d_last = r_ring - np.interp(psin_prev, pn_mid, rr)
    d1 = a.d1 if a.d1 else d_last
    print("outboard midplane: R_ring=%.3f R_sep=%.3f; last inner radial spacing %.3f cm; d1=%.3f" %
          (r_ring, r_sep, d_last, d1))

    def monotone_cut(p0, direc, length, sign):
        """psi_N along p0 + u*direc, u in [0, length], kept while inside the limiter and
        monotone (increasing for sign=+1, decreasing for -1). Returns (u, psin) or None."""
        u = np.linspace(0, length, 2000)
        P = p0[None, :] + u[:, None] * direc[None, :]
        pn = eq.psin(P[:, 0], P[:, 1])
        ins = lim_path.contains_points(P)
        ok = ins & np.concatenate([[True], sign * np.diff(pn) > 0])
        n = len(ok) if ok.all() else int(np.argmin(ok))
        if n < 5:
            return None
        return u[:n], pn[:n]

    cuts = {"sol": [], "pfr": []}  # (u, psin, required spacing as function of d_k, name)
    um = rr - r_ring
    msk = um >= 0
    cuts["sol"].append((um[msk], pn_mid[msk], lambda d: d, "outboard midplane"))
    c = monotone_cut(np.asarray(xpt, float), np.array([0.0, -1.0]), xpt[1] - eq.lim[:, 1].min(), -1)
    if c is not None:  # vertical line below the X-point, inside the limiter
        cuts["pfr"].append((c[0], c[1], lambda d: d, "below X-point"))
    for il, (pts, sh, eh) in enumerate(legs):
        sarc = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        # parametrise from the X-point end
        if np.linalg.norm(pts[0] - xpt) > np.linalg.norm(pts[-1] - xpt):
            pts, sarc = pts[::-1], sarc[-1] - sarc[::-1]
        for f in a.leg_cuts:
            sf = f * sarc[-1]
            p = np.array([np.interp(sf, sarc, pts[:, 0]), np.interp(sf, sarc, pts[:, 1])])
            j = min(max(np.searchsorted(sarc, sf), 1), len(pts) - 1)
            t = pts[j] - pts[j - 1]
            t = t / np.linalg.norm(t)
            nrm = np.array([-t[1], t[0]])
            hloc = float(hfun(p)[0])
            for side, sgn in (("sol", +1), ("pfr", -1)):
                # normal direction in which psi_N increases (sol) / decreases (pfr)
                test = float(eq.psin(*(p + 0.2 * nrm)))
                direc = nrm if (test - 1.0) * sgn > 0 else -nrm
                c = monotone_cut(p, direc, 40.0, sgn)
                if c is not None:
                    cuts[side].append((c[0], c[1], (lambda hh: (lambda d: max(d, a.leg_aspect * hh)))(hloc),
                                       "leg %d at %.0f%%" % (il, 100 * f)))

    def lim_dist(P):
        """distance of points P (n,2) to the limiter polygon"""
        a0, ab = lim[:-1], lim[1:] - lim[:-1]
        best = np.full(len(P), np.inf)
        for k in range(len(a0)):
            t = np.clip(((P - a0[k]) @ ab[k]) / (ab[k] @ ab[k]), 0, 1)
            best = np.minimum(best, np.linalg.norm(a0[k] + t[:, None] * ab[k] - P, axis=1))
        return best

    def level_points(lv):
        """points of the psi_N = lv contour inside the limiter, outside the core ring, and
        more than h away from the limiter (contours are trimmed there anyway)"""
        L = [ln for ln in gen.lines(lv) if len(ln) > 1]
        if not L:
            return np.zeros((0, 2))
        P = np.vstack(L)
        P = P[lim_path.contains_points(P) & ~ring_path.contains_points(P)]
        if len(P):
            P = P[lim_dist(P) > hfun(P)]
        return P

    from scipy.spatial import cKDTree

    def spacing_ok(lv, prev):
        """min over the contour of dist(contour lv, contour prev) / (gap_min * h) >= 1 ?"""
        A, B = level_points(lv), level_points(prev)
        if len(A) == 0 or len(B) == 0:
            return True
        dist, _ = cKDTree(B).query(A)
        return np.min(dist / (a.gap_min * hfun(A))) >= 1.0

    def choose_levels(side, start_level, stop):
        sgn = 1 if side == "sol" else -1
        levels, d, prev = [], d1, start_level
        while True:
            cand = []
            for u, pn, req, name in cuts[side]:
                if (prev - pn[-1]) * sgn >= 0:
                    continue  # previous level already beyond the end of this cut
                # previous level before the start of the cut (e.g. the core ring): measure from the start
                upos = 0.0 if (prev - pn[0]) * sgn <= 0 else np.interp(prev * sgn, pn * sgn, u)
                target = upos + req(d)
                if target > u[-1]:
                    continue
                cand.append(np.interp(target, u, pn))
            if not cand:
                break
            lv = float(max(cand) if sgn > 0 else min(cand))
            if os.environ.get("MESHDBG"): print(side, "prev %.5f d %.3f" % (prev, d), "cand", np.round(cand, 5))
            # widen the step until the new contour keeps its distance to the previous one
            # everywhere (not only on the cuts); the first level is checked against the separatrix
            ref = prev if (prev - 1.0) * sgn > 0 else 1.0
            for _ in range(40):
                if (lv - stop) * sgn > 0 or spacing_ok(lv, ref):
                    break
                lv = ref + (lv - ref) * 1.15
            if (lv - stop) * sgn > 0:
                break
            if abs(lv - 1.0) > 1e-4:
                levels.append(lv)
            prev = lv
            d = min(d * a.growth, a.dmax)
        return levels

    pmax = max(eq.psin(lim[:, 0], lim[:, 1]).max(), 1.0) + 1e-3
    pmin = eq.psin(lim[:, 0], lim[:, 1]).min()
    sol_levels = choose_levels("sol", float(psin_ring.mean()), pmax)
    pfr_levels = choose_levels("pfr", 1.0, pmin)
    print("cuts: SOL %d, PFR %d (leg cuts at %s of each of %d legs)" %
          (len(cuts["sol"]), len(cuts["pfr"]), a.leg_cuts, len(legs)))
    print("SOL levels (%d):" % len(sol_levels), np.round(sol_levels, 4))
    print("PFR levels (%d):" % len(pfr_levels), np.round(pfr_levels, 4))

    constraint_lines = []  # (points, start_hit, end_hit, label)

    def add_level(lv, label, keep=None):
        for line in gen.lines(lv):
            if len(line) < 3:
                continue
            # drop lines (or parts) inside the last ring (closed core contours)
            if ring_path.contains_points(line).any():
                continue
            # drop closed contours around the core (between the last ring and the separatrix)
            closed = np.linalg.norm(line[0] - line[-1]) < 1e-6
            if closed and Path(line).contains_point((eq.raxis, eq.zaxis)):
                continue
            for pts, sh, eh in clip_polyline(line, lim, lim_path):
                if keep is not None:
                    m = keep(pts)
                    if not m.all():
                        continue
                L = np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1))
                if L < 2 * a.hmin:
                    continue
                constraint_lines.append((pts, sh, eh, label))

    for lv in sol_levels:
        add_level(lv, "sol")
    for lv in pfr_levels:
        add_level(lv, "pfr")
    if not a.no_legs:
        for pts, sh, eh in legs:
            constraint_lines.append((pts, sh, eh, "leg"))

    # ---- assemble the PSLG
    if a.core == "native":
        # keep the grid_kind 2 triangles; the outer CDT starts at the last ring (no Steiner
        # points on it, so the two parts conform). Needs a fold-free core (theta_geom_flux=2).
        from fold_check import folds
        nf = len(folds(RZ_in, tri_in))
        if nf:
            raise SystemExit("grid_kind 2 core has %d folded edges in (R,Z); use theta_geom_flux = 2 "
                             "or --core cdt" % nf)
        verts = [p for p in ringRZ]
        segs = [(i, (i + 1) % len(ringRZ)) for i in range(len(ringRZ))]
        seg_marks = [2] * len(segs)
        nfix = len(ringRZ)
    else:
        # re-triangulate the core too: all grid_kind 2 vertices kept, every ring a constraint
        verts = [p for p in RZ_in]
        segs, seg_marks = [], []
        for sv_ in sv:
            idx = np.where(np.isclose(s_in, sv_))[0]
            for k in range(len(idx)):
                segs.append((idx[k], idx[(k + 1) % len(idx)]))
                seg_marks.append(2 if sv_ == sv[0] else 5)
        nfix = nin
    lim_extra = {k: [] for k in range(len(lim) - 1)}  # edge k -> list of (u, vertex index)

    def add_vertex(p):
        verts.append(np.asarray(p, float))
        return len(verts) - 1

    clearance = 0.35  # min distance (in units of local h) of interior contour points to the limiter

    def dist_to_lim(p):
        a0, b0 = lim[:-1], lim[1:]
        ab = b0 - a0
        t = np.clip(np.einsum("ij,ij->i", p - a0, ab) / np.einsum("ij,ij->i", ab, ab), 0, 1)
        return np.min(np.linalg.norm(a0 + t[:, None] * ab - p, axis=1))

    n_lines = 0
    for pts, sh, eh, label in constraint_lines:
        rs = resample(pts, hfun)
        # drop interior points too close to the limiter (except the end points that sit on it)
        keepm = np.ones(len(rs), bool)
        for j in range(1, len(rs) - 1):
            if dist_to_lim(rs[j:j + 1]) < clearance * hfun(rs[j])[0]:
                keepm[j] = False
        # also drop the second / second-to-last point if too close to an end on the limiter
        # vertices for kept points; segments only between neighbours that are both kept, so
        # a contour grazing the limiter is split instead of bridged by a long chord
        ids = [add_vertex(p) if keepm[j] else None for j, p in enumerate(rs)]
        if sh is not None:
            lim_extra[sh[0]].append((sh[1], ids[0]))
        if eh is not None:
            lim_extra[eh[0]].append((eh[1], ids[-1]))
        for j in range(len(ids) - 1):
            if ids[j] is not None and ids[j + 1] is not None:
                segs.append((ids[j], ids[j + 1]))
                seg_marks.append(3 if label == "leg" else 4)
        n_lines += 1
    print("constraint lines used: %d" % n_lines)

    # contour end points that land closer than snap*h to a limiter corner are moved onto the
    # corner (otherwise corner + end point form a needle triangle)
    nedge = len(lim) - 1
    snapped = set()
    for k in range(nedge):
        L = np.linalg.norm(lim[k + 1] - lim[k])
        moved = []
        for u, vid in list(lim_extra[k]):
            hh = float(hfun(verts[vid])[0])
            if u * L < a.snap * hh and k not in snapped:
                verts[vid] = lim[k].copy()
                snapped.add(k)
                moved.append(((u, vid), k, 0.0))
            elif (1 - u) * L < a.snap * hh and (k + 1) % nedge not in snapped:
                verts[vid] = lim[k + 1].copy()
                snapped.add((k + 1) % nedge)
                moved.append(((u, vid), (k + 1) % nedge, 0.0))
        for item, knew, unew in moved:
            lim_extra[k].remove(item)
            lim_extra[knew].append((unew, item[1]))
    if snapped:
        print("contour end points snapped to limiter corners: %d" % len(snapped))

    # limiter polygon, resampled, with contour end points inserted
    lim_ids = []
    for k in range(len(lim) - 1):
        pa, pb = lim[k], lim[k + 1]
        ex = sorted(lim_extra[k])
        stops = ([] if ex and ex[0][0] == 0.0 else [(0.0, None)]) + ex + [(1.0, None)]
        for (u0, i0), (u1, i1) in zip(stops[:-1], stops[1:]):
            # start vertex
            if i0 is None:
                if u0 == 0.0:
                    i0 = add_vertex(pa)
            if not lim_ids or lim_ids[-1] != i0:
                lim_ids.append(i0)
            q0 = pa + u0 * (pb - pa)
            q1 = pa + u1 * (pb - pa)
            sub = resample(np.vstack([q0, q1]), hfun)
            for p in sub[1:-1]:
                lim_ids.append(add_vertex(p))
            if i1 is not None:
                lim_ids.append(i1)
    # close: last stop of last edge is lim[-1] == lim[0], already the first vertex
    if np.linalg.norm(verts[lim_ids[-1]] - verts[lim_ids[0]]) < 1e-9 and lim_ids[-1] != lim_ids[0]:
        lim_ids.pop()
    # remove consecutive duplicates
    clean = [lim_ids[0]]
    for i in lim_ids[1:]:
        if i != clean[-1] and np.linalg.norm(verts[i] - verts[clean[-1]]) > 1e-9:
            clean.append(i)
    if clean[-1] == clean[0]:
        clean.pop()
    lim_ids = clean
    for j in range(len(lim_ids)):
        segs.append((lim_ids[j], lim_ids[(j + 1) % len(lim_ids)]))
        seg_marks.append(1)  # 1 = limiter

    V = np.array(verts)
    S = np.array(segs)
    seg_marks = np.array(seg_marks)
    # merge (nearly) duplicate vertices, keeping the first occurrence (ring vertices come first)
    tol = 1e-3
    remap = np.arange(len(V))
    order = np.lexsort((V[:, 1], V[:, 0]))
    for ii in range(len(order)):
        i = order[ii]
        jj = ii + 1
        while jj < len(order) and V[order[jj], 0] - V[i, 0] < tol:
            j = order[jj]
            if np.linalg.norm(V[j] - V[i]) < tol:
                lo, hi = min(remap[i], remap[j]), max(i, j)
                remap[hi] = lo
            jj += 1
    ndup = int((remap != np.arange(len(V))).sum())
    if ndup:
        print("merged %d duplicate vertices" % ndup)
    used = np.unique(remap)
    newidx = -np.ones(len(V), int)
    newidx[used] = np.arange(len(used))
    S = newidx[remap[S]]
    V = V[used]
    keepseg = S[:, 0] != S[:, 1]
    S, seg_marks = S[keepseg], seg_marks[keepseg]
    S_sorted = np.sort(S, axis=1)
    _, first = np.unique(S_sorted, axis=0, return_index=True)
    first = np.sort(first)
    S, seg_marks = S[first], seg_marks[first]
    nbad = count_crossings(V, S)
    print("PSLG: %d vertices, %d segments, %d segment crossings" % (len(V), len(S), nbad))
    if nbad:
        raise SystemExit("PSLG has crossing segments; aborting before Triangle")

    if a.xref_radius > 0:
        # X-point refinement: free (non-contour) points on circles around the X-point with
        # spacing h_x(r) = xref_h + hgrad * r, kept only where they are at least 0.5 h_x away
        # from every existing vertex and constraint segment (no slivers against contours),
        # inside the limiter and outside the core ring.
        A0, AB = V[S[:, 0]], V[S[:, 1]] - V[S[:, 0]]
        L2 = np.maximum(np.einsum("ij,ij->i", AB, AB), 1e-30)
        tree = cKDTree(V)
        new_pts = [] if ring_path.contains_point(xpt) else [xpt.copy()]
        r = a.xref_h
        while r < a.xref_radius:
            hx = a.xref_h + a.hgrad * r
            n = max(int(np.ceil(2 * np.pi * r / hx)), 6)
            ph = 2 * np.pi * (np.arange(n) + 0.5 * (len(new_pts) % 2)) / n
            for q in xpt + r * np.column_stack([np.cos(ph), np.sin(ph)]):
                if not lim_path.contains_point(q) or ring_path.contains_point(q):
                    continue
                if tree.query(q)[0] < 0.5 * hx:
                    continue
                t = np.clip(((q - A0) * AB).sum(1) / L2, 0, 1)
                if np.min(np.linalg.norm(A0 + t[:, None] * AB - q, axis=1)) < 0.5 * hx:
                    continue
                if new_pts and np.min(np.linalg.norm(np.array(new_pts) - q, axis=1)) < 0.7 * hx:
                    continue
                new_pts.append(q)
            r += hx
        if new_pts:
            V = np.vstack([V, new_pts])
        print("X-point refinement (r < %.1f cm, h = %.2f cm at the X-point): %d free points added"
              % (a.xref_radius, a.xref_h, len(new_pts)))

    pslg = dict(vertices=V, segments=S, segment_markers=seg_marks[:, None],
                holes=np.array([[eq.raxis, eq.zaxis]]))
    opts = "pY" + ("q%g" % a.minangle if a.minangle > 0 else "")
    out = triangle.triangulate(pslg, opts)
    Vo, To = out["vertices"], out["triangles"]
    fixed = ringRZ if a.core == "native" else RZ_in
    assert np.allclose(Vo[:nfix], fixed), "Triangle reordered the fixed input vertices"
    print("outer/CDT triangulation: %d verts (%d input), %d triangles" % (len(Vo), len(V), len(To)))
    if a.core == "native":
        mapping = np.empty(len(Vo), int)
        mapping[:nfix] = ring
        mapping[nfix:] = nin + np.arange(len(Vo) - nfix)
        RZ_all = np.vstack([RZ_in, Vo[nfix:]])
        T_all = np.vstack([tri_in, mapping[To]])
        region = np.concatenate([np.zeros(len(tri_in), int), np.ones(len(To), int)])
    else:
        RZ_all, T_all = Vo, To
        region = np.where(ring_path.contains_points(RZ_all[T_all].mean(axis=1)), 0, 1)
    hole_RZ = RZ_in[np.isclose(s_in, s_in.min())]
    if a.fill_axis:
        # close the hole around the magnetic axis (only meaningful in cylindrical coordinates):
        # innermost ring + axis point, triangulated with interior Steiner points only (Y)
        ids = np.where(np.isclose(s_in, s_in.min()))[0]
        P = RZ_in[ids]
        pa = 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])
        if pa < 0:
            ids, P = ids[::-1], P[::-1]
        n = len(P)
        hv = dict(vertices=np.vstack([P, [[eq.raxis, eq.zaxis]]]),
                  segments=np.array([(k, (k + 1) % n) for k in range(n)]))
        ho = triangle.triangulate(hv, "pYq25")
        assert np.allclose(ho["vertices"][:n], P)
        m = np.empty(len(ho["vertices"]), int)
        m[:n] = ids
        m[n:] = len(RZ_all) + np.arange(len(ho["vertices"]) - n)
        RZ_all = np.vstack([RZ_all, ho["vertices"][n:]])
        T_all = np.vstack([T_all, m[ho["triangles"]]])
        region = np.concatenate([region, np.full(len(ho["triangles"]), 2)])
        hole_RZ = None
        print("axis hole filled: %d triangles, %d new vertices (incl. the axis)"
              % (len(ho["triangles"]), len(ho["vertices"]) - n))
    nfold = len(__import__("fold_check").folds(RZ_all, T_all))
    print("folded edges in final mesh: %d" % nfold)
    if nfold:
        raise SystemExit("final mesh folds")
    report = check_mesh(RZ_all, T_all, eq, hole_RZ)
    # orient all triangles counter-clockwise
    T_all = orient_ccw(RZ_all, T_all)

    knots = os.path.join(a.outdir, "knots_aug_edge.dat")
    tris = os.path.join(a.outdir, "triangles_aug_edge.dat")
    with open(knots, "w") as f:
        f.write("%d 2\n" % len(RZ_all))
        for p in RZ_all:
            f.write("%.15e %.15e\n" % (p[0], p[1]))
    with open(tris, "w") as f:
        f.write("%d 3\n" % len(T_all))
        for t in T_all + 1:
            f.write("%d %d %d\n" % tuple(t))
    np.savetxt(os.path.join(a.outdir, "triangle_region.dat"), region, fmt="%d",
               header="0 = grid_kind 2 inner mesh, 1 = edge extension, 2 = axis fill")
    # flux coordinates of the core vertices (the first nin vertices are the grid_kind 2 vertices in
    # dump order) and the region of every triangle, for GORILLA (boole_extend_to_wall)
    st = -np.ones((len(RZ_all), 2))
    st[:nin, 0] = s_in
    st[:nin, 1] = theta_in
    with open(os.path.join(a.outdir, "edge_mesh_core.dat"), "w") as f:
        f.write("%d %d\n" % (len(RZ_all), len(T_all)))
        for v in st:
            f.write("%.16e %.16e\n" % (v[0], v[1]))
        for r in region:
            f.write("%d\n" % r)
    # settings record: core settings from the dump + the generator options that change the mesh
    rec = os.path.join(a.dumpdir, "core_settings.nml")
    if os.path.exists(rec):
        opts = " ".join(x for x in sys.argv[4:] if x != "--plot")
        lines = open(rec).read().rstrip().splitlines()
        assert lines[-1].strip() == "/"
        lines.insert(-1, " generator_options = '%s' ," % opts.replace("'", ""))
        with open(os.path.join(a.outdir, "edge_mesh.nml"), "w") as f:
            f.write("\n".join(lines) + "\n")
    else:
        print("note: no core_settings.nml in the dump dir (old dump tool?); edge_mesh.nml not written")
    with open(os.path.join(a.outdir, "mesh_report.txt"), "w") as f:
        f.write(" ".join(sys.argv) + "\n")
        for k, v in report.items():
            f.write("%s: %s\n" % (k, v))
    for k, v in report.items():
        print("  %s: %s" % (k, v))
    print("wrote", knots, tris)

    if a.plot:
        plot_mesh(RZ_all, T_all, region, eq, a.outdir, constraint_lines)


def count_crossings(V, S):
    """Number of pairs of segments that cross in their interiors (shared end points allowed)."""
    A, B = V[S[:, 0]], V[S[:, 1]]
    n = 0
    for k in range(len(S)):
        a, b = A[k], B[k]
        d1 = b - a
        d2 = B - A
        den = d1[0] * d2[:, 1] - d1[1] * d2[:, 0]
        w = A - a
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (w[:, 0] * d2[:, 1] - w[:, 1] * d2[:, 0]) / den
            u = (w[:, 0] * d1[1] - w[:, 1] * d1[0]) / den
        share = (S[:, 0] == S[k, 0]) | (S[:, 0] == S[k, 1]) | (S[:, 1] == S[k, 0]) | (S[:, 1] == S[k, 1])
        hit = (abs(den) > 1e-14) & (t > 1e-9) & (t < 1 - 1e-9) & (u > 1e-9) & (u < 1 - 1e-9) & ~share
        hit[: k + 1] = False
        if hit.any():
            for j in np.where(hit)[0][:3]:
                print("  crossing: seg", k, S[k], a, b, "with seg", j, S[j], A[j], B[j])
        n += int(hit.sum())
    return n


def orient_ccw(P, T):
    a = P[T[:, 0]]
    b = P[T[:, 1]]
    c = P[T[:, 2]]
    cr = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    T = T.copy()
    neg = cr < 0
    T[neg] = T[neg][:, [0, 2, 1]]
    return T


def check_mesh(P, T, eq, hole_RZ):
    from collections import Counter
    rep = {}
    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))
    rep["n_vertices"] = len(P)
    rep["n_triangles"] = len(T)
    rep["min_area_cm2"] = "%.3e" % area.min()
    # angles
    def ang(p, q, r):
        u, v = q - p, r - p
        return np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", u, v) /
                                            (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)), -1, 1)))
    angs = np.column_stack([ang(a, b, c), ang(b, c, a), ang(c, a, b)])
    rep["min_angle_deg"] = "%.2f" % angs.min()
    rep["max_angle_deg"] = "%.2f" % angs.max()
    cnt = Counter()
    for t in T:
        for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            cnt[tuple(sorted(e))] += 1
    rep["edges_shared_by_more_than_2"] = sum(1 for v in cnt.values() if v > 2)
    rep["boundary_edges"] = sum(1 for v in cnt.values() if v == 1)
    used = np.zeros(len(P), bool)
    used[T.ravel()] = True
    rep["unused_vertices"] = int((~used).sum())
    # Euler characteristic: V - E + F = 2 - 2g - (holes) ... for planar with h holes: V-E+F = 1 - h
    rep["euler_V-E+F"] = len(P) - len(cnt) + len(T)
    # boundary loops
    loops = boundary_loops(T)
    rep["boundary_loops"] = len(loops)
    rep["sum_area_cm2"] = "%.3f" % area.sum()
    # no overlaps: total area = area(limiter) - area(inner hole)
    lim = eq.lim
    A_lim = 0.5 * abs(np.sum(lim[:, 0] * np.roll(lim[:, 1], -1) - np.roll(lim[:, 0], -1) * lim[:, 1]))
    H = hole_RZ
    A_hole = 0.0 if H is None else 0.5 * abs(np.sum(H[:, 0] * np.roll(H[:, 1], -1) - np.roll(H[:, 0], -1) * H[:, 1]))
    rep["area_limiter_minus_hole_cm2"] = "%.3f" % (A_lim - A_hole)
    return rep


def plot_mesh(P, T, region, eq, outdir, lines):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    for name, box in (("full", None), ("xpoint", (120, 165, -125, -80)), ("outboard", (195, 225, -20, 40))):
        fig, ax = plt.subplots(figsize=(7, 10) if box is None else (8, 8))
        for reg, col in ((0, "tab:blue"), (1, "tab:red"), (2, "tab:orange")):
            m = region == reg
            pc = PolyCollection(P[T[m]], facecolors="none", edgecolors=col, linewidths=0.25)
            ax.add_collection(pc)
        ax.plot(*eq.lim.T, "k-", lw=1)
        for pts, sh, eh, label in lines:
            ax.plot(*pts.T, "-", color={"leg": "k", "sol": "tab:green", "pfr": "tab:purple"}.get(label, "y"),
                    lw=0.8 if label == "leg" else 0.5, alpha=0.8)
        ax.plot(*eq.xpt, "kx")
        ax.set_aspect("equal")
        if box:
            ax.set_xlim(box[0], box[1])
            ax.set_ylim(box[2], box[3])
        else:
            ax.autoscale_view()
        ax.set_xlabel("R [cm]")
        ax.set_ylabel("Z [cm]")
        ax.set_title("core (blue) + edge (red); contours: SOL green, PFR purple, legs black", fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "mesh_%s.png" % name), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()

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
    ap.add_argument("--no-legs", action="store_true", help="do not add the separatrix legs as constraints")
    ap.add_argument("--core", choices=["native", "cdt"], default="native",
                    help="native: keep grid_kind 2 triangles inside the last ring; cdt: re-triangulate")
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    eq = Equilibrium(a.gfile)
    RZ_in, s_in, tri_in = load_inner(a.dumpdir)
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

    # ---- contour levels: SOL measured along the outboard midplane from the last ring,
    #      PFR measured downward from the X-point.
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

    # SOL levels: until beyond the max psi_N inside the limiter
    pmax = max(eq.psin(lim[:, 0], lim[:, 1]).max(), 1.0) + 1e-3
    sol_levels, d, rpos = [], d1, r_ring
    # first SOL contour at distance d1 from the ring, then geometric
    while True:
        rpos = rpos + d
        lv = float(eq.psin(rpos, zmid))
        if lv > pmax:
            break
        if lv > 1.0 + 1e-4:
            sol_levels.append(lv)
        d = min(d * a.growth, a.dmax)
    # PFR levels: along a vertical line downward from the X-point
    zz = np.linspace(xpt[1], eq.lim[:, 1].min() - 5, 3000)
    pn_v = eq.psin(np.full_like(zz, xpt[0]), zz)
    pmin = eq.psin(lim[:, 0], lim[:, 1]).min()
    pfr_levels, d, zpos = [], d1, xpt[1]
    while True:
        zpos -= d
        if zpos < zz[-1]:
            break
        lv = float(np.interp(zpos, zz[::-1], pn_v[::-1]))
        if lv < pmin:
            break
        if lv < 1.0 - 1e-4:
            pfr_levels.append(lv)
        d = min(d * a.growth, a.dmax)
    print("SOL levels (%d):" % len(sol_levels), np.round(sol_levels, 4))
    print("PFR levels (%d):" % len(pfr_levels), np.round(pfr_levels, 4))

    # ---- contour lines on a fine grid
    Rg = np.linspace(lim[:, 0].min() - 2, lim[:, 0].max() + 2, 900)
    Zg = np.linspace(lim[:, 1].min() - 2, lim[:, 1].max() + 2, 1600)
    RR, ZZ = np.meshgrid(Rg, Zg)
    PN = eq.psin(RR, ZZ)
    gen = contourpy.contour_generator(Rg, Zg, PN, line_type="Separate")

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
        # separatrix legs: psi_N = 1 contour below the X-point, away from the X-point itself
        lv = 1.0
        for line in gen.lines(lv):
            if len(line) < 3:
                continue
            dx = np.linalg.norm(line - xpt, axis=1)
            below = (line[:, 1] < xpt[1] - 0.5 * a.hmin) & (dx > 3 * a.hmin)
            # split into runs where below is true
            idx = np.where(below)[0]
            if len(idx) == 0:
                continue
            breaks = np.where(np.diff(idx) > 1)[0]
            for part in np.split(idx, breaks + 1):
                sub = line[part]
                if len(sub) < 3:
                    continue
                for pts, sh, eh in clip_polyline(sub, lim, lim_path):
                    # keep only legs that reach the limiter (one end hit)
                    if sh is None and eh is None:
                        continue
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

    # limiter polygon, resampled, with contour end points inserted
    lim_ids = []
    for k in range(len(lim) - 1):
        pa, pb = lim[k], lim[k + 1]
        stops = [(0.0, None)] + sorted(lim_extra[k]) + [(1.0, None)]
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
    nfold = len(__import__("fold_check").folds(RZ_all, T_all))
    print("folded edges in final mesh: %d" % nfold)
    if nfold:
        raise SystemExit("final mesh folds")
    report = check_mesh(RZ_all, T_all, eq, RZ_in[np.isclose(s_in, s_in.min())])
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
               header="0 = grid_kind 2 inner mesh, 1 = edge extension")
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
    A_hole = 0.5 * abs(np.sum(H[:, 0] * np.roll(H[:, 1], -1) - np.roll(H[:, 0], -1) * H[:, 1]))
    rep["area_limiter_minus_hole_cm2"] = "%.3f" % (A_lim - A_hole)
    return rep


def plot_mesh(P, T, region, eq, outdir, lines):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    for name, box in (("full", None), ("xpoint", (120, 165, -125, -80)), ("outboard", (195, 225, -20, 40))):
        fig, ax = plt.subplots(figsize=(7, 10) if box is None else (8, 8))
        for reg, col in ((0, "tab:blue"), (1, "tab:red")):
            m = region == reg
            pc = PolyCollection(P[T[m]], facecolors="none", edgecolors=col, linewidths=0.25)
            ax.add_collection(pc)
        ax.plot(*eq.lim.T, "k-", lw=1)
        ax.plot(*eq.xpt, "kx")
        ax.set_aspect("equal")
        if box:
            ax.set_xlim(box[0], box[1])
            ax.set_ylim(box[2], box[3])
        else:
            ax.autoscale_view()
        ax.set_xlabel("R [cm]")
        ax.set_ylabel("Z [cm]")
        ax.set_title("AUG mesh: grid_kind 2 core (blue) + edge extension (red)")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "mesh_%s.png" % name), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()

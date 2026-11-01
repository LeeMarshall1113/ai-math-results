"""Vectorised cover13.solve: builds the SAME LP (bit-identical A, b, bounds; checked by test_equiv.py) with numpy
block assignment instead of a Python loop per facet, and solves it with the direct HiGHS call in hlp.py."""
import numpy as np
from cover13 import CUBE, EDGES, Ap, hrep
import hlp

_EU = [tuple(CUBE[b][i] - CUBE[a][i] for i in range(3)) for a, b in EDGES]

_CACHE = {}
def _dot3(HH, v):
    """bit-identical to the original per-row `n[:3] @ np.array(v, float)` (numpy's 1-D dot is an FMA chain, which a
    matrix-vector product does not reproduce), so it is evaluated row by row exactly as before -- but only once per
    (H-rep, point): rank() reuses the same H, Hm and points across up to 6000 structures."""
    key = (id(HH), v)
    r = _CACHE.get(key)
    if r is None:
        vv = np.array(v, float)
        r = np.array([n[:3] @ vv for n in HH]); _CACHE[key] = r
    return r

def clear_cache():
    _CACHE.clear()

def build(p, f0, v0, tau, H, Hm, facepts):
    Q = Ap(p)
    if H is None: H = hrep(Q)
    idx = {}; nv = 0
    for v in range(8): idx[('a', v)] = nv; nv += 3
    for f in range(6):
        if f != f0: idx[('b', f)] = nv; nv += 3
    tv = set(tau.values())
    for e in range(12):
        idx[('t', e)] = nv; nv += 1
        if e in tv: idx[('s', e)] = nv; nv += 1
    calls = []   # (const_pt, lin list, T)
    for v in range(8): calls.append((CUBE[v], [], idx[('a', v)]))
    for f in range(6):
        T = idx[('a', v0)] if f == f0 else idx[('b', f)]
        for pt in (facepts[f] if facepts is not None else [Q[f]]): calls.append((pt, [], T))
    tauinv = {e: f for f, e in tau.items() if e is not None}
    for e, (a, b) in enumerate(EDGES):
        u = _EU[e]
        if e not in tauinv:
            calls.append((CUBE[a], [(idx[('t', e)], u)], idx[('a', a)]))
            calls.append((CUBE[a], [(idx[('t', e)], u)], idx[('a', b)]))
        else:
            B = idx[('b', tauinv[e])]
            calls.append((CUBE[a], [(idx[('t', e)], u)], idx[('a', a)]))
            calls.append((CUBE[a], [(idx[('s', e)], u)], idx[('a', b)]))
            calls.append((CUBE[a], [(idx[('t', e)], u)], B))
            calls.append((CUBE[a], [(idx[('s', e)], u)], B))
    Tm = idx[('a', v0)]
    sizes = [len(Hm if (Hm is not None and T == Tm) else H) for _, _, T in calls]
    M = sum(sizes)
    A = np.zeros((M, nv + 1)); rhs = np.empty(M)
    r0 = 0
    for (pt, lin, T), sz in zip(calls, sizes):
        HH = Hm if (Hm is not None and T == Tm) else H
        Hn = HH[:, :3]
        rhs[r0:r0 + sz] = -HH[:, 3] - _dot3(HH, tuple(pt))
        for k, dvec in lin: A[r0:r0 + sz, k] += _dot3(HH, dvec)
        A[r0:r0 + sz, T:T + 3] -= Hn
        A[r0:r0 + sz, nv] = 1.0
        r0 += sz
    return A, rhs, nv

def solve(p, f0, v0, tau, H=None, Hm=None, facepts=None):
    A, rhs, nv = build(p, f0, v0, tau, H, Hm, facepts)
    c = np.zeros(nv + 1); c[-1] = -1
    lb = np.full(nv + 1, -hlp.INF); ub = np.full(nv + 1, hlp.INF); ub[-1] = 1.0
    ok, x = hlp.solve_ub(c, A, rhs, lb, ub)
    return x[-1] if ok else None

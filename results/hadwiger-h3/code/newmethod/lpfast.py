"""Fast float LP + vectorised exact integer check for reduction R1' (same mathematics as lpr1.py).

A structure gives a list of containment ITEMS.  Item = (base point b (3 Fractions), list of (var, direction u),
translate index T, merged flag).  Condition: z := b + sum var*u - a_T [- mu*lam*t*s if merged] lies in Q (or (1-lam)Q
if merged), Q = Q_P given by integer facets c0 + n.x >= 0.

Float stage: unit-normal reduced facets, maximise uniform slack.  Exact stage: round the solution to rationals with
denominator DEN, compute every z exactly (integers over a common denominator), and test it against the FULL integer
H-rep of Q_P with numpy object (python int) arithmetic.  The returned certificate (integer numerators over DEN) is
what certcheck.py verifies standalone."""
import itertools, math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog

CUBE = [tuple(v) for v in itertools.product((0, 1), repeat=3)]
EDGES = [(a, b) for a, b in itertools.combinations(range(8), 2)
         if sum(abs(x - y) for x, y in zip(CUBE[a], CUBE[b])) == 1]
FACES = [(i, j) for i in range(3) for j in (0, 1)]
def face_vertices(f): i, j = FACES[f]; return [k for k in range(8) if CUBE[k][i] == j]
def face_edges(f):
    fv = set(face_vertices(f)); return [e for e, (a, b) in enumerate(EDGES) if a in fv and b in fv]
def chord_dir(axis, sigma):
    s = [0, 0, 0]; j, k = [x for x in range(3) if x != axis]; s[j] = 1; s[k] = sigma; return s

LAM = F(7, 10); T = F(569, 1000); DEN = 10 ** 9

def items_for(R, f0, v0, combo, merged):
    """R: the six rectangles (lists of 4 Fraction points).  Returns (items, nvar, nmu).
    Variables: a_v (v=0..7) -> 3v..3v+2 ; b_f for f != f0 -> 24+3*rank ; then edge params ; then mu's."""
    others = [f for f in range(6) if f != f0]
    tau = dict(zip(others, combo))
    Tidx = {('a', v): 3 * v for v in range(8)}
    for r, f in enumerate(others): Tidx[('b', f)] = 24 + 3 * r
    nv = 39
    tauinv = {e: f for f, e in tau.items() if e is not None}
    evar = {}
    for e in range(12):
        evar[('t', e)] = nv; nv += 1
        if e in tauinv: evar[('s', e)] = nv; nv += 1
    items = []
    def add(base, lin, Tkey):
        items.append((tuple(F(x) for x in base), lin, Tidx[Tkey], merged and Tkey == ('a', v0)))
    for v in range(8): add(CUBE[v], [], ('a', v))
    for f in range(6):
        key = ('a', v0) if f == f0 else ('b', f)
        for pt in R[f]: add(pt, [], key)
    for e, (u0, u1) in enumerate(EDGES):
        u = tuple(CUBE[u1][i] - CUBE[u0][i] for i in range(3))
        if e not in tauinv:
            add(CUBE[u0], [(evar[('t', e)], u)], ('a', u0)); add(CUBE[u0], [(evar[('t', e)], u)], ('a', u1))
        else:
            B = ('b', tauinv[e])
            add(CUBE[u0], [(evar[('t', e)], u)], ('a', u0)); add(CUBE[u0], [(evar[('s', e)], u)], ('a', u1))
            add(CUBE[u0], [(evar[('t', e)], u)], B); add(CUBE[u0], [(evar[('s', e)], u)], B)
    nmu = sum(1 for it in items if it[3])
    return items, nv, nmu

def float_solve(items, nv, nmu, Nf, cf, s_vec):
    """Nf (m x 3) unit normals, cf (m) offsets: inside <=> cf + Nf.x >= 0.  Variables: nv + nmu + 1 (slack)."""
    m = len(cf); nt = nv + nmu + 1
    rows = []; rhs = []; mu = nv
    sv = np.array([float(x) for x in s_vec]) if s_vec is not None else None
    lamf = float(LAM)
    for base, lin, Ti, mg in items:
        A = np.zeros((m, nt))
        b = np.array([float(x) for x in base])
        # constraint: c' + N.(b + sum var u - a_T - mu*lam*t*s) >= slack
        cc = (cf * (1 - lamf) if mg else cf) + Nf @ b
        for k, u in lin: A[:, k] -= Nf @ np.array(u, float)
        A[:, Ti:Ti + 3] += Nf
        if mg:
            A[:, mu] += Nf @ (float(LAM * T) * sv); mu += 1
        A[:, -1] = 1.0
        rows.append(A); rhs.append(cc)
    A = np.vstack(rows); b = np.concatenate(rhs)
    bounds = [(None, None)] * nv + [(0, 1)] * nmu + [(None, 1.0)]
    c = np.zeros(nt); c[-1] = -1
    r = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
    if r.status != 0: return None, None
    return r.x[-1], r.x[:-1]

def round_solution(x, nv, nmu):
    q = [int(round(v * DEN)) for v in x]
    for k in range(nv, nv + nmu): q[k] = min(max(q[k], 0), DEN)
    return q

def exact_ok(items, nv, nmu, q, Hint, s_vec):
    """q: integer numerators over DEN (len nv+nmu).  Hint: list of integer (c0,n1,n2,n3) for the full Q_P.
    Merged items: need (1-lam) c0 + n.z >= 0, i.e. 3 c0 + 10 n.z >= 0 (lam = 7/10)."""
    # every z has denominator DEN * L where L = lcm of base denominators; compute numerators per item
    Zs = []; Ds = []; M = []
    mu = nv
    for base, lin, Ti, mg in items:
        L = 1
        for x in base: L = L * x.denominator // math.gcd(L, x.denominator)
        D = DEN * L
        z = []
        for i in range(3):
            v = base[i].numerator * (D // base[i].denominator)
            for k, u in lin: v += q[k] * u[i] * L
            v -= q[Ti + i] * L
            if mg:
                # mu*lam*t*s_i with lam*t = 3983/10000: write over D*10000 -> scale everything by 10000
                pass
            z.append(v)
        if mg:
            lt = LAM * T  # Fraction 3983/10000
            z = [zi * lt.denominator - q[mu] * L * lt.numerator * s_vec[i] for i, zi in enumerate(z)]
            D = D * lt.denominator; mu += 1
        Zs.append(z); Ds.append(D); M.append(mg)
    H = np.array(Hint, dtype=object)
    c0 = H[:, 0]; N = H[:, 1:]
    Z = np.array(Zs, dtype=object).T          # 3 x k
    Dv = np.array(Ds, dtype=object)
    vals_plain = c0[:, None] * Dv[None, :] + N.dot(Z)            # c0 D + n.Z
    Mv = np.array(M)
    if Mv.any():
        vals_m = 3 * c0[:, None] * Dv[None, :] + 10 * N.dot(Z)
        vals = np.where(Mv[None, :], vals_m, vals_plain)
    else:
        vals = vals_plain
    return bool((vals >= 0).all())

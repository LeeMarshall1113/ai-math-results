"""Standalone checker for the certify2.py certificate of H3 <= 13 (reduction R1').

Checks, using only the stored data plus hw.QP_hrep (the original exact Q_P code, independent of qpfast.py):
 1. coverage: the set of start boxes equals the set of 1/10-grid boxes containing a box of Prymak's U^D_20 (M=10),
    each appears exactly once, and inside each start box the leaves form an exact longest-side bisection tree
    (so they tile it with disjoint interiors);
 2. every leaf: its stored rational solution(s) satisfy every containment constraint EXACTLY against the full Q_P
    H-rep (plain: 13 translates of Q_P; chord: for both signs, 12 translates of Q_P and one of
    (1-lam)Q_P + lam t [0, e_j +- e_k]).
Usage: python certcheck.py OUTDIR [--sample K]   (K: check only every K-th leaf's LP part; coverage always full)"""
import sys, os, json, gzip, itertools, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
_v = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_vendor')
if os.path.isdir(_v) and _v not in sys.path: sys.path.insert(0, _v)
_u = r'C:\Users\hackathon\AppData\Roaming\Python\Python314\site-packages'
if os.path.isdir(_u) and _u not in sys.path: sys.path.append(_u)
from fractions import Fraction as F
import numpy as np
import hw

LAM = F(7, 10); T = F(569, 1000); DEN = 10 ** 9
CUBE = [tuple(v) for v in itertools.product((0, 1), repeat=3)]
EDGES = [(a, b) for a, b in itertools.combinations(range(8), 2)
         if sum(abs(x - y) for x, y in zip(CUBE[a], CUBE[b])) == 1]

def chord_dir(axis, sigma):
    s = [0, 0, 0]; j, k = [x for x in range(3) if x != axis]; s[j] = 1; s[k] = sigma; return s

def constraints(R, f0, v0, combo, merged):
    """independent re-implementation of the structure -> list of (base, lin, Tvar, merged)."""
    others = [f for f in range(6) if f != f0]
    assert len(combo) == 5
    es = [e for e in combo if e is not None]; assert len(set(es)) == len(es)
    tau = dict(zip(others, combo))
    Tv = {('a', v): 3 * v for v in range(8)}
    for r, f in enumerate(others): Tv[('b', f)] = 24 + 3 * r
    nv = 39; tauinv = {e: f for f, e in tau.items() if e is not None}; ev = {}
    for e in range(12):
        ev[('t', e)] = nv; nv += 1
        if e in tauinv: ev[('s', e)] = nv; nv += 1
    out = []
    def add(base, lin, key): out.append(([F(x) for x in base], lin, Tv[key], merged and key == ('a', v0)))
    for v in range(8): add(CUBE[v], [], ('a', v))
    for f in range(6):
        for pt in R[f]: add(pt, [], ('a', v0) if f == f0 else ('b', f))
    for e, (u0, u1) in enumerate(EDGES):
        u = [CUBE[u1][i] - CUBE[u0][i] for i in range(3)]
        if e not in tauinv:
            add(CUBE[u0], [(ev[('t', e)], u)], ('a', u0)); add(CUBE[u0], [(ev[('t', e)], u)], ('a', u1))
        else:
            B = ('b', tauinv[e])
            add(CUBE[u0], [(ev[('t', e)], u)], ('a', u0)); add(CUBE[u0], [(ev[('s', e)], u)], ('a', u1))
            add(CUBE[u0], [(ev[('t', e)], u)], B); add(CUBE[u0], [(ev[('s', e)], u)], B)
    nmu = sum(1 for c in out if c[3])
    return out, nv, nmu

def check_solution(cons, nv, nmu, q, H, s_vec):
    if len(q) != nv + nmu: return False
    qs = [F(x, DEN) for x in q]
    mus = qs[nv:]
    if any(m < 0 or m > 1 for m in mus): return False
    Z = []; mg = []; k = 0
    for base, lin, Ti, merged in cons:
        z = [base[i] - qs[Ti + i] + sum(qs[v] * u[i] for v, u in lin) for i in range(3)]
        if merged:
            z = [z[i] - mus[k] * LAM * T * s_vec[i] for i in range(3)]; k += 1
        Z.append(z); mg.append(merged)
    # exact: c0*w + n.z >= 0 with w = 1 or (1-lam); clear denominators per point
    Hc = np.array([h[0] for h in H], dtype=object); Hn = np.array([h[1:] for h in H], dtype=object)
    for z, merged in zip(Z, mg):
        L = 1
        for x in z: L = L * x.denominator // math.gcd(L, x.denominator)
        zi = np.array([int(x * L) for x in z], dtype=object)
        w = (1 - LAM) if merged else F(1)
        vals = Hc * (w.numerator * L) + Hn.dot(zi) * w.denominator
        if (vals < 0).any(): return False
    return True

def check_leaf(lf):
    a = [F(x) for x in lf['a']]; d = [F(x) for x in lf['d']]
    H = hw.QP_hrep(a, d); R = hw.rects(a, d)
    if lf['k'] == 'p':
        c = lf['c']; cons, nv, nmu = constraints(R, c['st'][0], c['st'][1], tuple(c['st'][2]), False)
        return check_solution(cons, nv, nmu, c['q'], H, None)
    if lf['k'] == 'c':
        for sg, c in ((1, lf['plus']), (-1, lf['minus'])):
            cons, nv, nmu = constraints(R, c['st'][0], c['st'][1], tuple(c['st'][2]), True)
            if not check_solution(cons, nv, nmu, c['q'], H, chord_dir(lf['ax'], sg)): return False
        return True
    return False

def tiles(a, d, leaves):
    """leaves: list of (a,d) tuples of Fractions; True iff they are exactly the leaves of a longest-side bisection of (a,d)."""
    if len(leaves) == 1 and leaves[0] == (tuple(a), tuple(d)): return True
    if any(l == (tuple(a), tuple(d)) for l in leaves) or not leaves: return False
    if any(any(l[1][k] > d[k] for k in range(6)) for l in leaves): return False
    i = max(range(6), key=lambda k: d[k])          # same tie rule as certify2 (first maximal index)
    d2 = list(d); d2[i] = d[i] / 2; a2 = list(a); a2[i] = a[i] + d2[i]
    L1 = [l for l in leaves if l[0][i] < a2[i]]; L2 = [l for l in leaves if l[0][i] >= a2[i]]
    for la, ld in L1:
        if not all(a[k] <= la[k] and la[k] + ld[k] <= (a[k] + d2[k] if k == i else a[k] + d[k]) for k in range(6)): return False
    return tiles(a, d2, L1) and tiles(a2, d2, L2)

def in_UD(k):
    k1, k2, k3, k4, k5, k6 = k
    return k2 + k4 <= k1 + k6 + 1 <= k3 + k5 + 2

if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    out = sys.argv[1]; samp = 1
    if '--sample' in sys.argv: samp = int(sys.argv[sys.argv.index('--sample') + 1])
    need = set(c for c in itertools.product(range(5), range(5), range(5), range(10), range(10), range(10))
               if any(in_UD([2 * c[i] + e[i] for i in range(6)]) for e in itertools.product((0, 1), repeat=6)))
    seen = set(); nleaf = 0; nchk = 0; bad = []
    for fn in sorted(os.listdir(out)):
        if not fn.endswith('.jsonl.gz'): continue
        for line in gzip.open(os.path.join(out, fn), 'rt'):
            r = json.loads(line); sb = tuple(r['sb'])
            if sb in seen: bad.append(('dup', sb)); continue
            seen.add(sb)
            a0 = [F(sb[i], 10) for i in range(6)]; d0 = [F(1, 10)] * 6
            L = [(tuple(F(x) for x in lf['a']), tuple(F(x) for x in lf['d'])) for lf in r['leaves']]
            if not tiles(a0, d0, L): bad.append(('tiling', sb)); continue
            for lf in r['leaves']:
                nleaf += 1
                if nleaf % samp == 0:
                    nchk += 1
                    if not check_leaf(lf): bad.append(('leaf', sb, lf['a'], lf['d']))
    missing = need - seen; extra = seen - need
    print('start boxes needed', len(need), 'present', len(seen), 'missing', len(missing), 'extra', len(extra))
    print('leaves', nleaf, 'LP-checked', nchk, 'bad', len(bad))
    if bad: print('first bad:', bad[:5])
    print('COMPLETE AND VALID' if not missing and not bad and not extra else 'NOT (YET) COMPLETE/VALID')

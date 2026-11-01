"""Unit equivalence tests for the fast helpers (run: python test_equiv.py).  Each must print 0 mismatches.
 1. hlp.solve_ub == scipy linprog(method='highs'): status and bit-identical x on random LPs
 2. qpfastf.QP_hrep_fast == qpfast.QP_hrep_fast (same sorted integer rows)
 3. vectorised tightness flags in certify2f.box_data == the per-row loop of certify2.box_data
 4. lpfastf.exact_ok == lpfast.exact_ok (fast float filter + exact fallback), including the forced-fallback path
 (cover13f.build == cover13.solve's matrix, bit for bit: test_cover13f.py)"""
import sys, os, random, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..')); sys.path.insert(0, os.path.join(here, '..')); sys.path.insert(0, here)
from fractions import Fraction as F
from scipy.optimize import linprog
import hlp, qpfast, qpfastf, certify2f as C, hw
LF = C.LF
rng = np.random.default_rng(0); bad = 0
for t in range(300):
    m, n = rng.integers(20, 400), rng.integers(5, 60)
    A = rng.normal(size=(m, n)); A[rng.random((m, n)) < 0.6] = 0; b = rng.random(m) + 0.1
    c = rng.normal(size=n); bnds = [(None, None)] * (n - 1) + [(None, 1.0)] if t % 2 == 0 else [(0, 1)] * n
    r = linprog(c, A_ub=A, b_ub=b, bounds=bnds, method='highs'); lb, ub = hlp.bounds_arrays(bnds)
    ok, x = hlp.solve_ub(c, A, b, lb, ub)
    bad += (r.status == 0) != ok or (ok and not np.array_equal(r.x, x))
print('1 hlp vs linprog mismatches', bad)
rnd = random.Random(0); bad = 0
def rbox(rnd):
    k = rnd.choice([10, 20, 40, 80, 160, 320, 640]); d = [F(1, k)] * 6
    for i in range(rnd.randint(0, 5)):
        if k < 640: d[i] = F(1, 2 * k)
    return [F(rnd.randrange(0, (5 if i < 3 else 10) * k // 10), k) for i in range(6)], d
for t in range(400):
    a, d = rbox(rnd); bad += qpfast.QP_hrep_fast(a, d) != qpfastf.QP_hrep_fast(a, d)
print('2 QP_hrep mismatches', bad)
bad = 0
for t in range(1000):
    a, d = rbox(rnd); bd = C.box_data(a, d)
    if bd is None: continue
    Hint, V = bd[0], bd[3]
    A = np.array([[float(h[1]), float(h[2]), float(h[3])] for h in Hint]); b = np.array([float(h[0]) for h in Hint])
    nr = np.linalg.norm(A, axis=1); A = A / nr[:, None]; b = b / nr
    ref = np.array([np.sum(np.abs(V @ row + bb) < 1e-9) >= 3 for row, bb in zip(A, b)])
    bad += not np.array_equal(ref, (np.abs(A @ V.T + b[:, None]) < 1e-9).sum(1) >= 3)
print('3 tightness mismatches', bad)
for rel in (LF.FILTER_REL, 10.0):
    LF.FILTER_REL = rel; bad = 0; npos = 0
    for t in range(60):
        a, d = rbox(random.Random(t))
        bd = C.box_data(a, d)
        if bd is None: continue
        Hint, Nf, cf, V = bd; R = hw.rects(a, d)
        for st in C.ALL[:3]:
            for merged in (False, True):
                s_vec = LF.chord_dir(t % 3, 1 - 2 * (t % 2)) if merged else None
                items, nv, nmu = LF.items_for(R, st[0], st[1], st[2], merged)
                m, x = LF.float_solve(items, nv, nmu, Nf, cf, s_vec)
                if m is None: continue
                q = LF.round_solution(x, nv, nmu)
                for qq in (q, [v + rnd.randrange(-3, 4) for v in q]):
                    r1 = LF._exact_ok_slow(items, nv, nmu, qq, Hint, s_vec); r2 = LF.exact_ok(items, nv, nmu, qq, Hint, s_vec)
                    bad += r1 != r2; npos += r1
    print('4 exact_ok mismatches (FILTER_REL=%g)' % rel, bad, '(true cases %d)' % npos)

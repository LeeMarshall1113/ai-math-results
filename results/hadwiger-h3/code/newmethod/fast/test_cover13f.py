import sys, os, itertools, numpy as np
_h = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(_h, '..')); sys.path.insert(0, _h)
import cover13 as C, cover13f as CF, lpfast as LF
import scipy.optimize
cap = {}
orig_lp = C.linprog
def spy(c, A_ub, b_ub, bounds, method):
    cap['A'] = A_ub; cap['b'] = b_ub; cap['bounds'] = bounds; return orig_lp(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method=method)
C.linprog = spy
rng = np.random.default_rng(1); bad = 0; n = 0
import time; t1 = t2 = 0
ALL = []
for f0 in range(6):
    for v0 in LF.face_vertices(f0):
        others = [f for f in range(6) if f != f0]
        for combo in itertools.product(*[LF.face_edges(f) + [None] for f in others]):
            es = [e for e in combo if e is not None]
            if len(set(es)) == len(es): ALL.append((f0, v0, combo))
for trial in range(300):
    p = list(rng.random(6) * [0.5, 0.5, 0.5, 1, 1, 1])
    Q = C.Ap(p); V = np.array(Q) + rng.normal(scale=0.01, size=(6, 3)); V = np.vstack([V, V + 0.02])
    H = C.hrep(V); Hm = C.hrep(0.3 * V) if trial % 3 == 0 else None
    fp = [[tuple(np.array(q) + rng.normal(scale=0.01, size=3)) for _ in range(4)] for q in Q]
    f0, v0, combo = ALL[rng.integers(len(ALL))]
    tau = dict(zip([f for f in range(6) if f != f0], combo))
    t = time.perf_counter(); m1 = C.solve(p, f0, v0, tau, H, Hm, fp); t1 += time.perf_counter() - t
    CF.clear_cache(); A, b, nv = CF.build(p, f0, v0, tau, H, Hm, fp)
    t = time.perf_counter(); m2 = CF.solve(p, f0, v0, tau, H, Hm, fp); t2 += time.perf_counter() - t
    n += 1
    if not (np.array_equal(A, cap['A']) and np.array_equal(b, cap['b'])) or not (m1 == m2 or (m1 is None and m2 is None)):
        bad += 1; print('diff', trial, m1, m2, np.abs(A - cap['A']).max(), np.abs(b - cap['b']).max())
print('tested', n, 'mismatch', bad, 'orig %.3fs fast %.3fs' % (t1, t2))

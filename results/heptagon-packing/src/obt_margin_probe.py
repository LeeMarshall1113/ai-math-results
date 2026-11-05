"""Probe the obtuse-mode slack  b - acrit - extra(T)  at random feasible canonical points in a
theta0 range, then locally refine the worst ones (Nelder-Mead on the slack, rejecting infeasible)."""
import sys, os, numpy as np, json, time
sys.path.insert(0, 'src')
import bnb2 as B2
B2.OBT = B2.ObtTable('runs/D3_table.npy', 'runs/D3_bins.json')
F = B2.FTrans('runs/f_w030.npy', 3, 2, mu=[-1.1684855e-12, -1.1684853e-12])
T0LO, T0HI, N = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(1); ACR = 1.5326754446782214
def slack(X):
    P = np.stack([X, X], 2)
    inf, _ = B2.infeasible(P)
    b = B2.b_A(F, [B2.R(X[:, q]) for q in range(6)], 1).v.lo
    ex = B2.obtuse_extra(P)
    return np.where(inf, np.inf, b - ACR - ex), b - ACR, ex
t0 = time.time(); pts = []
while sum(len(p) for p in pts) < N and time.time()-t0 < 150:
    m = 20000
    X = np.column_stack([rng.uniform(2.548, 4.0, m), rng.uniform(0, 2.0, m), rng.uniform(0.8, 2.2, m),
                         rng.uniform(T0LO, T0HI, m), rng.uniform(0, 0.8976, m), rng.uniform(0, 0.8976, m)])
    X[:, 1] = np.minimum(X[:, 1], X[:, 0]/2)
    s, bb, ex = slack(X); ok = np.isfinite(s)
    pts.append(np.column_stack([X[ok], s[ok], bb[ok], ex[ok]]))
P = np.concatenate(pts); print('feasible samples', len(P), f'[{time.time()-t0:.0f}s]')
o = np.argsort(P[:, 6]); print('min slack samples:', np.round(P[o[:5], 6], 5)); print(np.round(P[o[:5], :6], 4)); print('b-acrit', np.round(P[o[:5], 7], 4), 'extra', np.round(P[o[:5], 8], 4))
from scipy.optimize import minimize
best = []
for i in o[:6]:
    f = lambda x: float(slack(x[None])[0][0])
    r = minimize(f, P[i, :6], method='Nelder-Mead', options={'maxiter': 300, 'xatol': 1e-7, 'fatol': 1e-9})
    best.append((r.fun, r.x))
best.sort(key=lambda t: t[0])
print('refined min slack:', [round(v, 6) for v, _ in best[:6]])
x = best[0][1]; s, bb, ex = slack(x[None]); print('worst point', np.round(x, 5), 'b-acrit', bb[0], 'extra', ex[0])
d, X1, Y = x[0], x[1], x[2]; u = np.array([-X1, -Y]); v = np.array([d-X1, -Y]); print('long edge', d, 'cos gamma', u@v/np.linalg.norm(u)/np.linalg.norm(v), 'short edges', np.hypot(X1-d, Y), np.hypot(X1, Y))

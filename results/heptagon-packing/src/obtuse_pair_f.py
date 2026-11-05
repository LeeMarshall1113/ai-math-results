"""Obtuse triangle + neighbour across its longest edge, with the transfer f: min (b(T)+b(T'))/2."""
import numpy as np, sys, time
CP = sys.argv[1]; NST = int(sys.argv[2]); TLIM = float(sys.argv[3])
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'nonobtuse', '0.3', '0', '0.05']
sys.path.insert(0, 'src')
exec(open('src/lp_transfer3.py').read().split("if __name__")[0])
from min_dimer import cons as dcons, unpack
from scipy.optimize import minimize
c = np.load(CP)
def bpair(x):
    Cx, Tx = unpack(x)
    C2 = np.array([[Cx[0], Cx[1], Cx[2]], [Cx[0], Cx[3], Cx[1]]]); T2 = np.array([[Tx[0], Tx[1], Tx[2]], [Tx[0], Tx[3], Tx[1]]])
    a, M = rows(C2, T2); b = a + M @ c
    return b.sum(), b
def cons(x):
    Cx, Tx = unpack(x); g = list(dcons(x, False))
    g.append(-np.dot(Cx[0]-Cx[2], Cx[1]-Cx[2]))
    return np.array(g)
rng = np.random.default_rng(9); res = []; t0 = time.time()
for it in range(NST):
    if time.time()-t0 > TLIM: break
    x0 = np.array([rng.uniform(2.6, 3.9), rng.uniform(0.5, 2.5), rng.uniform(0.6, 1.6),
                   rng.uniform(0, 3), -rng.uniform(1.0, 3), *rng.uniform(0, 2*np.pi/7, 4)])
    x0[1] = min(x0[1], x0[0]-0.3)
    r = minimize(lambda x: bpair(x)[0], x0, method='SLSQP', constraints=[{'type': 'ineq', 'fun': cons}],
                 bounds=[(1.7, 4), (-4, 4), (0.01, 4), (-4, 4), (-4, -0.01)]+[(-1, 2)]*4, options={'maxiter': 500, 'ftol': 1e-12})
    if cons(r.x).min() > -1e-7: res.append((bpair(r.x)[0], r.x))
res.sort(key=lambda t: t[0])
print('starts', it+1, 'feasible', len(res), 'best avg/acrit', np.round([a/(2*ACRIT) for a, _ in res[:8]], 5))
a, x = res[0]; print('b per triangle', bpair(x)[1], 'x', np.round(x, 4))

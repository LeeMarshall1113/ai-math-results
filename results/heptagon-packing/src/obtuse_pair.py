"""Obtuse triangle + its neighbour across the longest edge: min total area (HK-style N rule test)."""
import numpy as np, sys, time
sys.path.insert(0, 'src')
from hept import *
from min_dimer import cons as dcons, area as darea, unpack, ACRIT
from scipy.optimize import minimize
def cons(x):
    C, T = unpack(x); g = list(dcons(x, False))
    g.append(-np.dot(C[0]-C[2], C[1]-C[2]))        # obtuse (or right) at c2
    return np.array(g)
rng = np.random.default_rng(9); res = []; t0 = time.time()
for it in range(400):
    if time.time()-t0 > 500: break
    x0 = np.array([rng.uniform(2.6, 3.9), rng.uniform(0.5, 2.5), rng.uniform(0.6, 1.6),
                   rng.uniform(0, 3), -rng.uniform(1.0, 3), *rng.uniform(0, 2*np.pi/7, 4)])
    x0[1] = min(x0[1], x0[0]-0.3)
    r = minimize(darea, x0, method='SLSQP', constraints=[{'type': 'ineq', 'fun': cons}],
                 bounds=[(1.7, 4), (-4, 4), (0.01, 4), (-4, 4), (-4, -0.01)]+[(-1, 2)]*4,
                 options={'maxiter': 500, 'ftol': 1e-12})
    if cons(r.x).min() > -1e-7: res.append((darea(r.x), r.x))
res.sort(key=lambda t: t[0])
print('starts', it+1, 'feasible', len(res), 'best totals', np.round([a for a, _ in res[:8]], 5))
a, x = res[0]; C, T = unpack(x)
print('best total', a, 'avg/acrit', a/(2*ACRIT))
for tri in ((0, 1, 2), (0, 3, 1)):
    P = C[list(tri)]; print('  tri', round(abs(tri_area(*P)), 5), 'edges', np.round([np.linalg.norm(P[(k+1)%3]-P[k]) for k in range(3)], 4), 'R', round(circumradius(*P), 4))

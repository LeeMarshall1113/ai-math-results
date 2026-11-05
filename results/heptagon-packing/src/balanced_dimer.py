"""Obstruction test for 'single triangle + antisymmetric pair transfer' certificates.
Dimer sharing a translation edge (theta0=theta1); the four outer edges are point-reflection
pairs (theta2=theta3=theta0+pi/7), on which any antisymmetric pair transfer vanishes.
So min total area of such a dimer < 2*acrit would kill that certificate class."""
import numpy as np, sys, time
sys.path.insert(0, 'src')
from hept import *
from min_dimer import cons as dcons, area as darea, ACRIT
from scipy.optimize import minimize

def full(y):
    # y = (d, x2, y2, x3, y3, th0)
    t = y[5]; return np.array([*y[:5], t, t, t+np.pi/7, t+np.pi/7])

NONOB = len(sys.argv) > 3 and sys.argv[3] == "nonobtuse"
if __name__ == "__main__":
    n = int(sys.argv[1]); tl = float(sys.argv[2]); rng = np.random.default_rng(5); res = []; t0 = time.time()
    for it in range(n):
        if time.time()-t0 > tl: break
        y0 = np.array([rng.uniform(1.8, 2.4), rng.uniform(-0.3, 2.2), rng.uniform(1.2, 2.2),
                       rng.uniform(-0.3, 2.2), -rng.uniform(1.2, 2.2), rng.uniform(0, 2*np.pi/7)])
        r = minimize(lambda y: darea(full(y)), y0, method='SLSQP',
                     constraints=[{'type': 'ineq', 'fun': lambda y: dcons(full(y), NONOB)}],
                     bounds=[(1.7, 4), (-4, 4), (0.01, 4), (-4, 4), (-4, -0.01), (-1, 2)],
                     options={'maxiter': 500, 'ftol': 1e-13})
        if dcons(full(r.x), NONOB).min() > -1e-7: res.append((darea(full(r.x)), r.x))
    res.sort(key=lambda t: t[0])
    print('starts', it+1, 'feasible', len(res))
    print('best totals', np.round([a for a, _ in res[:15]], 6))
    a, y = res[0]; x = full(y)
    print('best', a, '2acrit', 2*ACRIT, 'ratio', a/(2*ACRIT))
    C = np.array([[0, 0], [x[0], 0], [x[1], x[2]], [x[3], x[4]]])
    for tri in ((0, 1, 2), (0, 3, 1)):
        P = C[list(tri)]; print('  tri', abs(tri_area(*P)), 'edges', np.round([np.linalg.norm(P[(k+1)%3]-P[k]) for k in range(3)], 4), 'R', round(circumradius(*P), 4))
    np.save('runs/balanced_dimer_best_' + ('nonob' if NONOB else 'any') + '.npy', x)

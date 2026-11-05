"""Numerical (non-rigorous) minimum area of a Delaunay triangle of three non-overlapping
regular heptagons (circumradius 1), with circumradius <= 2; optionally nonobtuse."""
import numpy as np, sys, time
sys.path.insert(0, 'src')
from hept import *
from scipy.optimize import minimize

def unpack(x):
    C = np.array([[0, 0], [x[0], 0], [x[1], x[2]]]); T = x[3:6]
    return C, T

def cons(x, nonobtuse):
    C, T = unpack(x)
    H = [verts(C[i], T[i]) for i in range(3)]
    g = [sep(H[0], H[1]), sep(H[0], H[2]), sep(H[1], H[2]), 2 - circumradius(*C)]
    if nonobtuse:
        for k in range(3):
            a, b = C[(k+1) % 3]-C[k], C[(k+2) % 3]-C[k]
            g.append(np.dot(a, b))
    return np.array(g)

def area(x):
    return 0.5*x[0]*x[2]

def run(n, nonobtuse, seed, tlimit):
    rng = np.random.default_rng(seed); best = (np.inf, None); res_all = []
    t0 = time.time()
    for it in range(n):
        if time.time()-t0 > tlimit: break
        x0 = np.array([rng.uniform(1.8, 2.6), rng.uniform(-0.5, 2.5), rng.uniform(1.0, 2.4),
                       *rng.uniform(0, 2*np.pi/7, 3)])
        r = minimize(area, x0, method='SLSQP',
                     constraints=[{'type': 'ineq', 'fun': lambda x: cons(x, nonobtuse)}],
                     bounds=[(1.7, 4), (-4, 4), (0.01, 4), (-1, 2), (-1, 2), (-1, 2)],
                     options={'maxiter': 300, 'ftol': 1e-12})
        g = cons(r.x, nonobtuse)
        if g.min() > -1e-7:
            res_all.append((area(r.x), r.x))
            if area(r.x) < best[0]: best = (area(r.x), r.x)
    return best, res_all, it+1

if __name__ == '__main__':
    nonobtuse = sys.argv[1] == 'nonobtuse'; n = int(sys.argv[2]); tl = float(sys.argv[3])
    best, allr, its = run(n, nonobtuse, 1, tl)
    print('mode', sys.argv[1], 'starts', its, 'feasible', len(allr))
    vals = sorted(a for a, _ in allr)
    print('best areas', np.round(vals[:12], 6))
    a, x = best; C, T = unpack(x)
    print('best', a, 'x', np.round(x, 6))
    L = [np.linalg.norm(C[(k+1) % 3]-C[k]) for k in range(3)]
    print('edges', np.round(L, 5), 'R', circumradius(*C))
    print('acrit 1.5326754447; ratio', a/1.5326754447)

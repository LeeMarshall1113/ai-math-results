"""Numerical min of total area of two Delaunay triangles sharing an edge (4 heptagons)."""
import numpy as np, sys, time
sys.path.insert(0, 'src')
from hept import *
from scipy.optimize import minimize
ACRIT = 1.5326754446782214

def unpack(x):
    C = np.array([[0, 0], [x[0], 0], [x[1], x[2]], [x[3], x[4]]]); return C, x[5:9]

def ang(C, k, i, j):
    a, b = C[i]-C[k], C[j]-C[k]
    return np.arccos(np.clip(np.dot(a, b)/np.linalg.norm(a)/np.linalg.norm(b), -1, 1))

def cons(x, nonobtuse):
    C, T = unpack(x)
    H = [verts(C[i], T[i]) for i in range(4)]
    g = [sep(H[i], H[j]) for i in range(4) for j in range(i+1, 4)]
    g += [2-circumradius(C[0], C[1], C[2]), 2-circumradius(C[0], C[3], C[1])]
    g += [np.pi - ang(C, 2, 0, 1) - ang(C, 3, 0, 1)]          # local Delaunay
    if nonobtuse:
        for tri in ((0, 1, 2), (0, 3, 1)):
            for k in range(3):
                a, b = C[tri[(k+1) % 3]]-C[tri[k]], C[tri[(k+2) % 3]]-C[tri[k]]
                g.append(np.dot(a, b))
    return np.array(g)

def area(x):
    return 0.5*x[0]*(x[2]-x[4])

if __name__ == '__main__':
    nonobtuse = sys.argv[1] == 'nonobtuse'; n = int(sys.argv[2]); tl = float(sys.argv[3])
    rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 2); res = []; t0 = time.time()
    for it in range(n):
        if time.time()-t0 > tl: break
        x0 = np.array([rng.uniform(1.8, 2.4), rng.uniform(-0.3, 2.2), rng.uniform(1.2, 2.2),
                       rng.uniform(-0.3, 2.2), -rng.uniform(1.2, 2.2), *rng.uniform(0, 2*np.pi/7, 4)])
        r = minimize(area, x0, method='SLSQP', constraints=[{'type': 'ineq', 'fun': lambda x: cons(x, nonobtuse)}],
                     bounds=[(1.7, 4), (-4, 4), (0.01, 4)]+[(-4, 4), (-4, -0.01)]+[(-1, 2)]*4,
                     options={'maxiter': 500, 'ftol': 1e-12})
        if cons(r.x, nonobtuse).min() > -1e-7: res.append((area(r.x), r.x))
    res.sort(key=lambda t: t[0])
    print('mode', sys.argv[1], 'starts', it+1, 'feasible', len(res))
    print('best totals', np.round([a for a, _ in res[:10]], 6))
    a, x = res[0]; C, T = unpack(x)
    print('best', a, 'vs 2acrit', 2*ACRIT, 'ratio', a/(2*ACRIT))
    print('x', np.round(x, 6))
    for tri in ((0, 1, 2), (0, 3, 1)):
        P = C[list(tri)]; print('  tri area', abs(tri_area(*P)), 'edges', np.round([np.linalg.norm(P[(k+1)%3]-P[k]) for k in range(3)], 4))

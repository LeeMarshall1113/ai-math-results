"""Numerical min of the area of the Delaunay star of one heptagon (d neighbours).
Relaxation: only constraints internal to the star are imposed (valid lower bound model)."""
import numpy as np, sys, time
sys.path.insert(0, 'src')
from hept import *
from scipy.optimize import minimize
ACRIT = 1.5326754446782214

def unpack(x, d):
    C = np.vstack([[0, 0], x[:2*d].reshape(d, 2)]); T = np.concatenate([[0.0], x[2*d:3*d]])
    return C, T

def incirc(a, b, c, p):
    """>0 if p strictly inside circumcircle of ccw triangle abc"""
    m = np.array([[a[0]-p[0], a[1]-p[1], (a[0]-p[0])**2+(a[1]-p[1])**2],
                  [b[0]-p[0], b[1]-p[1], (b[0]-p[0])**2+(b[1]-p[1])**2],
                  [c[0]-p[0], c[1]-p[1], (c[0]-p[0])**2+(c[1]-p[1])**2]])
    return np.linalg.det(m)

def cons(x, d):
    C, T = unpack(x, d)
    H = [verts(C[i], T[i]) for i in range(d+1)]
    g = []
    for i in range(d+1):
        for j in range(i+1, d+1):
            if np.linalg.norm(C[i]-C[j]) < 2.3 or True:
                g.append(sep(H[i], H[j]))
    for i in range(1, d+1):
        a, b = C[i], C[i % d + 1]
        ta = tri_area(C[0], a, b)
        g.append(ta)                                   # ccw ordering
        g.append(2 - circumradius(C[0], a, b))
        for j in range(1, d+1):
            if j not in (i, i % d + 1):
                g.append(-incirc(C[0], a, b, C[j])/max(ta, 1e-3)**0)
    return np.array(g)

def area(x, d):
    C, _ = unpack(x, d)
    return sum(tri_area(C[0], C[i], C[i % d + 1]) for i in range(1, d+1))

def start(rng, d, local=None):
    if local is not None:
        return local + rng.normal(0, 0.03, len(local))
    r = rng.uniform(1.85, 2.1, d); ph = np.sort(rng.uniform(0, 2*np.pi, d)) if False else \
        2*np.pi*np.arange(d)/d + rng.normal(0, 0.15, d) + rng.uniform(0, 2*np.pi)
    P = np.stack([r*np.cos(ph), r*np.sin(ph)], 1)
    return np.concatenate([P.ravel(), rng.uniform(0, 2*np.pi/7, d)])

def dl_star():
    d = kk_double_lattice(); t1, t2 = d['t1'], d['t2']
    pts = []; rf = []
    for i in range(-2, 3):
        for j in range(-2, 3):
            pts.append(i*t1+j*t2); rf.append(False); pts.append(2*d['p2']+i*t1+j*t2); rf.append(True)
    pts = np.array(pts); rf = np.array(rf)
    dist = np.linalg.norm(pts, axis=1); idx = np.argsort(dist)[1:7]
    idx = idx[np.argsort(np.arctan2(pts[idx, 1], pts[idx, 0]))]
    # reflected heptagon = rotation by pi = rotation by pi - 6*(2pi/7)... pi mod 2pi/7 = pi/7
    th = np.where(rf[idx], np.pi/7, 0.0)
    return np.concatenate([pts[idx].ravel(), th])

if __name__ == '__main__':
    d = int(sys.argv[1]); n = int(sys.argv[2]); tl = float(sys.argv[3]); mode = sys.argv[4]
    rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv) > 5 else 3)
    x_dl = dl_star() if d == 6 else None
    if x_dl is not None:
        print('DL star area', area(x_dl, 6), '6*acrit', 6*ACRIT, 'min cons', cons(x_dl, 6).min())
    res = []; t0 = time.time()
    for it in range(n):
        if time.time()-t0 > tl: break
        x0 = start(rng, d, x_dl if mode == 'local' else None)
        r = minimize(lambda x: area(x, d), x0, method='SLSQP',
                     constraints=[{'type': 'ineq', 'fun': lambda x: cons(x, d)}],
                     options={'maxiter': 800, 'ftol': 1e-13})
        if cons(r.x, d).min() > -1e-7: res.append((area(r.x, d), r.x))
    res.sort(key=lambda t: t[0])
    print('d', d, 'mode', mode, 'starts', it+1, 'feasible', len(res))
    print('best', np.round([a for a, _ in res[:12]], 6))
    if res:
        a, x = res[0]; print('best/(d*acrit)', a/(d*ACRIT), ' best/(6*acrit)', a/(6*ACRIT))
        np.save(f'runs/star_d{d}_{mode}_best.npy', x)
        C, T = unpack(x, d)
        print('nbr dists', np.round(np.linalg.norm(C[1:], axis=1), 4))
        print('tri areas', np.round([tri_area(C[0], C[i], C[i % d + 1]) for i in range(1, d+1)], 4))
        print('angles mod 2pi/7 (deg)', np.round(np.degrees(np.mod(T, 2*np.pi/7)), 3))

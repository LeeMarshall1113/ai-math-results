"""Central triangle T0 = (c0,c1,c2) receiving K obtuse triangles across its edges (each obtuse
at its apex, so the shared edge is that triangle's longest edge). Min average b over the cluster."""
import numpy as np, sys, time
CP = sys.argv[1]; KK = int(sys.argv[2]); NST = int(sys.argv[3]); TLIM = float(sys.argv[4])
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'nonobtuse', '0.3', '0', '0.05']
sys.path.insert(0, 'src')
exec(open('src/lp_transfer3.py').read().split("if __name__")[0])
from hept import verts, sep, circumradius, tri_area
from scipy.optimize import minimize
c = np.load(CP)
EDG = [(0, 1), (1, 2)]; KK = 2
NH = 3 + KK
def unpack(x):
    C = np.zeros((NH, 2)); C[1, 0] = x[0]; C[2] = x[1:3]
    for k in range(KK): C[3+k] = x[3+2*k:5+2*k]
    T = x[3+2*KK:3+2*KK+NH]
    return C, T
def tris(C, T):
    tl = [(0, 1, 2)] + [(j, i, 3+k) for k, (i, j) in enumerate(EDG)]   # ccw: neighbour across (i,j) is (j,i,apex)
    return tl
def bvals(x):
    C, T = unpack(x); tl = tris(C, T)
    Ca = np.array([[C[a], C[b], C[cc]] for a, b, cc in tl]); Ta = np.array([[T[a], T[b], T[cc]] for a, b, cc in tl])
    A, M = rows(Ca, Ta); return A + M @ c
def cons(x):
    C, T = unpack(x); H = [verts(C[i], T[i]) for i in range(NH)]; g = []
    for i in range(NH):
        for j in range(i+1, NH):
            g.append(sep(H[i], H[j]))
    for a, b, cc in tris(C, T):
        g.append(tri_area(C[a], C[b], C[cc]))                 # ccw
        g.append(2 - circumradius(C[a], C[b], C[cc]))
    for k, (i, j) in enumerate(EDG):
        ap = C[3+k]
        if k == 0: g.append(-np.dot(C[i]-ap, C[j]-ap))      # T obtuse at its apex
        # local Delaunay across (i,j): angle at c_opp + angle at apex <= pi
        opp = C[3 - i - j] if (i, j) != (2, 0) else C[1]
        def ang(p, q, r):
            u, v = q-p, r-p; return np.arccos(np.clip(u@v/np.linalg.norm(u)/np.linalg.norm(v), -1, 1))
        g.append(np.pi - ang(opp, C[i], C[j]) - ang(ap, C[i], C[j]))
    g.append(-np.dot(C[1]-C[0], C[2]-C[0]))   # T0 obtuse at c0 (its longest edge is (1,2))
    return np.array(g)
rng = np.random.default_rng(4); res = []; t0 = time.time()
for it in range(NST):
    if time.time()-t0 > TLIM: break
    d = rng.uniform(1.85, 2.6); x0 = [d, rng.uniform(0.3, d-0.3), rng.uniform(1.4, 2.4)]
    C3 = np.array([[0, 0], [d, 0], x0[1:3]])
    for k, (i, j) in enumerate(EDG):
        m = (C3[i]+C3[j])/2; e = C3[j]-C3[i]; nrm = np.array([e[1], -e[0]])/np.linalg.norm(e)
        x0 += list(m + nrm*rng.uniform(0.4, 1.0))
    x0 += list(rng.uniform(0, P7, NH))
    r = minimize(lambda x: bvals(x).mean(), np.array(x0), method='SLSQP', constraints=[{'type': 'ineq', 'fun': cons}],
                 options={'maxiter': 600, 'ftol': 1e-12})
    if cons(r.x).min() > -1e-7: res.append((bvals(r.x).mean(), r.x))
res.sort(key=lambda t: t[0])
print('chain T->T0->T2:', 'starts', it+1, 'feasible', len(res), 'best avg/acrit', np.round([a/ACRIT for a, _ in res[:6]], 5))
if res: print('b per triangle / acrit', np.round(bvals(res[0][1])/ACRIT, 4))

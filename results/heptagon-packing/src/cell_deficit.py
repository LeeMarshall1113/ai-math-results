"""Per-cell max deficit of obtuse triangles: maximise  acrit - b(T) + inc(T)  over obtuse T with
long edge in [rb_i, rb_{i+1}) and cos(gamma) in [cb_j, cb_{j+1}) (float, multistart SLSQP).
inc(T) uses the current table (D3) for the short-edge incoming charges.
Usage: cell_deficit.py <i0> <i1>  (r-bin range), writes runs/celldef_<i0>_<i1>.npy"""
import numpy as np, sys, json, time
I0, I1 = int(sys.argv[1]), int(sys.argv[2])
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'any', '0.3', '0', '0.05']
sys.path.insert(0, 'src')
exec(open('src/lp_transfer3.py').read().split("if __name__")[0])
from scipy.optimize import minimize
from hept import verts, sep
c = np.load('runs/f_w030.npy'); bins = json.load(open('runs/D3_bins.json')); rb = np.array(bins['rb']); cb = np.array(bins['cb'])
D3 = np.load('runs/D3_table.npy'); REVD = np.maximum.accumulate(D3[:, ::-1], axis=1)[:, ::-1]
RTH = rb[0]; nr, nc = D3.shape
def dhat(r, ca):
    if r < RTH or ca <= 0: return 0.0
    i = min(max(np.searchsorted(rb, r, side='right')-1, 0), nr-1)
    jj = min(max(np.searchsorted(cb[1:], -ca, side='left'), 0), nc-1)
    return REVD[i, jj]
def geom(x):   # canonical: c0=0, c1=(d,0), c2=(x,y) apex; obtuse at c2
    d, X, Y = x[0], x[1], x[2]
    C0, C1, C2 = np.zeros(2), np.array([d, 0.0]), np.array([X, Y])
    u, v = C0-C2, C1-C2; cg = u@v/np.linalg.norm(u)/np.linalg.norm(v)
    r12, r20 = np.linalg.norm(C2-C1), np.linalg.norm(C2-C0)
    ca0 = (C1-C0)@(C2-C0)/d/r20; ca1 = (C0-C1)@(C2-C1)/d/r12
    return d, cg, r12, r20, ca0, ca1
def need(x):
    d, cg, r12, r20, ca0, ca1 = geom(x)
    bb = bval(x, c)
    return ACRIT - bb + dhat(r12, ca0) + dhat(r20, ca1)
def cons(x, i, j):
    d, cg, r12, r20, ca0, ca1 = geom(x)
    H = [verts(np.zeros(2), x[3]), verts(np.array([x[0], 0]), x[4]), verts(np.array([x[1], x[2]]), x[5])]
    area = 0.5*x[0]*x[2]; Rr = x[0]*r12*r20/(4*area)
    return np.array([sep(H[0], H[1]), sep(H[0], H[2]), sep(H[1], H[2]), 2-Rr, d-rb[i], rb[i+1]-d, cg-cb[j], cb[j+1]-cg,
                     d**2-r12**2, d**2-r20**2, x[2]-0.05, x[1], x[0]/2-x[1]])
rng = np.random.default_rng(I0); out = np.full((nr, nc), -np.inf); t0 = time.time()
for i in range(I0, I1):
    for j in range(nc):
        best = -np.inf
        for s in range(8):
            d = rng.uniform(rb[i], rb[i+1]); cgt = rng.uniform(cb[j], cb[j+1])
            # apex on the circle locus for angle gamma: choose isosceles-ish start
            gam = np.arccos(cgt); h = d/2/np.tan(gam/2)
            X = d/2 + rng.normal(0, 0.15)*(s > 3); Y = max(h*(1 + rng.normal(0, 0.1)*(s > 3)), 0.05)
            x0 = np.array([d, X, Y, *rng.uniform(0, P7, 3)])
            try:
                r = minimize(lambda x: -need(x), x0, method='SLSQP', constraints=[{'type': 'ineq', 'fun': lambda x: cons(x, i, j)}],
                             options={'maxiter': 300, 'ftol': 1e-12})
            except Exception as e:
                if s == 0 and j == 0:
                    import traceback; traceback.print_exc()
                continue
            if cons(r.x, i, j).min() > -1e-8: best = max(best, need(r.x))
        out[i, j] = best
    print(f'row {i} done [{time.time()-t0:.0f}s] max {np.max(out[i]):.4f}', flush=True)
np.save(f'runs/celldef_{I0}_{I1}.npy', out)

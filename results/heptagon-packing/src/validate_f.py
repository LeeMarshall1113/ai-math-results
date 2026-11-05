"""Independent check of a fitted transfer f: fresh samples (different seed) + local searches.
Reports min b(T) over nonobtuse Delaunay triangles found, relative to acrit."""
import numpy as np, sys, time
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'nonobtuse', '0.3'] + sys.argv[1:]
sys.path.insert(0, 'src')
src = open('src/lp_transfer.py').read().split("if __name__")[0]
exec(src)
rng = np.random.default_rng(int(sys.argv[7]) if len(sys.argv) > 7 else 777)
globals()['rng'] = rng
c = np.load('runs/lp_c_DR3_DF2_nonob.npy')
t0 = time.time()
C, T = sample(50000); A, M = rows(C, T); b = A + M @ c
print('fresh samples', len(b), 'min b/acrit', b.min()/ACRIT, 'min area/acrit', A.min()/ACRIT, flush=True)
idx = np.concatenate([np.argsort(b)[:60], rng.choice(len(b), 40, replace=False)])
best = np.inf
for w in idx:
    x0 = np.array([C[w, 1, 0], C[w, 2, 0], C[w, 2, 1], *T[w]])
    r = minimize(lambda x: bval(x, c), x0, method='SLSQP', constraints=[{'type': 'ineq', 'fun': bcons}],
                 options={'maxiter': 300, 'ftol': 1e-13})
    if bcons(r.x).min() > -1e-8:
        v = bval(r.x, c)
        if v < best: best = v; bx = r.x
print('local-search min b/acrit', best/ACRIT, 'at', np.round(bx, 5), f'[{time.time()-t0:.0f}s]')

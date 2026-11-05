import numpy as np, json, sys
S = np.load(sys.argv[1]); C = S.mean(2); b = json.load(open('runs/D3_bins.json')); rb = np.array(b['rb']); cb = np.array(b['cb'])
d, x, y = C[:, 0], C[:, 1], C[:, 2]
u = np.stack([-x, -y], 1); v = np.stack([d-x, -y], 1)
cg = (u*v).sum(1)/np.linalg.norm(u, axis=1)/np.linalg.norm(v, axis=1)
r12 = np.hypot(x-d, y); r20 = np.hypot(x, y)
def near(val, grid): return np.min(np.abs(val[:, None]-grid[None, :]), 1)
print('stuck boxes', len(S), 'widths', np.quantile((S[:, :, 1]-S[:, :, 0]).max(1), [0, .5, 1]))
print('dist of d to r-grid  q:', np.quantile(near(d, rb), [0, .5, .9]))
print('dist of cos(gamma) to c-grid q:', np.quantile(near(cg, cb), [0, .5, .9]))
print('dist of short edges to RTH q:', np.quantile(np.minimum(np.abs(r12-rb[0]), np.abs(r20-rb[0])), [0, .5, .9]))
print('cos gamma q', np.quantile(cg, [0, .5, 1]), ' d q', np.quantile(d, [0, .5, 1]))

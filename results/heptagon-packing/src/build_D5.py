"""Build node table D4 from per-cell deficit estimates (local search) and sampled needs.
cell value V = max(local-search max deficit, sampled need); target T = max(V, 0) + PAD
(+ EXTRA for cells within one cell of a near-tight cell, V > -0.01);
node value = max of T over its (up to 4) adjacent cells, but 0 if all adjacent V < -ZTOL."""
import numpy as np, json, glob, hashlib
PAD, EXTRA, ZTOL = 0.03, 0.02, 0.01
D2 = np.load('runs/D2_table.npy'); need_s = np.where(D2 > 0, D2 - 0.02, -np.inf)
nr, nc = D2.shape
V = np.full((nr, nc), -np.inf)
for f in glob.glob('runs/celldef_*_*.npy'):
    a = np.load(f); V = np.maximum(V, a)
V = np.maximum(V, need_s)
near = V > -0.01
from scipy.ndimage import maximum_filter
nearx = maximum_filter(near.astype(int), size=3, mode='nearest') > 0
T = np.where(np.isfinite(V), np.maximum(V, 0) + PAD + EXTRA*nearx, 0.0)
T = np.where(np.isfinite(V) & (V > -ZTOL), T, 0.0)     # comfortably negative or empty cells: 0
Dn = np.zeros((nr+1, nc+1))
for i in range(nr):
    for j in range(nc):
        for a in (i, i+1):
            for b in (j, j+1):
                Dn[a, b] = max(Dn[a, b], T[i, j])
np.save('runs/D5_nodes.npy', Dn)
b = json.load(open('runs/D3_bins.json')); json.dump(b, open('runs/D5_bins.json', 'w'))
h = hashlib.sha256(open('runs/D5_nodes.npy', 'rb').read()).hexdigest()
print('cells with V > -inf:', int(np.isfinite(V).sum()), ' V max', np.max(V[np.isfinite(V)]), ' near-tight cells', int(near.sum()))
print('D4 node max', Dn.max(), ' nonzero nodes', int((Dn > 0).sum()), ' sha256', h)

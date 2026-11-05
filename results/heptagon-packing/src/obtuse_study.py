"""Float study for the obtuse rule: deficits of obtuse triangles vs surpluses of receivers,
binned by the shared long edge length r and the angles opposite it (gamma on the obtuse side,
alpha on the receiver side; Delaunay requires alpha + gamma <= pi)."""
import numpy as np, sys, time
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'any', '0.3', '0', '0.05']
sys.path.insert(0, 'src')
exec(open('src/lp_transfer3.py').read().split("if __name__")[0])
c = np.load('runs/f_w030.npy')
rng = np.random.default_rng(5)
def samp(n):
    out_C, out_T = [], []
    while sum(len(x) for x in out_C) < n:
        m = 40000
        C = np.zeros((m, 3, 2)); C[:, 1, 0] = rng.uniform(2.4, 4.0, m)
        C[:, 2, 0] = C[:, 1, 0]*rng.uniform(0.1, 0.9, m); C[:, 2, 1] = rng.uniform(0.3, 2.2, m)
        T = rng.uniform(0, P7, (m, 3))
        cen = C.mean(1, keepdims=True); lo = np.full(m, 0.3); hi = np.full(m, 1.0)
        for _ in range(16):
            mid = (lo+hi)/2; a, R_, dm, s, E = tri_feats(cen + (C-cen)*mid[:, None, None], T)
            good = s >= 0; hi = np.where(good, mid, hi); lo = np.where(good, lo, mid)
        Cs = cen + (C-cen)*(hi*np.where(rng.random(m) < 0.6, 1.0, rng.uniform(1, 1.3, m)))[:, None, None]
        a, R_, dm, s, E = tri_feats(Cs, T)
        ok = (a > 0) & (R_ <= 2) & (s >= -1e-9)
        out_C.append(Cs[ok]); out_T.append(T[ok])
    return np.concatenate(out_C)[:n], np.concatenate(out_T)[:n]
t0 = time.time()
C, T = samp(int(sys.argv[1]) if False else 600000)
A, M = rows(C, T); b = A + M @ c
L = np.stack([np.linalg.norm(C[:, (k+1) % 3]-C[:, k], axis=1) for k in range(3)], 1)   # edge k: v_k -> v_{k+1}
# angle opposite edge k is at vertex k+2
def ang(k):
    p = C[:, (k+2) % 3]; u = C[:, k]-p; v = C[:, (k+1) % 3]-p
    return np.arccos(np.clip((u*v).sum(1)/np.linalg.norm(u, axis=1)/np.linalg.norm(v, axis=1), -1, 1))
ANG = np.stack([ang(k) for k in range(3)], 1)
rb = np.linspace(2.548, 4.0, 30); gb = np.linspace(np.pi/2, np.pi, 19)
deficit = np.full((len(rb)-1, len(gb)-1), -np.inf); surplus = np.full((len(rb)-1, len(gb)-1), np.inf)
for k in range(3):
    obt = ANG[:, k] > np.pi/2
    ir = np.digitize(L[:, k], rb)-1; ig = np.digitize(ANG[:, k], gb)-1
    sel = obt & (ir >= 0) & (ir < len(rb)-1) & (ig >= 0) & (ig < len(gb)-1)
    np.maximum.at(deficit, (ir[sel], ig[sel]), ACRIT - b[sel])
    # receivers: edge k long, opposite angle alpha < pi/2 ; bin by alpha mapped to gamma_max = pi - alpha
    rec = (ANG[:, k] < np.pi/2) & (L[:, k] >= 2.548)
    ga = np.pi - ANG[:, k]; ig2 = np.digitize(ga, gb)-1
    sel2 = rec & (ir >= 0) & (ir < len(rb)-1) & (ig2 >= 0) & (ig2 < len(gb)-1)
    np.minimum.at(surplus, (ir[sel2], ig2[sel2]), b[sel2] - ACRIT)
np.save('runs/obtuse_deficit.npy', deficit); np.save('runs/obtuse_surplus.npy', surplus)
print('samples', len(b), 'obtuse', int((ANG.max(1) > np.pi/2).sum()), f'[{time.time()-t0:.0f}s]')
print('max obtuse deficit overall', np.max(deficit[np.isfinite(deficit)]))
# receivers with alpha: deficit allowed for gamma <= pi - alpha : cumulative max over gamma bins <= ig2
cum = np.maximum.accumulate(np.where(np.isfinite(deficit), deficit, -1), axis=1)
need = surplus - cum
fin = np.isfinite(surplus)
print('min (surplus - cumulative deficit) over bins:', np.min(need[fin]) if fin.any() else None)
i = np.unravel_index(np.argmin(np.where(fin, need, np.inf)), need.shape)
print('worst bin r', rb[i[0]], 'gamma_max', gb[i[1]+1], 'surplus', surplus[i], 'cum deficit', cum[i])

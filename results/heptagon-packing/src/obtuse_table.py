"""Build the obtuse-transfer table D(r, gamma) (float, from samples), with chain recursion, and
check the single-triangle conditions on the samples:
  obtuse T:      acrit - b(T) + sum_in Dhat <= D(bin(T))
  non-obtuse T:  b(T) - sum_in Dhat >= acrit
Dhat(r, alpha) = max_{gamma <= pi - alpha} D(r, gamma) (Delaunay: alpha + gamma <= pi)."""
import numpy as np, sys, time, json
PAD = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'any', '0.3', '0', '0.05']
sys.path.insert(0, 'src')
exec(open('src/obtuse_study.py').read().split("t0 = time.time()")[0])
RTH = 2*0.9009688679024191*np.sqrt(2)
rb = np.linspace(RTH, 4.0, 30); gb = np.linspace(np.pi/2, np.pi, 19)
C, T = samp(800000); A, M = rows(C, T); b = A + M @ c
L = np.stack([np.linalg.norm(C[:, (k+1) % 3]-C[:, k], axis=1) for k in range(3)], 1)
def ang(k):
    p = C[:, (k+2) % 3]; u = C[:, k]-p; v = C[:, (k+1) % 3]-p
    return np.arccos(np.clip((u*v).sum(1)/np.linalg.norm(u, axis=1)/np.linalg.norm(v, axis=1), -1, 1))
ANG = np.stack([ang(k) for k in range(3)], 1)
kmax = np.argmax(ANG, 1); obt = ANG.max(1) > np.pi/2
nr, ng = len(rb)-1, len(gb)-1
D = np.zeros((nr, ng))
def Dhat(r, alpha):
    ir = np.clip(np.digitize(r, rb)-1, 0, nr-1); gmax = np.pi - alpha
    ig = np.clip(np.digitize(gmax, gb)-1, -1, ng-1)       # bins with gamma_lo <= gmax
    cum = np.maximum.accumulate(D, axis=1)
    return np.where((r >= RTH) & (alpha < np.pi/2) & (ig >= 0), cum[ir, np.maximum(ig, 0)], 0.0)
def incoming(skip_long):
    tot = np.zeros(len(b))
    for k in range(3):
        use = ~(skip_long & (kmax == k))
        tot += np.where(use, Dhat(L[:, k], ANG[:, k]), 0.0)
    return tot
# recursion by increasing r-bin (incoming edges of an obtuse triangle are shorter than its long edge)
rl = L[np.arange(len(b)), kmax]; gl = ANG.max(1)
irT = np.digitize(rl, rb)-1; igT = np.digitize(gl, gb)-1
for it in range(3):
    inc = incoming(obt)
    need = ACRIT - b + inc
    for i in range(nr):
        sel = obt & (irT == i) & (igT >= 0) & (igT < ng)
        if sel.any():
            Dn = np.full(ng, -np.inf); np.maximum.at(Dn, igT[sel], need[sel])
            D[i] = np.maximum(0.0, np.where(np.isfinite(Dn), Dn + PAD, 0.0))
    # smooth: monotone envelope over neighbouring bins
    D = np.maximum(D, np.maximum(np.roll(D, 1, 0)*0, 0))
inc = incoming(obt)
lhs_obt = (D[np.clip(irT, 0, nr-1), np.clip(igT, 0, ng-1)] - (ACRIT - b + inc))[obt]
lhs_non = (b - inc - ACRIT)[~obt]
print('samples', len(b), 'obtuse', int(obt.sum()), 'max D', D.max())
print('obtuse condition min slack', lhs_obt.min(), ' nonobtuse condition min slack', lhs_non.min())
print('nonobtuse with incoming>0: count', int((inc[~obt] > 0).sum()), 'min slack among them', (lhs_non[inc[~obt] > 0]).min())
np.save('runs/D_table.npy', D); json.dump(dict(rb=rb.tolist(), gb=gb.tolist(), pad=PAD), open('runs/D_table_bins.json', 'w'))
print(np.round(D[:, ::3], 3))

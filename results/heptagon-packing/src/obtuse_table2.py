"""Obtuse-transfer table D over (r, c) bins, c = cos(gamma) in [-1, 0) (gamma = obtuse angle,
r = long edge). Receiver across an edge (length r, opposite angle alpha): compatible obtuse
neighbours have cos gamma >= -cos alpha (Delaunay alpha + gamma <= pi).
  Dhat(r, ca) = max over c-bins with c_hi >= -ca of D(r-bin, c-bin)
Recursion over increasing r (incoming edges of an obtuse triangle are shorter than its long edge).
Float, from samples; the rigorous B&B then verifies the conditions."""
import numpy as np, sys, json
PAD = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 800000
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'any', '0.3', '0', '0.05']
sys.path.insert(0, 'src')
exec(open('src/obtuse_study.py').read().split("t0 = time.time()")[0])
RTH = 2*0.9009688679024191*np.sqrt(2)
rb = np.linspace(RTH, 4.0, 30); cb = np.linspace(-1.0, 0.0, 21)
nr, nc = len(rb)-1, len(cb)-1
C, T = samp(NS); A, M = rows(C, T); b = A + M @ c
L = np.stack([np.linalg.norm(C[:, (k+1) % 3]-C[:, k], axis=1) for k in range(3)], 1)
def cosang(k):
    p = C[:, (k+2) % 3]; u = C[:, k]-p; v = C[:, (k+1) % 3]-p
    return (u*v).sum(1)/np.linalg.norm(u, axis=1)/np.linalg.norm(v, axis=1)
CA = np.stack([cosang(k) for k in range(3)], 1)       # cos of angle opposite edge k
kmax = np.argmax(L, 1); cg = CA[np.arange(len(b)), kmax]; obt = cg < 0
D = np.zeros((nr, nc))
def Dhat(r, ca):
    ir = np.clip(np.digitize(r, rb)-1, 0, nr-1)
    # bins j with cb[j+1] >= -ca  -> j >= first index jj ; take max over j >= jj
    rev = np.maximum.accumulate(D[:, ::-1], axis=1)[:, ::-1]      # rev[i, j] = max_{j' >= j} D[i, j']
    jj = np.clip(np.searchsorted(cb[1:], -ca, side='left'), 0, nc-1)
    val = rev[ir, jj]
    return np.where((r >= RTH) & (ca > 0) & (-ca <= 0), val, 0.0)
rl = L[np.arange(len(b)), kmax]
irT = np.clip(np.digitize(rl, rb)-1, 0, nr-1); icT = np.clip(np.digitize(cg, cb)-1, 0, nc-1)
for it in range(4):
    inc = np.zeros(len(b))
    for k in range(3):
        inc += np.where(obt & (kmax == k), 0.0, Dhat(L[:, k], CA[:, k]))
    need = ACRIT - b + inc
    Dn = np.full((nr, nc), -np.inf); np.maximum.at(Dn, (irT[obt], icT[obt]), need[obt])
    D = np.maximum(0.0, np.where(np.isfinite(Dn), Dn + PAD, 0.0))
inc = np.zeros(len(b))
for k in range(3):
    inc += np.where(obt & (kmax == k), 0.0, Dhat(L[:, k], CA[:, k]))
slack_obt = (D[irT, icT] - (ACRIT - b + inc))[obt]; slack_non = (b - ACRIT - inc)[~obt]
print('samples', len(b), 'obtuse', int(obt.sum()), 'max D', round(D.max(), 4), 'nonzero bins', int((D > 0).sum()))
print('obtuse min slack', slack_obt.min(), ' non-obtuse min slack', slack_non.min(),
      ' non-obtuse with incoming>0 min slack', slack_non[inc[~obt] > 0].min() if (inc[~obt] > 0).any() else None)
np.save('runs/D2_table.npy', D); json.dump(dict(rb=rb.tolist(), cb=cb.tolist(), pad=PAD, RTH=RTH), open('runs/D2_bins.json', 'w'))

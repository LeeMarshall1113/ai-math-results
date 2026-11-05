import sys, numpy as np
sys.path.insert(0, 'src')
import bnb2 as B2
B2.OBT = B2.ObtTable4('runs/D4_nodes.npy', 'runs/D4_bins.json')
F = B2.FTrans('runs/f_w030.npy', 3, 2, mu=[-1.1684855e-12, -1.1684853e-12])
ACR = 1.5326754446782214
S = np.load(sys.argv[1]); idx = np.random.default_rng(0).choice(len(S), min(4000, len(S)), replace=False); S = S[idx]
C = S.mean(2); P = np.stack([C, C], 2); w = (S[:, :, 1]-S[:, :, 0])
inf, _ = B2.infeasible(P); ex = B2.obtuse_extra(P)
b = B2.b_A(F, [B2.R(C[:, q]) for q in range(6)], 1).v.lo
sl = b - ACR - ex; fe = ~inf
print('open sample', len(S), 'width q', np.quantile(w.max(1), [0, .5, .9]), 'feasible centres', fe.sum())
print('slack at feasible centres q', np.quantile(sl[fe], [0, .05, .5]))
d, x, y = C[:, 0], C[:, 1], C[:, 2]
u = np.stack([-x, -y], 1); v = np.stack([d-x, -y], 1); cg = (u*v).sum(1)/np.linalg.norm(u, axis=1)/np.linalg.norm(v, axis=1)
print('d q', np.round(np.quantile(d[fe], [0, .5, 1]), 3), 'cos g q', np.round(np.quantile(cg[fe], [0, .5, 1]), 3), 'y q', np.round(np.quantile(y[fe], [0, .5, 1]), 3))
print('b-acrit q', np.round(np.quantile((b-ACR)[fe], [0, .5]), 4), 'extra q', np.round(np.quantile(ex[fe], [0, .5, 1]), 4))
exB = B2.obtuse_extra(S); print('extra(box) - extra(centre) q', np.round(np.quantile((exB-ex)[fe], [0, .5, .9, 1]), 4))
lb1 = B2.lb_first(F, S[fe]); print('first-order lb - target q', np.round(np.quantile(lb1 - (ACR + exB[fe]), [0, .5, .9]), 4))

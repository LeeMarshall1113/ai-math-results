import sys, numpy as np
sys.path.insert(0, 'src')
import bnb2 as B2
B2.OBT = B2.ObtTable4('runs/D4_nodes.npy', 'runs/D4_bins.json')
F = B2.FTrans('runs/f_w030.npy', 3, 2, mu=[-1.1684855e-12, -1.1684853e-12])
ACR = 1.5326754446782214
S = np.load(sys.argv[1]); w = (S[:, :, 1]-S[:, :, 0]).max(1)
print('open', len(S), 'width q', np.quantile(w, [0, .1, .5]))
tiny = S[w < 1e-5]; print('tiny boxes (<1e-5)', len(tiny))
T = tiny[np.random.default_rng(0).choice(len(tiny), min(3000, len(tiny)), replace=False)] if len(tiny) else S[:3000]
C = T.mean(2); P = np.stack([C, C], 2)
inf, _ = B2.infeasible(P); ex = B2.obtuse_extra(P)
b = B2.b_A(F, [B2.R(C[:, q]) for q in range(6)], 1).v.lo
sl = b - ACR - ex; fe = ~inf
print('feasible centres', fe.sum(), 'slack q', np.quantile(sl[fe], [0, .1, .5]) if fe.any() else None)
d, x, y = C[:, 0], C[:, 1], C[:, 2]
u = np.stack([-x, -y], 1); v = np.stack([d-x, -y], 1); cg = (u*v).sum(1)/np.linalg.norm(u, axis=1)/np.linalg.norm(v, axis=1)
r12 = np.hypot(x-d, y); r20 = np.hypot(x, y)
for nm, a in (('d', d), ('cos g', cg), ('r12', r12), ('r20', r20), ('b-acrit', b-ACR), ('extra', ex)):
    print(nm, np.round(np.quantile(a[fe], [0, .5, 1]), 4))
j = np.argsort(sl[fe])[:3]; print('worst centres', np.round(C[fe][j], 5), sl[fe][j])

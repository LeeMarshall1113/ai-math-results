import sys, numpy as np, json
sys.path.insert(0, 'src')
import bnb2 as B2
B2.OBT = B2.ObtTable('runs/D3_table.npy', 'runs/D3_bins.json')
F = B2.FTrans('runs/f_w030.npy', 3, 2, mu=[-1.1684855e-12, -1.1684853e-12])
ACR = 1.5326754446782214
S = np.concatenate([np.load(f'runs/stuck913_k{k}.npy') for k in (9, 10, 11, 12, 13)])
print('stuck boxes sampled', len(S), 'max width q', np.quantile((S[:, :, 1]-S[:, :, 0]).max(1), [0, .5, 1]))
C = S.mean(2); P = np.stack([C, C], 2)
inf, _ = B2.infeasible(P); ex = B2.obtuse_extra(P)
b = B2.b_A(F, [B2.R(C[:, q]) for q in range(6)], 1).v.lo
sl = b - ACR - ex
print('centres feasible', (~inf).sum(), ' slack at feasible centres q[0,.1,.5]:', np.quantile(sl[~inf], [0, .1, .5]) if (~inf).any() else None)
print('centre ranges lo', np.round(C.min(0), 3), 'hi', np.round(C.max(0), 3))
d, x, y = C[:, 0], C[:, 1], C[:, 2]
u = np.stack([-x, -y], 1); v = np.stack([d-x, -y], 1)
cg = (u*v).sum(1)/np.linalg.norm(u, axis=1)/np.linalg.norm(v, axis=1)
r12 = np.hypot(x-d, y); r20 = np.hypot(x, y)
print('cos gamma q', np.round(np.quantile(cg, [0, .5, 1]), 4), ' d q', np.round(np.quantile(d, [0, .5, 1]), 4))
print('short edges q', np.round(np.quantile(r12, [0, .5, 1]), 4), np.round(np.quantile(r20, [0, .5, 1]), 4))
print('b - acrit q', np.round(np.quantile(b-ACR, [0, .5, 1]), 4), ' extra q', np.round(np.quantile(ex, [0, .5, 1]), 4))
# which interval term is undecidable: compare extra on the box vs at centre
exB = B2.obtuse_extra(S)
print('extra over box - extra at centre q', np.round(np.quantile(exB-ex, [0, .5, 1]), 4))
# R near 2?
area = 0.5*d*y; Rr = d*r12*r20/(4*area); print('circumradius q', np.round(np.quantile(Rr, [0, .5, 1]), 5))
# any configuration-level degeneracy: angles at c0, c1
ca0 = (x*d)/(r20*d); ca1 = ((d-x)*d)/(r12*d); print('cos angles at c0,c1 q', np.round(np.quantile(ca0, [0, .5, 1]), 4), np.round(np.quantile(ca1, [0, .5, 1]), 4))

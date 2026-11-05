"""Float check of the constant obtuse transfer: every obtuse triangle sends D to its long-edge
neighbour. Conditions (K = number of edges with length >= 2u sqrt2, counted as potential
incoming; an obtuse triangle's long edge is its outgoing edge):
  obtuse:     b - acrit + 2D - D*K >= 0      (K includes the long edge)
  non-obtuse: b - acrit - D*K >= 0"""
import numpy as np, sys
DV = [float(x) for x in sys.argv[1:]] or [0.15, 0.17, 0.2]
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'any', '0.3', '0', '0.05']
sys.path.insert(0, 'src')
exec(open('src/obtuse_study.py').read().split("t0 = time.time()")[0])
RTH = 2*0.9009688679024191*np.sqrt(2)
C, T = samp(800000); A, M = rows(C, T); b = A + M @ c
L = np.stack([np.linalg.norm(C[:, (k+1) % 3]-C[:, k], axis=1) for k in range(3)], 1)
d2 = np.sort(L, 1); obt = d2[:, 2]**2 > d2[:, 0]**2 + d2[:, 1]**2
K = (L >= RTH).sum(1)
for D in DV:
    lhs = np.where(obt, b - ACRIT + 2*D - D*K, b - ACRIT - D*K)
    w = np.argmin(lhs)
    print(f'D={D}: min slack {lhs.min():.4f} (obtuse min {lhs[obt].min():.4f}, nonobtuse K>=1 min {lhs[~obt & (K>=1)].min():.4f}); worst K={K[w]} obtuse={obt[w]} edges={np.round(L[w],3)}')

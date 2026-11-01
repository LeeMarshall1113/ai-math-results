"""Vectorised exact Q_P H-rep, equal (as a set of gcd-normalised integer rows) to hw.QP_hrep(a, d)."""
import itertools, math
from fractions import Fraction as F
import numpy as np

TRIPLES = list(itertools.combinations(range(6), 3))

def QP_hrep_fast(a, d):
    a = [F(x) for x in a]; d = [F(x) for x in d]
    S = 1
    for x in a + d: S = S * x.denominator // math.gcd(S, x.denominator)
    A = [int(x * S) for x in a]; D = [int(x * S) for x in d]
    # 64 vertices of the box P; for each, the 6 points of V_v (scaled by S)
    pts = []
    for e in itertools.product((0, 1), repeat=6):
        v = [A[i] + e[i] * D[i] for i in range(6)]
        p1, p2, p3, p4, p5, p6 = v
        pts.append([(0, p1, p2), (S, p1, p2), (p3, 0, p4), (p3, S, p4), (p5, p6, 0), (p5, p6, S)])
    P = np.array(pts, dtype=np.int64)                      # 64 x 6 x 3
    T = np.array(TRIPLES)
    X, Y, Z = P[:, T[:, 0]], P[:, T[:, 1]], P[:, T[:, 2]]  # 64 x 20 x 3
    U = Y - X; V = Z - X
    N = np.cross(U, V)                                     # 64 x 20 x 3
    c = -(N * X).sum(-1)                                   # 64 x 20
    vals = c[:, :, None] + np.einsum('vtk,vpk->vtp', N, P) # 64 x 20 x 6
    nz = np.any(N != 0, axis=-1)
    pos = (vals >= 0).all(-1); neg = (vals <= 0).all(-1)
    out = set()
    for vi, ti in zip(*np.nonzero(nz & (pos | neg))):
        n = [int(x) for x in N[vi, ti]]; cc = int(c[vi, ti])
        if not pos[vi, ti]: n = [-x for x in n]; cc = -cc
        # back to original coords: c + n.(S x) >= 0  ->  (c, S n)
        row = [cc, S * n[0], S * n[1], S * n[2]]
        g = 0
        for x in row: g = math.gcd(g, abs(x))
        out.add(tuple(x // g for x in row))
    out.add((0, 1, 0, 0))
    return sorted(out)

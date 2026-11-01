"""qpfast.QP_hrep_fast with the per-facet Python loop (sign flip, scale, gcd, set) done in numpy int64.
Returns exactly the same sorted list of gcd-normalised integer rows (checked in test_equiv.py).
int64 is safe: coordinates are integers <= S (S | 640 here, asserted <= 2**15), so |n| <= 2 S^2, |c| <= 6 S^3, and the rows
(c, S n) are < 2**52."""
import itertools, math
from fractions import Fraction as F
import numpy as np
from qpfast import TRIPLES

_T = np.array(TRIPLES)
_E = np.array(list(itertools.product((0, 1), repeat=6)), dtype=np.int64)   # same vertex order as qpfast

def QP_hrep_fast(a, d):
    a = [F(x) for x in a]; d = [F(x) for x in d]
    S = 1
    for x in a + d: S = S * x.denominator // math.gcd(S, x.denominator)
    assert S <= 2 ** 15
    A = np.array([int(x * S) for x in a], dtype=np.int64); D = np.array([int(x * S) for x in d], dtype=np.int64)
    v = A[None, :] + _E * D[None, :]                       # 64 x 6
    p1, p2, p3, p4, p5, p6 = v.T
    z = np.zeros(64, dtype=np.int64); s = np.full(64, S, dtype=np.int64)
    P = np.stack([np.stack([z, p1, p2], -1), np.stack([s, p1, p2], -1), np.stack([p3, z, p4], -1),
                  np.stack([p3, s, p4], -1), np.stack([p5, p6, z], -1), np.stack([p5, p6, s], -1)], 1)  # 64 x 6 x 3
    X, Y, Z = P[:, _T[:, 0]], P[:, _T[:, 1]], P[:, _T[:, 2]]
    N = np.cross(Y - X, Z - X)
    c = -(N * X).sum(-1)
    vals = c[:, :, None] + np.einsum('vtk,vpk->vtp', N, P)
    nz = np.any(N != 0, axis=-1)
    pos = (vals >= 0).all(-1); neg = (vals <= 0).all(-1)
    sel = nz & (pos | neg)
    n = N[sel]; cc = c[sel]; flip = ~pos[sel]
    n = np.where(flip[:, None], -n, n); cc = np.where(flip, -cc, cc)
    rows = np.concatenate([cc[:, None], S * n], axis=1)
    g = np.gcd.reduce(np.abs(rows), axis=1)
    rows = rows // g[:, None]
    rows = np.vstack([rows, np.array([[0, 1, 0, 0]], dtype=np.int64)])
    rows = np.unique(rows, axis=0)                          # lexicographic = sorted() of tuples
    return [tuple(int(x) for x in r) for r in rows]

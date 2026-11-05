"""Vectorised heptagon-triangle geometry for batches."""
import numpy as np
from hept import K, U
P7 = 2*np.pi/7
EN = None
def hverts(c, t):
    """c (n,2), t (n,) -> (n,7,2)"""
    a = t[:, None] + 2*np.pi*np.arange(7)[None, :]/7
    return np.stack([np.cos(a), np.sin(a)], -1) + c[:, None, :]

def bsep(A, B):
    """batched SAT separation, A,B (n,7,2)"""
    best = np.full(A.shape[0], -np.inf)
    for X, Y in ((A, B), (B, A)):
        E = np.roll(X, -1, 1) - X
        Nn = np.stack([E[..., 1], -E[..., 0]], -1); Nn /= np.linalg.norm(Nn, axis=-1, keepdims=True)
        px = np.einsum('nvd,ned->nve', X, Nn).max(1)
        qy = np.einsum('nvd,ned->nve', Y, Nn).min(1)
        best = np.maximum(best, (qy - px).max(1))
    return best

def tri_feats(C, T):
    """C (n,3,2), T (n,3). Returns area (ccw positive), circumradius, min dot (nonobtuse>=0),
    min sep, and pair configs (r,a,b) for oriented edges 0->1,1->2,2->0."""
    a = 0.5*((C[:, 1, 0]-C[:, 0, 0])*(C[:, 2, 1]-C[:, 0, 1]) - (C[:, 1, 1]-C[:, 0, 1])*(C[:, 2, 0]-C[:, 0, 0]))
    L = np.stack([np.linalg.norm(C[:, (k+1) % 3]-C[:, k], axis=1) for k in range(3)], 1)
    R = L.prod(1)/(4*np.abs(a)+1e-300)
    dots = np.stack([np.einsum('nd,nd->n', C[:, (k+1) % 3]-C[:, k], C[:, (k+2) % 3]-C[:, k]) for k in range(3)], 1)
    H = [hverts(C[:, k], T[:, k]) for k in range(3)]
    s = np.minimum(np.minimum(bsep(H[0], H[1]), bsep(H[1], H[2])), bsep(H[2], H[0]))
    E = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        v = C[:, j]-C[:, i]; phi = np.arctan2(v[:, 1], v[:, 0])
        E.append(np.stack([np.linalg.norm(v, axis=1), np.mod(T[:, i]-phi, P7), np.mod(T[:, j]-phi, P7)], 1))
    return a, R, dots.min(1), s, E

"""Numerical sanity probe for the SAT core GP+PL+AD: look for 10 points with the 10 faces of
{4,4}_(3,1) coplanar and, for every edge-adjacent pair (F,G), the two vertices of G off the
common edge strictly on the same side of plane(F).  The SAT result says this is impossible;
if this probe ever succeeds, the encoding has a bug.  usage: numtest_ad.py seed iters [map]"""
import os, sys, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import least_squares

def run(F, seed, iters, DEL=0.05):
    N = len({u for f in F for u in f}); M = len(F)
    Fi = np.array(F)
    adj = []
    for i in range(M):
        for j in range(M):
            s = set(F[i]) & set(F[j])
            if i != j and len(s) == 2:
                x, y = [u for u in F[j] if u not in s]; adj.append((i, x, y))
    adj = np.array(adj)
    def resid(z):
        X = z[:3*N].reshape(N, 3); P = z[3*N:].reshape(M, 4)
        n = P[:, :3]; c = P[:, 3]
        r1 = (np.einsum('fkd,fd->fk', X[Fi], n) - c[:, None]).ravel() * 10
        r2 = (n * n).sum(1) - 1
        sx = (X[adj[:, 1]] * n[adj[:, 0]]).sum(1) - c[adj[:, 0]]
        sy = (X[adj[:, 2]] * n[adj[:, 0]]).sum(1) - c[adj[:, 0]]
        r3 = np.minimum(sx * sy - DEL ** 2, 0) * 10
        r4 = np.minimum(1 - np.abs(X).ravel(), 0)
        return np.concatenate([r1, r2, r3, r4])
    rng = np.random.default_rng(seed); best = 1e9
    for it in range(iters):
        X = rng.uniform(-1, 1, size=(N, 3)); P = np.zeros((M, 4))
        for k, f in enumerate(F):
            Y = X[list(f)]; cc = Y.mean(0); nn = np.linalg.svd(Y - cc)[2][2]; P[k] = [*nn, nn @ cc]
        sol = least_squares(resid, np.concatenate([X.ravel(), P.ravel()]), max_nfev=4000)
        best = min(best, 2 * sol.cost)
        if 2 * sol.cost < 1e-16: return it, sol.x
    return None, best

if __name__ == "__main__":
    seed, iters = int(sys.argv[1]), int(sys.argv[2])
    which = sys.argv[3] if len(sys.argv) > 3 else "M10"
    if which == "M10": F = [tuple((i + d) % 10 for d in (0, 1, 8, 7)) for i in range(10)]
    elif which == "T9": F = [[0,1,4,3],[1,2,5,4],[2,0,3,5],[3,4,7,6],[4,5,8,7],[5,3,6,8],[6,7,1,0],[7,8,2,1],[8,6,0,2]]
    elif which == "Q12":  # {4,4}_(4,0)-like 4x3 grid, lattice <(4,0),(0,3)>
        idx = lambda x, y: (x % 4) + 4 * (y % 3)
        F = [(idx(x, y), idx(x+1, y), idx(x+1, y+1), idx(x, y+1)) for y in range(3) for x in range(4)]
    t = time.time(); it, res = run(F, seed, iters)
    print(which, "FOUND at iter" if it is not None else "not found, best cost", it if it is not None else res, f"{time.time()-t:.0f}s")

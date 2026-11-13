"""indep_homology.py -- a SECOND, independent homology computation for the Karp n=10 replay (Z5 fixed complexes).

Does not import replay_lib / replay_smith.  Faces are bitmasks over the k orbit indices; the empty face is 0 (so the homology is REDUCED).
Three routes to the same numbers, all different from replay_smith.homology_by_size (sparse F_p elimination, pivot = smallest column index):

  (1) Smith normal form over Z (unit-pivot elimination on sparse rows, then dense SNF of the remainder) -> the full integral homology
      (Betti numbers and torsion), and from it dim H_s(;F_p) by the universal coefficient theorem:
          dim H_s(C (x) F_p) = betti_s + #{torsion factors of H_s divisible by p} + #{torsion factors of H_{s-1} divisible by p};
  (2) rank over F_p of every boundary matrix by DENSE numpy Gaussian elimination (different data structure, pivot = first nonzero in the column of
      largest index: a different elimination order);
  (3) the Euler characteristic (alternating face count), a consistency check of both.

Relative homology H_*(K, L) of a pair is computed as the homology of the quotient chain complex (generators = faces of K not in L, boundary terms
falling into L dropped); the replay's "relative complex of the X-faces" is exactly this with L = the subcomplex of disconnected unions.
"""
from collections import defaultdict
import numpy as np


def popcount(x):
    return bin(x).count("1")


def boundary_rows(src, tgt_index):
    """sparse rows (one per face of `src`): {column index of target face: coefficient}; terms whose target is not a generator are dropped"""
    rows = []
    for t in src:
        r = {}
        pos = 0
        for b in range(32):
            if (t >> b) & 1:
                j = tgt_index.get(t ^ (1 << b))
                if j is not None:
                    r[j] = -1 if pos % 2 else 1
                pos += 1
        rows.append(r)
    return rows


def _dense_snf(M):
    """invariant factors of a small dense integer matrix (list of lists), by repeated minimal-pivot reduction"""
    M = [row[:] for row in M if any(row)]
    invs = []
    while M:
        # smallest nonzero |entry|
        best = None
        for i, row in enumerate(M):
            for j, v in enumerate(row):
                if v and (best is None or abs(v) < best[0]):
                    best = (abs(v), i, j)
        _, pi, pj = best
        while True:
            p = M[pi][pj]
            done = True
            for i in range(len(M)):
                if i != pi and M[i][pj]:
                    q = M[i][pj] // p
                    M[i] = [a - q * b for a, b in zip(M[i], M[pi])]
                    if M[i][pj]:
                        done = False
            for j in range(len(M[0])):
                if j != pj and M[pi][j]:
                    q = M[pi][j] // p
                    for row in M:
                        row[j] -= q * row[pj]
                    if M[pi][j]:
                        done = False
            if done:
                # make sure the pivot divides everything else
                bad = None
                for i in range(len(M)):
                    if i != pi:
                        for j in range(len(M[0])):
                            if j != pj and M[i][j] % p:
                                bad = (i, j)
                                break
                    if bad:
                        break
                if bad is None:
                    break
                M[pi] = [a + b for a, b in zip(M[pi], M[bad[0]])]       # bring a non-multiple into the pivot row and reduce again
            # new smaller pivot in row pi or column pj
            cand = [(abs(M[pi][j]), pi, j) for j in range(len(M[0])) if M[pi][j]] + [(abs(M[i][pj]), i, pj) for i in range(len(M)) if M[i][pj]]
            _, pi, pj = min(cand)
        invs.append(abs(M[pi][pj]))
        M = [[v for j, v in enumerate(row) if j != pj] for i, row in enumerate(M) if i != pi]
        M = [row for row in M if any(row)]
    return invs


def invariant_factors(rows):
    """nonzero invariant factors (with multiplicity) of the integer matrix given as sparse rows"""
    R = {i: dict(r) for i, r in enumerate(rows) if r}
    cols = defaultdict(set)
    for i, r in R.items():
        for c in r:
            cols[c].add(i)
    units = 0
    while True:
        found, best = None, None
        for i, r in R.items():
            for c, v in r.items():
                if v == 1 or v == -1:
                    cost = (len(r) - 1) * (len(cols[c]) - 1)
                    if best is None or cost < best:
                        best, found = cost, (i, c)
                        if cost == 0:
                            break
            if best == 0:
                break
        if found is None:
            break
        i, c = found
        pr = R.pop(i)
        pv = pr[c]
        for cc in pr:
            cols[cc].discard(i)
        for j in list(cols[c]):
            r = R[j]
            f = r[c] * pv                               # pv = +-1, so pv^{-1} = pv
            for cc, v in pr.items():
                nv = r.get(cc, 0) - f * v
                if nv:
                    r[cc] = nv
                    cols[cc].add(j)
                elif cc in r:
                    del r[cc]
                    cols[cc].discard(j)
            if not r:
                del R[j]
        units += 1
    invs = [1] * units
    if R:
        cs = sorted({c for r in R.values() for c in r})
        ci = {c: k for k, c in enumerate(cs)}
        dense = []
        for r in R.values():
            row = [0] * len(cs)
            for c, v in r.items():
                row[ci[c]] = v
            dense.append(row)
        invs += _dense_snf(dense)
    return invs


def rank_mod_p_dense(rows, ncols, p):
    """rank over F_p of the sparse-row matrix, by dense numpy elimination processing the columns from the LAST to the first"""
    if not rows or ncols == 0:
        return 0
    A = np.zeros((len(rows), ncols), dtype=np.int64)
    for i, r in enumerate(rows):
        for j, v in r.items():
            A[i, j] = v % p
    rank = 0
    r0 = 0
    for c in range(ncols - 1, -1, -1):
        nz = np.nonzero(A[r0:, c])[0]
        if len(nz) == 0:
            continue
        piv = r0 + nz[0]
        if piv != r0:
            A[[r0, piv]] = A[[piv, r0]]
        inv = pow(int(A[r0, c]), -1, p)
        A[r0] = (A[r0] * inv) % p
        others = np.nonzero(A[:, c])[0]
        for i in others:
            if i != r0:
                A[i] = (A[i] - A[i, c] * A[r0]) % p
        r0 += 1
        rank += 1
        if r0 == A.shape[0]:
            break
    return rank


def homology(K, L=frozenset(), primes=(2, 3, 5), want_z=True):
    """(reduced) homology of the pair (K, L): K, L are sets of bitmask faces (K contains 0, the empty face, if reduced homology is wanted).
    Returns dict with: closed (K is closed under taking faces; only meaningful for L empty), n[s], betti[s], torsion[s] (list), fp[p][s] (via SNF/UCT),
    fp_dense[p][s] (via dense rank mod p) and euler."""
    K = set(K)
    L = set(L)
    assert L <= K
    gens = defaultdict(list)
    for t in sorted(K - L):
        gens[popcount(t)].append(t)
    sizes = sorted(gens)
    smax = sizes[-1] if sizes else 0
    closed = True
    for t in K:
        for b in range(32):
            if (t >> b) & 1 and (t ^ (1 << b)) not in K:
                closed = False
    idx = {s: {t: i for i, t in enumerate(gens[s])} for s in sizes}
    rows = {}
    for s in sizes:
        if s >= 1 and (s - 1) in idx:
            rows[s] = boundary_rows(gens[s], idx[s - 1])
        else:
            rows[s] = []
    n = {s: len(gens[s]) for s in range(0, smax + 2)}
    out = dict(closed=closed, n=n, euler=sum((-1) ** s * n[s] for s in n))
    invs = {s: (invariant_factors(rows[s]) if want_z else []) for s in rows}
    for s in range(0, smax + 3):
        invs.setdefault(s, [])
    rk = {s: len(invs[s]) for s in range(0, smax + 3)}
    if want_z:
        out["betti"] = {s: n[s] - rk[s] - rk[s + 1] for s in range(0, smax + 1)}
        out["torsion"] = {s: [d for d in invs[s + 1] if d > 1] for s in range(0, smax + 1)}
        out["fp"] = {}
        for p in primes:
            out["fp"][p] = {s: out["betti"][s] + sum(1 for d in out["torsion"][s] if d % p == 0) + (sum(1 for d in out["torsion"][s - 1] if d % p == 0) if s >= 1 else 0)
                            for s in range(0, smax + 1)}
    out["fp_dense"] = {}
    for p in primes:
        rp = {s: (rank_mod_p_dense(rows[s], n[s - 1], p) if s in rows and s >= 1 else 0) for s in range(0, smax + 3)}
        out["fp_dense"][p] = {s: n[s] - rp[s] - rp.get(s + 1, 0) for s in range(0, smax + 1)}
    if want_z:
        for p in primes:
            assert out["fp"][p] == out["fp_dense"][p], ("SNF/UCT and dense F_p ranks disagree", p, out["fp"][p], out["fp_dense"][p])
    return out

"""Second implementation (referee), level one (written from the definitions; shares no code
with lr.py / tight6.py / tight7.py).

Objects.  k speeds, prime p > k+1, H = (p-1)/2, m = floor((p-1)/(k+1)).
A nonzero class c mod p (folded to 1..H by c ~ -c) COVERS a level-one time s in 1..H
iff dist(s*c mod p, 0) <= m  (equivalently (k+1)*dist <= p; equality is impossible
because (k+1) does not divide p, so open and closed covering coincide here).
A k-multiset M of classes is LEVEL-ONE IMPROPER iff every s in 1..H is covered.

Normalisation lemma (Lemma N in PROOF.md).  Let n_s = number of elements
of M (with multiplicity) covering s.  Sum_s n_s = k*m and every n_s >= 1, so the number
P of times with n_s = 1 satisfies P >= 2H - k*m.  Hence some element e of M is the
unique coverer of at least q = ceil((2H - k m)/k) times (q >= 1 always).  Multiplying
M by e^{-1} keeps improperness and maps e to the class 1.  So every unit orbit of
improper multisets has a member that contains the class 1 as a distinguished element
which is the unique coverer of >= q times.

Search.  DFS on the uncovered time with fewest live candidates, with sibling exclusion.
For every such normalised M there is a root-to-leaf path whose chosen set S is a
subset of M (as a set) and covers everything; the private times of 1 w.r.t. S contain
those w.r.t. M, so pruning "private(1) < q" is sound, and so is "uncovered > r*m".
Leaves S (|S| <= k) are completed with all multisets of k-|S| further classes
(any class, repetitions allowed) and filtered by private(1) >= q.
"""
import itertools
import numpy as np
from numba import njit


def tables(k, p):
    H = (p - 1) // 2
    m = (p - 1) // (k + 1)
    W = (H + 63) // 64
    cov = np.zeros((H + 1, W), dtype=np.uint64)
    covb = np.zeros((H + 1, H + 1), dtype=np.bool_)
    for c in range(1, H + 1):
        for s in range(1, H + 1):
            r = (s * c) % p
            if min(r, p - r) <= m:
                covb[c, s] = True
                cov[c, (s - 1) // 64] |= np.uint64(1) << np.uint64((s - 1) % 64)
    assert all(covb[c].sum() == m for c in range(1, H + 1))
    num = 2 * H - k * m
    q = -(-num // k)
    assert q >= 1
    return H, m, W, cov, covb, q


@njit(cache=True)
def _pc(x):
    c = 0
    while x:
        x &= x - np.uint64(1)
        c += 1
    return c


@njit(cache=True)
def dfs(k, H, m, W, cov, covl, q, out, maxout):
    """Returns (#leaves, #nodes). Leaves (chosen sets incl. class 1, 0-padded) -> out."""
    D = k + 1
    unc = np.zeros((D + 1, W), dtype=np.uint64)
    oth = np.zeros((D + 1, W), dtype=np.uint64)       # OR of covers of non-1 elements
    alive = np.zeros((D + 1, H + 1), dtype=np.bool_)
    cand = np.zeros((D + 1, H + 1), dtype=np.int64)
    ncand = np.full(D + 1, -1, dtype=np.int64)
    pos = np.zeros(D + 1, dtype=np.int64)
    chosen = np.zeros(D + 1, dtype=np.int64)
    okc = np.zeros(H + 1, dtype=np.bool_)
    topr = np.zeros(k + 1, dtype=np.int64)
    full = np.zeros(W, dtype=np.uint64)
    for s in range(1, H + 1):
        full[(s - 1) // 64] |= np.uint64(1) << np.uint64((s - 1) % 64)
    for w in range(W):
        unc[1, w] = full[w] & ~cov[1, w]
    for c in range(2, H + 1):
        alive[1, c] = True
    chosen[1] = 1
    d = 1
    nleaves = 0
    nnodes = 0
    while d >= 1:
        if ncand[d] == -1:
            nnodes += 1
            nu = 0
            for w in range(W):
                nu += _pc(unc[d, w])
            if nu == 0:
                if nleaves < maxout:
                    for i in range(k):
                        out[nleaves, i] = chosen[i + 1] if i + 1 <= d else 0
                nleaves += 1
                d -= 1
                continue
            r = k - d
            if r == 0 or nu > r * m:
                d -= 1
                continue
            # candidates must keep private(1) >= q: precompute ok[c] for alive c
            for c in range(2, H + 1):
                okc[c] = False
                if alive[d, c]:
                    pr = 0
                    for ww in range(W):
                        pr += _pc(cov[1, ww] & ~(oth[d, ww] | cov[c, ww]))
                    okc[c] = pr >= q
                    if okc[c]:
                        g = 0
                        for ww in range(W):
                            g += _pc(cov[c, ww] & unc[d, ww])
                        # keep the r largest gains in topr (descending)
                        if g > topr[r - 1]:
                            j = r - 1
                            while j > 0 and topr[j - 1] < g:
                                topr[j] = topr[j - 1]
                                j -= 1
                            topr[j] = g
            sg = 0
            for j in range(r):
                sg += topr[j]
                topr[j] = 0
            if sg < nu:          # even the r best remaining classes cannot cover U
                d -= 1
                continue
            # most constrained uncovered time
            best = -1
            bestc = 1 << 40
            for s in range(1, H + 1):
                w = (s - 1) // 64
                b = np.uint64(1) << np.uint64((s - 1) % 64)
                if unc[d, w] & b == 0:
                    continue
                cnt = 0
                for j in range(m):
                    if okc[covl[s, j]]:
                        cnt += 1
                if cnt < bestc:
                    bestc = cnt
                    best = s
                    if cnt == 0:
                        break
            if bestc == 0:
                d -= 1
                continue
            n = 0
            for j in range(m):
                c = covl[best, j]
                if okc[c]:
                    cand[d, n] = c
                    n += 1
            ncand[d] = n
            pos[d] = 0
        if pos[d] >= ncand[d]:
            ncand[d] = -1
            d -= 1
            continue
        c = cand[d, pos[d]]
        pos[d] += 1
        nd = d + 1
        for w in range(W):
            unc[nd, w] = unc[d, w] & ~cov[c, w]
            oth[nd, w] = oth[d, w] | cov[c, w]
        for cc in range(H + 1):
            alive[nd, cc] = alive[d, cc]
        alive[nd, c] = False
        alive[d, c] = False          # sibling exclusion
        chosen[nd] = c
        ncand[nd] = -1
        d = nd
    return nleaves, nnodes


def canon(row, p):
    """Orbit normal form under units, signs, permutations: min over ALL units u of
    sorted(fold(u*row)). Brute force over u (obviously orbit invariant)."""
    row = np.asarray(row, dtype=np.int64)
    u = np.arange(1, p, dtype=np.int64)[:, None]
    y = (u * row[None, :]) % p
    y = np.minimum(y, p - y)
    y.sort(axis=1)
    # lexicographic min
    idx = np.lexsort(y.T[::-1])
    return tuple(int(a) for a in y[idx[0]])


def private1(M, covb):
    cnt = np.zeros(covb.shape[1], dtype=np.int64)
    for c in M:
        cnt += covb[c]
    return int((covb[1] & (cnt == 1)).sum()) if list(M).count(1) == 1 else 0


def level_one_orbits(k, p, maxout=2_000_000):
    H, m, W, cov, covb, q = tables(k, p)
    out = np.zeros((maxout, k), dtype=np.int64)
    covl = np.zeros((H + 1, m), dtype=np.int64)
    for s in range(1, H + 1):
        cs = [c for c in range(1, H + 1) if covb[c, s]]
        assert len(cs) == m
        covl[s, :] = cs
    nl, nn = dfs(k, H, m, W, cov, covl, q, out, maxout)
    if nl > maxout:
        raise RuntimeError('maxout')
    orbits = set()
    for row in out[:nl]:
        S = [int(a) for a in row if a]
        e = k - len(S)
        for E in itertools.combinations_with_replacement(range(1, H + 1), e):
            M = S + list(E)
            if e and private1(M, covb) < q:
                continue
            orbits.add(canon(M, p))
    return orbits, dict(leaves=int(nl), nodes=int(nn), q=q, m=m, H=H)


def brute_level_one_orbits(k, p):
    """No pruning at all: every multiset of k classes containing 1, test full cover."""
    H = (p - 1) // 2
    m = (p - 1) // (k + 1)
    covb = np.zeros((H + 1, H), dtype=np.bool_)
    for c in range(1, H + 1):
        for s in range(1, H + 1):
            r = (s * c) % p
            covb[c, s - 1] = min(r, p - r) <= m
    orbits = set()
    combos = np.array(list(itertools.combinations_with_replacement(range(1, H + 1), k - 1)),
                      dtype=np.int64)
    acc = covb[1][None, :].copy() | covb[combos[:, 0]]
    for j in range(1, k - 1):
        acc |= covb[combos[:, j]]
    good = combos[acc.all(axis=1)]
    for rest in good:
        orbits.add(canon([1] + [int(a) for a in rest], p))
    return orbits


if __name__ == '__main__':
    import sys, time
    k = int(sys.argv[1])
    for p in map(int, sys.argv[2:]):
        t = time.time()
        o, st = level_one_orbits(k, p)
        print(k, p, len(o), st, round(time.time() - t, 2), flush=True)

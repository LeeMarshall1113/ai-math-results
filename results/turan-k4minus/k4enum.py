"""Enumerate K4^- -free 3-graphs on n<=7 vertices up to isomorphism.
K4^-: 4 vertices spanning >=3 edges.  Graph = bitmask over triples of range(n)
(using the global triple index for 7 vertices, so masks are comparable across n).
"""
import itertools, pickle, sys, time
import numpy as np

N = 7
TRI = list(itertools.combinations(range(N), 3))
TIDX = {t: i for i, t in enumerate(TRI)}

def tid(a, b, c):
    return TIDX[tuple(sorted((a, b, c)))]

def k4m_free(mask, verts):
    for q in itertools.combinations(verts, 4):
        c = 0
        for t in itertools.combinations(q, 3):
            if mask >> TIDX[t] & 1:
                c += 1
        if c >= 3:
            return False
    return True

PERM_MAPS = {}
def perm_maps(n):
    if n not in PERM_MAPS:
        tri_n = [t for t in TRI if max(t) < n]
        idx_n = [TIDX[t] for t in tri_n]
        perms = list(itertools.permutations(range(n)))
        M = np.zeros((len(perms), len(TRI)), dtype=np.int64)
        for pi, p in enumerate(perms):
            for t in tri_n:
                M[pi, TIDX[t]] = TIDX[tuple(sorted(p[x] for x in t))]
        PERM_MAPS[n] = (M, perms)
    return PERM_MAPS[n]

def canon(mask, n):
    M, _ = perm_maps(n)
    E = [i for i in range(len(TRI)) if mask >> i & 1]
    if not E:
        return 0
    vals = (np.left_shift(np.int64(1), M[:, E])).sum(axis=1)
    return int(vals.min())

def extend(classes, n):
    """classes: canonical masks on n-1 vertices -> all K4^- -free on n vertices."""
    v = n - 1
    pairs = list(itertools.combinations(range(v), 2))
    out = set()
    for g in classes:
        # candidate links: subsets of pairs; build incrementally with pruning
        def rec(i, link, mask):
            if i == len(pairs):
                out.add(canon(mask, n))
                return
            rec(i + 1, link, mask)
            a, b = pairs[i]
            m2 = mask | (1 << tid(a, b, v))
            # check 4-sets {v,a,b,c}
            ok = True
            for c in range(v):
                if c == a or c == b:
                    continue
                cnt = 0
                for t in ((a, b, v), (a, c, v), (b, c, v), (a, b, c)):
                    if m2 >> tid(*t) & 1:
                        cnt += 1
                if cnt >= 3:
                    ok = False
                    break
            if ok:
                rec(i + 1, link, m2)
        rec(0, None, g)
    return sorted(out)

if __name__ == "__main__":
    classes = {3: [0, 1 << tid(0, 1, 2)]}
    for n in range(4, N + 1):
        t0 = time.time()
        classes[n] = extend(classes[n - 1], n)
        print(n, len(classes[n]), f"{time.time()-t0:.1f}s", flush=True)
    classes[1] = [0]; classes[2] = [0]
    pickle.dump(classes, open(sys.argv[1] if len(sys.argv) > 1 else "classes.pkl", "wb"))

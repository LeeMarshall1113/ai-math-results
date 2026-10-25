"""Independent catalogue check: K4^- -free, pairwise non-isomorphic, orbit sum = labelled count.
Labelled counts: n<=6 by brute force over all 2^C(n,3) graphs; n=7 by summing, over every
6-vertex labelled K4^- -free G and every link L of vertex 6 (a graph on 6 vertices), the
validity test done directly in numpy (no superset-sum transform)."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import gzip, json, itertools
import numpy as np

def labelled(n):
    tri = list(itertools.combinations(range(n), 3)); ti = {t: i for i, t in enumerate(tri)}
    m = np.arange(1 << len(tri), dtype=np.int64)
    ok = np.ones(len(m), bool)
    for q in itertools.combinations(range(n), 4):
        c = sum(((m >> ti[t]) & 1) for t in itertools.combinations(q, 3))
        ok &= c < 3
    return m[ok], tri

for n in (4, 5, 6):
    print(n, len(labelled(n)[0]))
L6, tri6 = labelled(6)
ti6 = {t: i for i, t in enumerate(tri6)}
pairs = list(itertools.combinations(range(6), 2)); pi = {p: i for i, p in enumerate(pairs)}
# for each link L: forbidden 6-vertex triples (those with 2 link pairs); L invalid if a triple has 3
L = np.arange(1 << 15, dtype=np.int64)
valid = np.ones(len(L), bool); bad = np.zeros(len(L), np.int64)
for t in tri6:
    c = sum((L >> pi[p]) & 1 for p in itertools.combinations(t, 2))
    valid &= c < 3
    bad |= np.where(c == 2, 1 << ti6[t], 0)
badv = bad[valid]
# count pairs (G, L) with G & bad == 0, by histogramming bad masks then a direct loop over G chunks
ub, mult = np.unique(badv, return_counts=True)
total = 0
for i in range(0, len(ub), 64):
    blk = ub[i:i + 64]
    hit = (L6[:, None] & blk[None, :]) == 0
    total += int((hit.sum(axis=0) * mult[i:i + 64]).sum())
print("labelled n=7:", total)

cert = json.load(gzip.open(sys.argv[1], "rt"))
graphs = [int(g) for g in cert["graphs"]]
TRI = list(itertools.combinations(range(7), 3)); TI = {t: i for i, t in enumerate(TRI)}
P = np.array(list(itertools.permutations(range(7))))
T0 = np.zeros((len(P), 35), np.int64)
for t, (a, b, c) in enumerate(TRI):
    x = np.sort(P[:, [a, b, c]], axis=1)
    T0[:, t] = [TI[tuple(r)] for r in x.tolist()]
W = np.array([1 << t for t in range(35)], dtype=np.int64)
quads = [[TI[t] for t in itertools.combinations(q, 3)] for q in itertools.combinations(range(7), 4)]
seen = set(); orbit = 0
for H in graphs:
    hb = np.array([(H >> i) & 1 for i in range(35)], dtype=np.int64)
    assert all(hb[q].sum() < 3 for q in quads)
    codes = hb[T0] @ W
    can = int(codes.min()); aut = int((codes == H).sum())
    assert can not in seen; seen.add(can)
    orbit += 5040 // aut
print("classes", len(graphs), "orbit sum", orbit, "match:", orbit == total)

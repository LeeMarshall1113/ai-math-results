"""Exact flag-algebra data for pi(K4^-) on 7 vertices.

Writes data.npz with, all in exact integers:
  edges[h]                 number of edges of H_h (density = edges/35)
  for each block b (type sigma of size s in {1,3,5}, flags of size (7+s)/2):
     coefficient of H_h in <Q_b, [[F F^T]]_sigma> is  sum_{i,j} Q_ij * cnt_b[h,i,j] / DEN[s]
     with cnt = number of (ordered injective theta inducing sigma exactly,
     ordered split (A,B) of the rest) such that theta+A ~ F_i, theta+B ~ F_j,
     DEN[s] = (7)_s * #ordered splits   (so cnt/DEN is a probability)
  stat rows (degree stationarity), for each 1-flag F on 4 vertices:
     row_F[h] = P_{v,S,P}(F at (v,S), vPe in E) - P_{v,S,T}(F at (v,S), T in E)
     = num1/420 - num2/140, stored as integer 420*row = num1 - 3*num2.
Valid identity in an extremal limit: [[ (e1 - rho) F ]]_1 = 0 (all vertex degrees equal rho),
and rho*[[F]]_1 = e*[[F]]_1 is linear in phi.
"""
import itertools, pickle, sys, time
sys.path.append(r"C:/Users/hackathon/AppData/Roaming/Python/Python314/site-packages")
from math import comb, perm
import numpy as np
from k4enum import TRI, TIDX

N = int(sys.argv[1]) if len(sys.argv) > 1 else 7
SS = (1, 3, 5) if N == 7 else (0, 2, 4)
cl = pickle.load(open("classes.pkl", "rb"))
H7 = cl[N]; nH = len(H7)
TT = {}  # (a,b,c) unordered -> bit
for t in TRI:
    for p in itertools.permutations(t):
        TT[p] = 1 << TIDX[t]

def has(m, a, b, c):
    return m & TT[(a, b, c)]

def relabel(m, verts):
    out = 0
    for t in itertools.combinations(range(len(verts)), 3):
        if m & TT[(verts[t[0]], verts[t[1]], verts[t[2]])]:
            out |= 1 << TIDX[t]
    return out

def flag_key(m, typ, ext):
    best = None
    for p in itertools.permutations(ext):
        k = relabel(m, list(typ) + list(p))
        if best is None or k < best:
            best = k
    return best

# type representatives: labelled masks; every labelled K4^- -free graph on s vertices is a type
def labelled(s):
    tri = [t for t in TRI if max(t) < s]
    out = []
    for bits in range(1 << len(tri)):
        m = 0
        for i, t in enumerate(tri):
            if bits >> i & 1:
                m |= 1 << TIDX[t]
        ok = all(sum(1 for t in itertools.combinations(q, 3) if m & TT[t]) <= 2
                 for q in itertools.combinations(range(s), 4))
        if ok:
            out.append(m)
    return out

def canon_lab(m, s):
    return min(relabel(m, list(p)) for p in itertools.permutations(range(s)))

# one representative per isomorphism class of type (labelled representative = canonical form)
TYPES = {s: sorted(set(canon_lab(m, s) for m in labelled(s))) for s in SS}
print({s: len(v) for s, v in TYPES.items()}, flush=True)
DEN = {}
blocks = {}   # (s, typemask) -> {flagkey: idx}
trip = {}     # (s, typemask) -> list of (h, i, j, cnt)
t0 = time.time()
V = list(range(N))
for h, H in enumerate(H7):
    for s in SS:
        f = (N + s) // 2; m = f - s
        tset = set(TYPES[s])
        loc = {}
        nsplit = 0
        for theta in itertools.permutations(V, s):
            tm = relabel(H, list(theta))
            if tm not in tset:
                continue
            rest = [v for v in V if v not in theta]
            keys = {}
            for A in itertools.combinations(rest, m):
                keys[A] = flag_key(H, theta, A)
            for A in itertools.combinations(rest, m):
                B = tuple(v for v in rest if v not in A)
                d = blocks.setdefault((s, tm), {})
                ka, kb = keys[A], keys[B]
                ia = d.setdefault(ka, len(d)); ib = d.setdefault(kb, len(d))
                loc[(tm, ia, ib)] = loc.get((tm, ia, ib), 0) + 1
        for (tm, i, j), c in loc.items():
            trip.setdefault((s, tm), []).append((h, i, j, c))
    if h % 1000 == 0:
        print(h, f"{time.time()-t0:.0f}s", flush=True)
for s in SS:
    f = (N + s) // 2; m = f - s
    DEN[s] = perm(N, s) * comb(N - s, m)

edges = np.array([bin(H).count("1") for H in H7], dtype=np.int64)

# stationarity rows over 1-flags of size 4 (root first)
K = N - 3
fl4 = {}
for g in cl[K]:
    for root in range(K):
        fl4.setdefault(flag_key(g, [root], [v for v in range(K) if v != root]), None)
F4 = sorted(fl4)
stat = np.zeros((len(F4), nH), dtype=np.int64)
fidx = {F: i for i, F in enumerate(F4)}
for h, H in enumerate(H7):
    for v in V:
        R = [x for x in V if x != v]
        for S in itertools.combinations(R, K - 1):
            k = fidx[flag_key(H, [v], S)]
            R2 = [x for x in R if x not in S]
            n1 = sum(1 for P in itertools.combinations(R2, 2) if has(H, v, P[0], P[1]))
            n2 = 1 if has(H, *R2) else 0
            stat[k, h] += n1 - 3 * n2
print("stat rows", len(F4), "rank", np.linalg.matrix_rank(stat.astype(float)), flush=True)

keys = sorted(blocks)
out = {"N": N, "edges": edges, "stat420": stat, "F4": np.array(F4), "nblocks": len(keys)}
for b, k in enumerate(keys):
    s, tm = k
    T = np.array(trip[k], dtype=np.int64)
    out[f"b{b}_s"] = s; out[f"b{b}_type"] = tm; out[f"b{b}_den"] = DEN[s]
    fk = sorted(blocks[k].items(), key=lambda x: x[1])
    out[f"b{b}_flags"] = np.array([x[0] for x in fk], dtype=np.int64)
    out[f"b{b}_trip"] = T
    print(b, s, tm, "flags", len(fk), "nnz", len(T), flush=True)
np.savez_compressed(f"data{N}.npz", **out)
print("saved", f"{time.time()-t0:.0f}s")

"""PGL(3,3)-orbits on the Baer-invariant planar (2,5)-sets with empty fixed pattern (size 26)."""
import itertools, json, os
from ff import GF
from pg import PG
import plane_tables as T
F = GF(9); P2 = PG(3, 9, F)
fam, inst = T.invariant_family(P2, [[1,0,0],[0,1,0],[0,0,1]], 1)
std = set(i for i, v in enumerate(P2.points) if all(c in (0,1,2) for c in v))
empt = [X for X in fam[26] if not (set(X) & std)]
print("empty-pattern sets:", len(empt), flush=True)
# PGL(3,3) = GL(3,3) matrices acting on P2
perms = []
seenp = set()
for entries in itertools.product(range(3), repeat=9):
    M = [list(entries[0:3]), list(entries[3:6]), list(entries[6:9])]
    det = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])) % 3
    if det == 0: continue
    p = tuple(P2.perm_of_matrix(M, 0))
    if p in seenp: continue
    seenp.add(p); perms.append(p)
print("PGL(3,3) permutations:", len(perms), flush=True)
orbit_id = {}
reps = []
for X in empt:
    if X in orbit_id: continue
    idx = len(reps); reps.append(X)
    for p in perms:
        Y = tuple(sorted(p[a] for a in X))
        orbit_id[Y] = idx
sizes = {}
for X, i in orbit_id.items(): sizes[i] = sizes.get(i, 0) + 1
print("orbits on the empty-pattern sets:", len(reps), "sizes", sorted(sizes.values()), flush=True)
json.dump(dict(reps=reps), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "baer_empty_reps.json"), "w"))

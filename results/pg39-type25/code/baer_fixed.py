"""Possible sets S = K n PG(3,3) for a Baer-involution-invariant (2,5)-set K of PG(3,9):
every invariant plane (the 40 planes of PG(3,3)) meets K in a Baer-invariant planar (2,5)-set, so the
fixed-point pattern S n pi must be one of the patterns computed in baer_planar.py (transported to pi)."""
import itertools, json, os, time
from collections import Counter
from pysat.solvers import Cadical195
from ff import GF
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
F = GF(9)
P3 = PG(4, 9, F)       # PG(3,9)
P2 = PG(3, 9, F)       # PG(2,9)
pats = json.load(open(os.path.join(HERE, "..", "results", "planar_invariant_patterns.json")))["baer"]
allowed_std = set()
for k in ("26", "35"):
    for p in pats[k]:
        allowed_std.add(tuple(sorted(p)))
# standard plane PG(2,3): points of PG(2,9) with GF(3) coordinates
std = [i for i, v in enumerate(P2.points) if all(c in (0, 1, 2) for c in v)]
assert len(std) == 13
# PG(3,3)
B = [i for i, v in enumerate(P3.points) if all(c in (0, 1, 2) for c in v)]
assert len(B) == 40
Bset = set(B)
# planes of PG(3,3): plane normals with GF(3) coordinates
planes = []
for a in P3.points:
    if all(c in (0, 1, 2) for c in a):
        pl = [i for i in B if P3._dot(a, P3.points[i]) == 0]
        assert len(pl) == 13
        planes.append(pl)
assert len(planes) == 40
def mul3(c, v): return tuple(F.mul[c][x] for x in v)
def add(v, w): return tuple(F.add[x][y] for x, y in zip(v, w))
def transport(pl):
    """a collineation from the standard plane (13 GF(3)-points of PG(2,9)) to plane pl: dictionary std point -> PG(3,3) point"""
    vecs = [P3.points[i] for i in pl]
    # choose 3 independent GF(3)-vectors among vecs
    basis = []
    span = {(0, 0, 0, 0)}
    for v in vecs:
        if v not in span:
            basis.append(v)
            new = set()
            for s in span:
                for c in (0, 1, 2):
                    new.add(add(s, mul3(c, v)))
            span = new
        if len(basis) == 3: break
    assert len(basis) == 3
    mp = {}
    for i in std:
        y = P2.points[i]
        w = (0, 0, 0, 0)
        for c, b in zip(y, basis):
            w = add(w, mul3(c, b))
        mp[i] = P3.pidx[P3.normalize(w)]
    assert sorted(mp.values()) == sorted(pl)
    return mp
# SAT over the 40 fixed points
var = {p: j + 1 for j, p in enumerate(B)}
clauses = []
t0 = time.time()
for pl in planes:
    mp = transport(pl)
    allowed = set(tuple(sorted(mp[p] for p in pat)) for pat in allowed_std)
    pts = sorted(pl)
    for bits in itertools.product((0, 1), repeat=13):
        S = tuple(p for p, b in zip(pts, bits) if b)
        if S not in allowed:
            clauses.append([(-var[p] if b else var[p]) for p, b in zip(pts, bits)])
print("plane constraints: %d forbidden patterns (%.1fs)" % (len(clauses), time.time() - t0), flush=True)
# lines of PG(3,3): at most 3 points of S on a line (f <= 3)
Ls = set()
for L in P3.lines:
    inter = [p for p in L if p in Bset]
    if len(inter) == 4:
        Ls.add(tuple(inter))
assert len(Ls) == 130
for L in Ls:
    clauses.append([-var[p] for p in L])
s = Cadical195(bootstrap_with=clauses)
sols = []
t0 = time.time()
while s.solve():
    m = s.get_model(); val = {abs(l): l > 0 for l in m}
    S = tuple(p for p in B if val[var[p]])
    sols.append(S)
    s.add_clause([(-var[p] if val[var[p]] else var[p]) for p in B])
    if len(sols) % 20000 == 0: print("  ", len(sols), flush=True)
print("number of admissible sets S = K n PG(3,3): %d  sizes %s  (%.1fs)" % (len(sols), dict(sorted(Counter(len(x) for x in sols).items())), time.time() - t0), flush=True)
json.dump(sols, open(os.path.join(HERE, "..", "results", "baer_fixed_sets.json"), "w"))

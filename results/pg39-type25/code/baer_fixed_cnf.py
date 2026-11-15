"""CNF for 'S = K n PG(3,3) of a sigma-invariant (2,5)-set K': 40 variables (the points of PG(3,3)), for each of the 40 planes of PG(3,3)
the fixed pattern S n plane must be the fixed part of a sigma-invariant planar (2,5)-set (patterns from certs/table_baer.json), and a line of
PG(3,3) contains at most 3 points of S (f + 2t in {2,5} with t <= 3 forces f <= 3).  Plus one blocking clause per solution of
results/baer_fixed_sets.json => UNSAT (completeness of the list of 81 sets)."""
import itertools
import json
import os

import baer3d
from cnfenc import CNF
from ff import GF
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def build():
    F = GF(9)
    P3 = PG(4, 9, F)
    P2 = PG(3, 9, F)
    Z = {int(k): [tuple(x) for x in v] for k, v in json.load(open(os.path.join(CERTS, "table_baer.json"))).items()}
    std = set(i for i, v in enumerate(P2.points) if all(c in (0, 1, 2) for c in v))
    allowed_std = set()
    for k in (26, 35):
        for X in Z[k]:
            allowed_std.add(tuple(sorted(p for p in X if p in std)))
    B = [i for i, v in enumerate(P3.points) if all(c in (0, 1, 2) for c in v)]
    cnf = CNF()
    var = {p: cnf.new_var() for p in B}
    for (a, pts, fixed, basis) in baer3d.gf3_planes(P3):
        mp = baer3d.transport(P2, P3, basis)
        allowed = set(tuple(sorted(mp[p] for p in pat)) for pat in allowed_std)
        fx = sorted(fixed)
        for bits in itertools.product((0, 1), repeat=len(fx)):
            S = tuple(p for p, b in zip(fx, bits) if b)
            if S not in allowed:
                cnf.add(*[(-var[p] if b else var[p]) for p, b in zip(fx, bits)])
    Bset = set(B)
    Ls = set()
    for L in P3.lines:
        inter = tuple(p for p in L if p in Bset)
        if len(inter) == 4:
            Ls.add(inter)
    for L in sorted(Ls):
        cnf.add(*[-var[p] for p in L])
    sols = json.load(open(os.path.join(HERE, "..", "results", "baer_fixed_sets.json")))
    for S in sols:
        Ss = set(S)
        cnf.add(*[(-var[p] if p in Ss else var[p]) for p in B])
    return cnf, sols

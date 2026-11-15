"""Fixed-point patterns of Baer-invariant planar (2,5)-sets of PG(2,9), and the same for the planar homology."""
import itertools, json, os, sys, time
from collections import Counter
from pysat.solvers import Cadical195
from ff import GF
from pg import PG
from orbsat import OrbitInstance

F = GF(9)
P = PG(3, 9, F)
I3 = [[1,0,0],[0,1,0],[0,0,1]]
res = {}
for name, (M, f) in {"baer": (I3, 1), "homology": ([[1,0,0],[0,1,0],[0,0,F.neg[1]]], 0)}.items():
    perm = P.perm_of_matrix(M, f)
    inst = OrbitInstance(P, [perm])
    fixed = [i for i, o in enumerate(inst.orbs) if len(o) == 1]
    out = {}
    for k in (26, 35):
        cnf = inst.build_cnf(k, (2, 5), amo=True)
        s = Cadical195(bootstrap_with=cnf.clauses)
        pats = set()
        nsol = 0
        t = time.time()
        while s.solve():
            m = s.get_model(); val = {abs(l): l > 0 for l in m}
            pat = tuple(inst.orbs[i][0] for i in fixed if val[inst.z[i]])
            pats.add(pat)
            s.add_clause([(-inst.z[i] if val[inst.z[i]] else inst.z[i]) for i in fixed])
        out[k] = sorted(pats)
        print("%-8s k=%d: %d fixed-point patterns (fixed points %d), sizes of patterns %s  (%.1fs)" % (
            name, k, len(pats), len(fixed), dict(sorted(Counter(len(p) for p in pats).items())), time.time() - t), flush=True)
        s.delete()
    res[name] = {k: [list(p) for p in v] for k, v in out.items()}
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "planar_invariant_patterns.json"), "w"))

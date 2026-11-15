"""Number of Baer-invariant planar (2,5)-sets (sizes 26, 35) of PG(2,9) per fixed-point pattern."""
import json, os, time
from collections import Counter
from pysat.solvers import Cadical195
from ff import GF
from pg import PG
from orbsat import OrbitInstance
F = GF(9); P = PG(3, 9, F)
perm = P.perm_of_matrix([[1,0,0],[0,1,0],[0,0,1]], 1)
inst = OrbitInstance(P, [perm])
fixed = [i for i, o in enumerate(inst.orbs) if len(o) == 1]
pairs = [i for i, o in enumerate(inst.orbs) if len(o) == 2]
res = {}
for k in (26, 35):
    cnf = inst.build_cnf(k, (2, 5), amo=True)
    s = Cadical195(bootstrap_with=cnf.clauses)
    cnt = Counter(); n = 0; t = time.time()
    allsols = []
    while s.solve():
        m = s.get_model(); val = {abs(l): l > 0 for l in m}
        pat = tuple(inst.orbs[i][0] for i in fixed if val[inst.z[i]])
        pr = tuple(i for i in pairs if val[inst.z[i]])
        cnt[len(pat)] += 1; n += 1
        allsols.append((pat, pr))
        s.add_clause([(-inst.z[i] if val[inst.z[i]] else inst.z[i]) for i in fixed + pairs])
    print("k=%d: %d solutions in %.1fs; by |fixed pattern|: %s" % (k, n, time.time()-t, dict(sorted(cnt.items()))), flush=True)
    res[k] = allsols
    s.delete()
json.dump({str(k): v for k, v in res.items()}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "baer_planar_tables.json"), "w"))

"""Quick probes of the planar (2,5) SAT instance (speed of finding / refuting)."""
import sys
import time
from collections import Counter

from pysat.solvers import Cadical195

from cnfenc import CNF, sum_in_set, sum_equals
from pg import PG
import planar

P = PG(3, 9)
l0 = planar.line0_points(P)


def build(k, amo=True):
    cnf = CNF()
    x = [cnf.new_var() for _ in range(P.N)]
    for L in P.lines:
        sum_in_set(cnf, [x[p] for p in L], [1] * len(L), {2, 5}, amo=amo)
    sum_equals(cnf, x, [1] * P.N, k, amo=False)
    return cnf, x


def run(k, extra_units, budget=3000000, label=""):
    cnf, x = build(k)
    for u in extra_units(x):
        cnf.add(u)
    s = Cadical195(bootstrap_with=cnf.clauses)
    t = time.time()
    s.conf_budget(budget)
    r = s.solve_limited()
    print("k=%d %s: result %s in %.1fs (vars %d, clauses %d)" % (k, label, r, time.time() - t, cnf.nvars, len(cnf.clauses)), flush=True)
    if r:
        m = s.get_model()
        val = {abs(l): l > 0 for l in m}
        K = set(p for p in range(P.N) if val[x[p]])
        print("   |K| =", len(K), dict(Counter(sum(1 for p in L if p in K) for L in P.lines)), flush=True)
    s.delete()
    return r


if __name__ == "__main__":
    for k in (26, 35):
        run(k, lambda x: [x[0]], label="plain, point 0 in K")
    P1, reps = planar.pgammal2_orbits_on_5subsets()

    def lift(i):
        v = P1.points[i]
        return P.pidx[P.normalize((v[0], v[1], 0))]
    origin = P.pidx[(0, 0, 1)]
    for k in (26, 35):
        for ri, (R1, osz) in enumerate(reps):
            Rp = set(lift(i) for i in R1)
            def extra(x, Rp=Rp):
                u = [x[p] if p in Rp else -x[p] for p in l0]
                u.append(x[origin])
                return u
            run(k, extra, label="R#%d (orbit %d) + origin" % (ri, osz))

"""Classification of (2,5)-sets of size 26 and 35 in PG(2,9) up to PGammaL(3,9).

Every such K has a 5-secant l0 (b_P = 5 or 8 > 0).  Moving l0 to the line x3 = 0, K meets it
in a 5-subset R, and R can be moved to one representative per PGammaL(2,9)-orbit of 5-subsets
of the 10 points of l0 (there are two orbits, of sizes 180 and 72).  The stabiliser of (l0, R)
contains all translations of the affine plane AG(2,9) = PG(2,9) - l0, so we may also assume
that the affine origin (0,0,1) is in K.  What remains of the stabiliser is
    H_R = { diag(A,1) * frob^e : A in GL(2,9), e in {0,1}, A.frob^e stabilises R }
(order 64 for the 180-orbit, 160 for the 72-orbit) which acts on the 81 affine points; we add the
lex-leader constraint x <=_lex x.h for every h in H_R (sound: each H_R-orbit of solutions keeps its
lex-least member).

Stage 1 (this file, Windows python with PySAT/CaDiCaL): enumerate ALL sets K with
  K meets l0 exactly in R, origin in K, |K| = k, every line meeting K in 2 or 5 points,
  x <=_lex x.h for h in H_R;
solutions are saved to results/planar_k<k>_R<i>.json.
Stage 2 (planar_classify.py, WSL python with pynauty): canonical forms of the incidence structure.
"""
import itertools
import json
import os
import sys
import time

from pysat.solvers import Cadical195

from cnfenc import CNF, sum_in_set, sum_equals
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))


def line0_points(P):
    return [i for i, v in enumerate(P.points) if v[2] == 0]


def pgammal2_orbits_on_5subsets(sz=5):
    """orbits of PGammaL(2,9) (acting on the 10 points of PG(1,9)) on sz-subsets (default 5)."""
    P1 = PG(2, 9)      # PG(1,9): 10 points
    F = P1.F
    omega = F.prim
    gens = [
        P1.perm_of_matrix([[1, 1], [0, 1]]),
        P1.perm_of_matrix([[omega, 0], [0, 1]]),
        P1.perm_of_matrix([[0, 1], [1, 0]]),
        P1.perm_of_matrix([[1, 0], [0, 1]], frob=1),
    ]
    idp = tuple(range(10))
    seen = {idp}
    frontier = [idp]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                c = tuple(h[g[i]] for i in range(10))
                if c not in seen:
                    seen.add(c)
                    nxt.append(c)
        frontier = nxt
    assert len(seen) == 1440, len(seen)
    subsets = list(itertools.combinations(range(10), sz))
    orbit_of = {}
    reps = []
    for s in subsets:
        if s in orbit_of:
            continue
        idx = len(reps)
        orb = set()
        for g in seen:
            orb.add(tuple(sorted(g[i] for i in s)))
        for t in orb:
            orbit_of[t] = idx
        reps.append((s, len(orb)))
    return P1, reps


def stabilizer_perms(P, Rpts):
    """all diag(A,1)*frob^e (A in GL(2,9)) as permutations of the 91 points, that stabilise
    the set Rpts of points of l0 setwise.  Returns list of perms (including the identity)."""
    F = P.F
    out = []
    q = F.q
    for a, b, c, d in itertools.product(range(q), repeat=4):
        det = F.sub(F.mul[a][d], F.mul[b][c])
        if det == 0:
            continue
        M = [[a, b, 0], [c, d, 0], [0, 0, 1]]
        for e in (0, 1):
            perm = P.perm_of_matrix(M, e)
            if set(perm[p] for p in Rpts) == set(Rpts):
                out.append(perm)
    return out


def add_lex_leader(cnf, x, perm, order):
    """x <=_lex (x o perm) along the variable order `order` (list of point indices)."""
    prev = None       # variable meaning 'equal so far' (None = true)
    for p in order:
        a, b = x[p], x[perm[p]]
        if a == b:
            continue
        if prev is None:
            cnf.add(-a, b)
        else:
            cnf.add(-prev, -a, b)
        e = cnf.new_var()
        if prev is None:
            cnf.add(a, b, e)
            cnf.add(-a, -b, e)
        else:
            cnf.add(-prev, a, b, e)
            cnf.add(-prev, -a, -b, e)
        prev = e


def build_instance(P, k, Rpts, l0, origin, Hperms, amo=True):
    cnf = CNF()
    x = [cnf.new_var() for _ in range(P.N)]
    y = []        # y_l <-> |l n K| = 5
    for L in P.lines:
        fin = sum_in_set(cnf, [x[p] for p in L], [1] * len(L), {2, 5}, amo=amo)
        y.append(fin[5])
    sum_equals(cnf, x, [1] * P.N, k, amo=False)
    # redundant (implied) local counts: a point of K is on d+3 five-secants, a point off K on d
    d = (k - 20) // 3
    for p in range(P.N):
        thr = [y[j] for j in P.lines_through[p]]
        fin = sum_in_set(cnf, thr, [1] * len(thr), {d, d + 3}, amo=False)
        cnf.add(-x[p], fin[d + 3])
        cnf.add(x[p], -fin[d + 3])
    l0s = set(l0)
    for p in l0:
        cnf.add(x[p] if p in Rpts else -x[p])
    cnf.add(x[origin])
    order = [p for p in range(P.N) if p not in l0s]
    ident = list(range(P.N))
    for h in Hperms:
        if h != ident:
            add_lex_leader(cnf, x, h, order)
    return cnf, x, order


def enumerate_planar(k, Rpts, P, l0, origin, Hperms, limit=None, log=None, max_seconds=None, prior=()):
    """returns (solutions, finished).  `prior` = solutions found in earlier chunks (they are blocked again)."""
    cnf, x, order = build_instance(P, k, Rpts, l0, origin, Hperms)
    s = Cadical195(bootstrap_with=cnf.clauses)
    sols = [tuple(K) for K in prior]
    for K in sols:
        Ks = set(K)
        s.add_clause([(-x[p] if p in Ks else x[p]) for p in order])
    t0 = time.time()
    finished = False
    while True:
        if max_seconds is not None:
            # CaDiCaL (PySAT) cannot be interrupted portably: use a conflict budget per call
            s.conf_budget(2000000)
            r = s.solve_limited()
            if r is None:
                if time.time() - t0 > max_seconds:
                    break
                continue
        else:
            r = s.solve()
        if not r:
            finished = True
            break
        m = s.get_model()
        val = {abs(l): l > 0 for l in m}
        K = tuple(p for p in range(P.N) if val[x[p]])
        sols.append(K)
        s.add_clause([(-x[p] if val[x[p]] else x[p]) for p in order])
        if log and len(sols) % 500 == 0:
            log("   ... %d solutions, %.0fs" % (len(sols), time.time() - t0))
        if limit and len(sols) >= limit:
            break
        if max_seconds is not None and time.time() - t0 > max_seconds:
            break
    s.delete()
    return sols, finished


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, action="append")
    ap.add_argument("--ri", type=int, action="append")
    ap.add_argument("--chunk", type=int, default=None, help="stop after this many seconds and save a partial result")
    ap.add_argument("--sec", type=int, default=5, help="size of the secant line section that is fixed (5 or 2)")
    a = ap.parse_args()
    ks = a.k or [26, 35]
    P = PG(3, 9)
    P.verify()
    l0 = line0_points(P)
    assert len(l0) == 10
    P1, reps = pgammal2_orbits_on_5subsets(a.sec)
    print("PGammaL(2,9) orbits on %d-subsets of the 10 points of a line: %d  (sizes %s)" % (a.sec, len(reps), [r[1] for r in reps]), flush=True)

    def lift(i):
        v = P1.points[i]
        return P.pidx[P.normalize((v[0], v[1], 0))]
    origin = P.pidx[(0, 0, 1)]
    res = {}
    for ri, (R1, osz) in enumerate(reps):
        if a.ri is not None and ri not in a.ri:
            continue
        Rpts = set(lift(i) for i in R1)
        H = stabilizer_perms(P, Rpts)
        print("R#%d: orbit size %d, |H_R| = %d" % (ri, osz, len(H)), flush=True)
        for k in ks:
            path = os.path.join(HERE, "..", "results", ("planar_k%d_R%d.json" if a.sec == 5 else "planar2sec_k%d_R%d.json") % (k, ri))
            if os.path.exists(path):
                print("skip existing", path, flush=True)
                continue
            t = time.time()
            ppath = path + ".partial"
            prior = json.load(open(ppath))["solutions"] if os.path.exists(ppath) else []
            sols, finished = enumerate_planar(k, Rpts, P, l0, origin, H, log=lambda s: print(s, flush=True),
                                              max_seconds=a.chunk, prior=prior)
            print("k=%d R#%d (orbit size %d): %d solutions so far, finished=%s (%.1fs)" % (k, ri, osz, len(sols), finished, time.time() - t), flush=True)
            if finished:
                json.dump(dict(k=k, R=sorted(Rpts), orbit_size=osz, Hsize=len(H), solutions=sols), open(path, "w"))
                if os.path.exists(ppath):
                    os.remove(ppath)
            else:
                json.dump(dict(solutions=sols), open(ppath, "w"))
            res[(k, ri)] = len(sols)
    print(res)


if __name__ == "__main__":
    main()

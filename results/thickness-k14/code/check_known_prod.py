"""Soundness test of the PRODUCTION case CNFs (run_all.case_clauses: hull units, sweeps, summed
Euler, lex, colour precedence): a real drawing must satisfy the CNF of its own hull case for some
choice of p1 / mirror / colour permutation. usage: python check_known_prod.py witness.json"""
import sys, json, itertools, functools
from pysat.solvers import Solver
from run_all import case_clauses


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


d = json.load(open(sys.argv[1])); P = [tuple(p) for p in d["points"]]; n = len(P)
t = d.get("layers") or max(L for _, _, L in d["edges"])
col = {(min(u, v), max(u, v)): L - 1 for u, v, L in d["edges"]}
found = False
for mirror in (False, True):
    Pm = [(-x, y) for x, y in P] if mirror else P
    for h in range(n):
        others = [i for i in range(n) if i != h]
        if not any(all(orient(Pm[h], Pm[a], Pm[b]) > 0 for b in others if b != a) for a in others):
            continue
        order = sorted(others, key=functools.cmp_to_key(lambda a, b: -1 if orient(Pm[h], Pm[a], Pm[b]) > 0 else 1))
        lab = [h] + order
        hull = [i for i in range(n) if all(not (orient(Pm[a], Pm[b], Pm[i]) * orient(Pm[a], Pm[b], Pm[c]) > 0 and orient(Pm[b], Pm[c], Pm[i]) * orient(Pm[b], Pm[c], Pm[a]) > 0 and orient(Pm[c], Pm[a], Pm[i]) * orient(Pm[c], Pm[a], Pm[b]) > 0) for a, b, c in itertools.combinations([j for j in range(n) if j != i], 3))]
        S = tuple(sorted(lab.index(i) + 1 for i in hull if lab.index(i) + 1 not in (1, 2, n)))
        kind = {0: "h3", 1: "h4", 2: "h5"}.get(len(S), "ge6")
        if kind == "ge6": S = S[:3]
        nv, cl, e = case_clauses(n, t, kind, S)
        s = Solver(name="cadical195", bootstrap_with=cl)
        xs = [v if orient(Pm[lab[a-1]], Pm[lab[b-1]], Pm[lab[c-1]]) > 0 else -v for (a, b, c), v in e.X.items()]
        for perm in itertools.permutations(range(t)):
            cs = []
            for ((a, b), k), v in e.C.items():
                o = (min(lab[a-1], lab[b-1]), max(lab[a-1], lab[b-1]))
                cs.append(v if perm[col[o]] == k else -v)
            if s.solve(assumptions=xs + cs):
                print("consistent: mirror", mirror, "p1", h, "case", kind, S, "perm", perm); found = True; break
        if found: break
    if found: break
print("PRODUCTION ENCODING ACCEPTS WITNESS" if found else "REJECTED -> encoding unsound!")
sys.exit(0 if found else 1)

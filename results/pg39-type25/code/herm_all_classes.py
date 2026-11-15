"""Positive control for EVERY class: is there a set with line counts in {1,4,10} and |K| = 280 invariant under the class
representative?  (The Hermitian surface is invariant exactly under the classes that meet PGammaU(4,3); predictions in NOTES.)"""
import json, os, sys, time
import groups
from orbsat import OrbitInstance, run_solver, check_K
from controls import line_stats, plane_stats
HERE = os.path.dirname(os.path.abspath(__file__))
F, P, gens, reps = groups.build_all()
timeout = int(sys.argv[1]) if len(sys.argv) > 1 else 120
res = {}
for name, (M, f) in reps.items():
    inst = OrbitInstance(P, [P.perm_of_matrix(M, f)])
    cnf = inst.build_cnf(280, (1, 4, 10), amo=True)
    path = os.path.join(HERE, "..", "cnf", "herm_%s.cnf" % name.replace("+", "p"))
    cnf.to_dimacs(path)
    st, model, t, info = run_solver(path, "cadical", timeout)
    extra = ""
    if st == "SATISFIABLE":
        tr = set(m for m in model if m > 0)
        K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
        ok, bad = check_K(P, K, (1, 4, 10))
        ls, ps = dict(line_stats(P, K)), dict(plane_stats(P, K))
        extra = "|K|=%d exact=%s lines=%s planes=%s" % (len(K), ok, ls, ps)
    print("%-20s %-14s %6.1fs %s" % (name, st, t, extra), flush=True)
    res[name] = dict(status=st, seconds=round(t, 1), extra=extra)
json.dump(res, open(os.path.join(HERE, "..", "results", "hermitian_all_classes.json"), "w"), indent=1)

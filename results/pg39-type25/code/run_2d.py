"""Baer involution class 2d: the 3 admissible fixed-point sets S = K n PG(3,3) (baer_fixed.py):
 S empty (k = 260), S a point (k = 287), S = complement of a plane of PG(3,3) (k = 287)."""
import json, os, sys, time
import groups
from orbsat import OrbitInstance, run_solver, check_K, wsl_path
from pg import PG
HERE = os.path.dirname(os.path.abspath(__file__))
F, P, gens, reps = groups.build_all()
M, f = reps["2d_baer"]
inst = OrbitInstance(P, [P.perm_of_matrix(M, f)])
fixed = [i for i, o in enumerate(inst.orbs) if len(o) == 1]
fpts = {inst.orbs[i][0]: i for i in fixed}
sets = json.load(open(os.path.join(HERE, "..", "results", "baer_fixed_sets.json")))
by = {}
for S in sets:
    by.setdefault(len(S), []).append(S)
cases = [("S_empty", by[0][0], 260), ("S_point", by[1][0], 287), ("S_27", by[27][0], 287)]
timeout = int(sys.argv[1]) if len(sys.argv) > 1 else 600
for name, S, k in cases:
    Sset = set(S)
    units = [(fpts[p], 1 if p in Sset else 0) for p in fpts]
    cnf = inst.build_cnf(k, (2, 5), amo=True, units=units)
    path = os.path.join(HERE, "..", "cnf", "2d_%s_k%d.cnf" % (name, k))
    cnf.to_dimacs(path, comments=["class 2d_baer k %d case %s fixed set size %d" % (k, name, len(S))])
    st, model, t, info = run_solver(path, "cadical", timeout)
    print("2d %s k=%d : %s in %.1fs %s (vars %d clauses %d)" % (name, k, st, t, info, cnf.nvars, len(cnf.clauses)), flush=True)
    if st == "SATISFIABLE":
        tr = set(m for m in model if m > 0)
        K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
        ok, bad = check_K(P, K, (2, 5))
        print("  MODEL |K| = %d exact line check %s" % (len(K), ok), flush=True)

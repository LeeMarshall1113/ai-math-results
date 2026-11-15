"""Experiment: effect of redundant constraints (local counts, plane sizes) on a hard orbit instance."""
import sys, time, json, os
import groups, symbreak
from orbsat import OrbitInstance, run_solver
from pg import PG
cls, k = sys.argv[1], int(sys.argv[2])
variants = sys.argv[3].split(",")        # e.g. base,local,planes,both
budget = int(sys.argv[4]) if len(sys.argv) > 4 else 300
case = int(sys.argv[5]) if len(sys.argv) > 5 else 0
order = [int(x) for x in sys.argv[6].split(",")] if len(sys.argv) > 6 else None
F, P, gens, reps = groups.build_all()
N = symbreak.load_normalizers()
M, f = reps[cls]
inst = OrbitInstance(P, [P.perm_of_matrix(M, f)])
cs = symbreak.cases(inst.orbs, inst.oid, N[cls]["gens"], order)
members, rep, forb = cs[case]
units = [(rep, 1)] + [(o, 0) for o in forb]
for v in variants:
    t0 = time.time()
    cnf = inst.build_cnf(k, (2, 5), amo=True, units=units, local=(v in ("local", "both")), planes=(v in ("planes", "both")))
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cnf", "exp_%s_k%d_%s.cnf" % (cls, k, v))
    cnf.to_dimacs(path)
    tb = time.time() - t0
    st, model, t, info = run_solver(path, "cadical", budget)
    print("%s k=%d variant %-7s : %d vars %d clauses (build %.0fs) -> %s in %.1fs %s" % (cls, k, v, cnf.nvars, len(cnf.clauses), tb, st, t, info), flush=True)

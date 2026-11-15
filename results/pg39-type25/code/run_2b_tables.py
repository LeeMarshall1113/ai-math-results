"""2b hard case (|F1|,|F2|,k) = (2,2,260) with plane-table constraints on the 20 invariant planes."""
import json, os, sys, time
import cases as C, groups, plane_tables as T
from orbsat import OrbitInstance, run_solver, check_K
from pg import PG
from ff import GF
HERE = os.path.dirname(os.path.abspath(__file__))
F, P3, gens, reps = groups.build_all()
M, f = reps["2b_homology_2+2"]
inst3 = OrbitInstance(P3, [P3.perm_of_matrix(M, f)])
P2 = PG(3, 9, F)
t0 = time.time()
fam, inst2 = T.invariant_family(P2, [[1, 0, 0], [0, 1, 0], [0, 0, F.neg[1]]], 0)
print("homology-invariant planar sets: %s (%.1fs)" % ({k: len(v) for k, v in fam.items()}, time.time() - t0), flush=True)
fixed2 = set(i for i, v in enumerate(P2.points) if v[2] == 0) | {P2.pidx[(0, 0, 1)]}
E = {1: (1, 0, 0, 0), 2: (0, 1, 0, 0), 3: (0, 0, 1, 0), 4: (0, 0, 0, 1)}
l1 = [i for i, v in enumerate(P3.points) if v[2] == 0 and v[3] == 0]
l2 = [i for i, v in enumerate(P3.points) if v[0] == 0 and v[1] == 0]
cases = [c for c in C.cases_2b(P3)]
names = sys.argv[2].split(",") if len(sys.argv) > 2 else [c["name"] for c in cases if "F2_2_k260" in c["name"]]
timeout = int(sys.argv[1]) if len(sys.argv) > 1 else 600
opt = {p: o for p, o in enumerate(inst3.oid)}
for case in cases:
    if case["name"] not in names:
        continue
    ins = set(case["ins"])
    cnf = C.build_case_cnf(inst3, case, opt)
    planes_rows = []
    for p2 in l2:
        mp = T.linear_map_points(P2, P3, [E[1], E[2], P3.points[p2]])
        want = (ins & set(l1)) | ({p2} if p2 in ins else set())
        planes_rows.append((sorted(mp), T.allowed_rows(fam, mp, fixed2, want)))
    for p1 in l1:
        mp = T.linear_map_points(P2, P3, [E[3], E[4], P3.points[p1]])
        want = (ins & set(l2)) | ({p1} if p1 in ins else set())
        planes_rows.append((sorted(mp), T.allowed_rows(fam, mp, fixed2, want)))
    print(case["name"], "allowed rows per plane:", sorted(set(len(r) for pts, r in planes_rows)), flush=True)
    T.add_plane_tables(cnf, inst3, planes_rows)
    path = os.path.join(HERE, "..", "cnf", "2bT_%s.cnf" % case["name"])
    cnf.to_dimacs(path, comments=["class 2b case %s tables" % case["name"]])
    st, model, t, info = run_solver(path, "cadical", timeout)
    print("2bT %s : %s in %.1fs %s (vars %d clauses %d)" % (case["name"], st, t, info, cnf.nvars, len(cnf.clauses)), flush=True)
    if st == "SATISFIABLE":
        tr = set(m for m in model if m > 0)
        K = sorted(p for i, v in enumerate(inst3.z) if v in tr for p in inst3.orbs[i])
        ok, bad = check_K(P3, K, (2, 5))
        print("  MODEL |K| = %d exact line check %s" % (len(K), ok), flush=True)

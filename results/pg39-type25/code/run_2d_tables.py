"""2d with plane-table constraints, 3 fixed-set cases."""
import json, os, sys, time
import groups, baer3d
from orbsat import OrbitInstance, run_solver, check_K
from pg import PG
from ff import GF
HERE = os.path.dirname(os.path.abspath(__file__))
F, P3, gens, reps = groups.build_all()
M, f = reps["2d_baer"]
inst3 = OrbitInstance(P3, [P3.perm_of_matrix(M, f)])
t0 = time.time()
Z, inst2, P2 = baer3d.planar_Z(PG(3, 9, F))
print("planar sigma-invariant sets: %s (%.1fs)" % ({k: len(v) for k, v in Z.items()}, time.time() - t0), flush=True)
json.dump({str(k): v for k, v in Z.items()}, open(os.path.join(HERE, "..", "certs", "baer_Z.json"), "w"))
fixed3 = [i for i, o in enumerate(inst3.orbs) if len(o) == 1]
fpts = {inst3.orbs[i][0]: i for i in fixed3}
sets = json.load(open(os.path.join(HERE, "..", "results", "baer_fixed_sets.json")))
by = {}
for S in sets:
    by.setdefault(len(S), []).append(S)
cases = [("S_empty", by[0][0], 260), ("S_point", by[1][0], 287), ("S_27", by[27][0], 287)]
timeout = int(sys.argv[1]) if len(sys.argv) > 1 else 600
for name, S, k in cases:
    Sset = set(S)
    units = [(fpts[p], 1 if p in Sset else 0) for p in sorted(fpts)]
    cnf = inst3.build_cnf(k, (2, 5), amo=True, units=units)
    tables = baer3d.plane_tables(Z, inst2, P2, P3, S)
    print("%s: allowed sections per plane: min %d max %d" % (name, min(len(r) for a, r in tables), max(len(r) for a, r in tables)), flush=True)
    baer3d.add_tables_full(cnf, inst3, P3, tables)
    # symmetry breaking: reference plane with empty fixed pattern; its 432 sections form ONE orbit of PGL(3,3) (baer_orbits.py)
    planes = baer3d.gf3_planes(P3)
    ref = next(i for i, (a, pts, fixed, basis) in enumerate(planes) if not (set(fixed) & Sset))
    a_ref, rows_ref = tables[ref]
    assert len(rows_ref) == 432, len(rows_ref)
    row0 = set(rows_ref[0])
    pts_ref = planes[ref][1]
    for o in sorted(set(inst3.oid[p] for p in pts_ref)):
        cnf.add(inst3.z[o] if inst3.orbs[o][0] in row0 else -inst3.z[o])
    print("reference plane %d: section fixed (%d points)" % (ref, len(row0)), flush=True)
    path = os.path.join(HERE, "..", "cnf", "2dT_%s_k%d_ref.cnf" % (name, k))
    cnf.to_dimacs(path, comments=["class 2d_baer k %d case %s tables refsection" % (k, name)])
    st, model, t, info = run_solver(path, "cadical", timeout)
    print("2dT %s k=%d : %s in %.1fs %s (vars %d clauses %d)" % (name, k, st, t, info, cnf.nvars, len(cnf.clauses)), flush=True)
    if st == "SATISFIABLE":
        tr = set(m for m in model if m > 0)
        K = sorted(p for i, v in enumerate(inst3.z) if v in tr for p in inst3.orbs[i])
        ok, bad = check_K(P3, K, (2, 5))
        print("  MODEL |K| = %d exact line check %s" % (len(K), ok), flush=True)

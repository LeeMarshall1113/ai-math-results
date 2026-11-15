"""LRAT certificates for the table-based cases:

  python certify_tables.py 2b      : class 2b, case (|F1|,|F2|,k) = (2,2,260), orbit CNF + plane tables of the 20 invariant planes
  python certify_tables.py 2d      : class 2d, the three fixed-set cases, orbit CNF + S units + plane tables + reference section
  python certify_tables.py baerS   : completeness of the 81 fixed sets S (40 variables, forbidden plane patterns, <= 3 points per Baer line)

Each CNF is written to certs/<tag>.cnf, proved with CaDiCaL LRAT, checked with cake_lpr, and recorded in certs/manifest.jsonl.
The lists of invariant planar sets used by the tables are certs/table_hom.json / table_baer.json (completeness: table_cert.py).
"""
import itertools
import json
import os
import sys

import baer3d
import cases as C
import groups
import plane_tables as T
from certify import cake_check, sha256
from ff import GF
from orbsat import OrbitInstance, run_solver
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def record(tag, cnf, status, t, info, extra=None):
    cnfp, lratp = os.path.join(CERTS, tag + ".cnf"), os.path.join(CERTS, tag + ".lrat")
    rec = dict(cls=tag, tag=tag, status=status, solver_seconds=round(t, 1), vars=cnf.nvars, clauses=len(cnf.clauses),
               cnf_sha256=sha256(cnfp), **info)
    rec.update(extra or {})
    if status == "UNSATISFIABLE":
        rec["lrat_bytes"], rec["lrat_sha256"] = os.path.getsize(lratp), sha256(lratp)
        v, ct, rc = cake_check(cnfp, lratp)
        rec["cake_lpr"], rec["cake_seconds"], rec["cake_rc"] = v, round(ct, 1), rc
    print(json.dumps(rec), flush=True)
    open(os.path.join(CERTS, "manifest.jsonl"), "a").write(json.dumps(rec) + "\n")


def run_2b():
    F, P3, gens, reps = groups.build_all()
    M, f = reps["2b_homology_2+2"]
    inst3 = OrbitInstance(P3, [P3.perm_of_matrix(M, f)])
    P2 = PG(3, 9, F)
    fam = {int(k): [tuple(x) for x in v] for k, v in json.load(open(os.path.join(CERTS, "table_hom.json"))).items()}
    fixed2 = set(i for i, v in enumerate(P2.points) if v[2] == 0) | {P2.pidx[(0, 0, 1)]}
    E = {1: (1, 0, 0, 0), 2: (0, 1, 0, 0), 3: (0, 0, 1, 0), 4: (0, 0, 0, 1)}
    l1 = [i for i, v in enumerate(P3.points) if v[2] == 0 and v[3] == 0]
    l2 = [i for i, v in enumerate(P3.points) if v[0] == 0 and v[1] == 0]
    opt = {p: o for p, o in enumerate(inst3.oid)}
    for case in C.cases_2b(P3):
        if "F2_2_k260" not in case["name"]:
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
        T.add_plane_tables(cnf, inst3, planes_rows)
        tag = "2bT_" + case["name"]
        cnfp, lratp = os.path.join(CERTS, tag + ".cnf"), os.path.join(CERTS, tag + ".lrat")
        cnf.to_dimacs(cnfp, comments=["class 2b case %s plane tables" % case["name"]])
        st, model, t, info = run_solver(cnfp, "cadical", 3000, lratp)
        record(tag, cnf, st, t, info, dict(k=case["k"]))


def run_2d():
    F, P3, gens, reps = groups.build_all()
    M, f = reps["2d_baer"]
    inst3 = OrbitInstance(P3, [P3.perm_of_matrix(M, f)])
    P2 = PG(3, 9, F)
    Z = {int(k): [tuple(x) for x in v] for k, v in json.load(open(os.path.join(CERTS, "table_baer.json"))).items()}
    inst2 = OrbitInstance(P2, [P2.perm_of_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1)])
    fixed3 = [i for i, o in enumerate(inst3.orbs) if len(o) == 1]
    fpts = {inst3.orbs[i][0]: i for i in fixed3}
    sets = json.load(open(os.path.join(HERE, "..", "results", "baer_fixed_sets.json")))
    by = {}
    for S in sets:
        by.setdefault(len(S), []).append(S)
    for name, S, k in (("S_empty", by[0][0], 260), ("S_point", by[1][0], 287), ("S_27", by[27][0], 287)):
        Sset = set(S)
        units = [(fpts[p], 1 if p in Sset else 0) for p in sorted(fpts)]
        cnf = inst3.build_cnf(k, (2, 5), amo=True, units=units)
        tables = baer3d.plane_tables(Z, inst2, P2, P3, S)
        baer3d.add_tables_full(cnf, inst3, P3, tables)
        planes = baer3d.gf3_planes(P3)
        ref = next(i for i, (a, pts, fixed, basis) in enumerate(planes) if not (set(fixed) & Sset))
        row0 = set(tables[ref][1][0])
        for o in sorted(set(inst3.oid[p] for p in planes[ref][1])):
            cnf.add(inst3.z[o] if inst3.orbs[o][0] in row0 else -inst3.z[o])
        tag = "2dT_%s_k%d" % (name, k)
        cnfp, lratp = os.path.join(CERTS, tag + ".cnf"), os.path.join(CERTS, tag + ".lrat")
        cnf.to_dimacs(cnfp, comments=["class 2d_baer k %d case %s tables refsection" % (k, name)])
        st, model, t, info = run_solver(cnfp, "cadical", 3000, lratp)
        record(tag, cnf, st, t, info, dict(k=k))


def run_baerS():
    """completeness of the 81 sets S: same CNF as baer_fixed.py plus blocking clauses"""
    import baer_fixed_cnf
    cnf, sols = baer_fixed_cnf.build()
    tag = "baerS_complete"
    cnfp, lratp = os.path.join(CERTS, tag + ".cnf"), os.path.join(CERTS, tag + ".lrat")
    cnf.to_dimacs(cnfp, comments=["Baer fixed sets S: %d solutions blocked" % len(sols)])
    st, model, t, info = run_solver(cnfp, "cadical", 3000, lratp)
    record(tag, cnf, st, t, info, dict(solutions_blocked=len(sols)))


if __name__ == "__main__":
    {"2b": run_2b, "2d": run_2d, "baerS": run_baerS}[sys.argv[1]]()

"""Completeness certificates of the section lists of the star cover (level 1: 10 cases, level 2: the split cubes).

  python cover_cert.py level1 [A0 A1 A2 C0 ... C6]       (default: all 10 cases)
  python cover_cert.py level2 [ID ...]                   (default: A1_c04 A2_c01 A2_c04 A2_c05 A2_c17 C0_c15 C2_c00 C6_c08)

A list of sections of a plane through P0 is COMPLETE when   CNF(plane)  +  one blocking clause per listed section   is UNSAT.
  CNF(plane) = Setup.plane_cnf: the 91 points of the plane, every line of the plane meets the set in 2 or 5 points, |section| = 26 or 35 (determined by the
               pencil of Sigma), the local counts (implied by counting, NOTES 5.1), the line types through P0 and the point units, plus (level 2) the section of the
               first plane as unit clauses.
  level 1   : the blocked sets are ALL homothety images (8 homotheties diag(1,1,1,lambda)) of the cube representatives of certs/cubes_<case><cls>.json: UNSAT says
              that every section of the plane L1 is a cube representative up to a homothety (the orbit expansion is repeated by an independent stdlib script,
              replay/rebuild_all.py, which rebuilds the same CNF byte for byte from the representatives).
  level 2   : the blocked sets are the children of certs/split_<ID>.json.
Blocking clauses are sorted by the sorted tuple of member points, literals by point index (canonical, so that the stdlib replay reproduces the file exactly).
Writes certs/<tag>.cnf/.lrat (text LRAT), checks with cake_lpr and appends a record to certs/manifest.jsonl.
"""
import json
import os
import sys
import time

from certify import cake_check, sha256
from orbsat import run_solver
from trivial_tree import Setup

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")
RES = os.path.join(HERE, "..", "results")
DEFAULT_L2 = ["A1_c04", "A2_c01", "A2_c04", "A2_c05", "A2_c17", "C0_c15", "C1_c04", "C2_c00", "C6_c08"]


def finish(cnf, tag, comment, extra):
    cnfp = os.path.join(CERTS, tag + ".cnf")
    lratp = os.path.join(CERTS, tag + ".lrat")
    cnf.to_dimacs(cnfp, comments=[comment])
    status, model, t, info = run_solver(cnfp, "cadical", 3000, lratp)
    rec = dict(kind="cover_complete", tag=tag, status=status, solver_seconds=round(t, 2), vars=cnf.nvars, clauses=len(cnf.clauses), cnf_sha256=sha256(cnfp),
               time=time.strftime("%Y-%m-%d %H:%M:%S"), **info, **extra)
    if status == "UNSATISFIABLE":
        rec["lrat_bytes"] = os.path.getsize(lratp)
        rec["lrat_sha256"] = sha256(lratp)
        v, ct, rc = cake_check(cnfp, lratp, heap=800, stack=300)
        rec["cake_lpr"], rec["cake_seconds"], rec["cake_rc"] = v, round(ct, 1), rc
    else:
        rec["flag"] = "NOT UNSAT: the list is NOT complete (a missed section exists)"
    print(json.dumps(rec), flush=True)
    with open(os.path.join(CERTS, "manifest.jsonl"), "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return rec


def level1(name):
    case, ci = name[0], int(name[1:])
    S = Setup(case, ci)
    cj = json.load(open(os.path.join(CERTS, "cubes_%s%d.json" % (case, ci))))
    L = cj["plane"]
    homs = [S.P.perm_of_matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, lam]], 0) for lam in range(1, 9)]
    blocked = set()
    for cube in cj["cubes"]:
        key = sorted(p for p, v in cube["section"] if v)
        for h in homs:
            blocked.add(tuple(sorted(h[p] for p in key)))
    assert len(blocked) == cj["sections"], (len(blocked), cj["sections"])
    cnf, x, pts = S.plane_cnf(L, {})
    free = [p for p in pts if p not in S.pt_units]
    for b in sorted(blocked):
        bs = set(b)
        cnf.add(*[(-x[p] if p in bs else x[p]) for p in free])
    return finish(cnf, "cover1_%s%d" % (case, ci), "level-1 completeness case %s class %d plane %d: %d sections = homothety orbits of %d cube representatives" % (
        case, ci, L, len(blocked), len(cj["cubes"])), dict(level=1, case=case, cls=ci, plane=L, sections_blocked=len(blocked), representatives=len(cj["cubes"])))


def level2(jid):
    job = None
    for fn in ("a_slice_jobs.json", "c_slice_jobs.json"):
        job = next((j for j in json.load(open(os.path.join(RES, fn))) if j["id"] == jid), None)
        if job is not None:
            break
    S = Setup(job["case"], job["cls"])
    fixed = {abs(u) - 1: int(u > 0) for u in job["units"]}
    sp = json.load(open(os.path.join(CERTS, "split_%s.json" % jid)))
    L2 = sp["plane"]
    cnf, x, pts = S.plane_cnf(L2, fixed)
    free = [p for p in pts if p not in fixed and p not in S.pt_units]
    blocked = set()
    for ch in sp["children"]:
        assert sorted(abs(u) - 1 for u in ch) == free
        blocked.add(tuple(sorted(abs(u) - 1 for u in ch if u > 0)))
    assert len(blocked) == sp["n"] == len(sp["children"])
    for b in sorted(blocked):
        bs = set(b)
        cnf.add(*[(-x[p] if p in bs else x[p]) for p in free])
    return finish(cnf, "cover2_%s" % jid, "level-2 completeness cube %s plane %d: %d children" % (jid, L2, len(blocked)),
                  dict(level=2, cube=jid, plane2=L2, sections_blocked=len(blocked), free_points=len(free), **({"via": "local"} if jid == "C1_c04" else {})))


def main():
    mode = sys.argv[1]
    names = sys.argv[2:]
    if mode == "level1":
        for n in names or ["A0", "A1", "A2", "C0", "C1", "C2", "C3", "C4", "C5", "C6"]:
            level1(n)
    else:
        for n in names or DEFAULT_L2:
            level2(n)


if __name__ == "__main__":
    main()

"""Certify the case-split classes:  python certify_cases.py FAMILY [--only name1,name2] [--timeout S]
FAMILY in {2a, 2b, 2d}.  For each case: CNF -> certs/<name>.cnf, CaDiCaL LRAT -> certs/<name>.lrat, cake_lpr check,
record in certs/manifest.jsonl, and the case descriptor in certs/cases.json."""
import argparse
import json
import os
import sys
import time

import cases as C
import groups
from certify import sha256, cake_check
from orbsat import OrbitInstance, run_solver
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("family")
    ap.add_argument("--only", default="")
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--no-lrat", action="store_true")
    a = ap.parse_args()
    F, P, gens, reps = groups.build_all()
    cs = C.FAMILIES[a.family](P)
    only = set(a.only.split(",")) if a.only else None
    casesfile = os.path.join(CERTS, "cases.json")
    known = json.load(open(casesfile)) if os.path.exists(casesfile) else {}
    insts = {}
    out = open(os.path.join(CERTS, "manifest.jsonl"), "a")
    for case in cs:
        if only and case["name"] not in only:
            continue
        cls = case["cls"]
        if cls not in insts:
            M, f = reps[cls]
            inst = OrbitInstance(P, [P.perm_of_matrix(M, f)])
            insts[cls] = (inst, {p: o for p, o in enumerate(inst.oid)})
        inst, opt = insts[cls]
        cnf = C.build_case_cnf(inst, case, opt)
        tag = case["name"]
        cnfp = os.path.join(CERTS, tag + ".cnf")
        lratp = os.path.join(CERTS, tag + ".lrat")
        cnf.to_dimacs(cnfp, comments=["class %s k %d case %s amo=True" % (cls, case["k"], tag)])
        status, model, t, info = run_solver(cnfp, "cadical", a.timeout, None if a.no_lrat else lratp)
        rec = dict(cls=cls, k=case["k"], tag=tag, status=status, solver_seconds=round(t, 2), vars=cnf.nvars,
                   clauses=len(cnf.clauses), point_orbits=len(inst.orbs), cnf_sha256=sha256(cnfp), case=True, **info)
        if status == "UNSATISFIABLE" and not a.no_lrat:
            rec["lrat_bytes"] = os.path.getsize(lratp)
            rec["lrat_sha256"] = sha256(lratp)
            verdict, ct, rc = cake_check(cnfp, lratp)
            rec["cake_lpr"] = verdict
            rec["cake_seconds"] = round(ct, 1)
            rec["cake_rc"] = rc
        if status == "SATISFIABLE":
            tr = set(m for m in model if m > 0)
            K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
            from orbsat import check_K
            ok, bad = check_K(P, K, (2, 5))
            rec["model_size"] = len(K)
            rec["line_check"] = ok
            json.dump(K, open(os.path.join(CERTS, "MODEL_" + tag + ".json"), "w"))
        known[tag] = case
        json.dump(known, open(casesfile, "w"))
        print(json.dumps(rec), flush=True)
        out.write(json.dumps(rec) + "\n")
        out.flush()


if __name__ == "__main__":
    main()

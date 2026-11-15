"""Completeness certificates for the planar enumerations:  CNF(planar instance, frame R, origin, lex-leader) plus one
blocking clause per enumerated solution is UNSAT (LRAT, checked by cake_lpr).  Hence the JSON solution lists are complete.

  python planar_cert.py [k R ...]         (default: 26 0 26 1 35 0 35 1)
"""
import json
import os
import sys
import time

import planar
from certify import cake_check, sha256
from orbsat import run_solver
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def main():
    args = [int(x) for x in sys.argv[1:]] or [26, 0, 26, 1, 35, 0, 35, 1]
    P = PG(3, 9)
    l0 = planar.line0_points(P)
    P1, reps = planar.pgammal2_orbits_on_5subsets()

    def lift(i):
        v = P1.points[i]
        return P.pidx[P.normalize((v[0], v[1], 0))]
    origin = P.pidx[(0, 0, 1)]
    out = open(os.path.join(CERTS, "manifest.jsonl"), "a")
    for j in range(0, len(args), 2):
        k, ri = args[j], args[j + 1]
        d = json.load(open(os.path.join(HERE, "..", "results", "planar_k%d_R%d.json" % (k, ri))))
        R1, osz = reps[ri]
        Rpts = set(lift(i) for i in R1)
        assert sorted(Rpts) == d["R"]
        H = planar.stabilizer_perms(P, Rpts)
        cnf, x, order = planar.build_instance(P, k, Rpts, l0, origin, H)
        for K in d["solutions"]:
            Ks = set(K)
            cnf.add(*[(-x[p] if p in Ks else x[p]) for p in order])
        tag = "planar_complete_k%d_R%d" % (k, ri)
        cnfp = os.path.join(CERTS, tag + ".cnf")
        lratp = os.path.join(CERTS, tag + ".lrat")
        cnf.to_dimacs(cnfp, comments=["planar completeness k %d R#%d blocked %d solutions" % (k, ri, len(d["solutions"]))])
        status, model, t, info = run_solver(cnfp, "cadical", 3000, lratp)
        rec = dict(cls="planar", k=k, tag=tag, status=status, solver_seconds=round(t, 2), vars=cnf.nvars,
                   clauses=len(cnf.clauses), cnf_sha256=sha256(cnfp), solutions_blocked=len(d["solutions"]), **info)
        if status == "UNSATISFIABLE":
            rec["lrat_bytes"] = os.path.getsize(lratp)
            rec["lrat_sha256"] = sha256(lratp)
            v, ct, rc = cake_check(cnfp, lratp)
            rec["cake_lpr"], rec["cake_seconds"], rec["cake_rc"] = v, round(ct, 1), rc
        print(json.dumps(rec), flush=True)
        out.write(json.dumps(rec) + "\n")
        out.flush()


if __name__ == "__main__":
    main()

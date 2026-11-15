"""Completeness certificates for the lists of involution-invariant planar (2,5)-sets used as plane tables:

  tablelemma_<h>_k<k> : orbit CNF of the planar involution h in PG(2,9) (|S| = k) plus one blocking clause per listed set
                        => UNSAT (LRAT + cake_lpr).   h in {hom (diag(1,1,-1)), baer (Frobenius)}, k in {26,35}.
Also dumps the lists to certs/table_<h>.json (used by the replay).
"""
import json
import os
import sys

import plane_tables as T
from certify import cake_check, sha256
from ff import GF
from orbsat import OrbitInstance, run_solver
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def main():
    F = GF(9)
    P2 = PG(3, 9, F)
    hs = {"hom": ([[1, 0, 0], [0, 1, 0], [0, 0, F.neg[1]]], 0), "baer": ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1)}
    out = open(os.path.join(CERTS, "manifest.jsonl"), "a")
    for h, (M, f) in hs.items():
        fam, inst = T.invariant_family(P2, M, f)
        json.dump({str(k): v for k, v in fam.items()}, open(os.path.join(CERTS, "table_%s.json" % h), "w"))
        for k in (26, 35):
            cnf = inst.build_cnf(k, (2, 5), amo=True)
            for X in fam[k]:
                Xs = set(X)
                cnf.add(*[(-inst.z[i] if inst.orbs[i][0] in Xs else inst.z[i]) for i in range(len(inst.orbs))])
            tag = "tablelemma_%s_k%d" % (h, k)
            cnfp, lratp = os.path.join(CERTS, tag + ".cnf"), os.path.join(CERTS, tag + ".lrat")
            cnf.to_dimacs(cnfp, comments=["table lemma %s k %d blocked %d sets" % (h, k, len(fam[k]))])
            st, model, t, info = run_solver(cnfp, "cadical", 3000, lratp)
            rec = dict(cls="tablelemma_" + h, k=k, tag=tag, status=st, solver_seconds=round(t, 1), vars=cnf.nvars,
                       clauses=len(cnf.clauses), cnf_sha256=sha256(cnfp), sets=len(fam[k]), **info)
            if st == "UNSATISFIABLE":
                rec["lrat_bytes"], rec["lrat_sha256"] = os.path.getsize(lratp), sha256(lratp)
                v, ct, rc = cake_check(cnfp, lratp)
                rec["cake_lpr"], rec["cake_seconds"], rec["cake_rc"] = v, round(ct, 1), rc
            print(json.dumps(rec), flush=True)
            out.write(json.dumps(rec) + "\n")
            out.flush()


if __name__ == "__main__":
    main()

"""Run the orbit-collapsed SAT instance for each prime-order class and each k in {260, 287}.

  python run_classes.py [--timeout 300] [--classes a,b,c] [--ks 260,287] [--lrat]
Writes results/class_runs.jsonl (one line per run) and certs/<class>_k<k>.{cnf,lrat} for UNSAT runs
when --lrat is given.
"""
import argparse
import json
import os
import sys
import time

import groups
from orbsat import OrbitInstance, run_solver, check_K
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--classes", default="")
    ap.add_argument("--ks", default="260,287")
    ap.add_argument("--lrat", action="store_true")
    ap.add_argument("--noamo", action="store_true")
    args = ap.parse_args()
    F, P, gens, reps = groups.build_all()
    names = args.classes.split(",") if args.classes else list(reps)
    ks = [int(x) for x in args.ks.split(",")]
    out = open(os.path.join(HERE, "..", "results", "class_runs.jsonl"), "a")
    for name in names:
        M, f = reps[name]
        perm = P.perm_of_matrix(M, f)
        inst = OrbitInstance(P, [perm])
        summ = inst.summary()
        for k in ks:
            cnf = inst.build_cnf(k, (2, 5), amo=not args.noamo)
            tag = "%s_k%d" % (name.replace("(", "").replace(")", "").replace(",", "_").replace("+", "p"), k)
            cnfp = os.path.join(HERE, "..", "cnf", tag + ".cnf")
            cnf.to_dimacs(cnfp, comments=["class %s k %d allowed 2,5 amo=%s" % (name, k, not args.noamo)])
            lr = os.path.join(HERE, "..", "certs", tag + ".lrat") if args.lrat else None
            status, model, t, info = run_solver(cnfp, "cadical", args.timeout, lr)
            rec = dict(cls=name, k=k, status=status, seconds=round(t, 2), vars=cnf.nvars,
                       clauses=len(cnf.clauses), orbits=summ["point_orbits"], lorbits=summ["line_orbits"], **info)
            if status == "SATISFIABLE":
                tr = set(m for m in model if m > 0)
                K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
                ok, bad = check_K(P, K, (2, 5))
                rec["model_size"] = len(K)
                rec["line_check"] = ok
                rec["K"] = K
            print(json.dumps({k_: v for k_, v in rec.items() if k_ != "K"}), flush=True)
            out.write(json.dumps(rec) + "\n")
            out.flush()


if __name__ == "__main__":
    main()

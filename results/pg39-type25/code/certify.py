"""Produce and check UNSAT certificates (CaDiCaL 3.0.1 LRAT + cake_lpr) for orbit instances.

  python certify.py CLASS K [CLASS K ...]            (class names as in groups.build_all())

For each pair: writes certs/<tag>.cnf and certs/<tag>.lrat, runs cake_lpr (copying both files
to WSL-local disk first), and appends a record to certs/manifest.jsonl with the sha256 of the
CNF and of the LRAT file, the solver time, and the cake_lpr verdict.
"""
import hashlib
import json
import os
import subprocess
import sys
import time

import groups
from orbsat import OrbitInstance, run_solver, wsl_path, WSL_CAKE
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tag_of(name, k):
    return "%s_k%d" % (name.replace("(", "").replace(")", "").replace(",", "_").replace("+", "p"), k)


def cake_check(cnf, lrat, heap=1500, stack=500):
    """copy to WSL-local disk and run cake_lpr; returns (verdict_line, seconds)"""
    cmd = ("W=$HOME/pg39work/cake.$$; mkdir -p $W; cp '%s' $W/f.cnf && cp '%s' $W/f.lrat && "
           "%s --CML_HEAP_SIZE=%d --CML_STACK_SIZE=%d $W/f.cnf $W/f.lrat; rc=$?; rm -rf $W; exit $rc"
           % (wsl_path(cnf), wsl_path(lrat), WSL_CAKE, heap, stack))
    t = time.time()
    r = subprocess.run(["wsl", "-d", "Ubuntu-24.04", "-e", "bash", "-lc", cmd], capture_output=True, text=True)
    out = (r.stdout.replace("\x00", "") + r.stderr.replace("\x00", ""))
    lines = [l for l in out.splitlines() if l.startswith("s ") or "ERROR" in l or "error" in l]
    return (lines[-1] if lines else out.strip()[-200:]), time.time() - t, r.returncode


def certify(P, reps, name, k, timeout=900):
    M, f = reps[name]
    perm = P.perm_of_matrix(M, f)
    inst = OrbitInstance(P, [perm])
    cnf = inst.build_cnf(k, (2, 5), amo=True)
    tag = tag_of(name, k)
    cnfp = os.path.join(CERTS, tag + ".cnf")
    lratp = os.path.join(CERTS, tag + ".lrat")
    cnf.to_dimacs(cnfp, comments=["class %s k %d allowed 2,5 amo=True" % (name, k)])
    status, model, t, info = run_solver(cnfp, "cadical", timeout, lratp)
    rec = dict(cls=name, k=k, tag=tag, status=status, solver_seconds=round(t, 2), vars=cnf.nvars,
               clauses=len(cnf.clauses), point_orbits=len(inst.orbs), cnf_sha256=sha256(cnfp), **info)
    if status == "UNSATISFIABLE":
        rec["lrat_bytes"] = os.path.getsize(lratp)
        rec["lrat_sha256"] = sha256(lratp)
        verdict, ct, rc = cake_check(cnfp, lratp)
        rec["cake_lpr"] = verdict
        rec["cake_seconds"] = round(ct, 1)
        rec["cake_rc"] = rc
    return rec


def main():
    args = sys.argv[1:]
    F, P, gens, reps = groups.build_all()
    out = open(os.path.join(CERTS, "manifest.jsonl"), "a")
    os.makedirs(os.path.expanduser("~"), exist_ok=True)
    for i in range(0, len(args), 2):
        name, k = args[i], int(args[i + 1])
        rec = certify(P, reps, name, k)
        print(json.dumps(rec), flush=True)
        out.write(json.dumps(rec) + "\n")
        out.flush()


if __name__ == "__main__":
    main()

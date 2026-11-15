"""Level-1 cube-and-conquer driver for the star formulation (the proposed full run).

  python cube_run.py CASE CLASS [--only i,j,k] [--lrat] [--timeout S] [--keep-lrat]

CASE A (k = 287, P0 in K, 3 classes of size 26) or C (k = 260, P0 in K, 7 classes of size 35).
Cubes = the sections of the first plane of size 26 (L1, L2 or L3 of the normalised triple) up to the 8 homotheties
diag(1,1,1,lambda) (one representative per orbit).  For every cube: base CNF + the 91 point units of the section -> CaDiCaL;
with --lrat the proof is written to the WSL disk, checked by cake_lpr, and only its sha256/size are kept.
Results are appended to results/cube_run_<case><class>.jsonl (resumable).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

import full
from orbsat import wsl_path, WSL_CADICAL, WSL_CAKE
from trivial_tree import Setup, write_cube_cnf, run_win

HERE = os.path.dirname(os.path.abspath(__file__))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 22), b""):
            h.update(ch)
    return h.hexdigest()


def cubes_for(S, cap=20000):
    v1, v2, v3 = S.tri
    Ls = [S.quotient_line_through(v1, v2), S.quotient_line_through(v1, v3), S.quotient_line_through(v2, v3)]
    L = next(L for L in Ls if S.plane_size(L) == 26)
    secs = S.sections(L, {}, cap=cap)
    assert len(secs) < cap
    homs = [S.P.perm_of_matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, lam]], 0) for lam in range(1, 9)]
    reps, seen = [], set()
    for sec in secs:
        key = tuple(sorted(p for p, v in sec.items() if v))
        if key in seen:
            continue
        orbit = set(tuple(sorted(h[p] for p in key)) for h in homs)
        seen |= orbit
        reps.append((sec, len(orbit)))
    return L, reps, len(secs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("ci", type=int)
    ap.add_argument("--only", default="")
    ap.add_argument("--lrat", action="store_true")
    ap.add_argument("--timeout", type=int, default=7200)
    a = ap.parse_args()
    S = Setup(a.case, a.ci)
    L, reps, nsec = cubes_for(S)
    print("case %s class %d (|Aut| %d): plane %d size %d: %d sections, %d homothety orbits" % (
        a.case, a.ci, S.cl["aut_order"], L, S.plane_size(L), nsec, len(reps)), flush=True)
    base = os.path.join(HERE, "..", "cnf", "full_%s_c%d.cnf" % (a.case, a.ci))
    if not os.path.exists(base):
        full.make_case(a.case, a.ci, base)
    out_path = os.path.join(HERE, "..", "results", "cube_run_%s%d.jsonl" % (a.case, a.ci))
    done = set()
    if os.path.exists(out_path):
        done = {json.loads(l)["cube"] for l in open(out_path)}
    only = set(int(x) for x in a.only.split(",")) if a.only else None
    for i, (sec, orb) in enumerate(reps):
        if i in done or (only is not None and i not in only):
            continue
        lits = [((p + 1) if v else -(p + 1)) for p, v in {**S.pt_units, **sec}.items()]
        cnf = os.path.join(HERE, "..", "cnf", "cube_%s%d_%d.cnf" % (a.case, a.ci, i))
        write_cube_cnf(base, lits, cnf)
        rec = dict(case=a.case, cls=a.ci, cube=i, orbit_size=orb, plane=L)
        t = time.time()
        if not a.lrat:
            st, secs, conf = run_win(cnf, a.timeout)
            rec.update(status=st, seconds=round(secs, 1), conflicts=conf)
        else:
            W = "$HOME/pg39work"
            cmd = ("mkdir -p %s; cp '%s' %s/c.cnf && %s -t %d --lrat=true --binary=false --congruence=false %s/c.cnf %s/c.lrat > %s/c.out; rc=$?; "
                   "echo rc=$rc; if [ $rc -eq 20 ]; then sha256sum %s/c.cnf %s/c.lrat; ls -l %s/c.lrat | awk '{print $5}'; "
                   "%s --CML_HEAP_SIZE=3000 --CML_STACK_SIZE=1000 %s/c.cnf %s/c.lrat | tail -1; fi; rm -f %s/c.lrat %s/c.cnf"
                   % (W, wsl_path(cnf), W, WSL_CADICAL, a.timeout, W, W, W, W, W, W, WSL_CAKE, W, W, W, W))
            r = subprocess.run(["wsl", "-d", "Ubuntu-24.04", "-e", "bash", "-lc", cmd], capture_output=True, text=True)
            rec.update(status="UNSAT" if "rc=20" in r.stdout else ("SAT" if "rc=10" in r.stdout else "OTHER"), seconds=round(time.time() - t, 1),
                       log=r.stdout.replace(chr(0), "")[-600:])
        rec["cnf_sha256"] = sha256(cnf)
        os.remove(cnf)
        print(json.dumps({k: v for k, v in rec.items() if k != "log"}), flush=True)
        open(out_path, "a").write(json.dumps(rec) + "\n")
        if rec["status"] == "SAT":
            print("!!! SAT cube: a candidate set exists -- run replay/check_set.py on the model", flush=True)
            break


if __name__ == "__main__":
    main()

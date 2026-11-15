"""Speed-check driver (runs inside a lease; at most 4 solver threads).

  python bench_driver.py JOBS.json [--workers 4]

JOBS.json : list of {name, cnf (basename in pg39/cnf), solver: cad|kis, opts: [...], tlim: seconds,
                     proof: null|lrat|lratbin|drat, check: null|cake|drattrim}
cad = CaDiCaL 3.0.1 in WSL (CNF copied to WSL-local scratch), kis = kissat 4.0.4 (Windows build).
Records (results/bench.jsonl): solver CPU seconds (from the solver's own statistics), wall seconds, conflicts, status,
proof bytes, and for check=cake the cake_lpr verdict/time (one cake_lpr at a time).
Scratch is WSL-local ~/pg39work/bench ; proofs are deleted right after checking.
"""
import json
import os
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
CNFD = os.path.join(ROOT, "cnf")
OUT = os.path.join(ROOT, "results", "bench.jsonl")
CAD = "/home/user/ext/cadical/build/cadical"
CAKE = "/home/user/lrchk/bin/cake_lpr"
DRATTRIM = "/home/user/erdos580/drat-trim/drat-trim"
KIS = r"C:\Users\user\src\build-kissat\rel\kissat.exe"
WS = "/home/user/pg39work/bench"
cake_lock = threading.Lock()
rec_lock = threading.Lock()


def wsl(cmd, timeout=None):
    return subprocess.run(["wsl", "-d", "Ubuntu-24.04", "-e", "bash", "-lc", cmd], capture_output=True, text=True, timeout=timeout)


def to_wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def run(job):
    name, cnf, solver, opts, tlim = job["name"], job["cnf"], job["solver"], job.get("opts", []), job["tlim"]
    proof, check = job.get("proof"), job.get("check")
    t0 = time.time()
    rec = dict(name=name, cnf=cnf, solver=solver, opts=" ".join(opts), tlim=tlim, proof=proof, check=check, start=time.strftime("%H:%M:%S"))
    pf = "%s/%s.proof" % (WS, name)
    if solver == "cad":
        args = [CAD, "-t", str(tlim)] + opts
        if proof in ("lrat", "lratbin"):
            args += ["--lrat=true", "--binary=%s" % ("true" if proof == "lratbin" else "false")]
        args += ["%s/%s" % (WS, cnf)]
        if proof:
            args += [pf]
        cmd = " ".join(args) + " > %s/%s.out 2>&1; echo rc=$?" % (WS, name)
        r = wsl(cmd, timeout=tlim + 600)
        rc = int(re.search(r"rc=(\d+)", r.stdout).group(1)) if re.search(r"rc=(\d+)", r.stdout) else -1
        outtxt = wsl("cat %s/%s.out" % (WS, name)).stdout
        m = re.search(r"c total process time since initialization:\s+([\d.]+)", outtxt)
        cpu = float(m.group(1)) if m else None
        m2 = re.search(r"c conflicts:\s+(\d+)", outtxt)
    else:
        args = [KIS, "--time=%d" % tlim] + opts + [os.path.join(CNFD, cnf)]
        if proof == "drat":
            args += [os.path.join(CNFD, name + ".drat"), "-f"] if False else []
        p = subprocess.run(args, capture_output=True, text=True, timeout=tlim + 600)
        rc = p.returncode
        outtxt = p.stdout
        m = re.search(r"c process-time:\s+\S+\s+([\d.]+)\s+seconds", outtxt)
        cpu = float(m.group(1)) if m else None
        m2 = re.search(r"^c conflicts:\s+(\d+)", outtxt, re.M)
    wall = time.time() - t0
    rec.update(rc=rc, status={10: "SAT", 20: "UNSAT"}.get(rc, "LIMIT/OTHER"), wall_s=round(wall, 1), cpu_s=cpu, conflicts=int(m2.group(1)) if m2 else None)
    if solver == "cad" and proof:
        sz = wsl("stat -c %%s %s 2>/dev/null || echo 0" % pf).stdout.strip()
        rec["proof_bytes"] = int(sz or 0)
        if rc == 20 and check == "cake":
            with cake_lock:
                tc = time.time()
                r = wsl("%s --CML_HEAP_SIZE=3000 --CML_STACK_SIZE=1000 %s/%s %s | tail -1" % (CAKE, WS, cnf, pf))
                rec["cake_verdict"] = r.stdout.replace(chr(0), "").strip()
                rec["cake_s"] = round(time.time() - tc, 1)
        if rc == 20 and check == "drattrim":
            tc = time.time()
            r = wsl("%s %s/%s %s -L %s.lrat -t %d | tail -3" % (DRATTRIM, WS, cnf, pf, pf, job.get("check_tlim", 1200)), timeout=job.get("check_tlim", 1200) + 600)
            rec["drattrim_tail"] = r.stdout.replace(chr(0), "").strip()[-200:]
            rec["drattrim_s"] = round(time.time() - tc, 1)
        wsl("rm -f %s %s.lrat %s/%s.out" % (pf, pf, WS, name))
    else:
        wsl("rm -f %s/%s.out" % (WS, name))
    with rec_lock:
        open(OUT, "a").write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)
    return rec


def main():
    jobs = json.load(open(sys.argv[1]))
    nw = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4
    wsl("mkdir -p %s" % WS)
    for c in sorted(set(j["cnf"] for j in jobs)):
        wsl("[ -f %s/%s ] || cp '%s' %s/%s" % (WS, c, to_wsl(os.path.join(CNFD, c)), WS, c), timeout=600)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        list(ex.map(run, jobs))


if __name__ == "__main__":
    main()

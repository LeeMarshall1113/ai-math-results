"""Modal driver for the level-1 cubes of pg39 (NOT launched without the coordinator's go-ahead; see NOTES section 5.6).

One container per cube:  unpack the gzipped base CNF -> append the 92 cube units -> CaDiCaL 3.0.1 (binary LRAT) -> cake_lpr -> return ONLY
{verdict, sha256(cnf), sha256(lrat), lrat bytes, cpu/wall seconds, conflicts, cake seconds}.  The LRAT never leaves the container and is never downloaded.

  python -m modal run modal_cubes.py::calibrate                          one cube (CAL_b_A1s1, about $0.03): Modal core speed, bit-identical search check
  python -m modal run modal_cubes.py::main --ids A0_c00,A0_c01 --budget 1  named jobs
  python -m modal run modal_cubes.py::main --ids A --budget 5            all cubes of the A cases (prefix A, A0, C, C3 ...; or @file with one id per line)

Spend control (all of it local, before any container starts):
  * the metered cost of the month is read with `modal billing summary`; the run is refused when metered + --budget > --cap (default 48 dollars);
  * --budget is the most this invocation may add: results are consumed as they finish and when the running estimate (container seconds x rate)
    exceeds it the consumer stops and the app exits, which cancels the inputs that have not started; containers already running are
    bounded by the 3 h container timeout (at most 0.24 dollars each) and by PG39_MAXC (default 20 containers at once);
  * a ledger line per cube is appended to results/modal_ledger.jsonl (estimate per cube, measured seconds).
Price: 1 physical core 0.0473 $/h + 4 GiB x 0.008 $/GiB-h = 0.0793 $/container-hour (modal billing rates, 2026-11-08).  Disk: an explicit disk request is billed by
raising the memory request 20:1 and the API minimum is 512 GiB (25.6 GiB of memory, 3.6x the price), so NO ephemeral_disk is requested; the default quota of 512 GiB is free.
"""
import hashlib
import json
import os
import re
import subprocess
import sys
import time

import modal

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
STAGE = os.path.join(ROOT, "modal_stage")           # bin/cadical, bin/cake_lpr, bases/<name>.cnf.gz (see modal_stage_prep.py)

CPU_RATE = 0.0473        # $ / physical core / hour
MEM_RATE = 0.008         # $ / GiB / hour
CORES, MEM_MIB = 1.0, 4096
BASE_METERED = 16.95     # metered cost of the month before the first pg39 Modal run (modal billing summary, 2026-11-08 13:05); billing lags, so the ledger is added to it
RATE_PER_H = CORES * CPU_RATE + MEM_MIB / 1024.0 * MEM_RATE
MAXC = min(int(os.environ.get("PG39_MAXC", "20")), 90)   # Modal free plan: 100 containers per workspace, keep headroom

app = modal.App("pg39-cubes")
image = (modal.Image.from_registry("ubuntu:24.04", add_python="3.12")
         .add_local_dir(os.path.join(STAGE, "bin"), "/opt/bin", copy=True)
         .add_local_dir(os.path.join(STAGE, "bases"), "/opt/bases", copy=True))


RATE_BIG = CORES * CPU_RATE + 8 * MEM_RATE        # the 8 GiB variant for proofs whose cake_lpr check needs a bigger heap (C1_c04: "CakeML heap space exhausted" at 3000 MB)


@app.function(image=image, cpu=CORES, memory=MEM_MIB, timeout=3 * 3600, retries=0, max_containers=MAXC)      # NO ephemeral_disk: the API minimum is 512 GiB = a 25.6 GiB memory charge
def run_cube(job: dict) -> dict:
    return _run(job, job.get("cake_heap", 3000))


@app.function(image=image, cpu=CORES, memory=8192, timeout=4 * 3600, retries=0, max_containers=MAXC)
def run_cube_big(job: dict) -> dict:
    return _run(job, job.get("cake_heap", 6500))


def _run(job, heap):
    """job = {id, kind: cube|control, base (no suffix), units: [signed ints], tlim, lrat: 'bin'|'text'|'none'}"""
    import gzip
    import shutil
    t0 = time.time()
    for b in ("cadical", "cake_lpr"):
        os.chmod("/opt/bin/" + b, 0o755)               # files added from Windows carry no exec bit
    w = "/tmp/w"
    os.makedirs(w, exist_ok=True)
    base_gz = "/opt/bases/%s.cnf.gz" % job["base"]
    cnf, lrat = "%s/c.cnf" % w, "%s/c.lrat" % w
    with gzip.open(base_gz, "rt") as fh, open(cnf, "w", newline="\n") as o:
        line = fh.readline()
        while line.startswith("c"):
            o.write(line)
            line = fh.readline()
        parts = line.split()
        o.write("p cnf %s %d\n" % (parts[2], int(parts[3]) + len(job["units"])))
        shutil.copyfileobj(fh, o, 1 << 22)
        for u in job["units"]:
            o.write("%d 0\n" % u)

    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for ch in iter(lambda: f.read(1 << 24), b""):
                h.update(ch)
        return h.hexdigest()
    rec = dict(id=job["id"], cnf_sha256=sha(cnf), disk_free_gb=round(shutil.disk_usage(w).free / 2 ** 30, 1))
    use_lrat = job.get("lrat", "bin") in ("bin", "text")
    cmd = ["/opt/bin/cadical", "-t", str(job.get("tlim", 7200)), "--congruence=false"] + job.get("opts", [])
    if use_lrat:
        cmd += ["--lrat=true", "--binary=%s" % ("true" if job.get("lrat", "bin") == "bin" else "false")]
    cmd += [cnf] + ([lrat] if use_lrat else [])
    ts = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    rec["solver_wall"] = round(time.time() - ts, 1)
    m = re.search(r"c total process time since initialization:\s+([\d.]+)", p.stdout)
    rec["solver_cpu"] = float(m.group(1)) if m else None
    m = re.search(r"c conflicts:\s+(\d+)", p.stdout)
    rec["conflicts"] = int(m.group(1)) if m else None
    m = re.search(r"c maximum resident set size of process:\s+([\d.]+)", p.stdout)
    rec["solver_rss_mb"] = float(m.group(1)) if m else None
    rec["rc"] = p.returncode
    if p.returncode == 20 and use_lrat:
        rec["lrat_bytes"] = os.path.getsize(lrat)
        rec["lrat_sha256"] = sha(lrat)
        tc = time.time()
        q = subprocess.run(["/opt/bin/cake_lpr", "--CML_HEAP_SIZE=%d" % heap, "--CML_STACK_SIZE=1000", cnf, lrat], capture_output=True, text=True)
        rec["cake_heap_mb"] = heap
        rec["cake_wall"] = round(time.time() - tc, 1)
        lines = [l for l in (q.stdout.replace(chr(0), "") + q.stderr.replace(chr(0), "")).splitlines() if l.strip()]
        rec["verdict"] = lines[-1] if lines else "(no output)"
        rec["status"] = "UNSAT_VERIFIED" if rec["verdict"] == "s VERIFIED UNSAT" else "UNSAT_CAKE_FAILED"
    elif p.returncode == 10:
        rec["status"] = "SAT"
        pos = [int(v) for line in p.stdout.splitlines() if line.startswith("v ") for v in line[2:].split()]
        rec["model"] = pos[:2000]                       # the first 820 variables are the point variables (point p <-> variable p+1)
    elif p.returncode == 0:
        rec["status"] = "TIMEOUT"
    else:
        rec["status"] = "ERROR"
        rec["tail"] = (p.stdout + p.stderr)[-300:]
    rec["total_wall"] = round(time.time() - t0, 1)
    return rec


def container_cost(seconds, big=False):
    return seconds / 3600.0 * (RATE_BIG if big else RATE_PER_H)


def metered_cost():
    """metered cost of the current month from `modal billing summary` (dollars) or None"""
    try:
        out = subprocess.run([sys.executable, "-m", "modal", "billing", "summary"], capture_output=True, text=True, timeout=120).stdout
        m = re.search(r"Metered Cost:[\s|]*\$?([\d.]+)", out)        # plain and boxed table formats
        return float(m.group(1)) if m else None
    except Exception:
        return None


def ledger_spend():
    """sum of the per-cube cost estimates of all earlier pg39 Modal runs (results/modal_ledger.jsonl)"""
    tot = 0.0
    path = os.path.join(ROOT, "results", "modal_ledger.jsonl")
    if os.path.exists(path):
        for l in open(path):
            try:
                tot += json.loads(l).get("est_cost_usd", 0.0)
            except Exception:
                pass
    return tot


def select(js, ids):
    """ids: comma list of job ids / class prefixes (A, A0, C, C3), or @file"""
    if not ids:
        return []
    want = [l.strip() for l in open(ids[1:]) if l.strip()] if ids.startswith("@") else [x for x in ids.split(",") if x]
    return [j for j in js if any(j["id"] == w or re.fullmatch(re.escape(w) + r"\d?_c\d\d", j["id"]) for w in want)]


def go(js, budget, cap, note, big=False):
    m_meter = metered_cost()
    m_led = BASE_METERED + ledger_spend()
    m0 = None if m_meter is None else max(m_meter, m_led)         # billing lags: never trust a figure below base + own ledger
    print("metered cost this month: %s ; base + ledger %.2f ; used %s ; this run may add at most %.2f ; cap %.2f" % (m_meter, m_led, m0, budget, cap), flush=True)
    if m0 is None or m0 + budget > cap:
        print("REFUSED: metered cost unknown or max(metered, base + ledger) + budget > cap")
        return []
    ledger = open(os.path.join(ROOT, "results", "modal_ledger.jsonl"), "a")
    spent, res = 0.0, []
    print("%d jobs, rate %.4f $/container-hour, at most %d containers at once" % (len(js), RATE_PER_H, MAXC), flush=True)
    fn = run_cube_big if big else run_cube
    for r in fn.map(js, order_outputs=False, return_exceptions=True):
        if isinstance(r, Exception):
            print("EXC", repr(r)[:300], flush=True)
            continue
        c = container_cost(r.get("total_wall", 0), big)
        spent += c
        r["est_cost_usd"] = round(c, 4)
        r["note"] = note
        r["time"] = time.strftime("%Y-%m-%d %H:%M:%S")
        ledger.write(json.dumps({k: v for k, v in r.items() if k != "model"}) + "\n")
        ledger.flush()
        res.append(r)
        print(r["id"], r.get("status"), "cpu %s" % r.get("solver_cpu"), "confl %s" % r.get("conflicts"), "lrat %s" % r.get("lrat_bytes"),
              "cake %s" % r.get("cake_wall"), "cost $%.3f (run total $%.2f)" % (c, spent), flush=True)
        if r.get("status") == "SAT":
            print("!!! SAT: model saved, the exact check is done locally (replay/check_set.py)", flush=True)
            json.dump(r, open(os.path.join(ROOT, "results", "MODAL_SAT_%s.json" % r["id"]), "w"))
        if spent > budget:
            print("BUDGET EXCEEDED (%.2f > %.2f): the consumer stops, the app exit cancels what has not started" % (spent, budget), flush=True)
            break
    print("run total (estimate) $%.2f ; metered before %s ; check `modal billing summary` later (billing lags)" % (spent, m0), flush=True)
    return res


@app.local_entrypoint()
def calibrate():
    js = [j for j in json.load(open(os.path.join(ROOT, "results", "modal_jobs.json"))) if j["id"] == "CAL_b_A1s1"]
    res = go(js, budget=1.0, cap=48.0, note="calibration")
    for r in res:
        e = js[0]["expect"]
        print("expected vs got: conflicts %s / %s ; lrat bytes %s / %s ; lrat sha256 equal: %s ; cnf sha256 equal: %s ; local CPU 345.1 s vs Modal %s s" % (
            e["conflicts"], r.get("conflicts"), e["lrat_bytes"], r.get("lrat_bytes"), e["lrat_sha256"] == r.get("lrat_sha256"),
            e["cnf_sha256"] == r.get("cnf_sha256"), r.get("solver_cpu")))


@app.local_entrypoint()
def main(ids: str = "", budget: float = 5.0, cap: float = 48.0, note: str = "", tlim: int = 0, big: bool = False):
    js = select(json.load(open(os.path.join(ROOT, "results", "modal_jobs.json"))), ids)
    if not js:
        print("no jobs selected (give --ids)")
        return
    if tlim:                                   # override of the solver limit (probes): seconds
        for j in js:
            j["tlim"] = tlim
    go(js, budget=budget, cap=cap, note=note or ids, big=big)

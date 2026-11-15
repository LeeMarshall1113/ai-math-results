#!/usr/bin/env python3
"""WSL-side driver of the A slice (runs inside the broker lease, python3 stdlib only).

  python3 a_slice_driver.py JOBS.json [--workers 4] [--cap-gb 17] [--tlim 7200] [--only ID1,ID2] [--dry] [--out NAME]

For every job (results/a_slice_jobs.json, written by a_slice_prep.py):
  cube    : CNF = base + units  ->  CaDiCaL 3.0.1 --lrat=true --binary=true --congruence=false  (wall limit 2 h)
            UNSAT -> cake_lpr (one at a time, heap 2 GB, retried at 3.6 GB; a rejected LRAT <= 10 GB is kept in ~/pg39work/quarantine) -> sha256 + size + verdict recorded -> LRAT, CNF, log DELETED at once
            SAT   -> model decoded, exact check of all 7462 lines, model saved, whole batch stopped
            2 h   -> status TIMEOUT_2h (to be split to level 2 by split_level2.py, not run on)
  control : SAT control (Hermitian surface), model checked exactly
Scratch (LRAT + CNFs) is capped (default 17 GB, the coordinator's limit is 20 GB): when the sum of the running files exceeds the cap the
largest running LRAT is killed (job deferred, re-run alone at the end with the full cap; if it exceeds the cap alone it is SKIPPED_SCRATCH).
Records are appended to results/a_slice.jsonl and (cubes) certs/manifest.jsonl.  Resumable.  If the GAMING flag appears only one worker runs.
"""
import hashlib
import json
import os
import re
import subprocess
import sys
import threading
import time
from collections import deque

ROOT = "/mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39"
SCR = os.path.expanduser("~/pg39work/" + (sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else "a_slice"))
CADICAL = "/home/user/ext/cadical/build/cadical"
CAKE = "/home/user/lrchk/bin/cake_lpr"
GAMING = "/mnt/c/Users/user/.compute-broker/GAMING"
args = sys.argv[1:]
JOBS = args[0]
DRY = "--dry" in args           # test run: records go to results/a_slice_dry*.jsonl, never to the manifest
OUTNAME = args[args.index("--out") + 1] if "--out" in args else "a_slice"        # results/<OUTNAME>.jsonl (e.g. c_slice for the local C run)
OUT = ROOT + ("/results/%s_dry.jsonl" % OUTNAME if DRY else "/results/%s.jsonl" % OUTNAME)
MAN = ROOT + ("/results/a_slice_dry_manifest.jsonl" if DRY else "/certs/manifest.jsonl")
NW = int(args[args.index("--workers") + 1]) if "--workers" in args else 4
CAP = int(float(args[args.index("--cap-gb") + 1]) * 2 ** 30) if "--cap-gb" in args else 17 * 2 ** 30
TLIM = int(args[args.index("--tlim") + 1]) if "--tlim" in args else 7200
ONLY = set(args[args.index("--only") + 1].split(",")) if "--only" in args else None
CAKE_HEAP = 2000          # cake_lpr peak RSS 2.0 GB; heap 1600 MB also verified the calibration proof (speed check, NOTES 5.5); retried at CAKE_HEAP_RETRY if it fails
CAKE_HEAP_RETRY = 3600
CAKE_HEAP_RETRY2 = 5000   # only when a single worker runs (nothing else resident): solver 0.7 GB + 5.0 GB heap stays under the 6 GB WSL cap
# options of the full local regeneration (2026-11-11); all optional, the defaults reproduce the earlier runs
MINFREE = int(float(args[args.index("--min-free-gb") + 1]) * 2 ** 30) if "--min-free-gb" in args else 0         # stop launching when free space on C: is lower
PREEMPT = args[args.index("--preempt-file") + 1] if "--preempt-file" in args else None                           # lease-run.sh preempt flag: stop launching
STOP_AFTER = float(args[args.index("--stop-after-min") + 1]) * 60 if "--stop-after-min" in args else None        # slice length: stop launching new jobs after this
T_START = time.time()
stop_reason = []

lock = threading.Lock()
rec_lock = threading.Lock()
cake_lock = threading.Lock()
stop_flag = threading.Event()
done_evt = threading.Event()
running = {}
deferred = []
peak = [0]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 24), b""):
            h.update(ch)
    return h.hexdigest()


def mem_available_mb():
    with open("/proc/meminfo") as fh:
        for line in fh:
            if line.startswith("MemAvailable"):
                return int(line.split()[1]) // 1024
    return 0


def record(rec, manifest=False):
    with rec_lock:
        with open(OUT, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        if manifest:
            with open(MAN, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
    log("RECORD", json.dumps({k: v for k, v in rec.items() if k not in ("lrat_sha256", "cnf_sha256")}))


def write_cube_cnf(base, units, out):
    with open(base, "r") as fh, open(out, "w", newline="\n") as o:
        line = fh.readline()
        while line.startswith("c"):
            o.write(line)
            line = fh.readline()
        parts = line.split()
        assert parts[0] == "p"
        o.write("p cnf %s %d\n" % (parts[2], int(parts[3]) + len(units)))
        while True:
            chunk = fh.read(1 << 22)
            if not chunk:
                break
            o.write(chunk)
        for u in units:
            o.write("%d 0\n" % u)


def parse_out(path):
    d = {}
    model = []
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith("v "):
                    model.extend(int(x) for x in line[2:].split())
                elif line.startswith("c conflicts:"):
                    d["conflicts"] = int(line.split()[2])
                elif line.startswith("c total process time since initialization:"):
                    d["solver_cpu_seconds"] = float(line.split()[-2])
                elif line.startswith("c maximum resident set size of process:"):
                    d["max_rss_mb"] = float(line.split()[-2])
    except OSError:
        pass
    return d, [m for m in model if m != 0]


def exact_check(Kpts, allowed):
    sys.path.insert(0, ROOT + "/code")
    from pg import PG
    P = PG(4, 9)
    S = set(Kpts)
    bad = sum(1 for L in P.lines if sum(1 for p in L if p in S) not in allowed)
    return bad == 0, P


def monitor():
    while not done_evt.is_set():
        total, biggest = 0, None
        with lock:
            for jid, r in running.items():
                sz = 0
                lsz = 0
                for f in (r["cnf"],):
                    try:
                        sz += os.path.getsize(f)
                    except OSError:
                        pass
                if r["lrat"]:
                    try:
                        lsz = os.path.getsize(r["lrat"])
                    except OSError:
                        lsz = 0
                total += sz + lsz
                if r["phase"] == "solve" and r["proc"] is not None and (biggest is None or lsz > biggest[1]):
                    biggest = (jid, lsz)
            peak[0] = max(peak[0], total)
            if total > CAP and biggest and biggest[1] > 0:
                r = running[biggest[0]]
                r["killed"] = True
                try:
                    r["proc"].kill()
                except Exception:
                    pass
                log("SCRATCH %.1f GB > cap: killing %s (LRAT %.1f GB)" % (total / 2 ** 30, biggest[0], biggest[1] / 2 ** 30))
        time.sleep(2)


def run_job(job, solo):
    jid = job["id"]
    kind = job["kind"]
    cnf = "%s/%s.cnf" % (SCR, jid)
    lrat = "%s/%s.lrat" % (SCR, jid)
    out = "%s/%s.out" % (SCR, jid)
    t0 = time.time()
    tlim = int(job.get("local_tlim", TLIM))
    base_rec = dict(tag=jid, kind=kind, wall_start=time.strftime("%Y-%m-%d %H:%M:%S"))
    try:
        if kind == "control":
            with open(ROOT + "/cnf/" + job["cnf"], "rb") as fi, open(cnf, "wb") as fo:
                fo.write(fi.read())
            cmd = [CADICAL, "-t", "900", cnf]
            use_lrat = False
        else:
            write_cube_cnf("%s/%s" % (SCR, job["base"]), job["units"], cnf)
            cmd = [CADICAL, "-t", str(tlim), "--lrat=true", "--binary=true", "--congruence=false", cnf, lrat]      # binary LRAT: 2.1x smaller, cake_lpr 1.7x faster (speed check)
            use_lrat = True
        cnf_sha = sha256_file(cnf)
        exp = job.get("expect") or {}
        if exp.get("cnf_sha256") and exp["cnf_sha256"] != cnf_sha:
            record(dict(base_rec, status="DISAGREE_CNF", cnf_sha256=cnf_sha, expected_cnf_sha256=exp["cnf_sha256"]))
            stop_reason.append("DISAGREEMENT (CNF hash) on " + jid)
            stop_flag.set()
            return
        with open(out, "w") as fo:
            proc = subprocess.Popen(cmd, stdout=fo, stderr=subprocess.STDOUT)
        with lock:
            running[jid] = dict(proc=proc, cnf=cnf, lrat=lrat if use_lrat else None, phase="solve", killed=False)
        rc = proc.wait()
        wall = time.time() - t0
        info, model = parse_out(out)
        with lock:
            killed = running[jid]["killed"]
            running[jid]["phase"] = "cake" if rc == 20 else "done"
        rec = dict(base_rec, cnf_sha256=cnf_sha, solver_seconds=round(wall, 1), **info)
        if kind == "cube":
            rec.update(cls="star_%s%d" % (job["case"], job["cls"]), k=job["k"], case=job["case"], cube_class=job["cls"], cube=job["cube"],
                       orbit_size=job["orbit_size"], plane=job["plane"])
        if killed:
            rec["status"] = "SCRATCH_KILLED"
            if solo:
                rec["status"] = "SKIPPED_SCRATCH"
                rec["flag"] = "LRAT would exceed the %.0f GB scratch cap even when run alone; skipped" % (CAP / 2 ** 30)
                record(rec)
            else:
                deferred.append(job)
                log("DEFERRED", jid)
        elif rc == 20 and kind == "cube":
            with cake_lock:
                for _ in range(720):          # patient: up to 3 h (the solvers of the other workers end within 2 h each)
                    if mem_available_mb() > CAKE_HEAP + 300:
                        break
                    time.sleep(15)
                tc = time.time()
                r = subprocess.run([CAKE, "--CML_HEAP_SIZE=%d" % CAKE_HEAP, "--CML_STACK_SIZE=1000", cnf, lrat], capture_output=True, text=True)
                lines = [l for l in (r.stdout.replace(chr(0), "") + r.stderr.replace(chr(0), "")).splitlines() if l.strip()]
                verdict = lines[-1] if lines else "(no output)"
                if verdict != "s VERIFIED UNSAT":              # most likely heap exhaustion: one retry with a bigger heap, when memory allows
                    rec["cake_first_try"] = verdict[:200]
                    for _ in range(240):
                        if mem_available_mb() > CAKE_HEAP_RETRY + 300:
                            break
                        time.sleep(15)
                    r = subprocess.run([CAKE, "--CML_HEAP_SIZE=%d" % CAKE_HEAP_RETRY, "--CML_STACK_SIZE=1000", cnf, lrat], capture_output=True, text=True)
                    lines = [l for l in (r.stdout.replace(chr(0), "") + r.stderr.replace(chr(0), "")).splitlines() if l.strip()]
                    verdict = lines[-1] if lines else "(no output)"
                if verdict != "s VERIFIED UNSAT" and NW == 1 and "heap" in verdict.lower():        # alone on the machine: a last try with a 5.0 GB heap
                    rec["cake_second_try"] = verdict[:200]
                    for _ in range(240):
                        if mem_available_mb() > CAKE_HEAP_RETRY2 + 300:
                            break
                        time.sleep(15)
                    r = subprocess.run([CAKE, "--CML_HEAP_SIZE=%d" % CAKE_HEAP_RETRY2, "--CML_STACK_SIZE=1000", cnf, lrat], capture_output=True, text=True)
                    lines = [l for l in (r.stdout.replace(chr(0), "") + r.stderr.replace(chr(0), "")).splitlines() if l.strip()]
                    verdict = lines[-1] if lines else "(no output)"
                cake_s = time.time() - tc
            rec.update(status="UNSAT_VERIFIED" if verdict == "s VERIFIED UNSAT" else "UNSAT_CAKE_FAILED", cake_lpr=verdict,
                       cake_seconds=round(cake_s, 1), lrat_bytes=os.path.getsize(lrat), lrat_sha256=sha256_file(lrat))
            rec["total_cpu_hours"] = round((info.get("solver_cpu_seconds", wall) + cake_s) / 3600.0, 4)
            if rec["status"] != "UNSAT_VERIFIED" and rec["lrat_bytes"] < 10 * 2 ** 30:     # keep a rejected proof (at most 10 GB) for diagnosis
                qd = os.path.expanduser("~/pg39work/quarantine")
                os.makedirs(qd, exist_ok=True)
                os.replace(lrat, "%s/%s.lrat" % (qd, jid))
                rec["quarantine"] = "%s/%s.lrat" % (qd, jid)
            bad = []
            if exp and rec["status"] == "UNSAT_VERIFIED":
                bad = ["%s: local %s / modal %s" % (k, rec.get(k), exp[k]) for k in ("conflicts", "lrat_bytes", "lrat_sha256") if k in exp and rec.get(k) != exp[k]]
            if bad or rec["status"] != "UNSAT_VERIFIED":
                rec["flag"] = "DISAGREEMENT: " + "; ".join(bad) if bad else "cake_lpr did not verify"
                stop_reason.append("%s on %s" % (rec["flag"], jid))
                stop_flag.set()
                if "quarantine" not in rec:                      # keep the files of the item (CNF, LRAT, solver output) for diagnosis
                    qd = os.path.expanduser("~/pg39work/quarantine")
                    os.makedirs(qd, exist_ok=True)
                    for f in (cnf, lrat, out):
                        if os.path.exists(f):
                            os.replace(f, "%s/%s%s" % (qd, jid, os.path.splitext(f)[1]))
                    rec["quarantine"] = qd
            record(rec, manifest=True)
        elif rc == 10:
            pos = set(m for m in model if m > 0)
            K = sorted(p for p in range(820) if p + 1 in pos) if kind == "cube" else None
            if kind == "control":
                sys.path.insert(0, ROOT + "/code")
                from controls import hermitian_points
                from pg import PG
                P = PG(4, 9)
                S = set(p for p in range(820) if (p + 1) in pos)
                H = set(hermitian_points(P, 3))
                lines_ok = all(sum(1 for p in L if p in S) in (1, 4, 10) for L in P.lines)       # exact check of all 7462 lines
                ok = lines_ok and len(S) == 280                                                  # the model need not be the Hermitian surface itself
                rec.update(status="CONTROL_SAT_OK" if ok else "CONTROL_FAILED", size=len(S), lines_ok=lines_ok, equals_hermitian=(S == H), overlap_with_H=len(S & H))
                record(rec)
            else:
                ok, P = exact_check(K, {2, 5})
                rec.update(status="SAT_EXACT_PASS" if ok else "SAT_EXACT_FAIL", size=len(K))
                json.dump(K, open("%s/results/MODEL_%s.json" % (ROOT, jid), "w"))
                record(rec, manifest=True)
                log("!!!!!!!! SAT CUBE", jid, "exact check", "PASS" if ok else "FAIL", "- stopping the batch")
                stop_flag.set()
        elif rc == 0 and wall >= tlim - 30 and kind == "cube":
            rec["status"] = "TIMEOUT_LOCAL" if "local_tlim" in job else "TIMEOUT_2h"
            rec["flag"] = "not closed after %d s: recorded as a timeout, not retried" % tlim
            record(rec, manifest=True)
        else:
            rec["status"] = "ERROR"
            rec["rc"] = rc
            record(rec)
    except Exception as e:  # noqa
        record(dict(base_rec, status="ERROR", error=repr(e)))
    finally:
        with lock:
            running.pop(jid, None)
        for f in (cnf, lrat, out):
            try:
                os.remove(f)
            except OSError:
                pass


def worker(k, todo, solo):
    while not stop_flag.is_set():
        why = None
        if PREEMPT and os.path.exists(PREEMPT):
            why = "broker preempt flag"
        elif STOP_AFTER is not None and time.time() - T_START > STOP_AFTER:
            why = "slice length reached"
        elif MINFREE:
            st = os.statvfs("/mnt/c")
            if st.f_bavail * st.f_frsize < MINFREE:
                why = "free space on C: below %.0f GB" % (MINFREE / 2 ** 30)
        if why:
            if why not in stop_reason:
                stop_reason.append(why)
                log("STOP launching new jobs:", why)
            return
        if k >= (1 if os.path.exists(GAMING) else NW):
            time.sleep(30)
            if not todo:
                return
            continue
        with lock:
            job = todo.popleft() if todo else None
        if job is None:
            return
        run_job(job, solo)


def main():
    os.makedirs(SCR, exist_ok=True)
    jobs = json.load(open(JOBS))
    terminal = {"TIMEOUT_LOCAL", "DISAGREE_CNF", "UNSAT_VERIFIED", "UNSAT_CAKE_FAILED", "TIMEOUT_2h", "SKIPPED_SCRATCH", "SAT_EXACT_PASS", "SAT_EXACT_FAIL", "CONTROL_SAT_OK", "CONTROL_FAILED", "ERROR"}
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            try:
                r = json.loads(l)
                if r.get("status") in terminal:
                    done.add(r["tag"])
            except Exception:
                pass
    todo = deque(j for j in jobs if j["id"] not in done and (ONLY is None or j["id"] in ONLY))
    # base CNFs to WSL-local disk
    for b in sorted(set(j["base"] for j in todo if j["kind"] == "cube")):
        dst = "%s/%s" % (SCR, b)
        if not os.path.exists(dst):
            with open(ROOT + "/cnf/" + b, "rb") as fi, open(dst, "wb") as fo:
                while True:
                    ch = fi.read(1 << 24)
                    if not ch:
                        break
                    fo.write(ch)
    log("start: %d jobs to do (%d already done), workers %d, scratch cap %.1f GB, wall limit %d s per cube" % (len(todo), len(done), NW, CAP / 2 ** 30, TLIM))
    threading.Thread(target=monitor, daemon=True).start()
    ths = [threading.Thread(target=worker, args=(k, todo, False)) for k in range(NW)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    if deferred and not stop_flag.is_set():
        log("second pass: %d deferred cubes, one at a time" % len(deferred))
        todo2 = deque(deferred)
        deferred.clear()
        t = threading.Thread(target=worker, args=(0, todo2, True))
        t.start()
        t.join()
    done_evt.set()
    summary = dict(tag="a_slice_summary", status="BATCH_DONE", stop_reason=stop_reason, peak_scratch_gb=round(peak[0] / 2 ** 30, 2), stopped_on_sat=stop_flag.is_set(),
                   time=time.strftime("%Y-%m-%d %H:%M:%S"))
    record(summary)
    for b in os.listdir(SCR):
        try:
            os.remove(os.path.join(SCR, b))
        except OSError:
            pass


if __name__ == "__main__":
    try:
        main()
    finally:
        with lock:
            for r in running.values():
                try:
                    r["proc"].kill()
                except Exception:
                    pass

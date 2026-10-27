"""Driver: for each job "m:index" generate the eager CNF, solve with CaDiCaL (text LRAT),
check the LRAT with cake_lpr (CakeML-verified, in WSL) and gate_lrat_check (ours), record
hashes and verdicts in results/m{m}.jsonl, and delete the LRAT afterwards (KEEP_LRAT=1 keeps).

usage: python run_certs.py m:idx[@spec][,idx...] [m:idx ...]   spec = full | keep=p1.p2 (lks.variant_core)       e.g.  python run_certs.py 9:2,4,7,10,12
       python run_certs.py m:E      (all exceptional trees of classification_m{m}.json)
"""
import sys, os, json, subprocess, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CADICAL = os.environ.get("CADICAL", "cadical")
GATE = os.environ.get("GATE_LRAT_CHECK",
                      "gate_lrat_check")
CAKE = os.environ.get("CAKE_LPR", "~/lrchk/bin/cake_lpr")
WORK = os.environ.get("WORK", os.path.join(HERE, "work", "certs"))
HEAP = os.environ.get("CAKE_HEAP_MB", "2500")
TLIM = float(os.environ["CADICAL_TIMEOUT"]) if os.environ.get("CADICAL_TIMEOUT") else None


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def job(m, i, log, spec=""):
    os.makedirs(WORK, exist_ok=True)
    stem = os.path.join(WORK, f"m{m}_t{i:03d}" + (("_" + spec.replace("=", "").replace(".", "_")) if spec else ""))
    cnf, lrat = stem + ".cnf", stem + ".lrat"
    rec = dict(m=m, tree=i, mode=spec or "rooted-core")
    t = time.time()
    r = subprocess.run([sys.executable, os.path.join(HERE, "gen_eager.py"), str(m), str(i), cnf] + ([spec] if spec else []),
                       capture_output=True, text=True)
    rec["gen"] = r.stdout.strip(); rec["gen_s"] = round(time.time() - t, 1)
    if r.returncode:
        rec["error"] = r.stderr[-500:]; return rec
    rec["cnf_sha256"] = sha(cnf)
    t = time.time()
    r = subprocess.run([sys.executable, os.path.join(HERE, "audit_cnf.py"), str(m), str(i), cnf] + ([spec] if spec else []),
                       capture_output=True, text=True)
    rec["audit"] = (r.stdout + r.stderr).strip()[-200:]; rec["audit_s"] = round(time.time() - t, 1)
    if r.returncode:
        rec["error"] = "audit failed"; return rec
    t = time.time()
    try:
        r = subprocess.run([CADICAL, "-q", "--lrat=true", "--binary=false", cnf, lrat],
                           capture_output=True, text=True, timeout=TLIM)
        rec["cadical_rc"] = r.returncode
    except subprocess.TimeoutExpired:
        rec["cadical_rc"] = "timeout"
        if os.path.exists(lrat):
            os.remove(lrat)
    rec["cadical_s"] = round(time.time() - t, 1)
    if rec["cadical_rc"] != 20:
        rec["error"] = "not UNSAT (%s)" % rec["cadical_rc"]; return rec
    rec["lrat_bytes"] = os.path.getsize(lrat); rec["lrat_sha256"] = sha(lrat)
    t = time.time()
    r = subprocess.run([GATE, cnf, lrat], capture_output=True, text=True)
    rec["gate_rc"] = r.returncode; rec["gate_out"] = (r.stdout + r.stderr).strip()[-200:]
    rec["gate_s"] = round(time.time() - t, 1)
    t = time.time()
    r = subprocess.run(["wsl", "-d", "Ubuntu-24.04", "--", "bash", wsl(os.path.join(HERE, "cake_check.sh")),
                        wsl(cnf), wsl(lrat), HEAP], capture_output=True)
    out = (r.stdout + r.stderr).decode(errors="replace").replace("\0", "").strip()
    rec["cake_rc"] = r.returncode; rec["cake_out"] = out[-200:]
    rec["cake_s"] = round(time.time() - t, 1)
    rec["verified"] = ("s VERIFIED UNSAT" in out and r.returncode == 0 and rec["gate_rc"] == 0)
    if not os.environ.get("KEEP_LRAT"):
        os.remove(lrat)
    return rec


def main():
    jobs = []
    for a in sys.argv[1:]:
        m, idx = a.split(":")
        m = int(m)
        if idx == "E":
            rows = json.load(open(os.path.join(HERE, f"classification_m{m}.json")))
            jobs += [(m, r["index"], "") for r in rows if r["route"] == "E"]
        else:
            jobs += [(m, int(x.split("@")[0]), x.split("@")[1] if "@" in x else "") for x in idx.split(",")]
    RES = os.environ.get("RESULTS", os.path.join(HERE, "results"))
    os.makedirs(RES, exist_ok=True)
    for m, i, spec in jobs:
        out = os.path.join(RES, f"m{m}.jsonl")
        done = set()
        if os.path.exists(out):
            done = {(d["m"], d["tree"], "" if d.get("mode", "rooted-core") == "rooted-core" else d["mode"]) for d in map(json.loads, open(out)) if d.get("verified")}
        alts = spec.split("|")      # alternatives: try in turn until one verifies
        if any((m, i, sp) in done for sp in alts):
            continue
        for sp in alts:
            rec = job(m, i, None, sp)
            with open(out, "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(json.dumps({k: rec.get(k) for k in ("m", "tree", "mode", "audit_s", "cadical_s",
                  "cake_s", "lrat_bytes", "cake_out", "gate_rc", "verified", "error")}), flush=True)
            if rec.get("verified"):
                break


if __name__ == "__main__":
    main()

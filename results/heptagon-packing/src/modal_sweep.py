"""Modal driver for BOTH certificate sweeps (audit fixes 2026-11-04).

Jobs: (mode, slice tag, theta0 range). mode 'main' = nonobtuse sweep (no table, charts at the DL),
mode 'obt' = obtuse-mode sweep (table D5). Regions are widened outward (audit FLAW-2):
  t0: slice k covers [0.02k, 0.02(k+1)]; slice 0 starts at -1e-12, the last ends at 2pi/7 + 1e-12;
  t1, t2 in [-1e-12, 2pi/7 + 1e-12]; main-mode d starts at 1.8019377358048 (< 2u).
Every slice log starts with bnb2's provenance block (python, numpy, sha256 of every input file,
HEPT_* env). Checkpoints (every 10 min) live on the Volume heptagon-ckpt under /ckpt/<run>/ and are
KEPT after completion (renamed *.final.npy) until certification (audit UNCLEAR-2).
Status: DONE | TIMEOUT (dump returned) | STUCK (driver requests graceful stop, keeps collecting) | FAIL.
Usage (from heptagon-packing/):
    python -m modal run src/modal_sweep.py::main --mode main --run m1 [--only 39a,39b] [--resume 1]
    python -m modal run src/modal_sweep.py::request_stop
"""
import modal, os, json

P7 = 0.8975979010256552        # fl(2 pi / 7) (below 2pi/7 by 3.5e-17)
EPS = 1e-12
FILES = ["src/rig.py", "src/bnb2.py", "src/chart_excl.py", "runs/f_w030.npy", "runs/dl_ref.json",
         "runs/charts_w030c.json", "runs/D5_nodes.npy", "runs/D5_bins.json"]
image = modal.Image.debian_slim(python_version="3.13").pip_install("numpy==2.5.4")
for f in FILES:
    image = image.add_local_file(f, f"/root/hp/{f}")
app = modal.App("heptagon-sweeps")
ctl = modal.Dict.from_name("heptagon-ctl", create_if_missing=True)
vol = modal.Volume.from_name("heptagon-ckpt", create_if_missing=True)

def jobs(mode):
    out = []
    for k in range(45):
        lo = -EPS if k == 0 else round(0.02*k, 10)
        hi = P7 + EPS if k == 44 else round(0.02*(k+1), 10)
        if mode == "main" and k == 39:        # the DL slice: split 10 ways for parallelism
            for q in range(10):
                a = round(0.78 + 0.002*q, 10); b = round(0.78 + 0.002*(q+1), 10)
                out.append((f"39{chr(97+q)}", a, b))
        else:
            out.append((str(k), lo, hi))
    return out

def region(mode, lo, hi):
    if mode == "main":
        return {"lo": [1.8019377358048, 0, 1.2, lo, -EPS, -EPS], "hi": [4.0, 2.0, 3.5, hi, P7 + EPS, P7 + EPS]}
    return {"lo": [2.548, 0, 0.8, lo, -EPS, -EPS], "hi": [4.0, 2.0, 3.5, hi, P7 + EPS, P7 + EPS]}

@app.function(image=image, cpu=1.0, memory=1024, timeout=8*3600, volumes={"/ckpt": vol})
def run_slice(mode: str, tag: str, lo: float, hi: float, tl: float, run: str, resume: bytes = b"", nockpt: int = 0):
    import subprocess, sys, hashlib, time, signal
    os.chdir("/root/hp")
    hdr = [f"job mode={mode} tag={tag} t0=[{lo!r}, {hi!r}] run={run}"]
    reg = json.dumps(region(mode, lo, hi))
    os.makedirs(f"/ckpt/{run}", exist_ok=True)
    ck = f"/ckpt/{run}/open_{mode}_{tag}.npy"
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", HEPT_MU="[-1.1684855e-12, -1.1684853e-12]",
               HEPT_DUAL2="1", HEPT_CKPT=ck, HEPT_CKPT_EVERY="600")
    if nockpt:                     # from-scratch replay: no checkpoint; a preempted attempt restarts from scratch
        env.pop("HEPT_CKPT"); env.pop("HEPT_CKPT_EVERY")
        hdr.append("checkpointing disabled (from-scratch replay)")
    attempt = f"{int(time.time())}"
    vlog = f"/ckpt/{run}/log_{mode}_{tag}_{attempt}.txt"      # every attempt's log is kept on the Volume
    if mode == "main":
        env["HEPT_CHARTS"] = "runs/charts_w030c.json"
    else:
        env["HEPT_OBT4"] = "runs/D5_nodes.npy,runs/D5_bins.json"
    vol.reload()
    if os.path.exists(ck) and not nockpt:
        env["HEPT_RESUME"] = ck
        hdr.append(f"resumed from checkpoint {ck} sha256 {hashlib.sha256(open(ck, 'rb').read()).hexdigest()} (container restart)")
        keep = f"/ckpt/{run}/used_{mode}_{tag}_{int(time.time())}.npy"
        open(keep, "wb").write(open(ck, "rb").read())        # keep the checkpoint that was resumed from
    elif resume:
        open("/tmp/resume.npy", "wb").write(resume); env["HEPT_RESUME"] = "/tmp/resume.npy"
        hdr.append(f"resumed from dump sha256 {hashlib.sha256(resume).hexdigest()}")
    p = subprocess.Popen([sys.executable, "src/bnb2.py", "runs/f_w030.npy", "3", "2", "1.5326754446782214", str(tl), reg, "/tmp/open.npy"],
                         stdout=open("/tmp/out.txt", "w"), stderr=open("/tmp/err.txt", "w"), env=env)
    last_commit = time.time()
    while p.poll() is None:
        time.sleep(30)
        try:
            if time.time() - last_commit > 300:
                open(vlog, "w").write(chr(10).join(hdr) + chr(10) + open("/tmp/out.txt").read())
                vol.commit(); last_commit = time.time()
            if ctl.get("stop", False):
                p.send_signal(signal.SIGTERM); hdr.append("graceful stop requested"); p.wait(); break
        except Exception:
            pass
    out = open("/tmp/out.txt").read(); err = open("/tmp/err.txt").read()
    dump = open("/tmp/open.npy", "rb").read() if os.path.exists("/tmp/open.npy") else b""
    try:
        open(vlog, "w").write(chr(10).join(hdr) + chr(10) + out + err[-3000:])
        if os.path.exists(ck): os.replace(ck, ck.replace(".npy", ".final.npy"))
        vol.commit()
    except Exception:
        pass
    return mode, tag, "\n".join(hdr) + "\n" + out + err[-3000:], dump

@app.local_entrypoint()
def request_stop(flag: int = 1):
    ctl["stop"] = bool(flag); print("heptagon-ctl stop =", bool(flag))

def classify(log):
    last = [l for l in log.splitlines() if l.startswith("{")]
    if not last: return "FAIL", None
    r = json.loads(last[-1])
    if r.get("stuck", 0) > 0: return "STUCK", r
    return ("DONE" if r["done"] else "TIMEOUT"), r

@app.local_entrypoint()
def main(mode: str = "obt", run: str = "a1", only: str = "", tl: float = 13500.0, resume: int = 0, nockpt: int = 0):
    J = jobs(mode)
    base = {t: (a, b) for t, a, b in J}
    if only:
        sel = []
        for t in only.split(","):
            if t in base: sel.append((t, *base[t]))
            elif "s" in t and t.split("s")[0] in base:      # shard tag "<slice>s<i>": same region as <slice>
                sel.append((t, *base[t.split("s")[0]]))
        J = sel
    ctl["stop"] = False
    out_dir = f"runs/sweep_{mode}_{run}"; os.makedirs(out_dir, exist_ok=True)
    dumps = [open(f"{out_dir}/open_{t}.npy", "rb").read() if resume else b"" for t, _, _ in J]
    stopped = False
    args = ([mode]*len(J), [t for t, _, _ in J], [a for _, a, _ in J], [b for _, _, b in J], [tl]*len(J), [run]*len(J), dumps, [nockpt]*len(J))
    for md, tag, log, dump in run_slice.map(*args, order_outputs=False):
        sfx = "_resume" if resume else ""
        open(f"{out_dir}/k{tag}{sfx}.log", "w").write(log)
        if dump: open(f"{out_dir}/open_{tag}.npy", "wb").write(dump)
        st, r = classify(log)
        print(tag, st, (json.dumps({q: r[q] for q in ("processed", "stuck", "open", "secs")}) if r else ""), flush=True)
        if st == "STUCK" and not stopped:
            print("STOP RULE: slice", tag, "stuck; requesting graceful stop, collecting the rest", flush=True)
            ctl["stop"] = True; stopped = True

"""Run obtuse-mode B&B slices on Modal CPU containers (one slice per container).

Pins: numpy==2.5.4 in the image; the exact source files and the table are shipped and their
sha256 are printed into every slice log, together with numpy.__version__ / __file__.
- Checkpoints: bnb2 writes its open stack every 10 min to a Modal Volume (heptagon-ckpt,
  /ckpt/<run>/open_k{k}.npy). A preempted container restarts with the same input and resumes
  from that checkpoint (logged). Sound: a checkpoint is exactly the open set at that moment;
  every other box was already decided.
- Time cap: a slice reaching --tl dumps its open stack and is classified TIMEOUT; resume with
  --resume 1 (uses runs/sweepM/open_k{k}.npy).
- STUCK (stuck > 0): the driver sets the graceful-stop flag and KEEPS collecting, so every
  other slice returns its dump.
Status per slice: DONE (finished, stuck=0) | TIMEOUT | STUCK | FAIL.
Usage (from heptagon-packing/):
    python -m modal run src/modal_obt.py::main --ks 0,1,2 --tl 13500 --table D5 --run r1 [--resume 1]
    python -m modal run src/modal_obt.py::request_stop
"""
import modal, os, json

TABLE = os.environ.get("HEPT_TABLE", "D5")
FILES = ["src/rig.py", "src/bnb2.py", "runs/f_w030.npy", f"runs/{TABLE}_nodes.npy", f"runs/{TABLE}_bins.json"]
image = modal.Image.debian_slim(python_version="3.13").pip_install("numpy==2.5.4")
for f in FILES:
    image = image.add_local_file(f, f"/root/hp/{f}")
app = modal.App("heptagon-obtuse-sweep")
ctl = modal.Dict.from_name("heptagon-ctl", create_if_missing=True)   # ctl["stop"] = True -> graceful stop (dump)
vol = modal.Volume.from_name("heptagon-ckpt", create_if_missing=True)

@app.function(image=image, cpu=1.0, memory=1024, timeout=4*3600, volumes={"/ckpt": vol})
def run_slice(k: int, tl: float, resume: bytes = b"", run: str = "r", table: str = "D5"):
    import subprocess, sys, hashlib, numpy, time, signal
    os.chdir("/root/hp")
    hdr = [f"numpy {numpy.__version__} {numpy.__file__}", f"python {sys.version}", f"table {table}"]
    for f in FILES:
        hdr.append(f"sha256 {f} {hashlib.sha256(open(f, 'rb').read()).hexdigest()}")
    lo = round(0.02*k, 10); hi = min(round(0.02*(k+1), 10), 0.8975979010256552)
    reg = json.dumps({"lo": [2.548, 0, 0.8, lo, 0, 0], "hi": [4.0, 2.0, 3.5, hi, 0.8975979010256552, 0.8975979010256552]})
    os.makedirs(f"/ckpt/{run}", exist_ok=True)
    ck = f"/ckpt/{run}/open_k{k}.npy"
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", HEPT_MU="[-1.1684855e-12, -1.1684853e-12]",
               HEPT_DUAL2="1", HEPT_OBT4=f"runs/{table}_nodes.npy,runs/{table}_bins.json",
               HEPT_CKPT=ck, HEPT_CKPT_EVERY="600")
    vol.reload()
    if os.path.exists(ck):
        env["HEPT_RESUME"] = ck
        hdr.append(f"resumed from checkpoint {ck} sha256 {hashlib.sha256(open(ck, 'rb').read()).hexdigest()} (container restart)")
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
                vol.commit(); last_commit = time.time()
            if ctl.get("stop", False):
                p.send_signal(signal.SIGTERM); hdr.append("graceful stop requested"); p.wait(); break
        except Exception:
            pass
    out = open("/tmp/out.txt").read(); err = open("/tmp/err.txt").read()
    dump = open("/tmp/open.npy", "rb").read() if os.path.exists("/tmp/open.npy") else b""
    try:
        if os.path.exists(ck): os.remove(ck)        # finished attempt: clear checkpoint
        vol.commit()
    except Exception:
        pass
    return k, "\n".join(hdr) + "\n" + out + err[-3000:], dump

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
def main(ks: str = "0", tl: float = 13500.0, resume: int = 0, run: str = "r1", table: str = "D5"):
    assert table == TABLE, "set HEPT_TABLE env to the same table as --table"
    klist = [int(x) for x in ks.split(",")]
    ctl["stop"] = False
    out_dir = f"runs/sweepM_{run}"; os.makedirs(out_dir, exist_ok=True)
    dumps = [open(f"{out_dir}/open_k{k}.npy", "rb").read() if resume else b"" for k in klist]
    stopped = False
    for k, log, dump in run_slice.map(klist, [tl]*len(klist), dumps, [run]*len(klist), [table]*len(klist), order_outputs=False):
        tag = "_resume" if resume else ""
        open(f"{out_dir}/k{k}{tag}.log", "w").write(log)
        if dump: open(f"{out_dir}/open_k{k}.npy", "wb").write(dump)
        st, r = classify(log)
        print(k, st, (json.dumps({q: r[q] for q in ("processed", "stuck", "open", "secs")}) if r else ""), flush=True)
        if st == "STUCK" and not stopped:
            print("STOP RULE: slice", k, "has stuck boxes; requesting graceful stop, collecting remaining dumps", flush=True)
            ctl["stop"] = True; stopped = True

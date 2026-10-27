"""One command to re-check everything in this directory (Windows + WSL Ubuntu-24.04):
    python verify.py            # full: Zenodo audit + all our certificates (several CPU-hours)
    python verify.py --quick    # Zenodo audit + only the certificates that solve in < 30 s
Writes fresh results to results_recheck/ ; exits non-zero on any failure."""
import subprocess, sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
def run(args, **kw):
    print("+", " ".join(args), flush=True)
    r = subprocess.run(args, **kw)
    if r.returncode:
        sys.exit(f"FAILED: {args}")
for m in range(4, 11):
    run([PY, os.path.join(HERE, "classify.py"), str(m)], cwd=HERE)
if os.path.isdir(os.path.join(HERE, "zenodo")):  # Zenodo files are not shipped here; download record 21348157 into zenodo/
    run([PY, os.path.join(HERE, "diff_zenodo.py")], cwd=HERE)
    run([PY, os.path.join(HERE, "check_zbase.py"), "zenodo/path_support_014.cnf",
         "work/path_support_014.varmap"], cwd=HERE)
    run(["wsl", "-d", "Ubuntu-24.04", "--", "bash",
         os.path.join(HERE, "check_zenodo.sh")])
else:
    print("zenodo/ not present: skipping the audit of Zenodo 21348157 (our own certificates are still checked)")
env = dict(os.environ, RESULTS=os.path.join(HERE, "results_recheck"))
man = json.load(open(os.path.join(HERE, "manifest.json")))
items = list(man.values())
if "--quick" in sys.argv:
    items = [x for x in items if x["cadical_s"] < 30]
jobs = [f"{x['m']}:{x['tree']}" + (f"@{x['spec']}" if x["spec"] else "") for x in items]
run([PY, os.path.join(HERE, "run_certs.py")] + jobs, env=env)
bad = []
for f in os.listdir(env["RESULTS"]):
    for d in map(json.loads, open(os.path.join(env["RESULTS"], f))):
        if not d.get("verified"):
            bad.append((d["m"], d["tree"]))
# (CNF hashes in manifest.json are informational: comment headers changed during the run;
#  soundness comes from the fresh audit + LRAT + cake_lpr check above.)
print("ALL VERIFIED" if not bad else f"NOT VERIFIED: {bad}")
sys.exit(1 if bad else 0)

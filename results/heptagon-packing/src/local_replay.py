"""Local spot replay of sweep slices with exactly the Modal regions/env (modal_sweep.jobs/region)."""
import sys, os, json, subprocess
sys.path.insert(0, 'src')
P7 = 0.8975979010256552; EPS = 1e-12
def jobs(mode):
    out = []
    for k in range(45):
        lo = -EPS if k == 0 else round(0.02*k, 10)
        hi = P7 + EPS if k == 44 else round(0.02*(k+1), 10)
        if mode == "main" and k == 39:
            for q in range(10):
                out.append((f"39{chr(97+q)}", round(0.78 + 0.002*q, 10), round(0.78 + 0.002*(q+1), 10)))
        else:
            out.append((str(k), lo, hi))
    return out
def region(mode, lo, hi):
    if mode == "main":
        return {"lo": [1.8019377358048, 0, 1.2, lo, -EPS, -EPS], "hi": [4.0, 2.0, 3.5, hi, P7 + EPS, P7 + EPS]}
    return {"lo": [2.548, 0, 0.8, lo, -EPS, -EPS], "hi": [4.0, 2.0, 3.5, hi, P7 + EPS, P7 + EPS]}
mode, tag = sys.argv[1], sys.argv[2]
lo, hi = [(a, b) for t, a, b in jobs(mode) if t == tag][0]
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", HEPT_MU="[-1.1684855e-12, -1.1684853e-12]", HEPT_DUAL2="1")
if mode == "main": env["HEPT_CHARTS"] = "runs/charts_w030c.json"
else: env["HEPT_OBT4"] = "runs/D5_nodes.npy,runs/D5_bins.json"
os.makedirs("runs/replay_a1", exist_ok=True)
with open(f"runs/replay_a1/{mode}_k{tag}.log", "w") as fo:
    fo.write(f"job mode={mode} tag={tag} t0=[{lo!r}, {hi!r}] (local replay)\n"); fo.flush()
    subprocess.run([sys.executable, "src/bnb2.py", "runs/f_w030.npy", "3", "2", "1.5326754446782214", "40000", json.dumps(region(mode, lo, hi))],
                   stdout=fo, stderr=subprocess.STDOUT, env=env)

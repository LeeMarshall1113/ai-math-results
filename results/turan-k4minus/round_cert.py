"""Round a float SDP solution (sdp_ipm.py --out) to an exact integer certificate (JSON.gz).
usage: python round_cert.py data7.npz sol.npz cert.json.gz [--scale M] [--eig-tol r]"""
import sys, os, json, gzip, argparse
sys.path.append(r"C:/Users/hackathon/AppData/Roaming/Python/Python314/site-packages")
from fractions import Fraction
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("data"); ap.add_argument("sol"); ap.add_argument("out")
ap.add_argument("--scale", type=int, default=10**7)
ap.add_argument("--eig-tol", type=float, default=1e-11)
ap.add_argument("--muden", type=int, default=10**12)
ap.add_argument("--bound", default=None, help="claimed bound (default: exact max computed by numpy-int check)")
a = ap.parse_args()
D = np.load(a.data); S = np.load(a.sol)
M = a.scale
cert = {"N": 7, "graphs": [int(x) for x in np.load(a.data.replace("data", "graphs")) ]} if False else {"N": 7}
import pickle
cert["graphs"] = [int(x) for x in pickle.load(open("classes.pkl", "rb"))[7]]
cert["scale"] = M
blocks = []
for b in range(int(D["nblocks"])):
    Q = S[f"Q{b}"]; Q = 0.5 * (Q + Q.T)
    lam, U = np.linalg.eigh(Q)
    keep = lam > a.eig_tol * max(lam.max(), 1e-30)
    V = np.rint(M * (U[:, keep] * np.sqrt(lam[keep])).T).astype(np.int64)
    V = V[np.any(V != 0, axis=1)]
    blocks.append({"s": int(D[f"b{b}_s"]), "type": int(D[f"b{b}_type"]),
                   "flags": [int(x) for x in D[f"b{b}_flags"]], "factors": V.tolist()})
cert["blocks"] = blocks
w = S["w"]; pos = 0
assert not int(S["stat"]), "linear stat rows not supported in certificates"
if len(S["interval"]):
    t1, t2 = S["interval"]
    St = np.load(a.data.replace("data", "stat"))
    def rat(x):
        return str(Fraction(max(0, round(x * a.muden)), a.muden))
    itv = {"t1": str(Fraction(t1).limit_denominator(10**9)), "t2": str(Fraction(t2).limit_denominator(10**9))}
    if not int(S["nodeg"]):
        p, q = [], []
        for k in range(len(St["g"])):
            p.append(rat(w[pos])); q.append(rat(w[pos + 1])); pos += 2
        itv.update({"flags": [int(g) for g in St["g"]], "p": p, "q": q})
    if int(S["product"]):
        p, q = [], []
        for k in range(len(St["p0"])):
            p.append(rat(w[pos])); q.append(rat(w[pos + 1])); pos += 2
        itv.update({"graphs4": [int(g) for g in St["p0"]], "pp": p, "pq": q})
    cert["interval_stat"] = itv
assert pos == len(w), (pos, len(w))

# fast exact evaluation with the prep counts (independent re-check: check_cert.py)
m = len(D["edges"])
num = {}
vals = [Fraction(int(x), 35) for x in D["edges"]]
for b, blk in enumerate(blocks):
    V = [list(map(int, r)) for r in blk["factors"]]
    n = len(blk["flags"])
    G = [[0] * n for _ in range(n)]
    for v in V:
        nz = [(i, x) for i, x in enumerate(v) if x]
        for i, x in nz:
            for j, y in nz:
                G[i][j] += x * y
    T = D[f"b{b}_trip"]; den = int(D[f"b{b}_den"])
    acc = [0] * m
    for h, i, j, c in T.tolist():
        acc[h] += G[i][j] * c
    for h in range(m):
        vals[h] += Fraction(acc[h], den * M * M)
if "interval_stat" in cert:
    it = cert["interval_stat"]; t1, t2 = Fraction(it["t1"]), Fraction(it["t2"])
    for k in range(len(it.get("flags", []))):
        p, q = Fraction(it["p"][k]), Fraction(it["q"][k])
        cD, cG = St["cntD"][k].tolist(), St["cntG"][k].tolist()
        for h in range(m):
            if cG[h]:
                vals[h] += (p * (cD[h] - t1 * cG[h]) - q * (cD[h] - t2 * cG[h])) / int(St["den"])
    for k in range(len(it.get("graphs4", []))):
        p, q = Fraction(it["pp"][k]), Fraction(it["pq"][k])
        cE, cP = St["cntPE"][k].tolist(), St["cntP"][k].tolist()
        for h in range(m):
            if cP[h]:
                vals[h] += (p * (cE[h] - t1 * cP[h]) - q * (cE[h] - t2 * cP[h])) / int(St["pden"])
best = max(vals)
print("exact max_H val (prep counts):", best, float(best))
cert["bound"] = a.bound if a.bound else str(best)
with gzip.open(a.out, "wt") as fh:
    json.dump(cert, fh)
print("wrote", a.out, "factors", sum(len(b["factors"]) for b in blocks), "bound", float(Fraction(cert["bound"])))

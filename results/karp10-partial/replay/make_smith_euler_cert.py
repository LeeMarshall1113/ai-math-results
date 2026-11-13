"""make_smith_euler_cert.py -- certs/lemma_z5_euler.json: Z5 Smith patterns filtered by the Euler equations of the two groups of order 80 and 160
(classes 1014 and 1250 of gap/subgroups_s10.txt) that eliminate 8 of the 16 patterns passing the (necessary) relative-homology filter.  Claims = classes in / out of all surviving patterns
(computed in stage2/z5_z10_joint.py); replay_smith.py re-derives and checks everything."""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from k10common import edges_of  # noqa: E402

root = os.path.join(HERE, "..")
d = json.load(open(os.path.join(root, "stage2", "data2", "z5_joint.json")))
sub = {}
for line in open(os.path.join(root, "gap", "subgroups_s10.txt")):
    if line.startswith("R|"):
        _, idx, order, pp, ex, cg, gens = line.rstrip("\n").split("|")
        sub[int(idx)] = dict(order=int(order), gens=[[x - 1 for x in g] for g in eval(gens)])
wit = {}
for line in open(os.path.join(root, "gap", "witness_z5.out")):
    if line.startswith("W|") and "NOT_G1" not in line:
        parts = line.rstrip("\n").split("|")
        name, order, p, nord, ngens, g = parts[1:7]
        ng = [[int(x) - 1 for x in t.split(",")] for t in re.findall(r"\[([0-9,]+)\]", ngens)]
        wit[name] = dict(order=int(order), p=int(p), N_gens=ng, g=[int(x) - 1 for x in g.strip("[]").split(",")])
euler = []
for idx in (1014, 1250):
    w = wit[f"G{idx}"]
    euler.append(dict(kind="fix", gens=sub[idx]["gens"], order=w["order"], p=w["p"], N_gens=w["N_gens"], g=w["g"]))
cert = {
    "type": 5,
    "prime": 5,
    "claim": ("every non-evasive nontrivial monotone P of type 5 contains nine explicit graphs (C_10, the pentagon with pendant edges, Petersen, the pentagonal "
              "prism, K_5 with pendant edges, a 15-edge graph with degrees 4^5 2^5, a 20-edge 4-regular graph, (K_5 u C_5)+matching, K_5xK_2) and none of "
              "their complements (Smith theory for <(01234)(56789)> plus the Euler equations of two groups of order 80 and 160)"),
    "gens": [[1, 2, 3, 4, 0, 6, 7, 8, 9, 5]],
    "euler": euler,
    "claims_in": [[list(e) for e in edges_of(m)] for m in d["always"]],
    "claims_out": [[list(e) for e in edges_of(m)] for m in d["never"]],
    "expect": {"downsets": 2539, "relative_survivors": 16, "after_euler": [8, 8]},
}
json.dump(cert, open(os.path.join(root, "certs", "lemma_z5_euler.json"), "w"), indent=1)
print(len(cert["claims_in"]), "in,", len(cert["claims_out"]), "out; euler groups:", [e["order"] for e in euler])

#!/usr/bin/env python3
"""Job list for the LOCAL level-2 split of cube C1_c04 (stdlib only).

  python c1c04_prep.py

C1_c04 was refuted monolithically in the cloud, but the cake_lpr check of its 4.96 GB proof needs a 6500 MB heap, more than the 6 GB WSL VM of the local machine allows.
It is therefore split exactly like the seven cloud-timeout cubes (code/split_level2.py -> certs/split_C1_c04.json; completeness certificate certs/cover2_C1_c04 by
code/cover_cert.py level2 C1_c04) and its children are refuted and checked locally.
Child i = cube CNF (base full_C_c1.cnf + the 92 cube units of results/c_slice_jobs.json) + the 80 literals of the i-th section of the second plane as unit clauses;
id = C1_c04_p<plane>_k<i>.  Writes results/c1c04_jobs.json in the format of code/c1c04_driver.py (= a_slice_driver.py).  Idempotent."""
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
par = next(j for j in json.load(open(os.path.join(ROOT, "results", "c_slice_jobs.json"))) if j["id"] == "C1_c04")
sp = json.load(open(os.path.join(ROOT, "certs", "split_C1_c04.json")))
assert sp["id"] == "C1_c04" and sp["n"] == len(sp["children"]) and len(par["units"]) == 92
jobs = []
for i, ch in enumerate(sp["children"]):
    assert len(ch) == sp["assumption_size"] == 80
    jobs.append(dict(id="C1_c04_p%d_k%03d" % (sp["plane"], i), kind="cube", base=par["base"], units=par["units"] + ch, case=par["case"], cls=par["cls"], cube=par["cube"],
                     orbit_size=par["orbit_size"], plane=sp["plane"], k=par["k"], parent="C1_c04", child=i, plane2=sp["plane"], n_children=sp["n"]))
json.dump(jobs, open(os.path.join(ROOT, "results", "c1c04_jobs.json"), "w"))
print("c1c04_jobs.json: %d child jobs, plane %d" % (len(jobs), sp["plane"]))

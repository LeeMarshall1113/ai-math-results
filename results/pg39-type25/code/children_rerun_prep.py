#!/usr/bin/env python3
"""Seeded random sample of the 1,084 level-2 children for a LOCAL re-run (none had been re-run locally before).  Stdlib only.

  python children_rerun_prep.py [--seed 20261111] [--n 240]

Population: the 1,084 children (manifest tags <cube>_p<plane>_k<nnn>), each with a verified Modal record.  Sample: random.Random(seed).sample(sorted(tags), n).
Writes results/children_rerun_sample.json (seed, method, sample, expected CPU from the Modal records) and results/children_rerun_jobs.json
(job list in the format of code/a_slice_driver.py: id, kind, base (with .cnf), units).  The units are taken from results/modal_jobs.json, but their CNF hash is
checked against the manifest record by code/children_vs_local.py after the run (the driver records the sha256 of the CNF it really solved)."""
import json
import os
import random
import re
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 20261111
n = int(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 240
man = {}
for l in open(os.path.join(ROOT, "certs", "manifest.jsonl")):
    r = json.loads(l)
    if re.match(r"^[AC]\d_c\d\d_p\d+_k\d+$", r.get("tag", "")):
        man[r["tag"]] = r
assert len(man) == 1084
tags = sorted(man)
sample = random.Random(seed).sample(tags, n)          # kept in this random order: any prefix is itself a uniform random sample (so an early stop is unbiased)
jobs = {j["id"]: j for j in json.load(open(os.path.join(ROOT, "results", "modal_jobs.json")))}
out = []
for t in sample:
    j = jobs[t]
    assert j["kind"] == "cube" and len(j["units"]) == 172, (t, len(j["units"]))
    out.append(dict(id=t, kind="cube", base=j["base"] + ".cnf", units=j["units"], case=j["case"], cls=j["cls"], cube=j["cube"],
                    orbit_size=j["orbit_size"], plane=j["plane2"], k=j["k"]))
cpu = sum(man[t]["solver_cpu_seconds"] + man[t]["cake_seconds"] for t in sample) / 3600
json.dump(dict(seed=seed, n=n, population=len(tags), method="random.Random(seed).sample(sorted(manifest child tags), n)",
               expected_cpu_hours_from_modal_records=round(cpu, 2), sample=sample), open(os.path.join(ROOT, "results", "children_rerun_sample.json"), "w"), indent=1)
json.dump(out, open(os.path.join(ROOT, "results", "children_rerun_jobs.json"), "w"))
from collections import Counter
print("sample", n, "seed", seed, "expected CPU-hours (Modal solver + cake) %.2f" % cpu, dict(Counter(t[:6] for t in sample)))

#!/usr/bin/env python3
"""Compare the LOCAL re-run of the seeded sample of level-2 children (results/children_rerun_dry.jsonl, written by code/a_slice_driver.py --dry --out children_rerun)
with the Modal records of certs/manifest.jsonl, child by child.  Stdlib only; same pattern as code/modal_vs_local.py.

  python children_vs_local.py            exit code 1 on any disagreement; writes results/children_vs_local.json

Compared for every sampled child: status (UNSAT_VERIFIED), conflicts, lrat_bytes, lrat_sha256, cnf_sha256 (the driver hashes the CNF it really solved).
Sampled children without a local record yet are listed as 'not run' (not a disagreement)."""
import json
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sample = json.load(open(os.path.join(ROOT, "results", "children_rerun_sample.json")))
man = {}
for l in open(os.path.join(ROOT, "certs", "manifest.jsonl")):
    if l.strip():
        r = json.loads(l)
        man[r["tag"]] = r                                        # one Modal record per child tag (checked by replay/cover_asserts.py)
loc = {}
for l in open(os.path.join(ROOT, "results", "children_rerun_dry.jsonl")):
    try:
        r = json.loads(l)
    except Exception:
        continue
    if r.get("tag") in sample["sample"] and r.get("status") in ("UNSAT_VERIFIED", "UNSAT_CAKE_FAILED", "TIMEOUT_2h", "SKIPPED_SCRATCH", "SAT_EXACT_PASS", "SAT_EXACT_FAIL", "ERROR"):
        loc[r["tag"]] = r
agree, bad, notrun = [], [], []
cpu = 0.0
for t in sample["sample"]:
    r = loc.get(t)
    if r is None:
        notrun.append(t)
        continue
    cpu += (r.get("solver_cpu_seconds") or r.get("solver_seconds") or 0) + (r.get("cake_seconds") or 0)
    m = man[t]
    diffs = []
    if r.get("status") != "UNSAT_VERIFIED" or r.get("cake_lpr") != "s VERIFIED UNSAT":
        diffs.append("local status %s / %s" % (r.get("status"), r.get("cake_lpr")))
    for k in ("conflicts", "lrat_bytes", "lrat_sha256", "cnf_sha256"):
        if r.get(k) != m.get(k):
            diffs.append("%s: local %s / modal %s" % (k, r.get(k), m.get(k)))
    (bad if diffs else agree).append((t, diffs) if diffs else t)
for t, d in bad:
    print("DISAGREE", t, "; ".join(d))
print("SUMMARY: seed %d, sample %d of %d children: %d re-run locally and BIT-IDENTICAL to the Modal records (status, conflicts, LRAT size, LRAT sha256, CNF sha256), %d DISAGREEMENTS, %d not run; local CPU %.2f h (solver + cake_lpr)" % (
    sample["seed"], sample["n"], sample["population"], len(agree), len(bad), len(notrun), cpu / 3600))
from collections import Counter
print("per parent (identical):", dict(Counter(t[:6] for t in agree)))
json.dump(dict(seed=sample["seed"], n=sample["n"], identical=agree, disagree=bad, not_run=notrun, local_cpu_hours=round(cpu / 3600, 3)), open(os.path.join(ROOT, "results", "children_vs_local.json"), "w"), indent=1)
sys.exit(1 if bad else 0)

#!/usr/bin/env python3
"""Compare EVERY local record (results/a_slice.jsonl, c_slice.jsonl, children_rerun_dry.jsonl, full_rerun_dry.jsonl) with the Modal records, certificate by certificate.  Stdlib only.

  python full_vs_local.py     exit 1 on any disagreement; writes results/full_vs_local.json

Population: the 167 monolithic cubes and 1,084 children with a verified Modal record (best ledger record per id).  Compared: status UNSAT_VERIFIED + cake_lpr line, conflicts,
lrat_bytes, lrat_sha256, cnf_sha256 against the ledger AND (conflicts, sizes, hashes) against certs/manifest.jsonl.  Not run / timed out / cake failures are listed separately."""
import json
import os
import re
import sys
from collections import Counter

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = lambda *a: os.path.join(ROOT, *a)
led = {}
for l in open(R("results", "modal_ledger.jsonl")):
    try:
        r = json.loads(l)
    except Exception:
        continue
    if r.get("manual") or "id" not in r:
        continue
    o = led.get(r["id"])
    if o and o.get("status") == "UNSAT_VERIFIED" and r.get("status") != "UNSAT_VERIFIED":
        continue
    led[r["id"]] = r
man = {}
for l in open(R("certs", "manifest.jsonl")):
    if l.strip():
        r = json.loads(l)
        if r.get("status", "UNSAT_VERIFIED") != "TIMEOUT" and r.get("cake_lpr", "s VERIFIED UNSAT") == "s VERIFIED UNSAT" and r.get("lrat_bytes") is not None:
            man[r["tag"]] = r
pop = sorted(t for t, r in led.items() if r.get("status") == "UNSAT_VERIFIED" and re.match(r"^[AC]\d_c\d\d(_p\d+_k\d+)?$", t))
loc = {}
for f in ("a_slice.jsonl", "c_slice.jsonl", "children_rerun_dry.jsonl", "full_rerun_dry.jsonl"):
    p = R("results", f)
    if not os.path.exists(p):
        continue
    for l in open(p):
        try:
            r = json.loads(l)
        except Exception:
            continue
        t = r.get("tag")
        if t in led and r.get("status"):
            if t in loc and loc[t].get("status") == "UNSAT_VERIFIED" and r["status"] != "UNSAT_VERIFIED":
                continue
            loc[t] = r
agree, bad, notrun, other = [], [], [], []
cpu = 0.0
for t in pop:
    r = loc.get(t)
    if r is None:
        notrun.append(t)
        continue
    cpu += (r.get("solver_cpu_seconds") or r.get("solver_seconds") or 0) + (r.get("cake_seconds") or 0)
    if r["status"] != "UNSAT_VERIFIED":
        other.append((t, r["status"], r.get("flag", "")))
        continue
    m, mm = led[t], man.get(t)
    d = []
    if r.get("cake_lpr") != "s VERIFIED UNSAT":
        d.append("cake_lpr %s" % r.get("cake_lpr"))
    for k in ("conflicts", "lrat_bytes", "lrat_sha256", "cnf_sha256"):
        if r.get(k) != m.get(k):
            d.append("%s: local %s / ledger %s" % (k, r.get(k), m.get(k)))
        if mm is not None and r.get(k) != mm.get(k):
            d.append("%s: local %s / manifest %s" % (k, r.get(k), mm.get(k)))
    (bad if d else agree).append((t, d) if d else t)
for t, d in bad:
    print("DISAGREE", t, "; ".join(d))
for o in other:
    print("NOT VERIFIED LOCALLY", *o)
cubes = lambda L: sum(1 for x in L if "_p" not in (x[0] if isinstance(x, tuple) else x))
print("SUMMARY: %d certificates (%d cubes, %d children): %d re-run locally and bit-identical (%d cubes, %d children), %d DISAGREEMENTS, %d timed out / cake-failed / other, %d not run; local CPU %.1f h (solver + cake_lpr, over records of all runs)" % (
    len(pop), sum(1 for t in pop if "_p" not in t), sum(1 for t in pop if "_p" in t), len(agree), cubes(agree), len(agree) - cubes(agree), len(bad), len(other), len(notrun), cpu / 3600))
print("not run:", Counter("cube" if "_p" not in t else "child" for t in notrun))
json.dump(dict(identical=agree, disagree=bad, other=other, not_run=notrun, local_cpu_hours=round(cpu / 3600, 2)), open(R("results", "full_vs_local.json"), "w"), indent=1)
sys.exit(1 if bad or other else 0)

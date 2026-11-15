#!/usr/bin/env python3
"""Job lists for the FULL local regeneration: every verified certificate (167 monolithic cubes, 1,084 children) that has no identical local record yet.  Stdlib only.

  python full_rerun_prep.py

Population: best Modal ledger record per id (results/modal_ledger.jsonl, a verified record is never replaced by a later failure) with status UNSAT_VERIFIED, of jobs in results/modal_jobs.json
(monolithic cubes A*_cNN / C*_cNN and children *_pP_kNNN).  Already re-run locally (UNSAT_VERIFIED in results/a_slice.jsonl, c_slice.jsonl, children_rerun_dry.jsonl) are left out.
Each job carries `expect` (conflicts, lrat_bytes, lrat_sha256, cnf_sha256 of the Modal record: the driver stops the batch on any difference) and `local_tlim`
= max(1800, 4 x recorded Modal solver wall seconds).  Order: smallest recorded LRAT first.  Stages (written to results/full_rerun_stage*.json, same --out name full_rerun so completed tags are shared):
  1 : all children + cubes with LRAT <= 2.5 GB      (3 workers)
  2 : cubes with LRAT in (2.5, 4] GB                  (2 workers)
  3 : cubes with LRAT > 4 GB, one at a time           (1 worker; C1_c04 excluded: its cake_lpr check needs a 6500 MB heap, above the 6 GB WSL cap)"""
import json
import os
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = lambda *a: os.path.join(ROOT, "results", *a)
led = {}
for l in open(R("modal_ledger.jsonl")):
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
done = set()
for f in ("a_slice.jsonl", "c_slice.jsonl", "children_rerun_dry.jsonl", "full_rerun_dry.jsonl"):
    if os.path.exists(R(f)):
        for l in open(R(f)):
            try:
                r = json.loads(l)
            except Exception:
                continue
            if r.get("status") == "UNSAT_VERIFIED" and r.get("tag"):
                done.add(r["tag"])
jobs = {j["id"]: j for j in json.load(open(R("modal_jobs.json")))}
pop = [t for t, r in led.items() if r.get("status") == "UNSAT_VERIFIED" and t in jobs and re.match(r"^[AC]\d_c\d\d(_p\d+_k\d+)?$", t)]
cubes = [t for t in pop if "_p" not in t]
kids = [t for t in pop if "_p" in t]
print("population: %d cubes, %d children (expect 167, 1084)" % (len(cubes), len(kids)))
todo = sorted((t for t in pop if t not in done), key=lambda t: (led[t]["lrat_bytes"], t))
EXCL = {"C1_c04": "cake_lpr heap 6500 MB needed, above the 6 GB WSL cap"}
GB = 2 ** 30
stages = {1: [], 2: [], 3: []}
excluded = []
for t in todo:
    if t in EXCL:
        excluded.append((t, EXCL[t]))
        continue
    r, j = led[t], jobs[t]
    job = dict(id=t, kind="cube", base=j["base"] + ".cnf", units=j["units"], case=j["case"], cls=j["cls"], cube=j["cube"], orbit_size=j["orbit_size"],
               plane=j.get("plane2", j.get("plane", 0)), k=j["k"],
               local_tlim=int(max(1800, 4 * r["solver_wall"])),
               expect=dict(conflicts=r["conflicts"], lrat_bytes=r["lrat_bytes"], lrat_sha256=r["lrat_sha256"], cnf_sha256=r["cnf_sha256"]))
    if "plane" not in j and "plane2" not in j:
        job["plane"] = 0
    b = r["lrat_bytes"]
    stages[1 if ("_p" in t or b <= 2.5 * GB) else (2 if b <= 4 * GB else 3)].append(job)
cost = lambda t: (led[t].get("solver_cpu") or 0) + (led[t].get("cake_wall") or 0)
for s, L in stages.items():
    json.dump(L, open(R("full_rerun_stage%d.json" % s), "w"))
    print("stage %d: %d jobs (%d cubes), recorded cost %.1f CPU-h" % (s, len(L), sum(1 for x in L if "_p" not in x["id"]), sum(cost(x["id"]) for x in L) / 3600))
print("excluded:", excluded)
print("already done locally:", len(done & set(pop)), " to do:", len(todo))

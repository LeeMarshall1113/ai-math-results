#!/usr/bin/env python3
"""Close the LOCAL split of C1_c04 in the manifest (stdlib only; run after code/run_c1c04.sh has finished, with nothing else writing certs/manifest.jsonl).

  python c1c04_finish.py [--dry]

1. every child of certs/split_C1_c04.json must have a verified record in results/c1c04_manifest_local.jsonl (status UNSAT_VERIFIED, cake_lpr 's VERIFIED UNSAT'), the CNF hash
   of which is checked by replay/rebuild_all.py afterwards;  the records (via = local) are appended to certs/manifest.jsonl;
2. the original cloud record of the monolithic cube C1_c04 stays in the manifest, marked superseded_for_verification_by (its cake_lpr check needs a 6500 MB heap > the 6 GB WSL cap);
3. certs/PROOF-INDEX.tsv is rewritten (code/make_proof_index.py)."""
import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
MAN = os.path.join(ROOT, "certs", "manifest.jsonl")
DRY = "--dry" in sys.argv
sp = json.load(open(os.path.join(ROOT, "certs", "split_C1_c04.json")))
want = ["C1_c04_p%d_k%03d" % (sp["plane"], i) for i in range(sp["n"])]
loc = {}
for l in open(os.path.join(ROOT, "results", "c1c04_manifest_local.jsonl")):
    r = json.loads(l)
    if r.get("status") == "UNSAT_VERIFIED" and r.get("cake_lpr") == "s VERIFIED UNSAT" and r.get("via") == "local":
        loc[r["tag"]] = r
missing = [t for t in want if t not in loc]
print("children %d, verified local records %d, missing %d" % (len(want), len(loc), len(missing)))
if missing:
    print("MISSING:", missing[:10])
    sys.exit(1)
lines = open(MAN).read().splitlines()
have = {json.loads(l).get("tag") for l in lines if l.strip()}
out = []
for l in lines:
    r = json.loads(l) if l.strip() else None
    if r and r.get("tag") == "C1_c04" and r.get("via") == "modal":
        r["superseded_for_verification_by"] = "split C1_c04: cover2_C1_c04 + local children"
        r["superseded_note"] = ("cloud record kept as a second, independent refutation of the cube; its cake_lpr check needed a 6500 MB heap (over the 6 GB WSL cap), so the local proof "
                                "of the cube is the level-2 split (cover2_C1_c04 + the local children, all checked with cake_lpr heap <= 4500 MB)")
        l = json.dumps(r)
    out.append(l)
new = [json.dumps(loc[t]) for t in want if t not in have]
print("appending %d local child records" % len(new))
if not DRY:
    open(MAN, "w", newline="\n").write("\n".join(out + new) + "\n")
    subprocess.check_call([sys.executable, os.path.join(ROOT, "code", "make_proof_index.py")])

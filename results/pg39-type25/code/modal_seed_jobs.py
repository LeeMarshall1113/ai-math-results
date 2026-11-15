"""Add seed variants of cubes to results/modal_jobs.json (same units, extra option --seed=K, id = <cube>_s<K>).
  python modal_seed_jobs.py SEED ID1,ID2,...|@file [TLIM] [LRAT: bin|none]
A variant that closes with an LRAT checked by cake_lpr is a valid certificate for its cube (any seed gives a proof of the same CNF).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def main():
    seed = int(sys.argv[1])
    ids = sys.argv[2]
    want = [l.strip() for l in open(ids[1:]) if l.strip()] if ids.startswith("@") else ids.split(",")
    tlim = int(sys.argv[3]) if len(sys.argv) > 3 else 7200
    lrat = sys.argv[4] if len(sys.argv) > 4 else "bin"
    path = os.path.join(ROOT, "results", "modal_jobs.json")
    jobs = json.load(open(path))
    have = {j["id"] for j in jobs}
    by = {j["id"]: j for j in jobs}
    add = []
    for w in want:
        j = by[w]
        nid = "%s_s%d" % (w, seed)
        if nid in have:
            continue
        add.append(dict(j, id=nid, opts=["--seed=%d" % seed], tlim=tlim, lrat=lrat))
    json.dump(jobs + add, open(path, "w"))
    open(os.path.join(ROOT, "results", "modal_batch_%s_s%d.txt" % (os.path.basename(ids).replace("@", "").replace(".txt", ""), seed)), "w").write("\n".join(a["id"] for a in add) + "\n")
    print("added %d variants: %s" % (len(add), [a["id"] for a in add]))


if __name__ == "__main__":
    main()

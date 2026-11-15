"""Turn the level-2 splits of unclosed cubes (certs/split_<ID>.json, from split_level2.py) into Modal jobs.
  python modal_children_jobs.py ID1,ID2,... [TLIM]
Child i of cube ID = the cube CNF (base + the 92 cube units) plus the 80 literals of the i-th section of the second plane as unit clauses; id = <ID>_p<plane>_k<i>.
Appended to results/modal_jobs.json; the id list is written to results/modal_batch_children.txt.  A cube is closed when ALL its children are UNSAT_VERIFIED and the
list of sections is complete (the completeness LRAT of the planar enumeration is a separate certificate, NOTES 5.4 item 3).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def main():
    ids = sys.argv[1].split(",")
    tlim = int(sys.argv[2]) if len(sys.argv) > 2 else 7200
    path = os.path.join(ROOT, "results", "modal_jobs.json")
    jobs = json.load(open(path))
    by = {j["id"]: j for j in jobs}
    have = set(by)
    add = []
    for jid in ids:
        sp = os.path.join(ROOT, "certs", "split_%s.json" % jid)
        if not os.path.exists(sp):
            print("no split file for", jid)
            continue
        sp = json.load(open(sp))
        par = by[jid]
        for i, ch in enumerate(sp["children"]):
            nid = "%s_p%d_k%03d" % (jid, sp["plane"], i)
            if nid in have:
                continue
            add.append(dict(par, id=nid, units=par["units"] + ch, tlim=tlim, lrat="bin", parent=jid, child=i, plane2=sp["plane"], n_children=sp["n"]))
    json.dump(jobs + add, open(path, "w"))
    open(os.path.join(ROOT, "results", "modal_batch_children.txt"), "w").write("\n".join(a["id"] for a in add) + "\n")
    print("added %d child jobs for %s" % (len(add), ids))


if __name__ == "__main__":
    main()

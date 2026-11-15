"""Level-2 split of a cube that was not closed after the time limit (stop rule: split and log, do not run on).

  python split_level2.py JOB_ID [NPL] [CAP]

The cube JOB_ID (results/a_slice_jobs.json or c_slice_jobs.json) fixes the 91 points of one plane through P0 (a section).  A second plane through P0 meets it in a line;
its section is enumerated by SAT (Setup.sections, same planar encoder as the level-1 sections) among the first NPL planes (default 10), and the
plane with the fewest sections >= 2 (at most CAP, default 300) is kept.  The children are the section assignments of the free points of that plane
(signed literals, point p <-> variable p+1).  Result: certs/split_<JOB_ID>.json and one line in results/a_slice_splits.jsonl.  Nothing is solved here.
The children are a complete cover of the cube only if the enumeration was complete (cap not reached); a certificate for completeness is the planar CNF with
blocking clauses (as for level 1) and is produced when the children are actually run.
"""
import json
import os
import sys
import time

from trivial_tree import Setup

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def main():
    jid = sys.argv[1]
    npl = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 300
    job = None
    for fn in ("a_slice_jobs.json", "c_slice_jobs.json"):
        job = next((j for j in json.load(open(os.path.join(ROOT, "results", fn))) if j["id"] == jid), None)
        if job is not None:
            break
    assert job is not None, "unknown job %s" % jid
    S = Setup(job["case"], job["cls"])
    fixed = {abs(u) - 1: int(u > 0) for u in job["units"]}
    L1 = job["plane"]
    t0 = time.time()
    best, tried, sizes = None, 0, []
    for L in range(len(S.P2.lines)):
        if L == L1:
            continue
        secs = S.sections(L, fixed, cap=cap + 1)
        tried += 1
        sizes.append(len(secs))
        if 2 <= len(secs) <= cap and (best is None or len(secs) < len(best[1])):
            best = (L, secs)
        if tried >= npl or (best is not None and len(best[1]) <= 200):          # good enough: at most 200 children (the 140-section planes of size 26)
            break
    rec = dict(id=jid, planes_tried=tried, sections_per_plane=sizes, secs=round(time.time() - t0, 1), time=time.strftime("%Y-%m-%d %H:%M:%S"))
    if best is None:
        rec["status"] = "NO_SPLIT_PLANE_FOUND (all tried planes have 1 or more than %d sections)" % cap
    else:
        L, secs = best
        free_pts = sorted(set(secs[0]) - set(fixed) - set(S.pt_units))
        children = [[((p + 1) if sc[p] else -(p + 1)) for p in free_pts] for sc in secs]
        out = os.path.join(ROOT, "certs", "split_%s.json" % jid)
        json.dump(dict(id=jid, plane=L, n=len(children), assumption_size=len(free_pts), children=children), open(out, "w"))
        rec.update(status="SPLIT", plane=L, children=len(children), assumption_size=len(free_pts), file=os.path.basename(out))
    with open(os.path.join(ROOT, "results", "a_slice_splits.jsonl"), "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()

"""Level-2 children for the incremental test (speed check): sections of a second plane through P0, chosen as the plane
(among the first NPL planes tried) with the fewest sections >= 2.
  python bench_children.py CASE CLASS SEED_OFF NTH OUT.json [NPL] [CAP]
"""
import json
import sys
import time

from bench_make import cube_lits
from trivial_tree import Setup


def main():
    case, ci, off, nth, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    npl = int(sys.argv[6]) if len(sys.argv) > 6 else 8
    cap = int(sys.argv[7]) if len(sys.argv) > 7 else 300
    S = Setup(case, ci)
    lits, Ls, sec = cube_lits(S, ci, off, nth)
    best = None
    t0 = time.time()
    tried = 0
    for L in range(len(S.P2.lines)):
        if L == Ls[0]:
            continue
        secs = S.sections(L, sec, cap=cap + 1)
        tried += 1
        print("plane %d: %d sections%s (%.0fs)" % (L, len(secs), "+" if len(secs) > cap else "", time.time() - t0), flush=True)
        if 2 <= len(secs) <= cap and (best is None or len(secs) < len(best[1])):
            best = (L, secs)
        if tried >= npl:
            break
    L, secs = best
    fixed_pts = set(sec) | set(S.pt_units)
    children = [[((p + 1) if v else -(p + 1)) for p, v in sc.items() if p not in fixed_pts] for sc in secs]
    json.dump(dict(plane=L, n=len(secs), children=children), open(out, "w"))
    print("chosen plane %d: %d children, %d assumptions each (%.0fs)" % (L, len(secs), len(children[0]), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()

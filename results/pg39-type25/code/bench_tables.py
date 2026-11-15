"""Plane-table variant of a benchmark cube (speed check).

For every plane pi through P0 other than the cube plane, the section K n pi is a planar (2,5)-set containing P0 whose pencil at P0
is prescribed by Sigma, whose size is prescribed, and whose intersection with the cube plane (a line through P0) is already fixed by the cube.
All such sections are enumerated by SAT (planar CNF, same encoder) and added to the cube CNF as selector clauses: one selector per allowed
section, each selector implies the 91 point values, and at least one selector holds.  SOUND iff the enumerations are complete
(each would need its own LRAT completeness certificate; here we only measure the speed-up).
  python bench_tables.py CASE CLASS SEED_OFF NTH OUT.cnf [CAP]
"""
import json
import os
import sys
import time

from bench_make import cube_lits
from trivial_tree import Setup, write_cube_cnf

HERE = os.path.dirname(os.path.abspath(__file__))
CNF = os.path.join(HERE, "..", "cnf")


def main():
    case, ci, off, nth, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    cap = int(sys.argv[6]) if len(sys.argv) > 6 else 400
    S = Setup(case, ci)
    lits, Ls, sec = cube_lits(S, ci, off, nth)
    L1 = Ls[0]
    base = os.path.join(CNF, "full_%s_c%d.cnf" % (case, ci))
    tmp = out + ".tmp"
    write_cube_cnf(base, lits, tmp)
    with open(tmp) as fh:
        line = fh.readline()
        while line.startswith("c"):
            line = fh.readline()
        parts = line.split()
        nv, nc = int(parts[2]), int(parts[3])
    extra = []
    nsel = 0
    sizes = []
    best = None
    t0 = time.time()
    for L in range(len(S.P2.lines)):
        if L == L1:
            continue
        secs = S.sections(L, sec, cap=cap + 1)
        sizes.append(len(secs))
        if 2 <= len(secs) <= cap and (best is None or len(secs) < len(best[1])):
            best = (L, secs)
        if len(secs) > cap or not secs:
            continue                      # incomplete (or empty: then the cube is already refuted) tables are not used
        sels = []
        for sc in secs:
            nsel += 1
            s = nv + nsel
            sels.append(s)
            for p, v in sc.items():
                extra.append("%d %d 0\n" % (-s, (p + 1) if v else -(p + 1)))
        extra.append(" ".join(map(str, sels)) + " 0\n")
    with open(tmp) as fh, open(out, "w", newline="\n") as o:
        line = fh.readline()
        while line.startswith("c"):
            o.write(line)
            line = fh.readline()
        o.write("p cnf %d %d\n" % (nv + nsel, nc + len(extra)))
        while True:
            ch = fh.read(1 << 22)
            if not ch:
                break
            o.write(ch)
        o.writelines(extra)
    os.remove(tmp)
    if best is not None:
        L, secs = best
        fixed_pts = set(sec) | set(S.pt_units)
        children = [[((p + 1) if v else -(p + 1)) for p, v in sc.items() if p not in fixed_pts] for sc in secs]
        json.dump(dict(plane=L, n=len(secs), children=children), open(out + '.children.json', 'w'))
        print('level-2 split plane %d: %d children' % (L, len(secs)), flush=True)
    import collections
    print("planes tabled: %d of %d ; sections per plane: min %d median %d max %d ; %d selectors, %d extra clauses (%.0fs)" % (
        sum(1 for s in sizes if 0 < s <= cap), len(sizes), min(sizes), sorted(sizes)[len(sizes) // 2], max(sizes), nsel, len(extra), time.time() - t0), flush=True)
    print("zero-section planes (cube refuted by planes alone):", sum(1 for s in sizes if s == 0), flush=True)


if __name__ == "__main__":
    main()

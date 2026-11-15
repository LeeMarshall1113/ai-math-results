#!/usr/bin/env python3
"""Exact check of a candidate set K of PG(3,9) (stdlib only; independent of ../code).

  python check_set.py FILE.json [--allowed 2,5]

FILE.json : a list of point indices (our lexicographic indexing of normalised vectors) OR a list of 4-tuples of GF(9) elements
(pairs [a,b] meaning a+b*i, i^2=-1).  Two independent line enumerations are used:
  (1) row-reduced 2x4 matrices (replay.py),
  (2) brute force over all pairs of points (span of the pair) -- different algorithm, same field code.
For every one of the 7462 lines the number of points of K is counted and must lie in the allowed set.
"""
import itertools
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay as R


def lines_by_pairs(points):
    pidx = {p: i for i, p in enumerate(points)}
    seen = set()
    lines = []
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            u, v = points[i], points[j]
            members = {j}
            for b in range(9):
                members.add(pidx[R.normalize(tuple(R.e_add(x, R.e_mul(b, y)) for x, y in zip(u, v)))])
            L = tuple(sorted(members))
            if L not in seen:
                seen.add(L)
                lines.append(L)
    return lines


def main():
    path = sys.argv[1]
    allowed = {2, 5}
    if "--allowed" in sys.argv:
        allowed = set(int(x) for x in sys.argv[sys.argv.index("--allowed") + 1].split(","))
    data = json.load(open(path))
    points = R.build_points()
    pidx = {p: i for i, p in enumerate(points)}
    if data and isinstance(data[0], list):
        K = sorted(pidx[R.normalize(tuple(R.decode(e) for e in v))] for v in data)
    else:
        K = sorted(data)
    S = set(K)
    results = []
    for name, lines in (("rref", R.build_lines_rref(points)), ("pairs", None)):
        if lines is None:
            lines = lines_by_pairs(points)
        assert len(lines) == 7462
        bad = [L for L in lines if sum(1 for p in L if p in S) not in allowed]
        hist = {}
        for L in lines:
            c = sum(1 for p in L if p in S)
            hist[c] = hist.get(c, 0) + 1
        print("%-5s : |K| = %d ; lines %d ; line-count histogram %s ; violations %d" % (name, len(S), len(lines), dict(sorted(hist.items())), len(bad)))
        results.append(not bad)
    print("EXACT CHECK:", "PASS" if all(results) else "FAIL")
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main()

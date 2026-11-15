#!/usr/bin/env python3
"""Counting lemma for the two prime-order classes whose invariant lines partition PG(3,9)
(stdlib only; uses the independent geometry code in replay.py).

  Class 5a (g = multiplication by a 5th root of unity of GF(81) on GF(81)^2) and
  class 2c (g^2 = non-square scalar, i.e. multiplication by sqrt(lambda) in GF(81) on GF(81)^2):
  every GF(81)-line is g-invariant; these 82 lines (the Desarguesian spread) partition the 820 points
  and g acts on each of them without fixed points, with all orbits of one size s (s = 5 resp. 2).
  If K is g-invariant and meets every line in 2 or 5 points, then on a spread line |K n l| is a
  multiple of s lying in {2,5}: s = 5 forces 5, s = 2 forces 2.  Hence |K| = 82*5 = 410 resp. 82*2 = 164,
  which is neither 260 nor 287 (the only possible sizes of a (2,5)-set in PG(3,9), params_table.txt).
  The complement of a (5,8)-set is a (2,5)-set, so (5,8)-sets are covered as well.

This script verifies the geometric facts by exhaustive computation:
  1. the g-invariant lines are exactly 82, pairwise disjoint, covering all 820 points;
  2. g acts on each of them with all orbits of size s.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay as R

CERTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certs")


def check(name, s_expected, spread=True):
    classes = json.load(open(os.path.join(CERTS, "classes.json")))
    cls = classes[name]
    points = R.build_points()
    pidx = {p: i for i, p in enumerate(points)}
    lines = R.build_lines_rref(points)
    M = [[R.decode(e) for e in row] for row in cls["matrix"]]
    perm = R.point_perm(points, pidx, M, cls["frob"])
    lset = {L: i for i, L in enumerate(lines)}
    inv = [L for L in lines if tuple(sorted(perm[p] for p in L)) == L]
    cover = sorted(p for L in inv for p in L)
    print("%s: g-invariant lines: %d ; pairwise disjoint and covering all points: %s" % (
        name, len(inv), cover == list(range(len(points)))))
    assert len(inv) == 82 and cover == list(range(len(points)))
    for L in inv:
        orbs, _ = R.orbits_of(len(points), [perm])
        # orbits of g restricted to L
        seen = set()
        for p in L:
            if p in seen:
                continue
            o = [p]
            q = perm[p]
            while q != p:
                o.append(q)
                q = perm[q]
            assert all(x in L for x in o)
            seen.update(o)
            assert len(o) == s_expected, (name, len(o))
    allowed = {2, 5}
    mult = {t for t in range(0, 11, s_expected)} & allowed
    print("   all orbits on the invariant lines have size %d -> |K n l| in %s -> |K| = 82 * %s" % (
        s_expected, sorted(mult), sorted(mult)))
    assert len(mult) == 1
    size = 82 * next(iter(mult))
    print("   |K| would be %d, but only 260 or 287 are possible : NO invariant K" % size)
    assert size not in (260, 287)


if __name__ == "__main__":
    check("5a_scalar_GF81", 5)
    check("2c_fixedpointfree", 2)

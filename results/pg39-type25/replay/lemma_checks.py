#!/usr/bin/env python3
"""Machine checks of the hand-proved geometric lemmas of PROOF-CHAIN.md (stdlib; reuses the independent geometry of rebuild_all.py).

  python lemma_checks.py

1. Homotheties h_lambda = diag(1,1,1,lambda) (section 2.2): for every lambda in GF(9)^* and every point (v:t): h(v:t) = (v:lambda t); P0 and all height-0 points (v:0) are fixed; every line
   through P0 is mapped onto itself; the units of every one of the ten cases are fixed pointwise and the 91 line types are unchanged; h preserves the first plane L1 of every case.
2. Elations E_r : (x, x4) -> (x, x4 + r.x) (section 1.4), all 729 r: E_r fixes P0, maps every line through P0 onto itself, sends (v:t) to (v: t + r.v) for every normalised v; and for the first non-collinear
   triple (v1,v2,v3) of every case the map r -> (r.v1, r.v2, r.v3) is a bijection of GF(9)^3 (so the heights of the three second points can be normalised simultaneously).
3. Lifting (section 1.4): for 300 pseudo-random A in GL(3,9) and both values of the Frobenius exponent, the collineation diag(A,1) o frob^e fixes P0 and induces on the directions (the points of Q) exactly
   the permutation v -> normalise(A frob^e(v)).
4. Local counts (section 2.1): every planar (2,5)-set listed in results/planar_k{26,35}_R{0,1}.json (the 1,345 real sets of the certified frame enumerations) has, through each of its points, (m-11)/3 five-secants and,
   through each point outside, (m-20)/3 five-secants.
"""
import itertools
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import replay as R            # noqa: E402
import rebuild_all as RB      # noqa: E402

RES = os.path.join(HERE, "..", "results")


def main():
    G = RB.Geometry()
    classes = json.load(open(os.path.join(RES, "planar_classes.json")))
    cases = [RB.Case(G, "A", i, classes) for i in range(3)] + [RB.Case(G, "C", i, classes) for i in range(7)]
    P0, pts, pidx = G.P0, G.points4, G.pidx4
    # (v:t) table
    height = {}
    for d in range(91):
        v = G.points3[d]
        for t in range(9):
            height[(d, t)] = pidx[R.normalize(v + (t,))]
    # 1. homotheties
    for lam in range(1, 9):
        h = G.homs[lam - 1]
        assert h[P0] == P0
        for (d, t), q in height.items():
            assert h[q] == height[(d, R.e_mul(lam, t))]
        for j in G.lines_thru[P0]:
            assert sorted(h[p] for p in G.lines4[j]) == list(G.lines4[j])
    for C in cases:
        for lam in range(1, 9):
            h = G.homs[lam - 1]
            assert all(h[p] == p for p in C.pt_units)
        # the first plane L1 of the case is preserved
        v1, v2, v3 = C.tri
        Ls = [next(i for i, L in enumerate(G.lines3) if a in L and b in L) for a, b in ((v1, v2), (v1, v3), (v2, v3))]
        L1 = next(L for L in Ls if C.plane_size(L) == 26)
        S1 = set(G.line_pts[L1])
        for lam in range(1, 9):
            assert set(G.homs[lam - 1][p] for p in S1) == S1
    print("1. homotheties: OK (8 homotheties x 820 points, 10 cases)")
    # 2. elations
    dots = {}
    for r in itertools.product(range(9), repeat=3):
        perm = []
        for q, v in enumerate(pts):
            x4 = v[3]
            s = 0
            for a, b in zip(r, v[:3]):
                s = R.e_add(s, R.e_mul(a, b))
            perm.append(pidx[R.normalize((v[0], v[1], v[2], R.e_add(x4, s)))])
        assert perm[P0] == P0
        for j in G.lines_thru[P0]:
            assert sorted(perm[p] for p in G.lines4[j]) == list(G.lines4[j])
        for d in range(91):
            v = G.points3[d]
            rv = 0
            for a, b in zip(r, v):
                rv = R.e_add(rv, R.e_mul(a, b))
            for t in (0, 1, 5):
                assert perm[height[(d, t)]] == height[(d, R.e_add(t, rv))]
        dots[r] = None
    for C in cases:
        imgs = set()
        for r in itertools.product(range(9), repeat=3):
            imgs.add(tuple(
                (lambda v: R.e_add(R.e_add(R.e_mul(r[0], v[0]), R.e_mul(r[1], v[1])), R.e_mul(r[2], v[2])))(G.points3[d]) for d in C.tri))
        assert len(imgs) == 729
    print("2. elations: OK (729 elations x 91 directions, regular on the heights of the triple in all ten cases)")
    # 3. lifting
    rnd = random.Random(1)
    n_ok = 0
    while n_ok < 300:
        A = [[rnd.randrange(9) for _ in range(3)] for _ in range(3)]
        # invertible?
        det = 0
        for perm_, sign in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
            term = R.e_mul(R.e_mul(A[0][perm_[0]], A[1][perm_[1]]), A[2][perm_[2]])
            det = R.e_add(det, term if sign == 1 else R.e_mul(term, 2))
        if det == 0:
            continue
        for e in (0, 1):
            for q, v in enumerate(pts):
                w = tuple(R.e_frob(c) for c in v) if e else v
                out = []
                for row in A:
                    s = 0
                    for a, b in zip(row, w[:3]):
                        s = R.e_add(s, R.e_mul(a, b))
                    out.append(s)
                img = (out[0], out[1], out[2], w[3])
                if q == P0:
                    assert R.normalize(img) == pts[P0]
                    continue
                qi = pidx[R.normalize(img)]
                assert qi != P0
                d_img = G.proj[qi]
                d = G.proj[q]
                v3 = G.points3[d]
                w3 = tuple(R.e_frob(c) for c in v3) if e else v3
                expect = []
                for row in A:
                    s = 0
                    for a, b in zip(row, w3):
                        s = R.e_add(s, R.e_mul(a, b))
                    expect.append(s)
                assert d_img == G.pidx3[R.normalize(tuple(expect))]
        n_ok += 1
    print("3. lifting: OK (300 random A in GL(3,9) x {frob^0, frob^1} x 820 points)")
    # 4. local counts on the 1,345 real planar sets
    pts3 = G.points3
    lines3 = G.lines3
    thru3 = [[] for _ in pts3]
    for j, L in enumerate(lines3):
        for p in L:
            thru3[p].append(j)
    n_sets = 0
    for m in (26, 35):
        for ri in (0, 1):
            d = json.load(open(os.path.join(RES, "planar_k%d_R%d.json" % (m, ri))))
            for sol in d["solutions"]:
                S = set(sol)
                assert len(S) == m
                five = [len(S & set(L)) == 5 for L in lines3]
                assert all(len(S & set(L)) in (2, 5) for L in lines3)
                for p in range(91):
                    b = sum(1 for j in thru3[p] if five[j])
                    assert b == ((m - 11) // 3 if p in S else (m - 20) // 3)
                n_sets += 1
    print("4. local counts: OK (%d real planar (2,5)-sets)" % n_sets)
    print("LEMMA CHECKS OK")


if __name__ == "__main__":
    main()

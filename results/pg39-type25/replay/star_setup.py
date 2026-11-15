#!/usr/bin/env python3
"""Independent (stdlib, replay.py geometry) derivation of the star set-up of section 5.1 of NOTES.md for the 10 top cases
A0,A1,A2 (k = 287, |Sigma| = 26) and C0..C6 (k = 260, |Sigma| = 35):

  python star_setup.py            : prints, per case, the base units and verifies the structural facts used
  python star_setup.py --compare  : additionally compares with the main code (code/trivial_tree.Setup): units, line types, triple, plane sizes

Verified facts (per case): Sigma (class representative, certs/planar_classes.json) is a (2,5)-set of the quotient plane PG(2,9) of the right size;
the plane sizes (26 iff |Sigma n L| = 5) and their number (26 resp. 56); the triple of normalised directions is non-collinear; the point units are
P0 and the three 'height 0' points (v:0); the line units are 0 (two-secant) on the lines of Sigma and 1 on the other 65 (resp. 56).
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay as R

CERTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def build_points3():
    pts = []
    import itertools
    for v in itertools.product(range(9), repeat=3):
        nz = next((c for c in v if c), None)
        if nz == 1:
            pts.append(v)
    return pts


def lines3(points3):
    """lines of PG(2,9): {x : a.x = 0} for the 91 normalised a"""
    out = []
    for a in points3:
        L = []
        for i, v in enumerate(points3):
            s = 0
            for x, y in zip(a, v):
                s = R.e_add(s, R.e_mul(x, y))
            if s == 0:
                L.append(i)
        out.append(tuple(L))
    return out


def noncollinear_triple(S, lines):
    S = sorted(S)
    for a in S:
        for b in S:
            if b <= a:
                continue
            Lab = next(L for L in lines if a in L and b in L)
            for c in S:
                if c > b and c not in Lab:
                    return (a, b, c)
    raise RuntimeError


def setup(points4, pidx4, points3, pidx3, lines_q, rep, case):
    P0 = pidx4[(0, 0, 0, 1)]
    inK = True
    Sigma = set(rep)
    # (2,5)-set check
    cnt = {}
    for L in lines_q:
        c = len(Sigma & set(L))
        cnt[c] = cnt.get(c, 0) + 1
    assert set(cnt) == {2, 5}, cnt
    # lines through P0 : direction of the line = normalised (x1,x2,x3) of any other point
    dir_of_line = {}
    for q, v in enumerate(points4):
        if q == P0:
            continue
        d = pidx3[R.normalize(v[:3])]
        dir_of_line.setdefault(d, []).append(q)
    assert len(dir_of_line) == 91 and all(len(v) == 9 for v in dir_of_line.values())
    tri = noncollinear_triple(Sigma, lines_q)
    zero = {d: pidx4[R.normalize(points3[d] + (0,))] for d in range(91)}
    point_units = {P0: 1}
    for d in tri:
        point_units[zero[d]] = 1
    line_y = {d: (0 if d in Sigma else 1) for d in range(91)}        # y = 1 : five-secant
    sizes = {}
    for i, L in enumerate(lines_q):
        s5 = len(Sigma & set(L))
        sizes[i] = 26 if s5 == 5 else 35
    n26 = sum(1 for v in sizes.values() if v == 26)
    return dict(P0=P0, tri=tri, point_units=point_units, line_y=line_y, plane_sizes=sizes, n26=n26, cnt=cnt)


def main():
    classes = json.load(open(os.path.join(CERTS, "planar_classes.json")))
    points4 = R.build_points()
    pidx4 = {p: i for i, p in enumerate(points4)}
    points3 = build_points3()
    pidx3 = {p: i for i, p in enumerate(points3)}
    lq = lines3(points3)
    assert len(lq) == 91 and all(len(L) == 10 for L in lq)
    cases = [("A", i, c) for i, c in enumerate([c for c in classes if c["k"] == 26])] + \
            [("C", i, c) for i, c in enumerate([c for c in classes if c["k"] == 35])]
    results = {}
    for case, i, c in cases:
        st = setup(points4, pidx4, points3, pidx3, lq, c["rep"], case)
        exp_n26 = 26 if case == "A" else 56
        assert st["n26"] == exp_n26, (case, i, st["n26"])
        assert len(c["rep"]) == (26 if case == "A" else 35)
        print("%s%d  |Aut Sigma|=%-3d |Sigma|=%d  triple %s  P0=%d  zero points %s  size-26 planes through P0: %d  line units: %d x '2-secant', %d x '5-secant'" % (
            case, i, c["aut_order"], len(c["rep"]), st["tri"], st["P0"],
            sorted(v for v in st["point_units"] if v != st["P0"]), st["n26"],
            sum(1 for v in st["line_y"].values() if v == 0), sum(1 for v in st["line_y"].values() if v == 1)))
        results[(case, i)] = st
    if "--compare" in sys.argv:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
        from trivial_tree import Setup
        for (case, i), st in results.items():
            S = Setup(case, i)
            main_pu = {p: v for p, v in S.pt_units.items()}
            ok_pu = main_pu == st["point_units"]
            ok_tri = tuple(S.tri) == tuple(st["tri"])
            # line units: main keyed by line index of P (pair-covering order); compare through the direction
            ok_y = all(S.line_y[S.line_of[d]] == st["line_y"][d] for d in range(91))
            lq_index = {tuple(L): j for j, L in enumerate(lq)}
            ok_sz = all(S.plane_size(L) == st["plane_sizes"][lq_index[tuple(S.P2.lines[L])]] for L in range(91))
            print("compare %s%d: point units %s, triple %s, line types %s, plane sizes %s" % (
                case, i, "OK" if ok_pu else "DIFF", "OK" if ok_tri else "DIFF", "OK" if ok_y else "DIFF", "OK" if ok_sz else "DIFF"))
            assert ok_pu and ok_tri and ok_y and ok_sz
        print("ALL CASES AGREE WITH THE MAIN CODE")


if __name__ == "__main__":
    main()

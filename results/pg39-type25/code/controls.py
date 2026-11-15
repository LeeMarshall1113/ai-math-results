"""Brute-force controls for the parameter derivation (params.py).

 1. trivial two-intersection sets of PG(3,q), q = 3,4,5 : point, plane,
    complement of a plane, complement of a point.  Line/plane statistics are
    counted on the explicit set and compared with the table rows.
 2. PG(2,9): Hermitian curve (unital, 28 points, type (1,4)) and a Baer
    subplane (13 points, type (1,4)) : compared with the planar table.
 3. PG(3,9): the Hermitian surface x1^4+x2^4+x3^4+x4^4 = 0 : 280 points,
    line counts {1,4,10}, plane counts {28,37}.
"""
from collections import Counter
from pg import PG
import params


def line_stats(P, K):
    S = set(K)
    return Counter(sum(1 for p in L if p in S) for L in P.lines)


def plane_stats(P, K):
    S = set(K)
    return Counter(sum(1 for p in pl if p in S) for pl in P.planes)


def row_for(d, q, alpha, beta, k):
    for r in params.survey(d, [q]):
        if (r["alpha"], r["beta"], r["k"]) == (alpha, beta, k):
            return r
    return None


def check_trivial(q):
    P = PG(4, q)
    N = P.N
    pt = [0]
    plane = list(P.planes[0])
    comp_plane = [p for p in range(N) if p not in set(plane)]
    comp_pt = [p for p in range(N) if p != 0]
    out = []
    for name, K in (("point", pt), ("plane", plane), ("complement of plane", comp_plane),
                    ("complement of point", comp_pt)):
        ls = line_stats(P, K)
        ps = plane_stats(P, K)
        vals = sorted(ls)
        assert len(vals) == 2, (name, ls)
        alpha, beta = vals
        row = row_for(3, q, alpha, beta, len(K))
        assert row is not None, (name, alpha, beta, len(K))
        assert (row["xa"], row["xb"]) == (ls[alpha], ls[beta]), (name, row, ls)
        assert dict(ps) == row["plane_dist"][0], (name, dict(ps), row["plane_dist"])
        out.append("q=%d %-20s k=%-4d type (%d,%d) lines %s planes %s : matches table row" %
                   (q, name, len(K), alpha, beta, dict(ls), dict(ps)))
    return out


def hermitian_points(P, q0):
    """points of x1^(q0+1)+...+xn^(q0+1) = 0 in PG(n-1, q0^2)"""
    F = P.F
    e = q0 + 1
    res = []
    for i, v in enumerate(P.points):
        s = 0
        for c in v:
            s = F.add[s][F.pw(c, e)]
        if s == 0:
            res.append(i)
    return res


def baer_subplane_points(P2):
    """points of PG(2,9) with all coordinates in the subfield GF(3) (normalised
    coordinates in {0,1,2} up to scalar: use the fixed points of the Frobenius)"""
    F = P2.F
    res = []
    for i, v in enumerate(P2.points):
        w = tuple(F.frob[c] for c in v)
        if P2.normalize(w) == v:
            res.append(i)
    return res


def main():
    for q in (3, 4, 5):
        for line in check_trivial(q):
            print(line)
    P2 = PG(3, 9)
    P2.verify()
    unital = hermitian_points(P2, 3)
    baer = baer_subplane_points(P2)
    for name, K in (("Hermitian curve (unital)", unital), ("Baer subplane", baer)):
        ls = line_stats(P2, K)
        vals = sorted(ls)
        assert len(vals) == 2
        row = None
        for r in params.survey(2, [9]):
            if (r["alpha"], r["beta"], r["k"]) == (vals[0], vals[1], len(K)):
                row = r
        assert row is not None and (row["xa"], row["xb"]) == (ls[vals[0]], ls[vals[1]]), (name, ls)
        print("PG(2,9) %-26s k=%d type (%d,%d) lines %s : matches planar table (b_P=%d d_P=%s)" %
              (name, len(K), vals[0], vals[1], dict(ls), row["bP"], row["dP"]))
    P3 = PG(4, 9)
    H = hermitian_points(P3, 3)
    ls, ps = line_stats(P3, H), plane_stats(P3, H)
    print("PG(3,9) Hermitian surface: |H| = %d, line counts %s, plane counts %s" % (len(H), dict(ls), dict(ps)))
    assert len(H) == 280 and set(ls) == {1, 4, 10} and set(ps) == {28, 37}
    assert ps[37] == 280 and ps[28] == 540
    # a (2,5) check of the arithmetic on a would-be example: none exists, print the target statistics
    print("target (2,5) sets in PG(3,9): see results/params_table.txt")


if __name__ == "__main__":
    main()

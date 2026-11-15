#!/usr/bin/env python3
"""Counting identity that closes the classification gap of PROOF-CHAIN 1.3 (stdlib; uses rebuild_all.planar_items for the frame data).

The frame enumerations (results/planar_k<m>_R<i>.json, certified complete by certs/planar_complete_*) list, for each of the two PGammaL(2,9)-orbit representatives R#i of 5-subsets of the line l0 = {z = 0},
one lex-leader solution per orbit of H_R = {diag(A,1) frob^e preserving R} on the sets K with K n l0 = R and origin in K.  Hence

   |Omega0(R#i)| = sum over solutions S of |H_R| / |Stab_{H_R}(S)|          (sets with K n l0 = R, origin in K)
   N(R#i)        = 81 |Omega0(R#i)| / (m - 5)                                (translations of AG(2,9) are regular on the affine points; double count of pairs (K, affine point of K))
   N_total(m)    = 91 (180 N(R#0) + 72 N(R#1)) / x5(m),   x5 = (10 m - 182)/3  (pairs (K, 5-secant); PGammaL(3,9) is transitive on lines, 252 = 180 + 72 five-subsets of l0)

N_total(m) is the number of (2,5)-sets of size m in PG(2,9), obtained WITHOUT nauty.  The class list of results/planar_classes.json (nauty) predicts sum_i |PGammaL(3,9)| / |Aut_i|.  The two numbers must be equal:
a missing class makes the second smaller, a wrong |Aut| changes it.  The |Aut_i| are recomputed without nauty, and the representatives are shown to be pairwise inequivalent, by planar_aut.py.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rebuild_all as RB   # noqa: E402

RES = os.path.join(HERE, "..", "results")


def main():
    pts, pidx, lines, thru, l0, origin, results = RB.planar_items()
    q = 9
    G = 2 * ((q ** 3 - 1) * (q ** 3 - q) * (q ** 3 - q * q) // (q - 1))          # |PGammaL(3,9)| = 2 |PGL(3,9)|
    classes = json.load(open(os.path.join(RES, "planar_classes.json")))
    ok = True
    for m in (26, 35):
        x5 = (10 * m - 182) // 3
        assert 3 * x5 == 10 * m - 182
        tot = 0
        for ri, osz, Rpts, H in results:
            d = json.load(open(os.path.join(RES, "planar_k%d_R%d.json" % (m, ri))))
            om0 = 0
            for sol in d["solutions"]:
                S = set(sol)
                stab = sum(1 for h in H if set(h[p] for p in S) == S)
                assert len(H) % stab == 0
                om0 += len(H) // stab
            assert (81 * om0) % (m - 5) == 0, (m, ri, om0)
            N = 81 * om0 // (m - 5)
            print("m=%d R#%d (orbit %d, |H_R| = %d): %d lex-leader solutions, |Omega0| = %d, N(R) = %d" % (m, ri, osz, len(H), len(d["solutions"]), om0, N))
            tot += osz * N
        assert (91 * tot) % x5 == 0
        n_total = 91 * tot // x5
        cl = [c for c in classes if c["k"] == m]
        n_classes = sum(G // c["aut_order"] for c in cl)
        print("m=%d: number of (2,5)-sets by the frame count (no nauty): %d ; by the class list sum |G|/|Aut_i| (|G| = %d, Aut = %s): %d ; %s" % (
            m, n_total, G, [c["aut_order"] for c in cl], n_classes, "EQUAL" if n_total == n_classes else "DIFFERENT"))
        ok &= (n_total == n_classes)
    print("COUNTING IDENTITY HOLDS FOR m = 26 AND m = 35" if ok else "COUNTING IDENTITY FAILS")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

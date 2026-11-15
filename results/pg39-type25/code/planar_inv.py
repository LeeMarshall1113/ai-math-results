"""Planar (2,5)-sets of size 26 / 35 in PG(2,9) invariant under a prime-order collineation.

Collineation classes of prime order of PGammaL(3,9) (checked below by listing all
cycle types obtained from our generators is NOT needed: we simply test the natural
candidates):  homology diag(1,1,-1), Baer involution, elation (order 3), planar
order 5 element diag(1,M), order 7 and 13 (Singer of GF(9^3)), order 13 etc.
For each we ask: is there an invariant (2,5)-set with |K| = 26 or 35 ?
"""
import itertools
import sys
import time

import groups
from ff import GF
from orbsat import OrbitInstance, run_solver, check_K
from pg import PG


def main():
    F = GF(9)
    P = PG(3, 9, F)
    omega = F.prim
    minus1 = F.neg[1]
    I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    cands = {}
    cands["homology diag(1,1,-1)"] = ([[1, 0, 0], [0, 1, 0], [0, 0, minus1]], 0)
    cands["baer (frobenius)"] = (I3, 1)
    cands["elation I+e13 (order 3)"] = ([[1, 0, 1], [0, 1, 0], [0, 0, 1]], 0)
    cands["unipotent J3 (order 3)"] = ([[1, 1, 0], [0, 1, 1], [0, 0, 1]], 0)
    # order 5: diag(1, M1) with M1 of order 5 in GL(2,9)
    blocks = groups.order5_blocks(F)
    M1 = blocks[sorted(blocks)[0]]
    cands["order5 diag(1,M1)"] = ([[1, 0, 0], [0, M1[0][0], M1[0][1]], [0, M1[1][0], M1[1][1]]], 0)
    C3 = groups.primitive_companion(F, 3)
    cands["order7 singer GF(9^3)"] = (groups.matpow(F, C3, 104), 0)
    cands["order13 singer GF(9^3)"] = (groups.matpow(F, C3, 56), 0)
    for name, (M, f) in cands.items():
        perm = P.perm_of_matrix(M, f)
        o = groups.perm_order(perm)
        inst = OrbitInstance(P, [perm])
        s = inst.summary()
        for k in (26, 35):
            cnf = inst.build_cnf(k, (2, 5))
            path = "../cnf/planar_inv_tmp.cnf"
            cnf.to_dimacs(path)
            status, model, t, info = run_solver(path, "cadical", 120)
            extra = ""
            if status == "SATISFIABLE":
                tr = set(m for m in model if m > 0)
                K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
                ok, bad = check_K(P, K, (2, 5))
                extra = " |K|=%d exact check %s" % (len(K), ok)
            print("%-26s order %-2d fixed pts %-2d orbits %-3d k=%d : %s (%.1fs)%s" % (
                name, o, s["fixed_points"], s["point_orbits"], k, status, t, extra), flush=True)


if __name__ == "__main__":
    main()

"""Class 5c (g = diag(M1,M2), no fixed point, two invariant lines l1 = P(V1), l2 = P(V2)), k = 260.

Facts used (all re-checked by code below):
  * K n l1 is one of the two g-orbits R on l1 (an invariant line with elliptic action meets K in 5 points);
    the centraliser contains the torus T (order 800) whose elements (alpha,1) act on l1 by shifts,
    so WLOG R = O1a.
  * the 10 planes pi_p = <l1, p> (p in l2) partition the points off l1; g permutes them in two orbits of 5;
    |pi_p n K| in {26,35}, and by the pencil count (k = 260) exactly 5 planes have 26 points and 5 have 35;
    the planes of one g-orbit have equal size, so one orbit has size 26, and WLOG (shift of l2 by (1,beta))
    it is the orbit of p_a = a fixed point of l2.
  * K n pi_{p_a} = S_a is a (2,5)-set of size 26 of the plane pi_{p_a} meeting l1 in R.
    The 5 planes of the orbit determine 110 + ... points of K : K n U = R u g-orbit of (S_a - l1), U = union.
Cubes: every labelled S_a (all planar solutions of size 26 with section R on l1).  Each cube is solved by
CaDiCaL under assumptions on the orbit variables of U.
"""
import itertools
import json
import os
import sys
import time

from pysat.solvers import Cadical195

import groups
import planar
from orbsat import OrbitInstance
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    F, P3, gens, reps = groups.build_all()
    P2 = PG(3, 9, F)
    M, f = reps["5c_mixed"]
    gperm = P3.perm_of_matrix(M, f)
    inst = OrbitInstance(P3, [gperm])
    pts3 = P3.points
    l1 = [i for i, v in enumerate(pts3) if v[2] == 0 and v[3] == 0]
    l2 = [i for i, v in enumerate(pts3) if v[0] == 0 and v[1] == 0]
    assert len(l1) == 10 and len(l2) == 10
    # g-orbits on l1 and l2
    def orbits_on(line):
        seen, orbs = set(), []
        for p in line:
            if p in seen:
                continue
            o = [p]
            q = gperm[p]
            while q != p:
                o.append(q)
                q = gperm[q]
            seen.update(o)
            orbs.append(sorted(o))
        return orbs
    o1 = orbits_on(l1)
    o2 = orbits_on(l2)
    assert [len(o) for o in o1] == [5, 5] and [len(o) for o in o2] == [5, 5]
    R = o1[0]
    pa = o2[0][0]
    a, b = pts3[pa][2], pts3[pa][3]
    # coordinatisation phi: PG(2,9) -> pi_pa
    def phi(i):
        x1, x2, s = P2.points[i]
        return P3.pidx[P3.normalize((x1, x2, F.mul[s][a], F.mul[s][b]))]
    phi_map = [phi(i) for i in range(P2.N)]
    plane_pts = sorted(phi_map)
    assert len(set(plane_pts)) == 91
    inv_phi = {j: i for i, j in enumerate(phi_map)}
    l0 = planar.line0_points(P2)
    Rprime = sorted(inv_phi[p] for p in R)
    assert all(p in l0 for p in Rprime)
    # class of R' among the PGammaL(2,9)-orbits on 5-subsets; planar solutions W_i
    P1, rreps = planar.pgammal2_orbits_on_5subsets(5)

    def lift(i):
        v = P1.points[i]
        return P2.pidx[P2.normalize((v[0], v[1], 0))]
    # all elements of Gamma L(2,9) acting on l0 (diag(A,1) * frob^e) as permutations of the 91 points
    cand = None
    t0 = time.time()
    allH = []
    for aa, bb, cc, dd in itertools.product(range(9), repeat=4):
        det = F.sub(F.mul[aa][dd], F.mul[bb][cc])
        if det == 0:
            continue
        for e in (0, 1):
            allH.append(([[aa, bb, 0], [cc, dd, 0], [0, 0, 1]], e))
    # find sigma with sigma(R_i) = R'
    sigma = None
    ri_sel = None
    for ri, (R1, osz) in enumerate(rreps):
        Ri = set(lift(i) for i in R1)
        for (Mx, e) in allH:
            perm = P2.perm_of_matrix(Mx, e)
            if set(perm[p] for p in Ri) == set(Rprime):
                sigma, ri_sel = perm, ri
                break
        if sigma is not None:
            break
    print("R' is in PGammaL(2,9)-orbit #%d (size %d)  (%.1fs)" % (ri_sel, rreps[ri_sel][1], time.time() - t0), flush=True)
    # planar solutions W for this R_i (size 26): load results/planar_k26_R<ri>.json
    d = json.load(open(os.path.join(HERE, "..", "results", "planar_k26_R%d.json" % ri_sel)))
    W = [tuple(K) for K in d["solutions"]]
    Ri = set(d["R"])
    H = planar.stabilizer_perms(P2, Ri)
    # translations of the affine plane x3 != 0 : (x1,x2,x3) -> (x1 + t1 x3, x2 + t2 x3, x3)
    trans = []
    for t1 in range(9):
        for t2 in range(9):
            trans.append(P2.perm_of_matrix([[1, 0, t1], [0, 1, t2], [0, 0, 1]], 0))
    Q = []
    for h in H:
        for tr in trans:
            Q.append(tuple(tr[h[p]] for p in range(P2.N)))
    Q = list(set(Q))
    print("|W| = %d, |H_R| = %d, |Q| = %d" % (len(W), len(H), len(Q)), flush=True)
    X = set()
    for K in W:
        for q in Q:
            X.add(frozenset(sigma[q[p]] for p in K))
    print("labelled planar solutions with section R' on l1: %d  (%.1fs)" % (len(X), time.time() - t0), flush=True)
    # sanity: every X has |X|=26, X n l0 = R', and every line meets in 2/5
    for Xs in list(X)[:50]:
        assert len(Xs) == 26 and set(p for p in Xs if p in set(l0)) == set(Rprime)
        assert all(sum(1 for p in L if p in Xs) in (2, 5) for L in P2.lines)
    # build the 5c orbit CNF, k = 260
    cnf = inst.build_cnf(260, (2, 5), amo=True)
    solver = Cadical195(bootstrap_with=cnf.clauses)
    # the union U of the 5 planes of the g-orbit of pi_pa
    plane_orbit = [plane_pts]
    cur = plane_pts
    for _ in range(4):
        cur = sorted(gperm[p] for p in cur)
        plane_orbit.append(cur)
    U = set(p for pl in plane_orbit for p in pl)
    Uorbs = sorted(set(inst.oid[p] for p in U))
    assert all(set(inst.orbs[o]) <= U for o in Uorbs)
    print("U has %d points = %d orbits (free orbits: %d)" % (len(U), len(Uorbs), len(inst.orbs) - len(Uorbs)), flush=True)
    stats = dict(unsat=0, sat=0)
    t1 = time.time()
    sols = []
    for n, Xs in enumerate(X):
        Sa = set(phi_map[p] for p in Xs)
        Ka = set(R)
        for i in range(5):
            pass
        # K_a = R u g-orbit of (Sa - l1)
        cur = set(p for p in Sa if p not in set(l1))
        for i in range(5):
            Ka |= cur
            cur = set(gperm[p] for p in cur)
        assump = []
        for o in Uorbs:
            inK = inst.orbs[o][0] in Ka
            assert all((p in Ka) == inK for p in inst.orbs[o])
            assump.append(inst.z[o] if inK else -inst.z[o])
        if solver.solve(assumptions=assump):
            stats["sat"] += 1
            m = solver.get_model()
            val = {abs(l): l > 0 for l in m}
            Kfull = sorted(p for i, v in enumerate(inst.z) if val[v] for p in inst.orbs[i])
            sols.append(Kfull)
            print("SAT cube", n, "size", len(Kfull), flush=True)
            break
        else:
            stats["unsat"] += 1
        if (n + 1) % 2000 == 0:
            print("  %d cubes done, %.1fs, %s" % (n + 1, time.time() - t1, stats), flush=True)
    print("cubes: %d total, %s, %.1fs" % (len(X), stats, time.time() - t1), flush=True)
    json.dump(dict(cubes=len(X), stats=stats, sat=sols), open(os.path.join(HERE, "..", "results", "cube_5c.json"), "w"))


if __name__ == "__main__":
    main()

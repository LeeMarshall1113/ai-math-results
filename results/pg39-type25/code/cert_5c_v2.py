"""Class 5c, k = 260: certificate of the cube cover.

Sections S_a of the plane pi_a = <l1,p_a> (|S_a| = 26, S_a n l1 = R) are produced by GROUP EXPANSION from the planar
classification: X = { sigma q K : K in W, q in Q } where W = lex-leader solutions of planar_k26_R<i>.json (certified complete by
planar_complete_k26_R<i>, see planar_cert.py), Q = translations x stabiliser H_R, sigma a collineation moving R_i to R'.
The list is cross-checked against an independent SAT enumeration (certs/5c_cubes.json, 9882 sections).

Output : certs/5c_cover_k260.{cnf,lrat} (orbit CNF of 5c with |K| = 260 and one selector per cube) and a manifest record.
"""
import itertools
import json
import os
import sys
import time

import groups
import planar
from certify import cake_check, sha256
from orbsat import OrbitInstance, run_solver
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def main():
    t0 = time.time()
    F, P3, gens, reps = groups.build_all()
    P2 = PG(3, 9, F)
    M, f = reps["5c_mixed"]
    gperm = P3.perm_of_matrix(M, f)
    inst = OrbitInstance(P3, [gperm])
    pts3 = P3.points
    l1 = [i for i, v in enumerate(pts3) if v[2] == 0 and v[3] == 0]
    l2 = [i for i, v in enumerate(pts3) if v[0] == 0 and v[1] == 0]

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
    o1, o2 = orbits_on(l1), orbits_on(l2)
    R = o1[0]
    pa = o2[0][0]
    a, b = pts3[pa][2], pts3[pa][3]
    phi_map = [P3.pidx[P3.normalize((P2.points[i][0], P2.points[i][1], F.mul[P2.points[i][2]][a], F.mul[P2.points[i][2]][b]))] for i in range(P2.N)]
    plane = sorted(phi_map)
    inv_phi = {j: i for i, j in enumerate(phi_map)}
    l0 = planar.line0_points(P2)
    Rprime = sorted(inv_phi[p] for p in R)
    P1, rreps = planar.pgammal2_orbits_on_5subsets(5)

    def lift(i):
        v = P1.points[i]
        return P2.pidx[P2.normalize((v[0], v[1], 0))]
    sigma, ri_sel = None, None
    for ri, (R1, osz) in enumerate(rreps):
        Ri = set(lift(i) for i in R1)
        for aa, bb, cc, dd in itertools.product(range(9), repeat=4):
            if F.sub(F.mul[aa][dd], F.mul[bb][cc]) == 0:
                continue
            for e in (0, 1):
                perm = P2.perm_of_matrix([[aa, bb, 0], [cc, dd, 0], [0, 0, 1]], e)
                if set(perm[p] for p in Ri) == set(Rprime):
                    sigma, ri_sel = perm, ri
                    break
            if sigma:
                break
        if sigma:
            break
    d = json.load(open(os.path.join(HERE, "..", "results", "planar_k26_R%d.json" % ri_sel)))
    W = [tuple(K) for K in d["solutions"]]
    Ri = set(d["R"])
    H = planar.stabilizer_perms(P2, Ri)
    trans = [P2.perm_of_matrix([[1, 0, t1], [0, 1, t2], [0, 0, 1]], 0) for t1 in range(9) for t2 in range(9)]
    Q = list(set(tuple(tr[h[p]] for p in range(P2.N)) for h in H for tr in trans))
    X = set()
    for K in W:
        for q in Q:
            X.add(frozenset(phi_map[sigma[q[p]]] for p in K))
    print("group expansion: R' in orbit #%d, |W|=%d, |Q|=%d, sections %d (%.0fs)" % (ri_sel, len(W), len(Q), len(X), time.time() - t0), flush=True)
    cubes = json.load(open(os.path.join(CERTS, "5c_cubes.json")))
    sat_set = set(frozenset(S) for S in cubes["sections"])
    print("independent SAT enumeration: %d sections ; sets equal: %s" % (len(sat_set), sat_set == X), flush=True)
    assert sat_set == X
    sols = sorted(tuple(sorted(S)) for S in X)
    # cover CNF
    base = inst.build_cnf(260, (2, 5), amo=True)
    plane_orbit = [plane]
    cur = plane
    for _ in range(4):
        cur = sorted(gperm[p] for p in cur)
        plane_orbit.append(cur)
    U = set(p for pl in plane_orbit for p in pl)
    Uorbs = sorted(set(inst.oid[p] for p in U))
    l1s = set(l1)
    sel = []
    for S in sols:
        Ka = set(R)
        cur = set(p for p in S if p not in l1s)
        for _ in range(5):
            Ka |= cur
            cur = set(gperm[p] for p in cur)
        c = base.new_var()
        sel.append(c)
        for o in Uorbs:
            base.add(-c, inst.z[o] if inst.orbs[o][0] in Ka else -inst.z[o])
    base.add(*sel)
    tag = "5c_cover_k260"
    cnfp, lratp = os.path.join(CERTS, tag + ".cnf"), os.path.join(CERTS, tag + ".lrat")
    base.to_dimacs(cnfp, comments=["5c cover k 260: %d cubes" % len(sols)])
    st, model, t, info = run_solver(cnfp, "cadical", 3000, lratp)
    rec = dict(cls="5c_mixed", k=260, tag=tag, status=st, solver_seconds=round(t, 1), vars=base.nvars, clauses=len(base.clauses),
               cnf_sha256=sha256(cnfp), cubes=len(sols), **info)
    if st == "UNSATISFIABLE":
        rec["lrat_bytes"], rec["lrat_sha256"] = os.path.getsize(lratp), sha256(lratp)
        v, ct, rc = cake_check(cnfp, lratp)
        rec["cake_lpr"], rec["cake_seconds"], rec["cake_rc"] = v, round(ct, 1), rc
    print(json.dumps(rec), flush=True)
    open(os.path.join(CERTS, "manifest.jsonl"), "a").write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()

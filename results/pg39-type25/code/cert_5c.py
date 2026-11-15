"""Certificate for class 5c (k = 260):   two LRAT-checked CNFs.

  (1) certs/5c_plane_complete.cnf : the planar problem in the plane pi_a = <l1, p_a>  (|S| = 26, every line of pi_a meets S in 2 or 5
      points, S n l1 = R) plus one blocking clause per section in the list  =>  UNSAT, i.e. the list of sections is complete.
  (2) certs/5c_cover_k260.cnf     : the orbit CNF of 5c (|K| = 260) plus, for every listed section S, a selector c_S that forces
      K n pi_a-orbit-union to be the g-orbit-image of S, plus the clause  OR c_S   =>  UNSAT.

WLOG facts used (NOTES.md section 4.3): K n l1 = R is a g-orbit; the planes <l1,p> (p in l2) are permuted by g in two orbits of 5;
exactly one orbit has section size 26 (pencil count for k = 260); by the centraliser torus we may take R = O1a and the 26-plane to be pi_{p_a}.
"""
import json
import os
import sys
import time

from pysat.solvers import Cadical195

import groups
from certify import cake_check, sha256
from cnfenc import CNF, sum_in_set, sum_equals
from orbsat import OrbitInstance, run_solver
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def geometry():
    F, P3, gens, reps = groups.build_all()
    M, f = reps["5c_mixed"]
    gperm = P3.perm_of_matrix(M, f)
    pts = P3.points
    l1 = [i for i, v in enumerate(pts) if v[2] == 0 and v[3] == 0]
    l2 = [i for i, v in enumerate(pts) if v[0] == 0 and v[1] == 0]

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
    a, b = pts[pa][2], pts[pa][3]
    # points of pi_a = <l1, p_a>: vectors (x1,x2, s a, s b)
    pl = set()
    for x1 in range(9):
        for x2 in range(9):
            for s in range(9):
                if x1 == x2 == s == 0:
                    continue
                pl.add(P3.pidx[P3.normalize((x1, x2, F.mul[s][a], F.mul[s][b]))])
    plane = sorted(pl)
    assert len(plane) == 91
    plane_lines = [L for L in P3.lines if all(p in pl for p in L)]
    assert len(plane_lines) == 91
    return F, P3, reps, gperm, l1, l2, R, pa, plane, plane_lines


def plane_cnf(plane, plane_lines, l1, R, k=26):
    """planar problem on the 91 points of pi_a (variables 1..91 in the order of `plane`)"""
    cnf = CNF()
    x = {p: cnf.new_var() for p in plane}
    y = []
    for L in plane_lines:
        fin = sum_in_set(cnf, [x[p] for p in L], [1] * 10, {2, 5}, amo=True)
        y.append(fin[5])
    sum_equals(cnf, [x[p] for p in plane], [1] * 91, k, amo=False)
    d = (k - 20) // 3
    through = {p: [j for j, L in enumerate(plane_lines) if p in L] for p in plane}
    for p in plane:
        thr = [y[j] for j in through[p]]
        fin = sum_in_set(cnf, thr, [1] * len(thr), {d, d + 3}, amo=False)
        cnf.add(-x[p], fin[d + 3])
        cnf.add(x[p], -fin[d + 3])
    Rs = set(R)
    for p in l1:
        cnf.add(x[p] if p in Rs else -x[p])
    return cnf, x


def main():
    t0 = time.time()
    F, P3, reps, gperm, l1, l2, R, pa, plane, plane_lines = geometry()
    cnf, x = plane_cnf(plane, plane_lines, l1, R)
    free = [p for p in plane if p not in set(l1)]
    s = Cadical195(bootstrap_with=cnf.clauses)
    sols = []
    while s.solve():
        m = s.get_model()
        val = {abs(l): l > 0 for l in m}
        S = tuple(p for p in plane if val[x[p]])
        sols.append(S)
        s.add_clause([(-x[p] if val[x[p]] else x[p]) for p in free])
        if len(sols) % 1000 == 0:
            print("  enumerated", len(sols), "%.0fs" % (time.time() - t0), flush=True)
    s.delete()
    print("sections S (|S|=26, S n l1 = R): %d  (%.0fs)" % (len(sols), time.time() - t0), flush=True)
    json.dump(dict(R=R, pa=pa, plane=plane, sections=sols), open(os.path.join(CERTS, "5c_cubes.json"), "w"))
    out = open(os.path.join(CERTS, "manifest.jsonl"), "a")
    # (1) completeness CNF
    cnf2, x2 = plane_cnf(plane, plane_lines, l1, R)
    for S in sols:
        Ss = set(S)
        cnf2.add(*[(-x2[p] if p in Ss else x2[p]) for p in free])
    tag = "5c_plane_complete"
    cnfp, lratp = os.path.join(CERTS, tag + ".cnf"), os.path.join(CERTS, tag + ".lrat")
    cnf2.to_dimacs(cnfp, comments=["5c plane completeness: %d sections blocked" % len(sols)])
    st, model, t, info = run_solver(cnfp, "cadical", 3000, lratp)
    rec = dict(cls="5c_plane", k=26, tag=tag, status=st, solver_seconds=round(t, 1), vars=cnf2.nvars, clauses=len(cnf2.clauses),
               cnf_sha256=sha256(cnfp), sections=len(sols), **info)
    if st == "UNSATISFIABLE":
        rec["lrat_bytes"], rec["lrat_sha256"] = os.path.getsize(lratp), sha256(lratp)
        rec["cake_lpr"], rec["cake_seconds"], rec["cake_rc"] = cake_check(cnfp, lratp)[0], 0, 0
    print(json.dumps(rec), flush=True)
    out.write(json.dumps(rec) + "\n")
    out.flush()
    # (2) cover CNF
    M, f = reps["5c_mixed"]
    inst = OrbitInstance(P3, [gperm])
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
            inK = inst.orbs[o][0] in Ka
            base.add(-c, inst.z[o] if inK else -inst.z[o])
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
    out.write(json.dumps(rec) + "\n")
    out.flush()


if __name__ == "__main__":
    main()

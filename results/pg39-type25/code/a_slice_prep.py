"""Preparation of the A slice (level-1 cubes of the star cases A0, A1, A2 = 52 cubes) and of the SAT control.

Runs under the lease (see run_a_slice.sh).  Writes
  cnf/full_A_c<ci>.cnf             base CNF of the case (k = 287, units for P0, the three height-0 points and the 91 line types through P0)
  certs/cubes_A<ci>.json           the cube list: one homothety-orbit representative per orbit of the sections of the first plane of size 26
  cnf/control_hermitian.cnf        SAT control: allowed counts {1,4,10}, |K| = 280, units = membership of H(3,9) on 3 planes through a point of H
  results/a_slice_jobs.json        job list for the WSL driver (a_slice_driver.py)
"""
import json
import os
import sys
import time

import cube_run
import full
from cnfenc import sum_equals, sum_in_set
from controls import hermitian_points
from ff import GF
from pg import PG
from trivial_tree import Setup, write_cube_cnf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def main():
    t0 = time.time()
    jobs = []
    # ---------------------------------------------------------------- SAT control (Hermitian surface through the same cube machinery)
    F = GF(9)
    P = PG(4, 9, F)
    H = set(hermitian_points(P, 3))
    base = os.path.join(ROOT, "cnf", "control_hermitian_base.cnf")
    cnf = full.CNFFile(base)
    x = [cnf.new_var() for _ in range(P.N)]
    for L in P.lines:
        sum_in_set(cnf, [x[p] for p in L], [1] * 10, {1, 4, 10}, amo=True)
    sum_equals(cnf, x, [1] * P.N, 280, amo=False)
    cnf.close(["control: allowed 1,4,10 k=280 (Hermitian surface)"])
    P0 = min(H)
    pts = set()
    for pl in [pl for pl in P.planes if P0 in pl][:3]:
        pts |= set(pl)
    lits = [((p + 1) if p in H else -(p + 1)) for p in sorted(pts)]
    ctrl = os.path.join(ROOT, "cnf", "control_hermitian.cnf")
    write_cube_cnf(base, lits, ctrl)
    os.remove(base)
    jobs.append(dict(id="control_hermitian", kind="control", cnf="control_hermitian.cnf", allowed=[1, 4, 10], k=280,
                     note="SAT expected: Hermitian surface from its membership on %d points of 3 planes through a point" % len(pts)))
    # ---------------------------------------------------------------- the 52 cubes
    for ci in range(3):
        S = Setup("A", ci)
        L, reps, nsec = cube_run.cubes_for(S)
        bpath = os.path.join(ROOT, "cnf", "full_A_c%d.cnf" % ci)
        if not os.path.exists(bpath):
            full.make_case("A", ci, bpath)
        cubes = []
        for i, (sec, orb) in enumerate(reps):
            lits = [((p + 1) if v else -(p + 1)) for p, v in {**S.pt_units, **sec}.items()]
            jobs.append(dict(id="A%d_c%02d" % (ci, i), kind="cube", base="full_A_c%d.cnf" % ci, units=lits, case="A", cls=ci, cube=i,
                             orbit_size=orb, plane=L, k=287))
            cubes.append(dict(index=i, orbit_size=orb, plane=L, section=[[p, v] for p, v in sorted(sec.items())]))
        json.dump(dict(case="A", cls=ci, aut=S.cl["aut_order"], plane=L, sections=nsec, triple=list(S.tri), base_units={str(p): v for p, v in S.pt_units.items()},
                       cubes=cubes), open(os.path.join(ROOT, "certs", "cubes_A%d.json" % ci), "w"))
        print("A%d: plane %d, %d sections, %d cubes (%.0fs)" % (ci, L, nsec, len(reps), time.time() - t0), flush=True)
    json.dump(jobs, open(os.path.join(ROOT, "results", "a_slice_jobs.json"), "w"))
    print("jobs:", len(jobs), "(1 control + %d cubes)" % (len(jobs) - 1), flush=True)


if __name__ == "__main__":
    main()

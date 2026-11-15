"""Preparation of the C slice (level-1 cubes of the star cases C0..C6, k = 260, 119 cubes) for the Modal run / a later local run.

Runs under a lease (python + PySAT, about 1-2 GB).  Writes
  cnf/full_C_c<ci>.cnf        base CNF (if missing)
  certs/cubes_C<ci>.json      the cube list (one homothety-orbit representative per orbit of the sections of the first plane of size 26)
  results/c_slice_jobs.json   job list (same format as results/a_slice_jobs.json)
"""
import json
import os
import time

import cube_run
import full
from trivial_tree import Setup

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def main():
    t0 = time.time()
    jobs = []
    for ci in range(7):
        S = Setup("C", ci)
        L, reps, nsec = cube_run.cubes_for(S)
        bpath = os.path.join(ROOT, "cnf", "full_C_c%d.cnf" % ci)
        if not os.path.exists(bpath):
            full.make_case("C", ci, bpath)
        cubes = []
        for i, (sec, orb) in enumerate(reps):
            lits = [((p + 1) if v else -(p + 1)) for p, v in {**S.pt_units, **sec}.items()]
            jobs.append(dict(id="C%d_c%02d" % (ci, i), kind="cube", base="full_C_c%d.cnf" % ci, units=lits, case="C", cls=ci, cube=i,
                             orbit_size=orb, plane=L, k=260))
            cubes.append(dict(index=i, orbit_size=orb, plane=L, section=[[p, v] for p, v in sorted(sec.items())]))
        json.dump(dict(case="C", cls=ci, aut=S.cl["aut_order"], plane=L, sections=nsec, triple=list(S.tri), base_units={str(p): v for p, v in S.pt_units.items()},
                       cubes=cubes), open(os.path.join(ROOT, "certs", "cubes_C%d.json" % ci), "w"))
        print("C%d: plane %d, %d sections, %d cubes (%.0fs)" % (ci, L, nsec, len(reps), time.time() - t0), flush=True)
    json.dump(jobs, open(os.path.join(ROOT, "results", "c_slice_jobs.json"), "w"))
    print("C jobs:", len(jobs), flush=True)


if __name__ == "__main__":
    main()

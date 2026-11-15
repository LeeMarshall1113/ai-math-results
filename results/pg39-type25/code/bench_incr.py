"""Incremental level-2 test: ONE CaDiCaL instance (PySAT) loaded with the cube CNF; every level-2 child (section of a second plane) is solved
under assumptions, so learnt clauses are reused across the children.  Compare with the monolithic run of the parent cube.
  python bench_incr.py CUBE.cnf CHILDREN.json [TLIM_TOTAL]
"""
import json
import sys
import time

from pysat.solvers import Cadical195


def main():
    cnf, ch, tlim = sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 3000
    kids = json.load(open(ch))["children"]
    s = Cadical195()
    t = time.process_time()
    with open(cnf) as fh:
        for line in fh:
            if line[0] in "cp":
                continue
            s.add_clause([int(x) for x in line.split()[:-1]])
    print("loaded %.0fs cpu" % (time.process_time() - t), flush=True)
    t0 = time.process_time()
    per = []
    for i, a in enumerate(kids):
        t1 = time.process_time()
        r = s.solve(assumptions=a)
        per.append(time.process_time() - t1)
        if r:
            print("SAT child", i, flush=True)
            break
        if time.process_time() - t0 > tlim:
            print("time limit after %d children" % (i + 1), flush=True)
            break
    tot = time.process_time() - t0
    print(json.dumps(dict(children=len(kids), done=len(per), total_cpu_s=round(tot, 1), mean_child_s=round(tot / max(1, len(per)), 2),
                          max_child_s=round(max(per), 1), first5=[round(x, 1) for x in per[:5]])), flush=True)


if __name__ == "__main__":
    main()

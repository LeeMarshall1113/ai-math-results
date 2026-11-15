"""Plumbing test of the trivial-group cube pipeline on a problem WITH a solution:
allowed line counts {1,4,10}, |K| = 280 (Hermitian surface H(3,9)); the cube fixes the membership of the 244 points of three
planes through P0 in H (as the star sampler does for (2,5)-sets).  Expected: SAT, and the model passes the exact line check."""
import os
import subprocess
import time

import full
from cnfenc import sum_in_set, sum_equals
from controls import hermitian_points
from ff import GF
from pg import PG
from trivial_tree import write_cube_cnf, run_win

HERE = os.path.dirname(os.path.abspath(__file__))
F = GF(9)
P = PG(4, 9, F)
H = set(hermitian_points(P, 3))
base = os.path.join(HERE, "..", "cnf", "plumbing_base.cnf")
cnf = full.CNFFile(base)
x = [cnf.new_var() for _ in range(P.N)]
for L in P.lines:
    sum_in_set(cnf, [x[p] for p in L], [1] * 10, {1, 4, 10}, amo=True)
sum_equals(cnf, x, [1] * P.N, 280, amo=False)
cnf.close(["plumbing test: allowed 1,4,10 k=280"])
P0 = min(H)
# three planes through P0 and the lines through P0: take the planes spanned by P0 and three lines of the quotient forming a triangle
P0v = P.points[P0]
planes = [pl for pl in P.planes if P0 in pl][:3]
pts = set()
for pl in planes:
    pts |= set(pl)
lits = [((p + 1) if p in H else -(p + 1)) for p in sorted(pts)]
out = os.path.join(HERE, "..", "cnf", "plumbing_cube.cnf")
write_cube_cnf(base, lits, out)
t = time.time()
r = subprocess.run([full.__dict__.get("CADICAL_WIN", r"C:\Users\user\src\build-cadical\cadical.exe"), "--congruence=false", out],
                   capture_output=True, text=True, timeout=900)
print("fixed %d points of H on %d planes; cadical rc=%d (10 = SAT) in %.1fs" % (len(pts), len(planes), r.returncode, time.time() - t))
if r.returncode == 10:
    model = []
    for line in r.stdout.splitlines():
        if line.startswith("v "):
            model += [int(v) for v in line[2:].split()]
    val = {abs(v): v > 0 for v in model if v != 0}
    K = [p for p in range(P.N) if val[p + 1]]
    cnt = {}
    for L in P.lines:
        c = sum(1 for p in L if p in set(K))
        cnt[c] = cnt.get(c, 0) + 1
    print("model |K| =", len(K), "line counts", cnt, "equals H:", set(K) == H)
    assert set(c for c in cnt) <= {1, 4, 10} and len(K) == 280
    print("PLUMBING OK")

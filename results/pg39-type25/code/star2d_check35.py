"""Independent frame for the k = 35 planar classification: star at P0 in K.  A point of K lies on 8 five-secants and 2 two-secants;
Sigma = the 2 two-secant directions (one PGammaL(2,9)-orbit of 2-subsets, 45 of them).  Heights of both Sigma-directions are normalised
to 0 by the elation group (order 81, simply transitive on the pair of heights).
Predicted number of normalised solutions = (sum_C |G|/|Aut C|) * 35 / (91 * 45 * 81)  (every point of K has the same pencil type)."""
import json, os
from fractions import Fraction
from pysat.solvers import Cadical195
from cnfenc import CNF, sum_in_set, sum_equals
from ff import GF
from pg import PG
HERE = os.path.dirname(os.path.abspath(__file__))
F = GF(9); P = PG(3, 9, F); P1 = PG(2, 9, F)
GSIZE = 84913920
classes = [c for c in json.load(open(os.path.join(HERE, "..", "results", "planar_classes.json"))) if c["k"] == 35]
total = sum(Fraction(GSIZE, c["aut_order"]) for c in classes)
pred = total * 35 / (91 * 45 * 81)
print("classes:", [c["aut_order"] for c in classes], "labelled sets:", total, "predicted star-normalised count:", pred, flush=True)
P0 = P.pidx[(0, 0, 1)]
line_dir = {}
for j in P.lines_through[P0]:
    q = next(p for p in P.lines[j] if p != P0)
    line_dir[j] = P1.pidx[P1.normalize(P.points[q][:2])]
sigma = {0, 1}                    # two directions
cnf = CNF()
x = [cnf.new_var() for _ in range(P.N)]
y = []
for L in P.lines:
    fin = sum_in_set(cnf, [x[p] for p in L], [1] * 10, {2, 5}, amo=True)
    y.append(fin[5])
sum_equals(cnf, x, [1] * P.N, 35, amo=False)
d = (35 - 20) // 3
for p in range(P.N):
    thr = [y[j] for j in P.lines_through[p]]
    fin = sum_in_set(cnf, thr, [1] * len(thr), {d, d + 3}, amo=False)
    cnf.add(-x[p], fin[d + 3])
    cnf.add(x[p], -fin[d + 3])
cnf.add(x[P0])
for j, d in line_dir.items():
    cnf.add(-y[j] if d in sigma else y[j])
for d in sigma:
    vec = P1.points[d]
    cnf.add(x[P.pidx[P.normalize((vec[0], vec[1], 0))]])
sol = Cadical195(bootstrap_with=cnf.clauses)
n = 0
while sol.solve():
    m = sol.get_model(); val = {abs(l): l > 0 for l in m}
    n += 1
    sol.add_clause([(-x[p] if val[x[p]] else x[p]) for p in range(P.N)])
    if n % 1000 == 0: print("  ", n, flush=True)
print("observed %d, predicted %s : %s" % (n, pred, "MATCH" if Fraction(n) == pred else "MISMATCH"), flush=True)

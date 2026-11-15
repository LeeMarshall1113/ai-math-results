"""The trivial-group instance: all 820 points free, 7462 lines, |K| = k, symmetry broken at a point.

Star formulation at P0 = (0,0,0,1): lines through P0 <-> points v of the quotient plane PG(2,9) = {(x1,x2,x3)};
the secant type of each line through P0 is a planar (2,5)-set Sigma of the quotient (size 26 or 35):

   case A (k = 287, P0 in K):     Sigma = the 26 lines through P0 meeting K in exactly 2 points;
   case B (k = 260, P0 not in K): Sigma = the 26 lines through P0 meeting K in exactly 5 points.

(The other combinations lead to planar sets of size 35.)  By PGammaL(3,9) acting on the quotient (induced by the
stabiliser of P0) we may take Sigma = a class representative (3 classes of size 26, results/planar_classes.json).
Elations with centre P0 act on heights and let us normalise three non-collinear directions:
   case A: the second K-point on the lines l_v1, l_v2, l_v3 (v_i in Sigma) is the point (v_i : 0);
   case B: (v_i : 0) in K for three non-collinear v_i NOT in Sigma (a 2-secant has two K-points, translate one to 0).

CNF is streamed to disk (no clause list in memory).
"""
import json
import os
import sys
import time

from cnfenc import sum_in_set, sum_equals
from pg import PG
from ff import GF

HERE = os.path.dirname(os.path.abspath(__file__))


class CNFFile:
    def __init__(self, path):
        self.path = path
        self.tmp = path + ".body"
        self.fh = open(self.tmp, "w", newline="\n")
        self.nvars = 0
        self.nclauses = 0

    def new_var(self):
        self.nvars += 1
        return self.nvars

    def add(self, *lits):
        self.fh.write(" ".join(map(str, lits)) + " 0\n" if lits else "0\n")
        self.nclauses += 1

    def close(self, comments=()):
        self.fh.close()
        with open(self.path, "w", newline="\n") as out:
            for c in comments:
                out.write("c %s\n" % c)
            out.write("p cnf %d %d\n" % (self.nvars, self.nclauses))
            with open(self.tmp, "r") as body:
                while True:
                    chunk = body.read(1 << 22)
                    if not chunk:
                        break
                    out.write(chunk)
        os.remove(self.tmp)


def build_full(path, P, k, units_lits=None, point_units=(), line_units=(), local_points=(), comments=()):
    """point_units: iterable of (point, 0/1);  line_units: iterable of (line index, 0/1) meaning |l n K| = 5 (1) or 2 (0).
    local_points: points at which the (implied) local count 'number of 5-secants = base + 30 x_p' is added.
    Returns (nvars, nclauses, xvars, yvars)."""
    cnf = CNFFile(path)
    x = [cnf.new_var() for _ in range(P.N)]
    y = []
    for L in P.lines:
        fin = sum_in_set(cnf, [x[p] for p in L], [1] * len(L), {2, 5}, amo=True)
        y.append(fin[5])
    sum_equals(cnf, x, [1] * P.N, k, amo=False)
    base = (k - 182) // 3
    for p in local_points:
        lits = [y[j] for j in P.lines_through[p]]
        fin = sum_in_set(cnf, lits, [1] * len(lits), {base, base + 30}, amo=False)
        cnf.add(-x[p], fin[base + 30])
        cnf.add(x[p], -fin[base + 30])
        cnf.add(x[p], fin[base])
    for p, v in point_units:
        cnf.add(x[p] if v else -x[p])
    for j, v in line_units:
        cnf.add(y[j] if v else -y[j])
    nv, nc = cnf.nvars, cnf.nclauses
    cnf.close(comments)
    return nv, nc, x, y


def star_data(P, F=None):
    """P0 and for every point of the quotient plane PG(2,9) (index in the n=3 geometry) the line through P0,
    plus the 'height-0' point (v:0)."""
    F = F or P.F
    P0 = P.pidx[(0, 0, 0, 1)]
    P2 = PG(3, 9, F)
    line_of = {}
    for j in P.lines_through[P0]:
        L = P.lines[j]
        q = next(p for p in L if p != P0)
        v = P.points[q]
        base = P2.pidx[P2.normalize(v[:3])]
        line_of[base] = j
    assert len(line_of) == 91
    zero_pt = {i: P.pidx[P.normalize((v[0], v[1], v[2], 0))] for i, v in enumerate(P2.points)}
    return P0, P2, line_of, zero_pt


def noncollinear_triple(P2, S):
    S = sorted(S)
    for a in S:
        for b in S:
            if b <= a:
                continue
            Lab = next(L for L in P2.lines if a in L and b in L)
            for c in S:
                if c > b and c not in Lab:
                    return (a, b, c)
    raise RuntimeError


def make_case(case, class_index, path, local=False):
    """case 'A' (k=287, P0 in K) or 'B' (k=260, P0 not in K); class_index into planar_classes.json (k=26 classes)"""
    F = GF(9)
    P = PG(4, 9, F)
    classes = [c for c in json.load(open(os.path.join(HERE, "..", "results", "planar_classes.json"))) if c["k"] == (35 if case == "C" else 26)]
    cl = classes[class_index]
    P0, P2, line_of, zero_pt = star_data(P, F)
    Sigma = set(cl["rep"])
    k = 287 if case == "A" else 260
    inK = case in ("A", "C")
    point_units = [(P0, 1 if inK else 0)]
    line_units = []
    for v, j in line_of.items():
        in_sigma = v in Sigma
        # A: Sigma = 2-secants (y = 0) ; B: Sigma = 5-secants (y = 1)
        if inK:
            line_units.append((j, 0 if in_sigma else 1))
        else:
            line_units.append((j, 1 if in_sigma else 0))
    pool = Sigma if inK else set(range(P2.N)) - Sigma
    tri = noncollinear_triple(P2, pool)
    for v in tri:
        point_units.append((zero_pt[v], 1))
    lp = [P0] if local else []
    nv, nc, x, y = build_full(path, P, k, point_units=point_units, line_units=line_units, local_points=lp,
                              comments=["trivial group case %s class %d (aut %d) k %d" % (case, class_index, cl["aut_order"], k)])
    return dict(case=case, cls=class_index, aut=cl["aut_order"], k=k, vars=nv, clauses=nc, path=path)


if __name__ == "__main__":
    case, ci = sys.argv[1], int(sys.argv[2])
    path = os.path.join(HERE, "..", "cnf", "full_%s_c%d.cnf" % (case, ci))
    t = time.time()
    info = make_case(case, ci, path)
    print(json.dumps(info), "%.1fs" % (time.time() - t))

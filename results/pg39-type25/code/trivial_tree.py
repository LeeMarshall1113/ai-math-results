"""Probing tree for the trivial-group instance (case A or B, planar class ci of the star at P0).

Level 0 : (case, class, normalised heights)  -- full.py
Level r : sections of r planes through P0 (quotient lines L_1.. L_r), each a planar (2,5)-set containing the fixed data.

  python trivial_tree.py CASE CLASS count [CAP]              : number of sections of each of the first three planes
  python trivial_tree.py CASE CLASS sample N R T [CAP]       : N random cubes of depth R (R planes fixed), CaDiCaL for T seconds each
"""
import json
import os
import random
import subprocess
import sys
import time

from pysat.solvers import Cadical195

import full
from cnfenc import CNF, sum_in_set, sum_equals
from ff import GF
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
CADICAL_WIN = r"C:\Users\user\src\build-cadical\cadical.exe"


class Setup:
    def __init__(self, case, ci):
        self.case, self.ci = case, ci
        F = GF(9)
        self.F = F
        self.P = P = PG(4, 9, F)
        classes = [c for c in json.load(open(os.path.join(HERE, "..", "results", "planar_classes.json"))) if c["k"] == (35 if case == "C" else 26)]
        self.cl = classes[ci]
        self.P0, self.P2, self.line_of, self.zero_pt = full.star_data(P, F)
        self.Sigma = set(self.cl["rep"])
        self.k = 287 if case == "A" else 260
        inK = case in ("A", "C")      # A: k=287, P0 in K, |Sigma|=26 ; C: k=260, P0 in K, |Sigma|=35 ; B: k=260, P0 not in K, |Sigma|=26
        pool = self.Sigma if inK else set(range(self.P2.N)) - self.Sigma
        self.tri = full.noncollinear_triple(self.P2, pool)
        # base point units (point -> 0/1) and line types
        self.pt_units = {self.P0: 1 if inK else 0}
        for v in self.tri:
            self.pt_units[self.zero_pt[v]] = 1
        self.line_y = {}
        for v, j in self.line_of.items():
            ins = v in self.Sigma
            self.line_y[j] = (0 if ins else 1) if inK else (1 if ins else 0)
        # projection of every point of P to the quotient
        self.proj = {}
        for q, v in enumerate(P.points):
            if q == self.P0:
                continue
            self.proj[q] = self.P2.pidx[self.P2.normalize(v[:3])]
        self.line_pts = {}
        for L2 in range(len(self.P2.lines)):
            Ls = set(self.P2.lines[L2])
            self.line_pts[L2] = sorted([self.P0] + [q for q, v in self.proj.items() if v in Ls])

    def quotient_line_through(self, a, b):
        return next(i for i, L in enumerate(self.P2.lines) if a in L and b in L)

    def plane_size(self, L2):
        s5 = len(self.Sigma & set(self.P2.lines[L2]))      # |Sigma n L| in {2,5}
        if self.case in ("A", "C"):      # Sigma = 2-secants through P0 ; P0 in K ; b' five-secants = 10 - s5
            return 26 if s5 == 5 else 35
        # case B : P0 not in K ; Sigma = 5-secants ; d' = s5 five-secants through P0 in the plane: d'=2 -> 26, d'=5 -> 35
        return 26 if s5 == 2 else 35

    def plane_cnf(self, L2, fixed):
        """planar CNF on the 91 points of the plane over L2; fixed: dict point->0/1 (global points)"""
        P = self.P
        pts = self.line_pts[L2]
        ps = set(pts)
        cnf = CNF()
        x = {p: cnf.new_var() for p in pts}
        lines = [j for j, L in enumerate(P.lines) if all(p in ps for p in L)]
        assert len(lines) == 91
        y = {}
        for j in lines:
            fin = sum_in_set(cnf, [x[p] for p in P.lines[j]], [1] * 10, {2, 5}, amo=True)
            y[j] = fin[5]
        m = self.plane_size(L2)
        sum_equals(cnf, [x[p] for p in pts], [1] * 91, m, amo=False)
        d = (m - 20) // 3
        for p in pts:
            thr = [y[j] for j in lines if p in P.lines[j]]
            fin = sum_in_set(cnf, thr, [1] * len(thr), {d, d + 3}, amo=False)
            cnf.add(-x[p], fin[d + 3])
            cnf.add(x[p], -fin[d + 3])
        for j in lines:
            if j in self.line_y:
                cnf.add(y[j] if self.line_y[j] else -y[j])
        for p, v in {**self.pt_units, **fixed}.items():
            if p in x:
                cnf.add(x[p] if v else -x[p])
        return cnf, x, pts

    def sections(self, L2, fixed, cap=None):
        cnf, x, pts = self.plane_cnf(L2, fixed)
        s = Cadical195(bootstrap_with=cnf.clauses)
        free = [p for p in pts if p not in fixed and p not in self.pt_units]
        sols = []
        while s.solve():
            m = s.get_model()
            val = {abs(l): l > 0 for l in m}
            sols.append({p: (1 if val[x[p]] else 0) for p in pts})
            s.add_clause([(-x[p] if val[x[p]] else x[p]) for p in free])
            if cap and len(sols) >= cap:
                break
        s.delete()
        return sols


def write_cube_cnf(base_path, units, out_path):
    """copy the base CNF and append unit clauses (signed literals); header clause count is patched"""
    with open(base_path, "r") as fh, open(out_path, "w", newline="\n") as out:
        line = fh.readline()
        while line.startswith("c"):
            out.write(line)
            line = fh.readline()
        parts = line.split()
        assert parts[0] == "p"
        out.write("p cnf %s %d\n" % (parts[2], int(parts[3]) + len(units)))
        while True:
            chunk = fh.read(1 << 22)
            if not chunk:
                break
            out.write(chunk)
        for u in units:
            out.write("%d 0\n" % u)


def run_win(cnf, timeout):
    t = time.time()
    try:
        r = subprocess.run([CADICAL_WIN, "--congruence=false", cnf], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "TIMEOUT", timeout, None
    st = {10: "SAT", 20: "UNSAT"}.get(r.returncode, "ERR%d" % r.returncode)
    conf = None
    for line in r.stdout.splitlines():
        if line.startswith("c conflicts:"):
            conf = int(line.split()[2])
    return st, time.time() - t, conf


def main():
    case, ci, mode = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    S = Setup(case, ci)
    v1, v2, v3 = S.tri
    Ls = [S.quotient_line_through(v1, v2), S.quotient_line_through(v1, v3), S.quotient_line_through(v2, v3)]
    Ls = [L for L in Ls if S.plane_size(L) == 26] + [L for L in Ls if S.plane_size(L) != 26]
    print("case %s class %d (aut %d): triple %s ; plane sizes %s" % (case, ci, S.cl["aut_order"], S.tri, [S.plane_size(L) for L in Ls]), flush=True)
    if mode == "count":
        cap = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
        t = time.time()
        for name, L in zip(("L1", "L2", "L3"), Ls):
            secs = S.sections(L, {}, cap=cap)
            print("first plane %s: %d sections%s (%.1fs)" % (name, len(secs), " (cap)" if len(secs) >= cap else "", time.time() - t), flush=True)
    elif mode == "count26":
        # first-plane section counts for the planes of size 26 among L1,L2,L3, and orbits under the homothety group diag(1,1,1,lam)
        cap = int(sys.argv[4]) if len(sys.argv) > 4 else 30000
        F = S.F
        homs = [S.P.perm_of_matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, lam]], 0) for lam in range(1, 9)]
        for name, L in zip(("L1", "L2", "L3"), Ls):
            if S.plane_size(L) != 26:
                print("plane %s has size %d (skipped)" % (name, S.plane_size(L)), flush=True)
                continue
            secs = S.sections(L, {}, cap=cap)
            key = lambda sec: tuple(sorted(p for p, v in sec.items() if v))
            seen, orbs = set(), 0
            for sec in secs:
                k0 = key(sec)
                if k0 in seen:
                    continue
                orbs += 1
                for h in homs:
                    seen.add(tuple(sorted(h[p] for p in k0)))
            print("plane %s (size 26): %d sections, %d orbits under the 8 homotheties%s" % (name, len(secs), orbs, " (cap)" if len(secs) >= cap else ""), flush=True)
    elif mode == "sample":
        N, R, T = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
        cap = int(sys.argv[7]) if len(sys.argv) > 7 else 20000
        base = os.path.join(HERE, "..", "cnf", "full_%s_c%d.cnf" % (case, ci))
        if not os.path.exists(base):
            info = full.make_case(case, ci, base)
            print("base CNF:", info["vars"], "vars", info["clauses"], "clauses", flush=True)
        rng = random.Random(1000 + ci + 100 * (int(sys.argv[8]) if len(sys.argv) > 8 else 0))
        for n in range(N):
            fixed, bs, dead = {}, [], False
            for r in range(R):
                secs = S.sections(Ls[r], fixed, cap=cap)
                bs.append(len(secs) if len(secs) < cap else ">=%d" % cap)
                if not secs:
                    dead = True
                    break
                fixed.update(rng.choice(secs))
            if dead:
                print("sample %d: branching %s -> dead end (no section)" % (n, bs), flush=True)
                continue
            lits = [((p + 1) if v else -(p + 1)) for p, v in {**S.pt_units, **fixed}.items()]
            out = os.path.join(HERE, "..", "cnf", "cube_%s_c%d_R%d_n%d_s%s.cnf" % (case, ci, R, n, sys.argv[8] if len(sys.argv) > 8 else "0"))
            write_cube_cnf(base, lits, out)
            st, secs_, conf = run_win(out, T)
            print("sample %d depth %d: branching %s ; fixed %d points ; %s in %.1fs conflicts %s" % (n, R, bs, len(fixed), st, secs_, conf), flush=True)
            os.remove(out)


if __name__ == "__main__":
    main()

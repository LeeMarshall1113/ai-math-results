"""Validation of the star normalisation (pencil type Sigma + elation normalisation of heights) in the PLANE, where the answer is
independently known from the classification:

 predicted number of star-normalised solutions for a given 5-subset s of the pencil at P0 (k = 26, P0 in K, Sigma = 2-secants):
     c(s) = sum over classes C of  [ (|PGammaL(3,9)|/|Aut C|) * n_C(tau) / (91 * |orbit_tau|) ] / 81
 where n_C(tau) = number of points p of C whose pencil of 2-secants has PGammaL(2,9)-orbit type tau, and 81 = |elation group|
 (it acts freely on the solutions with the given Sigma because Sigma has >= 2 directions).
 observed : SAT enumeration of all sets K with  P0 in K, 2-secant directions = s, heights of two chosen directions of s equal to 0.
"""
import itertools
import json
import os
import sys
from fractions import Fraction

from pysat.solvers import Cadical195

import planar
from cnfenc import CNF, sum_in_set, sum_equals
from ff import GF
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))
F = GF(9)
P = PG(3, 9, F)                    # PG(2,9)
P1, rreps = planar.pgammal2_orbits_on_5subsets(5)     # PG(1,9) and the two orbit reps (+ orbit sizes)
GSIZE = 84913920
P0 = P.pidx[(0, 0, 1)]
# pencil: lines through P0 <-> directions (x1:x2) ; quotient point index in P1
line_dir = {}
for j in P.lines_through[P0]:
    L = P.lines[j]
    q = next(p for p in L if p != P0)
    v = P.points[q]
    line_dir[j] = P1.pidx[P1.normalize(v[:2])]
dir_line = {d: j for j, d in line_dir.items()}
assert len(dir_line) == 10


def orbit_type_of_subset(sub):
    """index of the PGammaL(2,9)-orbit (in rreps) of a 5-subset of PG(1,9) given as set of P1 indices"""
    # brute force over the group generated like in planar.py
    omega = F.prim
    gens = [P1.perm_of_matrix([[1, 1], [0, 1]]), P1.perm_of_matrix([[omega, 0], [0, 1]]), P1.perm_of_matrix([[0, 1], [1, 0]]),
            P1.perm_of_matrix([[1, 0], [0, 1]], frob=1)]
    seen = {tuple(sorted(sub))}
    fr = [tuple(sorted(sub))]
    while fr:
        nx = []
        for t in fr:
            for g in gens:
                u = tuple(sorted(g[i] for i in t))
                if u not in seen:
                    seen.add(u)
                    nx.append(u)
        fr = nx
    for ri, (R1, osz) in enumerate(rreps):
        if tuple(sorted(R1)) in seen:
            return ri
    raise RuntimeError


# ---------------------------------------------------------------- prediction from the classification
classes = [c for c in json.load(open(os.path.join(HERE, "..", "results", "planar_classes.json"))) if c["k"] == 26]
def det3(M):
    a = F.mul
    s = F.add
    def sub(x, y): return F.sub(x, y)
    t1 = a[M[0][0]][sub(a[M[1][1]][M[2][2]], a[M[1][2]][M[2][1]])]
    t2 = a[M[0][1]][sub(a[M[1][0]][M[2][2]], a[M[1][2]][M[2][0]])]
    t3 = a[M[0][2]][sub(a[M[1][0]][M[2][1]], a[M[1][1]][M[2][0]])]
    return s[sub(t1, t2)][t3]


def apply_mat_dir(M, q):
    v = P.points[q]
    img = P.normalize(tuple(__import__("functools").reduce(lambda acc, x: F.add[acc][x], [F.mul[a][b] for a, b in zip(row, v)], 0) for row in M))
    return img


def direction_matrix(p):
    """deterministic: M with M*p = P0 built by column operations"""
    v = list(P.points[p])
    # find invertible M with M v = (0,0,1)^T : take any basis b1,b2 of the kernel-complement; M^-1 has columns b1,b2,v
    # brute force: choose b1,b2 among unit vectors such that [b1 b2 v] invertible
    for b1, b2 in itertools.permutations([(1, 0, 0), (0, 1, 0), (0, 0, 1)], 2):
        Minv = [[b1[0], b2[0], v[0]], [b1[1], b2[1], v[1]], [b1[2], b2[2], v[2]]]
        if det3(Minv) != 0:
            # invert 3x3 over GF(9) by adjugate
            return inv3(Minv)
    raise RuntimeError


def inv3(M):
    a, s, sub = F.mul, F.add, F.sub
    def cof(i, j):
        rows = [r for r in range(3) if r != i]
        cols = [c for c in range(3) if c != j]
        m = [[M[r][c] for c in cols] for r in rows]
        d = sub(a[m[0][0]][m[1][1]], a[m[0][1]][m[1][0]])
        return d if (i + j) % 2 == 0 else F.neg[d]
    d = det3(M)
    di = F.inv[d]
    return [[a[di][cof(j, i)] for j in range(3)] for i in range(3)]


def pencil_type_at(Ks, p):
    """type tau of the 2-secant pencil of the set Ks at its point p"""
    M = direction_matrix(p)
    # lines through p: for each other point q of the plane, direction of line pq after moving p to P0
    dirs_2sec = set()
    for L in P.lines:
        if p in L:
            c = sum(1 for x in L if x in Ks)
            if c == 2:
                q = next(x for x in L if x != p)
                img = apply_mat_dir(M, q)
                dirs_2sec.add(P1.pidx[P1.normalize(img[:2])])
    assert len(dirs_2sec) == 5, dirs_2sec
    return orbit_type_of_subset(dirs_2sec)


pred = {0: Fraction(0), 1: Fraction(0)}
for c in classes:
    Ks = set(c["rep"])
    ntau = {0: 0, 1: 0}
    for p in Ks:
        ntau[pencil_type_at(Ks, p)] += 1
    for tau in (0, 1):
        N_s = Fraction(GSIZE, c["aut_order"]) * ntau[tau] / (91 * rreps[tau][1])
        pred[tau] += N_s / 81
    print("class |Aut|=%d : points by pencil type %s" % (c["aut_order"], ntau), flush=True)
print("predicted star-normalised counts:", {t: str(v) for t, v in pred.items()}, flush=True)

# ---------------------------------------------------------------- observed
def count_star(tau):
    s = rreps[tau][0]            # a 5-subset (tuple of P1 indices) of directions
    sigma = set(s)
    cnf = CNF()
    x = [cnf.new_var() for _ in range(P.N)]
    y = []
    for L in P.lines:
        fin = sum_in_set(cnf, [x[p] for p in L], [1] * 10, {2, 5}, amo=True)
        y.append(fin[5])
    sum_equals(cnf, x, [1] * P.N, 26, amo=False)
    cnf.add(x[P0])
    for j, d in line_dir.items():
        cnf.add(-y[j] if d in sigma else y[j])
    v1, v2 = sorted(sigma)[:2]
    for d in (v1, v2):
        vec = P1.points[d]
        zero = P.pidx[P.normalize((vec[0], vec[1], 0))]       # the height-0 point (v:0)
        cnf.add(x[zero])
    sol = Cadical195(bootstrap_with=cnf.clauses)
    n = 0
    while sol.solve():
        m = sol.get_model()
        val = {abs(l): l > 0 for l in m}
        n += 1
        sol.add_clause([(-x[p] if val[x[p]] else x[p]) for p in range(P.N)])
    sol.delete()
    return n


for tau in (0, 1):
    n = count_star(tau)
    print("type %d: observed %d, predicted %s  %s" % (tau, n, pred[tau], "MATCH" if Fraction(n) == pred[tau] else "MISMATCH"), flush=True)

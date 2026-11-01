"""Exact check that Prymak's Lemma 2.4(ii) fails at a box corner (arXiv 2112.10698, Lemma 2.4).

Lemma 2.4(ii) says: for a box P and p in P, some translate of Q_P lies in int O_p.
Box P = prod [a_i, a_i + 1/160] with a = (0, 0, 63/160, 5/8, 0, 57/80); p = the corner a + (1,1,0,0,1,0)/160.

The script proves, in exact rational arithmetic:
 (1) for several facets {x : c0 + n.x >= 0} of O_p, some point of Q_P has slack exactly 0
     (the point is an exact vertex of Q_P: three tight rows of Q_P's H-rep, solved exactly, then checked
     against every row);
 (2) the inner normals n of these touched facets satisfy sum_i lam_i n_i = 0 with lam_i >= 0, not all 0.
Then for every translation u, sum_i lam_i (n_i . u) = 0, so n_i . u <= 0 for some touched facet i, and the
translate Q_P + u still has a point with slack <= 0 on facet i: it is not inside int O_p.  (Gordan's lemma.)

Uses hw.py (exact Q_P H-rep, Prymak's definition: intersection of O_v over the 64 vertices v of P).
scipy is used only to locate candidate vertices and multipliers; every claim is re-verified exactly.
Usage: python corner_example.py      (about 1 s)"""
import itertools, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import linprog
import hw


def solve3(rows):
    """exact solution x of c0 + n.x = 0 for three rows (c0, n1, n2, n3); None if singular."""
    M = [[F(r[1]), F(r[2]), F(r[3]), F(-r[0])] for r in rows]
    for c in range(3):
        piv = next((r for r in range(c, 3) if M[r][c] != 0), None)
        if piv is None: return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(3):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(4)]
    return [M[i][3] / M[i][i] for i in range(3)]


def slack(row, x): return row[0] + sum(row[k + 1] * x[k] for k in range(3))


def main():
    a = [F(0), F(0), F(63, 160), F(5, 8), F(0), F(57, 80)]; d = [F(1, 160)] * 6
    corner = (1, 1, 0, 0, 1, 0)
    p = [a[i] + corner[i] * d[i] for i in range(6)]
    H = hw.QP_hrep(a, d)
    V = [[0, p[0], p[1]], [1, p[0], p[1]], [p[2], 0, p[3]], [p[2], 1, p[3]], [p[4], p[5], 0], [p[4], p[5], 1]]
    O = sorted(hw.hrep_points([[F(x) for x in q] for q in V]))
    A = np.array([[-float(h[1]), -float(h[2]), -float(h[3])] for h in H]); b = np.array([float(h[0]) for h in H])
    touched = []
    for f in O:
        r = linprog(np.array([float(x) for x in f[1:]]), A_ub=A, b_ub=b, bounds=[(None, None)] * 3, method='highs')
        x0 = r.x
        tight = sorted(H, key=lambda h: abs(slack([float(t) for t in h], x0)))[:12]
        found = None
        for tri in itertools.combinations(tight, 3):
            x = solve3(tri)
            if x is None: continue
            if all(slack(h, x) >= 0 for h in H) and slack(f, x) == 0:
                found = x; break
        if found is not None:
            touched.append((f, found))
    print('p =', [str(x) for x in p])
    print('Q_P rows:', len(H), ' O_p facets:', len(O), ' facets of O_p touched by Q_P (exact):', len(touched))
    for f, x in touched:
        print('  facet', f, ' point of Q_P on it:', [str(t) for t in x])
    # Gordan multipliers: find lam >= 0, sum lam = 1, sum lam n = 0 (float), then make exact on a support of size <= 4
    N = [f[1:] for f, _ in touched]
    k = len(N)
    Aeq = np.array([[float(n[j]) for n in N] for j in range(3)] + [[1.0] * k]); beq = np.array([0, 0, 0, 1.0])
    r = linprog(np.zeros(k), A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * k, method='highs')
    ok = False
    if r.status == 0:
        supp = [i for i in range(k) if r.x[i] > 1e-12]
        for size in range(2, min(4, len(supp)) + 1):
            for S in itertools.combinations(supp, size):
                # exact null vector of the 3 x size matrix [n_i], i in S
                cols = [[F(N[i][j]) for j in range(3)] for i in S]
                # solve sum_{i in S} lam_i n_i = 0 with lam_{S[-1]} = 1
                m = size - 1
                if m == 0: continue
                rows = [[cols[c][j] for c in range(m)] + [-cols[-1][j]] for j in range(3)]
                # least squares not allowed: pick m independent rows and solve, then check all 3
                for R in itertools.combinations(range(3), m):
                    Msub = [rows[j][:] for j in R]
                    try:
                        # Gaussian elimination
                        Mx = [r_[:] for r_ in Msub]
                        for c in range(m):
                            piv = next(rr for rr in range(c, m) if Mx[rr][c] != 0)
                            Mx[c], Mx[piv] = Mx[piv], Mx[c]
                            for rr in range(m):
                                if rr != c and Mx[rr][c] != 0:
                                    fct = Mx[rr][c] / Mx[c][c]
                                    Mx[rr] = [Mx[rr][t] - fct * Mx[c][t] for t in range(m + 1)]
                        lam = [Mx[i][m] / Mx[i][i] for i in range(m)] + [F(1)]
                    except StopIteration:
                        continue
                    if all(l >= 0 for l in lam) and all(sum(lam[c] * cols[c][j] for c in range(size)) == 0 for j in range(3)):
                        print('exact Gordan certificate: sum lam_i n_i = 0 with')
                        for c, i in enumerate(S): print('   lam =', lam[c], ' n =', N[i])
                        ok = True; break
                if ok: break
            if ok: break
    print('RESULT: no translate of Q_P lies in int O_p (Lemma 2.4(ii) fails here)' if ok else 'RESULT: not established')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

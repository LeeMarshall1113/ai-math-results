"""Parameter table for two-intersection sets (type (alpha, beta)) with respect to
lines in PG(d, q), d = 2 or 3.  Exact integer arithmetic only.

For K with |K| = k and every line meeting K in alpha or beta points
(alpha < beta, both occurring):

  V = #points, r = #lines through a point, L = #lines, x_a + x_b = L
  alpha x_a + beta x_b = k r                       (count incident pairs)
  alpha(alpha-1) x_a + beta(beta-1) x_b = k(k-1)   (count ordered point pairs)
  => k^2 - ((alpha+beta-1) r + 1) k + alpha beta L = 0

Local counts: a point of K lies on b_P lines meeting K in beta points,
  (alpha-1)(r - b_P) + (beta-1) b_P = k - 1;
a point off K lies on d_P such lines, alpha (r - d_P) + beta d_P = k.

In PG(3,q) every plane section is a two-intersection set in PG(2,q) with the
same (alpha, beta); its size m must satisfy the planar local conditions, and
the plane sizes obey three moment equations (planes: V of them, each point on
r_pl = q^2+q+1 planes, each point pair on q+1 planes).
"""
import sys
from math import isqrt


def gauss_binom(n, k, q):
    num = den = 1
    for i in range(k):
        num *= q ** (n - i) - 1
        den *= q ** (i + 1) - 1
    assert num % den == 0
    return num // den


def geometry(d, q):
    V = (q ** (d + 1) - 1) // (q - 1)
    r = (q ** d - 1) // (q - 1)
    L = V * r // (q + 1)
    assert V * r % (q + 1) == 0 and L == gauss_binom(d + 1, 2, q)
    return V, r, L


def candidates(d, q, alpha, beta):
    """all integer k solving the line moment equations (with x_a, x_b >= 1)"""
    V, r, L = geometry(d, q)
    B = (alpha + beta - 1) * r + 1
    C = alpha * beta * L
    disc = B * B - 4 * C
    if disc < 0:
        return []
    s = isqrt(disc)
    if s * s != disc:
        return []
    out = []
    for sign in (+1, -1):
        num = B + sign * s
        if num % 2:
            continue
        k = num // 2
        if not (0 <= k <= V):
            continue
        xb_num = k * r - alpha * L
        if xb_num % (beta - alpha):
            continue
        xb = xb_num // (beta - alpha)
        xa = L - xb
        if xa < 1 or xb < 1:
            continue
        # second moment holds by construction of the quadratic; assert anyway
        assert alpha * (alpha - 1) * xa + beta * (beta - 1) * xb == k * (k - 1)
        out.append((k, xa, xb))
    return sorted(set(out))


def local_counts(d, q, alpha, beta, k):
    """(b_P or None, d_P or None) as integers; None if not an integer in [0, r]
    (and the corresponding point type exists)."""
    V, r, L = geometry(d, q)
    ok = True
    bP = dP = None
    if k >= 1:
        num = k - 1 - (alpha - 1) * r
        if num % (beta - alpha):
            ok = False
        else:
            bP = num // (beta - alpha)
            ok = ok and 0 <= bP <= r
    if k <= V - 1:
        num = k - alpha * r
        if num % (beta - alpha):
            ok = False
        else:
            dP = num // (beta - alpha)
            ok = ok and 0 <= dP <= r
    return ok, bP, dP


def plane_sizes_3d(q, alpha, beta, k):
    """feasible sizes m of a plane section of K in PG(3,q): the section is a
    two-intersection set of PG(2,q) (or the plane is empty / full - excluded
    below since then lines of the plane meet in 0 / q+1 which is not alpha,beta
    unless alpha = 0 / beta = q+1), so local planar conditions must hold."""
    V2, r2, L2 = geometry(2, q)
    sizes = []
    for m in range(0, V2 + 1):
        # planar counts: need to be consistent with some (alpha,beta)-set of PG(2,q)
        if m == 0:
            if alpha == 0:
                sizes.append(m)
            continue
        if m == V2:
            if beta == q + 1:
                sizes.append(m)
            continue
        c = [x for x in candidates(2, q, alpha, beta) if x[0] == m]
        if not c:
            continue
        ok, _, _ = local_counts(2, q, alpha, beta, m)
        if ok:
            sizes.append(m)
    return sizes


def plane_distributions(q, alpha, beta, k, sizes):
    """non-negative integer solutions x_m (m in sizes) of the three moment
    equations over the planes of PG(3,q)."""
    V = geometry(3, q)[0]           # number of planes = number of points
    r_pl = q * q + q + 1             # planes through a point
    eq0, eq1, eq2 = V, k * r_pl, k * (k - 1) * (q + 1)
    sols = []
    if len(sizes) > 3:
        return None
    # brute force over the first len-2 unknowns is overkill; sizes is tiny in practice
    import itertools
    n = len(sizes)
    if n == 0:
        return []
    if n == 1:
        m = sizes[0]
        return [{m: V}] if (m * V == eq1 and m * (m - 1) * V == eq2) else []
    if n == 2:
        m1, m2 = sizes
        # x1 + x2 = V ; m1 x1 + m2 x2 = eq1
        num = eq1 - m2 * V
        den = m1 - m2
        if num % den:
            return []
        x1 = num // den
        x2 = V - x1
        if x1 < 0 or x2 < 0:
            return []
        if m1 * (m1 - 1) * x1 + m2 * (m2 - 1) * x2 != eq2:
            return []
        return [{m1: x1, m2: x2}]
    # n == 3 : one free parameter
    m1, m2, m3 = sizes
    for x3 in range(V + 1):
        # solve for x1, x2
        rhs0, rhs1 = V - x3, eq1 - m3 * x3
        num = rhs1 - m2 * rhs0
        den = m1 - m2
        if num % den:
            continue
        x1 = num // den
        x2 = rhs0 - x1
        if x1 < 0 or x2 < 0:
            continue
        if sum(m * (m - 1) * x for m, x in ((m1, x1), (m2, x2), (m3, x3))) == eq2:
            sols.append({m1: x1, m2: x2, m3: x3})
    return sols


def survey(d, qs, nontrivial_only=False):
    rows = []
    for q in qs:
        for alpha in range(0, q + 1):
            for beta in range(alpha + 1, q + 2):
                for (k, xa, xb) in candidates(d, q, alpha, beta):
                    ok, bP, dP = local_counts(d, q, alpha, beta, k)
                    if not ok:
                        continue
                    row = dict(d=d, q=q, alpha=alpha, beta=beta, k=k, xa=xa, xb=xb, bP=bP, dP=dP)
                    if d == 3:
                        sizes = plane_sizes_3d(q, alpha, beta, k)
                        row["plane_sizes"] = sizes
                        row["plane_dist"] = plane_distributions(q, alpha, beta, k, sizes)
                        row["planes_ok"] = bool(row["plane_dist"])
                    rows.append(row)
    return rows


def fmt(row):
    s = "q=%-3d (%d,%d) k=%-5d x_a=%-5d x_b=%-5d b_P=%s d_P=%s" % (
        row["q"], row["alpha"], row["beta"], row["k"], row["xa"], row["xb"], row["bP"], row["dP"])
    if "plane_sizes" in row:
        s += "  planes %s -> %s" % (row["plane_sizes"], row["plane_dist"])
    return s


def is_trivial(row):
    q = row["q"]
    return not (2 <= row["alpha"] and row["beta"] <= q - 1)


if __name__ == "__main__":
    qs3 = [3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 64, 81]
    print("=== PG(3,q): all (alpha,beta) with integral line moments and integral local counts ===")
    rows = survey(3, qs3)
    for row in rows:
        if row["planes_ok"]:
            print(("TRIVIAL " if is_trivial(row) else "NONTRIV ") + fmt(row))
    print()
    print("=== PG(3,q): non-trivial range 2 <= alpha < beta <= q-1: survivors after ALL tests ===")
    for q in qs3:
        sv = [r for r in rows if r["q"] == q and not is_trivial(r) and r["planes_ok"]]
        lines_only = [r for r in rows if r["q"] == q and not is_trivial(r)]
        print("q=%-3d  after line moments + local counts: %d   after plane tests too: %d" %
              (q, len(lines_only), len(sv)))
        for r in sv:
            print("      " + fmt(r))
    print()
    print("=== PG(2,q): planar check (2,5)-sets in PG(2,9) ===")
    for row in survey(2, [9]):
        if row["alpha"] == 2 and row["beta"] == 5:
            print(fmt(row))

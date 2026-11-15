#!/usr/bin/env python3
"""Exact integer arithmetic of the level-0 steps of PROOF-CHAIN.md (sections 1.1, 1.2, 1.3 step 1, 2.1); stdlib, runs in a fraction of a second.

  python level0_arith.py
"""
q = 9
LINES = (q * q + 1) * (q * q + q + 1)               # 7462 lines of PG(3,9)
PTS_PER_LINE = q + 1
LINES_THROUGH_POINT = q * q + q + 1                  # 91


def isqrt_exact(n):
    r = int(n ** 0.5)
    while r * r < n:
        r += 1
    while r * r > n:
        r -= 1
    assert r * r == n, n
    return r


# 1.1  k^2 - ((a+b-1)(q^2+q+1)+1) k + a b (q^2+1)(q^2+q+1) = 0   for (a,b) = (2,5)
a, b = 2, 5
B = (a + b - 1) * LINES_THROUGH_POINT + 1
C = a * b * LINES
disc = B * B - 4 * C
assert (B, C, disc) == (547, 74620, 729)
ks = sorted(((B - isqrt_exact(disc)) // 2, (B + isqrt_exact(disc)) // 2))
assert ks == [260, 287]
for k in ks:
    # line counts: x2 + x5 = LINES, 2 x2 + 5 x5 = 91 k  and the second moment 2 x2 + 20 x5 = k(k-1)
    x5 = (LINES_THROUGH_POINT * k - 2 * LINES) // 3
    x2 = LINES - x5
    assert 3 * x5 == LINES_THROUGH_POINT * k - 2 * LINES
    assert 2 * x2 + 20 * x5 == k * (k - 1)
    # through a point of K: a two-secants, b five-secants: a + b = 91, a + 4 b = k - 1
    bb = (k - 1 - LINES_THROUGH_POINT) // 3
    assert 3 * bb == k - 1 - LINES_THROUGH_POINT and 91 - bb == {260: 35, 287: 26}[k] and bb == {260: 56, 287: 65}[k]
    print("k = %d: %d two-secants, %d five-secants; through a point of K: %d two-secants (|Sigma|), %d five-secants" % (k, x2, x5, 91 - bb, bb))

# 1.2  plane sections: m^2 - 61 m + 910 = 0  (91 lines of PG(2,9), 10 points per line)
Lp = q * q + q + 1
Bp = (a + b - 1) * (q + 1) + 1
Cp = a * b * Lp
discp = Bp * Bp - 4 * Cp
assert (Bp, Cp, discp) == (61, 910, 81)
ms = sorted(((Bp - isqrt_exact(discp)) // 2, (Bp + isqrt_exact(discp)) // 2))
assert ms == [26, 35]
for m in ms:
    x5 = (10 * m - 2 * Lp) // 3                       # five-secants of a planar (2,5)-set of size m
    assert 3 * x5 == 10 * m - 2 * Lp
    bP = (m - 11) // 3                                # five-secants through P0 (in K) inside the plane
    aP = 10 - bP                                      # two-secants through P0 inside the plane
    assert 3 * bP == m - 11 and aP + 4 * bP == m - 1
    bo = (m - 20) // 3                                # five-secants through a point of the plane outside the section
    assert 3 * bo == m - 20 and 2 * (10 - bo) + 5 * bo == m
    print("plane section m = %d: %d five-secants (planes of size 26 through P0 for |Sigma| = 26: x5 = 26, for |Sigma| = 35: x5 = 56), through P0: %d two-secants, %d five-secants; "
          "local counts: member %d, non-member %d five-secants" % (m, x5, aP, bP, bP, bo))
assert {26: 5, 35: 2} == {m: 10 - (m - 11) // 3 for m in ms}              # |Sigma n ell| = 5 iff the plane has 26 points
# number of planes of size 26 through P0 = number of lines of Q meeting Sigma in 5 points = x5 of the planar set Sigma
assert (10 * 26 - 2 * Lp) // 3 == 26 and (10 * 35 - 2 * Lp) // 3 == 56

# 1.3 step 1: PGammaL(2,9) orbits on 5-subsets: 180 + 72 = C(10,5); |PGammaL(3,9)|
from math import comb
assert comb(10, 5) == 252 == 180 + 72
G = 2 * ((q ** 3 - 1) * (q ** 3 - q) * (q ** 3 - q * q) // (q - 1))
assert G == 84913920
# number of translations of AG(2,9) = 81, affine points of a section: m - 5
assert q * q == 81 and 26 - 5 == 21 and 35 - 5 == 30
print("|PGammaL(3,9)| = %d ;  elations with centre P0: %d = 9^3 ; homotheties diag(1,1,1,lambda): %d" % (G, q ** 3, q - 1))
print("LEVEL-0 ARITHMETIC OK")

"""Referee: exact ML = 1/(k+1) checks by rational interval covering (no breakpoint lemma).

ML(v) <= d  iff the CLOSED sets {t : ||t v_i|| <= d} cover [0,1].
ML(v) >= d  iff the OPEN sets {t : ||t v_i|| < d} do NOT cover [0,1] (a point remains,
             at which every ||t v_i|| >= d).
Both are decided exactly by sorting rational intervals and sweeping.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce


def intervals(v, d):
    out = []
    for a in v:
        for j in range(0, a + 1):
            out.append((Fr(j, a) - d / a, Fr(j, a) + d / a))
    return sorted(out)


def closed_cover(v, d):
    reach = Fr(0)
    for lo, hi in intervals(v, d):
        if lo > reach:
            return False
        reach = max(reach, hi)
    return reach >= 1


def open_cover(v, d):
    """True iff the union of open intervals covers [0,1] (a gap point, even a single
    shared endpoint, means not covered)."""
    reach = Fr(-1)       # sup of union so far, union open: points < reach covered (if contiguous)
    started = False
    for lo, hi in intervals(v, d):
        if not started:
            if lo >= 0:
                return False
            started = True
            reach = hi
            continue
        if lo >= reach:      # point `reach` uncovered (open intervals)
            return reach > 1
        reach = max(reach, hi)
    return reach > 1


def ml_equals(v, d):
    return closed_cover(v, d) and not open_cover(v, d)


def witness(v, d):
    """An explicit time t with min ||t v_i|| >= d, checked exactly."""
    reach = None
    for lo, hi in intervals(v, d):
        if reach is not None and lo >= reach and 0 <= reach <= 1:
            t = reach
            break
        reach = hi if reach is None else max(reach, hi)
    nd = lambda x: min(x - (x.numerator // x.denominator), 1 - (x - (x.numerator // x.denominator)))
    assert min(nd(t * a) for a in v) >= d
    return t


if __name__ == '__main__':
    cases = [((1, 2, 3, 4, 5, 6), Fr(1, 7), True), ((3, 6, 9, 12, 15, 18), Fr(1, 7), True),
             ((1, 2, 3, 4, 5, 7), Fr(1, 7), False),
             ((1, 2, 3, 4, 5, 6, 7), Fr(1, 8), True), ((1, 2, 3, 4, 5, 7, 12), Fr(1, 8), True),
             ((1, 4, 5, 6, 7, 11, 13), Fr(1, 8), True), ((2, 8, 10, 12, 14, 22, 26), Fr(1, 8), True),
             ((1, 2, 3, 4, 5, 6, 8), Fr(1, 8), False), ((1, 3, 4, 5, 9), Fr(1, 6), True)]
    ok = True
    for v, d, expect in cases:
        r = ml_equals(v, d)
        w = witness(v, d) if not open_cover(v, d) else None
        print(v, 'ML ==', d, ':', r, ' closed-cover', closed_cover(v, d), ' witness t =', w)
        ok &= (r == expect)
    print('ML CHECKS', 'OK' if ok else 'FAILED')

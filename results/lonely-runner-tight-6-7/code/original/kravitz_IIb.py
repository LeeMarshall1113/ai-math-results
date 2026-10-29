"""Exact check of condition (II.b) of Kravitz's Theorem 10.2 for a tight set T (ML = 1/n,
n = |T| + 1), following Lemma 10.1 literally (exact Fractions)."""
import sys
from fractions import Fraction as Fr
from math import lcm


def dist(x):
    x = x - (x.numerator // x.denominator)
    return min(x, 1 - x)


def equality_times(T, L):
    cands = set()
    for a in T:
        for b in T:
            for d in (a + b, abs(a - b)):
                if d:
                    for m in range(d):
                        cands.add(Fr(m, d))
    return sorted(t for t in cands if min(dist(t * v) for v in T) == L)


def check(T):
    T = sorted(T); n = len(T) + 1; L = Fr(1, n)
    ts = equality_times(T, L)
    assert ts, 'not tight'
    D = lcm(*[t.denominator for t in ts])
    report = []
    ok = True
    for Q in range(-((D - 1) // 2), D // 2 + 1):
        if not all(dist(t * Q) < L for t in ts):
            continue
        best = None
        for t in ts:
            u = t * Q; u = u - round(u)          # in (-1/2, 1/2]
            if u == Fr(-1, 2):
                u = Fr(1, 2)
            frac = lambda x: x - (x.numerator // x.denominator)
            al = max([v for v in T if frac(t * v) == L], default=None)
            be = max([v for v in T if frac(t * v) == 1 - L], default=None)
            for kind, v in (('a', al), ('b', be)):
                if v is None:
                    continue
                val = v * (L - u) if kind == 'a' else v * (L + u)
                # tie-breaks: alpha -> larger u; beta -> smaller u; alpha vs beta -> larger u
                key = (val, -u if kind == 'a' else u, 0 if kind == 'a' else 1)
                if best is None or key < best[0]:
                    best = (key, kind, v, u, t)
        _, kind, v, u, t = best
        if kind == 'a':
            c1 = (Q == -n * v * u); den = n * v * (1 - n * u)
        else:
            c1 = (Q == n * v * u); den = n * v * (1 + n * u)
        c2 = (Fr(D) / den).denominator == 1
        report.append(dict(Q=Q, kappa=kind, v=v, u=str(u), t=str(t), eqQ=c1, integer=c2))
        ok &= c1 and c2
    return ok, D, [str(t) for t in ts], report


if __name__ == '__main__':
    for T in [(1, 2, 3), (1, 2, 3, 4), (1, 3, 4, 7), (1, 2, 3, 4, 5), (1, 3, 4, 5, 9),
              (1, 2, 3, 4, 5, 6), (1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]:
        ok, D, ts, rep = check(T)
        print(T, 'II.b holds:', ok, ' D =', D, ' equality times', ts, ' admissible Q:', [r['Q'] for r in rep])
        for r in rep:
            if not (r['eqQ'] and r['integer']):
                print('   FAIL', r)

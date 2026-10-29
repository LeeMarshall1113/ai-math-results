"""Kravitz (Comb. Theory 2021) Conjecture 8.3, 'one very fast runner' form of the
Loneliness Spectrum Conjecture: for v_1 < ... < v_n with v_n huge,
ML = s/(ns+1) for some s, or ML >= 1/n.
By his Theorem 10.3 it reduces to the tight sets of n-1 speeds.
This script tests each tight (n-1)-set T numerically: ML(T + {x}) for every x in a window
[X0, X0 + W), exact rational arithmetic (numba integer kernel), and reports any value
below 1/n that is not of the form s/(ns+1).
"""
import os, sys, math
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numba import njit


@njit(cache=True)
def ml_num_den(v):
    """Exact ML of integer tuple v: returns (num, den) with ML = num/den (Cordella Lemma 2.2)."""
    k = v.shape[0]
    bn, bd = 0, 1
    for i in range(k):
        for j in range(i, k):
            for sgn in (1, -1):
                L = v[i] + sgn * v[j]
                if L < 0:
                    L = -L
                if L == 0:
                    continue
                for m in range(L):
                    mn = L
                    for a in range(k):
                        r = (m * v[a]) % L
                        if L - r < r:
                            r = L - r
                        if r < mn:
                            mn = r
                            if mn * bd <= bn * L:
                                break
                    if mn * bd > bn * L:
                        bn, bd = mn, L
    return bn, bd


def scan(T, X0, W):
    n = len(T) + 1
    bad, forms = [], {}
    for x in range(X0, X0 + W):
        v = np.array(list(T) + [x], dtype=np.int64)
        a, b = ml_num_den(v)
        f = Fraction(int(a), int(b))
        if f >= Fraction(1, n):
            forms['>=1/n'] = forms.get('>=1/n', 0) + 1
            continue
        # f = s/(ns+1) ?
        p, q = f.numerator, f.denominator
        k = q - n * p
        forms[f'k={k}'] = forms.get(f'k={k}', 0) + 1
        if k != 1:
            bad.append((x, str(f), k))
    return forms, bad


if __name__ == '__main__':
    tests = [
        ((1, 2, 3, 4, 5, 6), 2000, 400),          # n = 7 (Prop 10.4 predicts k = 1)
        ((1, 2, 3, 4, 5, 6, 7), 2000, 300),       # n = 8
        ((1, 2, 3, 4, 5, 7, 12), 2000, 300),      # n = 8, Wills
        ((1, 4, 5, 6, 7, 11, 13), 2000, 300),     # n = 8, Goddyn-Wong
    ]
    for T, X0, W in tests:
        forms, bad = scan(T, X0, W)
        print(T, f'x in [{X0},{X0+W})', forms, 'non-(k=1) below 1/n:', bad[:8], flush=True)

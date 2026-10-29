"""Bernstein nonnegativity test on [0,1]^6, integers only, with de Casteljau subdivision.

bern_nonneg(P) returns True only if it has exhibited a finite partition of [0,1]^6 into
dyadic boxes on each of which every Bernstein coefficient of P is >= 0 (which implies
P >= 0 there, since the Bernstein basis is nonnegative and sums to 1).

Method (deliberately different from the other two checkers, which re-expand the composed
polynomial from scratch after every bisection):
  * Bernstein coefficients are computed once, axis by axis, on a dense integer tensor:
    b_i = sum_{j<=i} C(i,j)/C(d,j) a_j, with every a_j pre-multiplied by
    lcm_j C(d,j) / C(d,j) so everything stays integral (a uniform positive scale).
  * A bisection at t_l = 1/2 is done directly on the Bernstein tensor with the
    de Casteljau triangle along axis l (unnormalised sums, uniform scale 2^d).
  * A box whose corner coefficient is negative (= sign of P at that corner) fails at once.
"""
from math import comb, lcm
from operator import sub as _sub
from poly import degs, W, FM, NV

STATS = {'calls': 0, 'boxes': 0}


class Budget(Exception):
    pass


def _tensor(P):
    d = degs(P)
    act = [l for l in range(NV) if d[l] > 0]
    dims = [d[l] + 1 for l in act]
    strides = [1] * len(act)
    for i in range(len(act) - 2, -1, -1):
        strides[i] = strides[i + 1] * dims[i + 1]
    N = 1
    for x in dims:
        N *= x
    wts = []
    for l in act:
        dl = d[l]
        L = lcm(*[comb(dl, j) for j in range(dl + 1)])
        wts.append([L // comb(dl, j) for j in range(dl + 1)])
    arr = [0] * N
    for k, c in P.items():
        idx = 0
        v = c
        for a, l in enumerate(act):
            e = (k >> (W * l)) & FM
            idx += e * strides[a]
            v *= wts[a][e]
        arr[idx] += v
    # power -> Bernstein, axis by axis: b_i = sum_j C(i,j) a'_j  (Pascal passes)
    for a in range(len(act)):
        s = strides[a]
        dl = dims[a] - 1
        S = s * (dl + 1)
        for o in range(0, N, S):
            for r in range(s):
                base = o + r
                f = arr[base:base + S:s]
                for kk in range(dl):
                    for i in range(dl, kk, -1):
                        f[i] += f[i - 1]
                arr[base:base + S:s] = f
    return arr, dims, strides, N


def _corners(dims, strides):
    idx = [0]
    for dm, s in zip(dims, strides):
        idx = [x + y for x in idx for y in (0, (dm - 1) * s)]
    return idx


def _split(arr, dims, strides, N, a):
    s = strides[a]
    dl = dims[a] - 1
    S = s * (dl + 1)
    Lf = [0] * N
    Rt = [0] * N
    for o in range(0, N, S):
        for r in range(s):
            base = o + r
            row = arr[base:base + S:s]
            left = [0] * (dl + 1)
            right = [0] * (dl + 1)
            left[0] = row[0] << dl
            right[dl] = row[dl] << dl
            for m in range(1, dl + 1):
                row = [row[i] + row[i + 1] for i in range(dl - m + 1)]
                left[m] = row[0] << (dl - m)
                right[dl - m] = row[-1] << (dl - m)
            Lf[base:base + S:s] = left
            Rt[base:base + S:s] = right
    return Lf, Rt


def _spread(arr, dims, strides, N, a):
    s = strides[a]
    S = s * dims[a]
    best = 0
    for o in range(0, N, S):
        blk = arr[o:o + S]
        m = max(map(abs, map(_sub, blk[s:], blk[:-s])))
        if m > best:
            best = m
    return best


def bern_nonneg(P, maxdepth=22, budget=40000):
    STATS['calls'] += 1
    if not P:
        return True
    arr, dims, strides, N = _tensor(P)
    if not dims:
        return arr[0] >= 0
    corners = _corners(dims, strides)
    cnt = [0]

    def rec(arr, depth):
        cnt[0] += 1
        STATS['boxes'] += 1
        if cnt[0] > budget:
            raise Budget
        if min(arr) >= 0:
            return True
        if any(arr[c] < 0 for c in corners):
            return False
        if depth >= maxdepth:
            return False
        a = max(range(len(dims)), key=lambda a: _spread(arr, dims, strides, N, a))
        L, R = _split(arr, dims, strides, N, a)
        del arr
        return rec(L, depth + 1) and rec(R, depth + 1)
    try:
        return rec(arr, 0)
    except Budget:
        return False

"""Second implementation (referee), the closed lift to modulus N = (k+1) p.

For a level-one class tuple c (mod p) and a digit vector d (mod g = k+1), x_i is the CRT
lift (x_i = c_i mod p, x_i = d_i mod g).  x is CLOSED-IMPROPER iff every t in 1..N-1 has
some i with g * dist(t x_i mod N, 0) <= N, i.e. ||t x_i / N|| <= 1/g.  (ALL times are
tested, including multiples of g and of p -- no reliance on level one.)

Algorithm (different from tight6/tight7, which loop over digit vectors and times):
for each coordinate i and digit d the set of covered times is a bitset C[i][d];
x is improper iff OR_i C[i][d_i] = all times.  Meet in the middle: split the
coordinates into A and B, tabulate the ORs for every partial digit vector, and for
each a find the b with C_B[b] containing the complement of C_A[a].

Lemma G filter (Lemma G in PROOF.md): a lift whose digit vector has at most one
coordinate not divisible by the Lemma-G modulus (7 for k = 6; 2 for k = 7) is ignored.
"""
import itertools
import numpy as np


def crt(c, d, p, g):
    # x = c mod p, x = d mod g, 0 <= x < g p
    return (c + p * (((d - c) * pow(p, -1, g)) % g)) % (g * p)


def covered_bits(x, N, g):
    t = np.arange(1, N, dtype=np.int64)
    r = (t * x) % N
    dist = np.minimum(r, N - r)
    b = (g * dist) <= N
    nbits = ((N - 1 + 63) // 64) * 64
    pad = np.zeros(nbits, dtype=np.bool_)
    pad[:N - 1] = b
    return np.packbits(pad, bitorder='little').view(np.uint64)


def all_bits(N):
    nbits = ((N - 1 + 63) // 64) * 64
    pad = np.zeros(nbits, dtype=np.bool_)
    pad[:N - 1] = True
    return np.packbits(pad, bitorder='little').view(np.uint64)


def lemma_g_mod(k):
    return {6: 7, 7: 2}[k]


def improper_lifts(cls, p, k, want_all=False, limit=10**9):
    """Return list of digit vectors d (tuples mod g) whose lift is closed-improper and
    NOT excluded by Lemma G (stops after `limit` hits unless want_all)."""
    g = k + 1
    N = g * p
    C = [[covered_bits(crt(c, d, p, g), N, g) for d in range(g)] for c in cls]
    ALL = all_bits(N)
    a = k // 2
    A, B = list(range(a)), list(range(a, k))

    def table(idx):
        digs = list(itertools.product(range(g), repeat=len(idx)))
        acc = np.zeros((len(digs), ALL.shape[0]), dtype=np.uint64)
        for j, dv in enumerate(digs):
            v = np.zeros(ALL.shape[0], dtype=np.uint64)
            for i, d in zip(idx, dv):
                v |= C[i][d]
            acc[j] = v
        return digs, acc

    dA, accA = table(A)
    dB, accB = table(B)
    dB_arr = np.array(dB, dtype=np.int64)
    LG = lemma_g_mod(k)
    hits = []
    for ia, da in enumerate(dA):
        need = ALL & ~accA[ia]
        ok = ((accB & need) == need).all(axis=1)
        if not ok.any():
            continue
        nzA = sum(1 for d in da if d % LG)
        nzB = (dB_arr[ok] % LG != 0).sum(axis=1)
        good = np.nonzero(nzA + nzB >= 2)[0]
        for j in good:
            hits.append(tuple(da) + tuple(int(z) for z in dB_arr[ok][j]))
            if not want_all and len(hits) >= limit:
                return hits
    return hits


def naive_improper_count(cls, p, k, chunk=200_000):
    """Direct check for validation: every digit vector, every time 1..N-1 (vectorised over
    digit vectors, loop over times).  Counts non-Lemma-G closed-improper lifts."""
    g = k + 1
    N = g * p
    LG = lemma_g_mod(k)
    total = g ** k
    inv = pow(p, -1, g)
    c = np.array(cls, dtype=np.int64)
    n = 0
    for lo in range(0, total, chunk):
        codes = np.arange(lo, min(total, lo + chunk), dtype=np.int64)
        D = np.zeros((len(codes), k), dtype=np.int64)
        cc = codes.copy()
        for i in range(k):
            D[:, i] = cc % g
            cc //= g
        X = (c[None, :] + p * (((D - c[None, :]) * inv) % g)) % N
        imp = np.ones(len(codes), dtype=np.bool_)
        for t in range(1, N):
            r = (t * X) % N
            dist = np.minimum(r, N - r)
            imp &= ((g * dist) <= N).any(axis=1)
        keep = (D % LG != 0).sum(axis=1) >= 2
        n += int((imp & keep).sum())
    return n

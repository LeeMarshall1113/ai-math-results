"""Exact sparse integer polynomials in up to 6 variables, Python standard library only.

A polynomial is a dict {key: int} with no zero values.  A monomial t0^e0 ... t5^e5 is
packed into one integer key = sum e_l << (W*l), W = 10 bits per exponent.  Monomial
multiplication is key addition.  Exponents are kept below 2^(W-1) = 512 (checked), so the
top bit of every field is a guard bit and fields never carry into each other.  The integer
order of keys is then a monomial order (lex, t5 most significant), which the exact division
and square-root routines rely on.

Many routines below work "up to a positive constant": the proof only ever needs the sign
of a polynomial on a box, so P may be replaced by q*P for any rational q > 0.  Those
routines say so; they return primitive integer polynomials.
"""
from math import comb, gcd, isqrt
import heapq

W = 10
NV = 6
FM = (1 << W) - 1
MAXE = (1 << (W - 1)) - 1
GUARD = sum(1 << (W * l + W - 1) for l in range(NV))
ONE = {0: 1}


def mono(exps):
    k = 0
    for l, e in enumerate(exps):
        assert 0 <= e <= MAXE
        k += e << (W * l)
    return k


def ex(k, l):
    return (k >> (W * l)) & FM


def exps(k):
    return [(k >> (W * l)) & FM for l in range(NV)]


def var(l):
    return {1 << (W * l): 1}


def const(c):
    return {0: c} if c else {}


def degs(P):
    d = [0] * NV
    for k in P:
        for l in range(NV):
            e = (k >> (W * l)) & FM
            if e > d[l]:
                d[l] = e
    return d


def _chk(d):
    assert max(d) <= MAXE, ('exponent overflow', d)


def add(P, Q):
    R = dict(P)
    for k, c in Q.items():
        v = R.get(k, 0) + c
        if v:
            R[k] = v
        else:
            R.pop(k, None)
    return R


def scale(P, c):
    return {k: c * v for k, v in P.items()} if c else {}


def neg(P):
    return {k: -v for k, v in P.items()}


def sub(P, Q):
    return add(P, neg(Q))


def mul(P, Q):
    if not P or not Q:
        return {}
    dP, dQ = degs(P), degs(Q)
    _chk([a + b for a, b in zip(dP, dQ)])
    R = {}
    get = R.get
    for k1, c1 in P.items():
        for k2, c2 in Q.items():
            k = k1 + k2
            R[k] = get(k, 0) + c1 * c2
    return {k: v for k, v in R.items() if v}


def pw(P, n):
    R = ONE
    for _ in range(n):
        R = mul(R, P)
    return R


def deriv(P, l):
    s = W * l
    one = 1 << s
    R = {}
    for k, c in P.items():
        e = (k >> s) & FM
        if e:
            R[k - one] = c * e
    return R


def subs01(P, l, v):
    """P with t_l = v, v in {0, 1}."""
    assert v in (0, 1)
    s = W * l
    R = {}
    for k, c in P.items():
        e = (k >> s) & FM
        if v == 0:
            if e == 0:
                R[k] = R.get(k, 0) + c
        else:
            k2 = k - (e << s)
            R[k2] = R.get(k2, 0) + c
    return {k: v for k, v in R.items() if v}


def prim(P):
    """P divided by the positive gcd of its coefficients (a positive multiple)."""
    if not P:
        return P
    g = gcd(*P.values())
    if g == 1:
        return P
    return {k: v // g for k, v in P.items()}


def affine(P, l, a, b, g):
    """g^d * P(t_l -> (a + b t_l)/g), d = deg_l P, g > 0: a positive multiple of
    P composed with the affine map.  Result is made primitive."""
    assert g > 0 and b != 0
    s = W * l
    d = max((k >> s) & FM for k in P) if P else 0
    tab = []
    for e in range(d + 1):
        ge = g ** (d - e)
        tab.append([ge * comb(e, j) * a ** (e - j) * b ** j for j in range(e + 1)])
    R = {}
    get = R.get
    for k, c in P.items():
        e = (k >> s) & FM
        base = k - (e << s)
        row = tab[e]
        for j in range(e + 1):
            cf = row[j]
            if cf:
                kk = base + (j << s)
                R[kk] = get(kk, 0) + c * cf
    return prim({k: v for k, v in R.items() if v})


def reflect(P, l):
    """P(t_l -> 1 - t_l) (exact, no scaling)."""
    return affine(P, l, 1, -1, 1) if P else P


def blowup_map(P, m, S):
    """P(t_j -> t_m t_j for j in S, j != m; t_m unchanged).  Exact (a monomial map)."""
    others = [j for j in S if j != m]
    sm = W * m
    R = {}
    for k, c in P.items():
        tot = 0
        for j in others:
            tot += (k >> (W * j)) & FM
        R[k + (tot << sm)] = c
    _chk(degs(R)) if R else None
    return R


def subst(P, mp):
    """Simultaneous substitution t_l -> mp[l] (polynomials) for l in mp; other variables
    are kept.  Exact."""
    keys = sorted(mp)
    groups = {}
    for k, c in P.items():
        es = tuple((k >> (W * l)) & FM for l in keys)
        base = k
        for l, e in zip(keys, es):
            base -= e << (W * l)
        groups.setdefault(es, {})[base] = c
    cache = {}

    def power(l, e):
        if (l, e) not in cache:
            cache[(l, e)] = ONE if e == 0 else mul(power(l, e - 1), mp[l])
        return cache[(l, e)]
    R = {}
    for es, rest in groups.items():
        X = ONE
        for l, e in zip(keys, es):
            if e:
                X = mul(X, power(l, e))
        R = add(R, mul(X, rest))
    return R


# ---------- exact division and square root (used to check 'factor' certificates) -----

def _divisible(k, kf):
    return ((k | GUARD) - kf) & GUARD == GUARD


def exact_div(R, f):
    """Q with R == Q*f exactly over Z, or None if f does not divide R in Z[t]."""
    if not f:
        return None
    kf = max(f)
    cf = f[kf]
    R = dict(R)
    Q = {}
    heap = [-k for k in R]
    heapq.heapify(heap)
    while R:
        while -heap[0] not in R:
            heapq.heappop(heap)
        k = -heapq.heappop(heap)
        c = R[k]
        if not _divisible(k, kf):
            return None
        q, r = divmod(c, cf)
        if r:
            return None
        dk = k - kf
        Q[dk] = q
        for kk, cc in f.items():
            nk = kk + dk
            v = R.get(nk, 0) - q * cc
            if v:
                if nk not in R:
                    heapq.heappush(heap, -nk)
                R[nk] = v
            else:
                R.pop(nk, None)
    return Q


def sqrt_poly(R):
    """S with S*S == R exactly (R integer, positive leading coefficient), or None."""
    if not R:
        return {}
    dR = degs(R)
    k0 = max(R)
    c0 = R[k0]
    if c0 <= 0:
        return None
    r0 = isqrt(c0)
    if r0 * r0 != c0:
        return None
    e0 = exps(k0)
    if any(e % 2 for e in e0):
        return None
    s0 = mono([e // 2 for e in e0])
    S = {s0: r0}
    res = dict(R)
    del res[k0]
    heap = [-k for k in res]
    heapq.heapify(heap)

    def upd(k, v):
        nv = res.get(k, 0) + v
        if nv:
            if k not in res:
                heapq.heappush(heap, -k)
            res[k] = nv
        else:
            res.pop(k, None)
    while res:
        while -heap[0] not in res:
            heapq.heappop(heap)
        k = -heapq.heappop(heap)
        c = res[k]
        if not _divisible(k, s0):
            return None
        nk = k - s0
        if nk >= s0 or any(2 * e > d for e, d in zip(exps(nk), dR)):
            return None
        q, r = divmod(c, 2 * r0)
        if r:
            return None
        # S <- S + q t^nk ; residual -= 2 q t^nk S + q^2 t^(2nk)
        for ks, cs in list(S.items()):
            upd(ks + nk, -2 * q * cs)
        upd(2 * nk, -q * q)
        S[nk] = q
    return S


def parse(s, names=('t0', 't1', 't2', 't3', 't4', 't5')):
    """Parse an integer polynomial written as 'c*t0^2*t3 - t1 + 4' (python-flint's
    str() format, as stored in the certificate JSON)."""
    idx = {n: i for i, n in enumerate(names)}
    s = s.strip().replace(' - ', ' + -')
    P = {}
    for term in s.split(' + '):
        term = term.strip()
        sg = 1
        if term.startswith('-'):
            sg = -1
            term = term[1:].strip()
        c = 1
        e = [0] * NV
        for fac in term.split('*'):
            if fac[0].isdigit():
                c *= int(fac)
            else:
                if '^' in fac:
                    v, p = fac.split('^')
                    e[idx[v]] += int(p)
                else:
                    e[idx[fac]] += 1
        k = mono(e)
        P[k] = P.get(k, 0) + sg * c
    return {k: v for k, v in P.items() if v}


def evalq(P, pt):
    """Exact value at a point (list of Fractions/ints)."""
    from fractions import Fraction
    tot = Fraction(0)
    for k, c in P.items():
        v = Fraction(c)
        for l in range(NV):
            e = (k >> (W * l)) & FM
            if e:
                v *= Fraction(pt[l]) ** e
        tot += v
    return tot

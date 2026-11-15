"""Finite fields GF(q), q = p^e, built from scratch (stdlib only).

Elements are the integers 0..q-1.  For q = p prime the arithmetic is mod p.
For q = p^e, e > 1, an element is the polynomial sum_i d_i x^i with
d_i the base-p digits of the integer; the field is GF(p)[x]/(f) for the
lexicographically first monic irreducible f of degree e (found by brute
force: no root / no factor of degree <= e/2).

All tables are checked against the field axioms by brute force in
`check_axioms` (q <= 16 is instantaneous).
"""
import itertools


def prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            e, r = 0, q
            while r % p == 0:
                r //= p
                e += 1
            if r != 1:
                raise ValueError("not a prime power: %d" % q)
            return p, e
    raise ValueError("bad q")


def _poly_mulmod(a, b, f, p):
    """a, b, f: coefficient lists, low degree first; f monic of degree e."""
    e = len(f) - 1
    prod = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                prod[i + j] = (prod[i + j] + x * y) % p
    # reduce modulo f
    for d in range(len(prod) - 1, e - 1, -1):
        c = prod[d]
        if c:
            for k in range(e + 1):
                prod[d - e + k] = (prod[d - e + k] - c * f[k]) % p
    prod = prod[:e] + [0] * (e - len(prod[:e]))
    return prod


def _is_irreducible(f, p):
    """f monic, coefficients low->high.  Brute force over all monic factors of degree 1..e//2."""
    e = len(f) - 1
    for d in range(1, e // 2 + 1):
        for coeffs in itertools.product(range(p), repeat=d):
            g = list(coeffs) + [1]
            # polynomial remainder f mod g
            r = list(f)
            for deg in range(len(r) - 1, d - 1, -1):
                c = r[deg]
                if c:
                    for k in range(d + 1):
                        r[deg - d + k] = (r[deg - d + k] - c * g[k]) % p
            if all(x == 0 for x in r[:d]):
                return False
    return True


class GF:
    def __init__(self, q):
        self.q = q
        self.p, self.e = prime_power(q)
        p, e = self.p, self.e
        if e == 1:
            self.modulus = None
            self.add = [[(a + b) % p for b in range(q)] for a in range(q)]
            self.mul = [[(a * b) % p for b in range(q)] for a in range(q)]
        else:
            f = None
            for coeffs in itertools.product(range(p), repeat=e):
                cand = list(coeffs) + [1]
                if _is_irreducible(cand, p):
                    f = cand
                    break
            assert f is not None
            self.modulus = f

            def to_poly(a):
                return [(a // p ** i) % p for i in range(e)]

            def from_poly(c):
                return sum(c[i] * p ** i for i in range(e))

            self.add = [[from_poly([(x + y) % p for x, y in zip(to_poly(a), to_poly(b))])
                         for b in range(q)] for a in range(q)]
            self.mul = [[from_poly(_poly_mulmod(to_poly(a), to_poly(b), f, p))
                         for b in range(q)] for a in range(q)]
        self.neg = [next(b for b in range(q) if self.add[a][b] == 0) for a in range(q)]
        self.inv = [0] + [next(b for b in range(1, q) if self.mul[a][b] == 1) for a in range(1, q)]
        # Frobenius x -> x^p
        self.frob = [self._pow(a, p) for a in range(q)]
        self.check_axioms()
        # a primitive element (generator of the multiplicative group)
        self.prim = next(g for g in range(2, q) if self._order(g) == q - 1) if q > 2 else 1

    def _pow(self, a, n):
        r = 1
        for _ in range(n):
            r = self.mul[r][a]
        return r

    def _order(self, a):
        r, n = a, 1
        while r != 1:
            r = self.mul[r][a]
            n += 1
        return n

    def check_axioms(self):
        q, A, M = self.q, self.add, self.mul
        rng = range(q)
        for a in rng:
            assert A[a][0] == a and M[a][1] == a and M[a][0] == 0
            assert A[a][self.neg[a]] == 0
            if a:
                assert M[a][self.inv[a]] == 1
        for a in rng:
            for b in rng:
                assert A[a][b] == A[b][a] and M[a][b] == M[b][a]
                for c in rng:
                    assert A[A[a][b]][c] == A[a][A[b][c]]
                    assert M[M[a][b]][c] == M[a][M[b][c]]
                    assert M[a][A[b][c]] == A[M[a][b]][M[a][c]]
        # no zero divisors
        for a in range(1, q):
            for b in range(1, q):
                assert M[a][b] != 0
        # Frobenius is a field automorphism
        for a in rng:
            for b in rng:
                assert self.frob[A[a][b]] == A[self.frob[a]][self.frob[b]]
                assert self.frob[M[a][b]] == M[self.frob[a]][self.frob[b]]

    def sub(self, a, b):
        return self.add[a][self.neg[b]]

    def pw(self, a, n):
        """a^n for n >= 0 by square and multiply."""
        r = 1
        while n:
            if n & 1:
                r = self.mul[r][a]
            a = self.mul[a][a]
            n >>= 1
        return r


if __name__ == "__main__":
    for q in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25):
        F = GF(q)
        print(q, "p=%d e=%d" % (F.p, F.e), "modulus", F.modulus, "prim", F.prim)

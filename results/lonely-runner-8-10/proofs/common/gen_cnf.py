"""Gate CNF generator (TRUSTED: part of the stated trusted base of every proof package).

    python gen_cnf.py K L P [out.cnf]        (stdout if no file given)

Writes the DIMACS CNF  G(K, L, P)  for a prime P and a prime-power level L with P not
dividing L.  Pure Python 3, no dependencies, deterministic (byte-identical output).

THEOREM (proof in ../../proofs/ENCODING.md, section 2):
    if G(K, L, P) is unsatisfiable then the improper set I(K, L, P) is empty,
where a K-tuple x of residues mod N = L*P, none divisible by P, is PROPER iff
    (a) gcd(x_j : j != i, N) > 1 for some i, or
    (b) some t in Z/N has (K+1) * dist_N(t*x_i) >= N for every i,
and improper otherwise (dist_N(r) = min(r mod N, N - r mod N)).

Only that direction is used, and its proof is a model construction: from an improper
tuple build the assignment described in each block's comment below and check that it
satisfies every clause.  No clause has to be "complete"; extra (valid) clauses only help
the solver.  Normal form used by the model construction (ENCODING.md 2.1):
    L = q^e.  An improper x has >= 2 coordinates prime to q (else (a) holds).  One of them
    is a unit mod N; multiplying by its inverse and folding signs keeps (im)properness.
    So WLOG x contains 1, and S := {fold_N(x_j)} minus {1} is a subset of
    D = {2..floor(N/2)} minus PZ with |S| <= K-1, and if |S| = K-1 then some s in S is
    prime to q.  Every t has some coordinate with (K+1) dist_N(t x_j) < N ("bad").
"""
import sys


def fold(r, n):
    r %= n
    return min(r, n - r)


def prime_power_base(L):
    q = next(d for d in range(2, L + 1) if L % d == 0)
    m = L
    while m % q == 0:
        m //= q
    if m != 1:
        raise ValueError('level must be a prime power')
    return q


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


class Builder:
    def __init__(self):
        self.nv = 0
        self.clauses = []

    def new(self, n):
        a = self.nv + 1
        self.nv += n
        return list(range(a, a + n))

    def add(self, c):
        self.clauses.append(list(c))

    def at_most(self, lits, K):
        """Sinz sequential counter, at-most-K over lits (K >= 1).  Aux r[i][j-1] means
        '>= j of lits[0..i] are true' (j = 1..min(i+1, K)).  MODEL: set every aux to that
        exact meaning; all clauses hold whenever at most K lits are true.  Returns the last
        row (r[n-1][j-1] = '>= j of all lits')."""
        prev = None
        for i, x in enumerate(lits):
            row = self.new(min(i + 1, K))
            self.add([-x, row[0]])                                   # x -> >=1
            if prev is not None:
                for j in range(len(prev)):
                    self.add([-prev[j], row[j]])                     # monotone in i
                for j in range(1, len(row)):
                    self.add([-x, -prev[j - 1], row[j]])             # x and >=j -> >=j+1
                if len(prev) == K:
                    self.add([-x, -prev[K - 1]])                     # never K+1
            prev = row
        return prev


def generate(K, L, P):
    if not (is_prime(P) and L % P and K >= 2 and P > 2):
        raise ValueError('need prime P > 2 not dividing L, K >= 2')
    q = prime_power_base(L)
    N = L * P
    H = (P - 1) // 2
    D = [v for v in range(2, N // 2 + 1) if v % P]
    cls = {v: fold(v, P) for v in D}            # class of v mod P, folded into 1..H
    b = Builder()
    nclauses_Y = None

    # ---- Block A (variables y_c and their counter; clauses mention only these) ----------
    # y_c, c = 2..H : MODEL y_c = [some s in S has fold_P(s) = c].
    Y = b.new(H - 1)
    y = {c: Y[c - 2] for c in range(2, H + 1)}
    # A1. level-1 time clauses.  For s in 1..H with (K+1)s >= P:  OR_{c bad at s} y_c,
    #     c bad at s iff (K+1) dist_P(s c) < P.  MODEL: at the level-L time t = L s,
    #     dist_N(t v) = L dist_P(s v), so v is bad at t iff its class is bad at s; the
    #     coordinates of class 1 are good at t because (K+1)s >= P; x is improper, so some
    #     s' in S of class c >= 2 is bad, and y_c is true.
    for s in range(1, H + 1):
        if (K + 1) * s >= P:
            b.add([y[c] for c in range(2, H + 1) if (K + 1) * fold(s * c, P) < P])
    # A2. at most K-1 of the y_c.  MODEL: the classes used by S number <= |S| <= K-1.
    b.at_most(Y, K - 1)
    nclauses_Y = len(b.clauses)

    # ---- Block B -----------------------------------------------------------------------
    # x_v, v in D (ascending) : MODEL x_v = [v in S].
    X = b.new(len(D))
    x = dict(zip(D, X))
    # B1. time clauses.  For t in [ceil(N/(K+1)), floor(N/2)]:  OR_{v in D bad at t} x_v.
    #     MODEL: coordinate 1 is good at t (dist_N(t) = t, (K+1)t >= N); every other
    #     coordinate equals 1 or lies in S (up to sign), so some s in S is bad at t.
    #     Times t in (N/2, N) mirror N - t; t < N/(K+1) need no clause.
    for t in range(-(-N // (K + 1)), N // 2 + 1):
        b.add([x[v] for v in D if (K + 1) * fold(t * v, N) < N])
    # B2. x_v -> y_cls(v) for cls(v) >= 2;  B3. y_c -> OR_{cls(v) = c} x_v.
    #     MODEL: both hold by the definitions of x and y.
    members = {c: [] for c in range(1, H + 1)}
    for v in D:
        members[cls[v]].append(x[v])
    for v in D:
        if cls[v] >= 2:
            b.add([-x[v], y[cls[v]]])
    for c in range(2, H + 1):
        b.add([-y[c]] + members[c])
    # B4. |S| <= K-1, counted as  #(classes >= 2 used) + #(class-1 members of S)
    #     + sum_c (#members of class c in S - 1)_+ .  Per class c >= 2 a counter over its
    #     members (at most K-1; MODEL: |S| <= K-1) whose outputs '>= j' for j >= 2 are the
    #     'extras'; then one counter over Z = Y ++ class-1 x's ++ extras, at most K-1.
    #     MODEL: with exact aux values the number of true literals of Z equals |S|.
    extras = []
    for c in range(2, H + 1):
        if len(members[c]) >= 2:
            last = b.at_most(members[c], K - 1)
            extras.extend(last[1:])                     # '>= 2', '>= 3', ...
    Z = Y + members[1] + extras
    last = b.at_most(Z, K - 1)
    # B5. |S| = K-1 -> some s in S prime to q.  MODEL: if |S| = K-1 the tuple is (1, S)
    #     with no padding; (a) fails, so >= 2 coordinates are prime to q (1 is one).
    if len(last) == K - 1:
        b.add([-last[K - 2]] + [x[v] for v in D if v % q])
    meta = dict(K=K, L=L, P=P, N=N, q=q, nY=len(Y), nclauses_Y=nclauses_Y, x_first=X[0] if X else None,
                nv=b.nv, nclauses=len(b.clauses))
    return b, meta


def dimacs(b, meta):
    head = ('c lonely-runner gate CNF v1  K=%(K)d L=%(L)d P=%(P)d  N=%(N)d q=%(q)d\n'
            'c y vars 1..%(nY)d ; clauses 1..%(nclauses_Y)d mention only block-A vars\n' % meta)
    body = ''.join(' '.join(map(str, c)) + ' 0\n' for c in b.clauses)
    return head + 'p cnf %d %d\n' % (b.nv, len(b.clauses)) + body


if __name__ == '__main__':
    K, L, P = map(int, sys.argv[1:4])
    b, meta = generate(K, L, P)
    s = dimacs(b, meta)
    if len(sys.argv) > 4:
        with open(sys.argv[4], 'w', newline='\n') as f:
            f.write(s)
    else:
        sys.stdout.write(s)

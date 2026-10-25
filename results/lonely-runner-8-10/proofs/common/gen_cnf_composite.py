"""Gate CNF for a level L with SEVERAL prime factors (e.g. k = 9, L = 10).  TRUSTED for such k
only; prime-power levels use gen_cnf.py.

    python gen_cnf_composite.py K L P [out.cnf]

THEOREM (ENCODING.md section 5): if this CNF is UNSAT then I(K, L, P) is empty.
Model construction from an improper tuple x mod N = L*P:
  * (a) fails  <=>  for every prime q | L at least two coordinates are prime to q
    (P divides no coordinate, so a common factor of all-but-one would be some q | L).
  * Anchor: let Q = q*^e be the largest prime-power part of L and M = L / Q.  Some x_i is prime
    to q*, hence a unit mod Q*P.  Take the unit u of Z/N with u = x_i^-1 (mod Q*P), u = 1
    (mod M).  Multiplying by u and folding signs keeps (im)properness, and now
    u*x_i = 1 (mod Q*P): the anchor lies in A = fold_N{a : a = 1 mod Q*P}.
  * S := set of the folded values of ALL coordinates (anchor included); S is a subset of
    D' = {1..floor(N/2)} minus PZ, |S| <= K, S meets A.
  * For each q | L: some s in S is prime to q, and if |S| = K (no repeated value) then at
    least two are.
  MODEL: x_v = [v in S], y_c = [S has an element of class c >= 2], every auxiliary at its
  intended meaning (counter rows 'at least j of the first i', o_i = 'one of the first i',
  t_i = 'two of the first i').
"""
import math
import sys

import gen_cnf
from gen_cnf import Builder, fold, is_prime


def prime_factors(n):
    out, d = [], 2
    while n > 1:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    return out


def generate(K, L, P):
    if not (is_prime(P) and L % P and K >= 2 and P > 2 and len(prime_factors(L)) >= 2):
        raise ValueError('need prime P > 2 not dividing L, K >= 2, L with >= 2 prime factors')
    N = L * P
    H = (P - 1) // 2
    qs = prime_factors(L)
    Qbest = max(q ** max(e for e in range(1, 64) if L % q ** e == 0) for q in qs)
    M = L // Qbest
    anchors = sorted({fold(1 + j * Qbest * P, N) for j in range(M)})
    D = [v for v in range(1, N // 2 + 1) if v % P]
    cls = {v: fold(v, P) for v in D}
    b = Builder()
    # Block A: identical in meaning to gen_cnf (anchor has class 1, so classes >= 2 used by S
    # number at most |S| - 1 <= K - 1).
    Y = b.new(H - 1)
    y = {c: Y[c - 2] for c in range(2, H + 1)}
    for s in range(1, H + 1):
        if (K + 1) * s >= P:      # at time L*s class-1 values (the anchor among them) are good
            b.add([y[c] for c in range(2, H + 1) if (K + 1) * fold(s * c, P) < P])
    b.at_most(Y, K - 1)
    nA = len(b.clauses)
    X = b.new(len(D))
    x = dict(zip(D, X))
    # anchor present
    b.add([x[a] for a in anchors])
    # time clauses, every t in 1..floor(N/2) (t = 0: everything bad; N - t mirrors t)
    for t in range(1, N // 2 + 1):
        b.add([x[v] for v in D if (K + 1) * fold(t * v, N) < N])
    members = {c: [] for c in range(1, H + 1)}
    for v in D:
        members[cls[v]].append(x[v])
    for v in D:
        if cls[v] >= 2:
            b.add([-x[v], y[cls[v]]])
    for c in range(2, H + 1):
        b.add([-y[c]] + members[c])
    # |S| <= K, counted as #classes>=2 used + #class-1 values + extras
    extras = []
    for c in range(2, H + 1):
        if len(members[c]) >= 2:
            extras.extend(b.at_most(members[c], K)[1:])
    Z = Y + members[1] + extras
    last = b.at_most(Z, K)
    full = last[K - 1] if len(last) == K else None        # '>= K of Z' (forced upward)
    # gcd conditions per prime q | L
    for q in qs:
        cop = [x[v] for v in D if v % q]
        b.add(cop)                                           # some s prime to q
        if full is None:
            continue
        # t_n = 'two of cop are true', with downward-sound definitions:
        #   o_i -> o_{i-1} or cop_i ;  t_i -> t_{i-1} or cop_i ;  t_i -> t_{i-1} or o_{i-1}
        o_prev = t_prev = None
        for i, lit in enumerate(cop):
            o = b.new(1)[0]
            b.add([-o, lit] + ([o_prev] if o_prev else []))
            if i >= 1:
                t = b.new(1)[0]
                b.add([-t, lit] + ([t_prev] if t_prev else []))
                b.add([-t, o_prev] + ([t_prev] if t_prev else []))
                t_prev = t
            o_prev = o
        if t_prev is None:
            b.add([-full])
        else:
            b.add([-full, t_prev])                           # |S| = K -> two prime to q
    meta = dict(K=K, L=L, P=P, N=N, q=qs, anchors=anchors, nY=len(Y), nclauses_Y=nA,
                x_first=X[0], nv=b.nv, nclauses=len(b.clauses))
    return b, meta


def dimacs(b, meta):
    head = ('c lonely-runner gate CNF composite-v1  K=%d L=%d P=%d  N=%d anchors=%s\n'
            % (meta['K'], meta['L'], meta['P'], meta['N'], meta['anchors']))
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

"""Validation of gen_cnf.py against brute force (NOT part of the trusted base; evidence for it).

Two kinds of check, for small (K, L, P):

FULL  : enumerate ALL K-tuples of B(L,P) = {0..N-1} minus PZ straight from the definition
        (improper_def, no normal form).  For each improper tuple, map it to its support S
        (unit-normalise, fold, drop 1s) and check that the CNF has a model extending
        x_v = [v in S], y_c = [S has class c]   -- this is exactly the model-construction
        argument of ENCODING.md, tested.  Then enumerate every model of the CNF projected on
        the x variables and check that each decodes to an improper tuple (so the CNF is exact,
        not just sound), and that the two support sets coincide.
NF    : for larger cases enumerate supports S directly (all subsets of D with |S| <= K-1)
        with an independent bitmask coverage test, and compare the set of improper supports
        with the model projection, as above.

usage: python test_encoding.py            (all cases, ~minutes)
"""
import itertools
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '.pydeps'))
from pysat.solvers import Solver  # noqa: E402
import gen_cnf  # noqa: E402


def improper_def(x, K, L, P):
    N = L * P
    for i in range(K):
        g = N
        for j in range(K):
            if j != i:
                g = math.gcd(g, x[j])
        if g > 1:
            return False
    for t in range(N):
        if all((K + 1) * min(t * v % N, N - t * v % N) >= N for v in x):
            return False
    return True


def support(x, L, P):
    N = L * P
    u = next(v for v in x if math.gcd(v, N) == 1)
    ui = pow(u, -1, N)
    return frozenset(f for f in (gen_cnf.fold(v * ui, N) for v in x) if f != 1)


def var_maps(K, L, P):
    N = L * P
    H = (P - 1) // 2
    D = [v for v in range(2, N // 2 + 1) if v % P]
    y = {c: c - 1 for c in range(2, H + 1)}
    b, meta = gen_cnf.generate(K, L, P)
    # x vars start right after block A's vars: recover from the first B1 clause's literals
    x = {v: meta['x_first'] + i for i, v in enumerate(D)}
    return b, meta, D, x, y


def model_supports(b, D, x, limit=10 ** 6):
    out = set()
    xs = [x[v] for v in D]
    inv = {x[v]: v for v in D}
    with Solver(name='cadical195', bootstrap_with=b.clauses) as s:
        while s.solve():
            m = s.get_model()
            S = frozenset(inv[l] for l in xs if m[l - 1] > 0)
            out.add(S)
            s.add_clause([-x[v] if v in S else x[v] for v in D])
            assert len(out) < limit
    return out


def extends(b, D, x, y, P, S):
    a = [x[v] if v in S else -x[v] for v in D]
    used = {gen_cnf.fold(v, P) for v in S}
    a += [y[c] if c in used else -y[c] for c in y]
    with Solver(name='cadical195', bootstrap_with=b.clauses) as s:
        return s.solve(assumptions=a)


def full_case(K, L, P):
    N = L * P
    B = [v for v in range(N) if v % P]
    b, meta, D, x, y = var_maps(K, L, P)
    imp = set()
    n = 0
    for tup in itertools.product(B, repeat=K):
        if improper_def(tup, K, L, P):
            n += 1
            imp.add(support(tup, L, P))
    for S in imp:
        assert extends(b, D, x, y, P, S), ('model construction fails', K, L, P, sorted(S))
    ms = model_supports(b, D, x)
    for S in ms:
        tup = [1] + sorted(S) + [1] * (K - 1 - len(S))
        assert improper_def(tup, K, L, P), ('model not improper', sorted(S))
    assert ms == imp, (K, L, P, len(ms), len(imp))
    return n, len(imp)


def nf_case(K, L, P):
    N = L * P
    q = gen_cnf.prime_power_base(L)
    b, meta, D, x, y = var_maps(K, L, P)
    T = range(-(-N // (K + 1)), N // 2 + 1)
    full = 0
    bad = {}
    for v in D:
        m = 0
        for i, t in enumerate(T):
            if (K + 1) * min(t * v % N, N - t * v % N) < N:
                m |= 1 << i
        bad[v] = m
    full = (1 << len(T)) - 1
    imp = set()
    for r in range(0, K):
        for S in itertools.combinations(D, r):
            m = 0
            for v in S:
                m |= bad[v]
            if m == full and (r < K - 1 or any(v % q for v in S)):
                imp.add(frozenset(S))
    # spot-check the bitmask test against the definition
    for S in list(imp)[:200]:
        assert improper_def([1] + sorted(S) + [1] * (K - 1 - len(S)), K, L, P)
    for S in imp:
        assert extends(b, D, x, y, P, S), ('model construction fails', K, L, P, sorted(S))
    ms = model_supports(b, D, x)
    assert ms == imp, (K, L, P, len(ms), len(imp))
    return len(imp)


FULL = [(2, 3, 5), (2, 3, 7), (2, 3, 11), (2, 3, 13), (2, 4, 5), (2, 2, 7),
        (3, 4, 5), (3, 4, 7), (3, 4, 11), (3, 4, 13), (3, 2, 7), (3, 8, 5), (3, 9, 5),
        (4, 5, 3), (4, 5, 7), (4, 4, 7), (4, 8, 3)]
NF = [(4, 5, 11), (4, 5, 13), (4, 5, 17), (6, 7, 5), (6, 7, 11), (6, 7, 13),
      (7, 8, 3), (7, 8, 5), (7, 8, 7), (8, 9, 5), (8, 9, 7), (6, 8, 11), (5, 8, 11), (5, 9, 7),
      (7, 8, 11), (7, 8, 13), (6, 7, 17), (4, 5, 19), (3, 4, 17), (8, 9, 5)]

if __name__ == '__main__' and '--composite' not in sys.argv:
    t0 = time.time()
    for c in FULL:
        n, s = full_case(*c)
        print('FULL K=%d L=%d P=%d: %d improper tuples, %d supports, CNF exact' % (c + (n, s)), flush=True)
    for c in NF:
        s = nf_case(*c)
        print('NF   K=%d L=%d P=%d: %d improper supports, CNF exact' % (c + (s,)), flush=True)
    print('all %d cases pass, %.1f s' % (len(FULL) + len(NF), time.time() - t0))


# ---------------------------------------------------------------- composite levels (gen_cnf_composite)
def composite_case(K, L, P):
    """Soundness direction in full: every improper tuple of B(L,P)^K (brute force from the
    definition) yields, after the anchor normalisation of gen_cnf_composite, an assignment
    that extends to a model.  Also counts models whose support admits no improper padding
    (allowed: the CNF is only claimed sound, not exact)."""
    import gen_cnf_composite as gc
    N = L * P
    b, meta = gc.generate(K, L, P)
    D = [v for v in range(1, N // 2 + 1) if v % P]
    x = {v: meta['x_first'] + i for i, v in enumerate(D)}
    y = {c: c - 1 for c in range(2, (P - 1) // 2 + 1)}
    qs = meta['q']
    Q = max(q ** max(e for e in range(1, 64) if L % q ** e == 0) for q in qs)
    qstar = next(q for q in qs if Q % q == 0)
    B = [v for v in range(N) if v % P]
    sup = set()
    n = 0
    for tup in itertools.product(B, repeat=K):
        if not improper_def(tup, K, L, P):
            continue
        n += 1
        xi = next(v for v in tup if v % qstar)
        inv = pow(xi, -1, Q * P)
        u = next(w for w in range(1, N) if w % (Q * P) == inv and w % (L // Q) == 1 % (L // Q))
        sup.add(frozenset(gen_cnf.fold(u * v, N) for v in tup))
    for S in sup:
        assert set(S) & set(meta['anchors']), S
        assert extends(b, D, x, y, P, S), ('composite model construction fails', K, L, P, sorted(S))
    ms = model_supports(b, D, x)
    assert sup <= ms
    spurious = 0
    for S in ms - sup:
        S = sorted(S)
        ok = any(improper_def(S + list(pad), K, L, P)
                 for pad in itertools.combinations_with_replacement(S, K - len(S)))
        spurious += not ok
    return n, len(sup), len(ms), spurious


COMPOSITE = [(3, 6, 5), (3, 6, 7), (4, 6, 5), (4, 6, 7), (5, 6, 7), (3, 10, 3), (3, 10, 7), (2, 6, 5),
             (4, 10, 3), (5, 6, 5)]

if __name__ == '__main__' and '--composite' in sys.argv:
    for c in COMPOSITE:
        n, s, m, sp = composite_case(*c)
        print('COMP K=%d L=%d P=%d: %d improper tuples, %d supports, %d model supports (%d spurious)'
              % (c + (n, s, m, sp)), flush=True)

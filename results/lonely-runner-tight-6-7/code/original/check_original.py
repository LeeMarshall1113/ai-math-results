"""Checker for the first implementation (lr.py + tight6.py + tight7.py).

Usage: python check_original.py OUT_K6.json OUT_K7.json [GATES_DIR]

Reads the sieve outputs written by `tight6.py sieve 997 OUT_K6.json` and
`tight7.py 11 1153 OUT_K7.json`, then
  1. lists the gates and non-gates;
  2. checks the bookkeeping of PROOF.md Section 6 in exact integer arithmetic:
        prod_{gates} p  >  (S0/k)^k * prod_T max_j U_j(T),   S0 = m^(k-1),
        U_j(T) = C(k,j) (S0^2/k)^j E_1(T)^j + S0^(2j) E_j(T),  E_j(T) = e_j(T^2);
  3. if GATES_DIR is given, compares every per-prime verdict and level-one orbit count
     with the published gate lists gates_k6.json / gates_k7.json.
Exit status 0 iff everything holds.
"""
import sys, os, json, math, itertools
from fractions import Fraction

KNOWN = {6: [(1, 2, 3, 4, 5, 6)],
         7: [(1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]}


def e(xs, j):
    return sum(math.prod(c) for c in itertools.combinations(xs, j))


def needed_bound(k):
    m = math.comb(k + 1, 2) + 1
    S0 = m ** (k - 1)
    B = Fraction(S0, k) ** k
    for T in KNOWN[k]:
        w = [a * a for a in T]
        E1 = e(w, 1)
        B *= max(math.comb(k, j) * Fraction(S0 * S0, k) ** j * E1 ** j + S0 ** (2 * j) * e(w, j)
                 for j in range(2, k + 1))
    return m, S0, B


def lnfrac(x):
    return math.log(x.numerator) - math.log(x.denominator)


def bookkeeping(k, rows, cutoffs):
    m, S0, B = needed_bound(k)
    ok = True
    print(f'k={k}: m={m}, S0=m^{k-1}={S0}, needed gate mass log(B) = {lnfrac(B):.6f} nats')
    for pmax in cutoffs:
        gates = [r['p'] for r in rows if r['gate'] and r['p'] <= pmax]
        G = math.prod(gates)
        holds = G * B.denominator > B.numerator
        print(f'  primes 11..{pmax}: {len(gates)} gates, mass {math.log(G):.6f} nats, '
              f'margin {math.log(G) - lnfrac(B):.6f} nats, exact inequality holds: {holds}')
        if pmax == cutoffs[-1]:
            ok &= holds
    return ok


def compare(k, rows, gates_file):
    pub = {r['p']: r for r in json.load(open(gates_file))['primes']}
    mine = {r['p']: r for r in rows}
    bad = [p for p in pub if p not in mine or mine[p]['gate'] != pub[p]['gate']
           or mine[p]['level1_orbits'] != pub[p]['level1_orbits']]
    extra = [p for p in mine if p not in pub]
    print(f'  k={k}: {len(pub)} published primes, verdict/orbit-count mismatches: {bad}, '
          f'unpublished primes: {extra}')
    return not bad and not extra


if __name__ == '__main__':
    r6 = json.load(open(sys.argv[1]))
    r7 = json.load(open(sys.argv[2]))
    ok = True
    for k, rows in ((6, r6), (7, r7)):
        ps = [r['p'] for r in rows]
        print(f'k={k}: primes {ps[0]}..{ps[-1]} ({len(ps)}), non-gates '
              f'{[r["p"] for r in rows if not r["gate"]]}')
    ok &= bookkeeping(6, r6, [997])
    ok &= bookkeeping(7, r7, [1129, 1153])
    if len(sys.argv) > 3:
        ok &= compare(6, r6, os.path.join(sys.argv[3], 'gates_k6.json'))
        ok &= compare(7, r7, os.path.join(sys.argv[3], 'gates_k7.json'))
    print('ORIGINAL CHECK', 'OK' if ok else 'FAILED')
    sys.exit(0 if ok else 1)

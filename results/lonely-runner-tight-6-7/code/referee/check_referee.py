"""Checker for the second implementation (ref_level1.py + ref_lift.py + ref_gates.py).

Usage: python check_referee.py OUT_K6.jsonl OUT_K7.jsonl [GATES_DIR]

Uses ref_bounds.bound (exact Fractions) for the inequality of PROOF.md Section 6 and,
if GATES_DIR is given, compares every per-prime verdict and level-one orbit count with
the published gate lists.  Exit status 0 iff everything holds.
"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_bounds as RB


def load(fn):
    return [json.loads(l) for l in open(fn) if l.strip()]


def check(k, rows, cutoffs):
    m = math.comb(k + 1, 2) + 1
    S0 = m ** (k - 1)
    B, _, _ = RB.bound(k, S0)
    lb = math.log(B.numerator) - math.log(B.denominator)
    ok = all(r['exempt_sane'] for r in rows)
    print(f'k={k}: m={m}, S0={S0}, needed log(B) = {lb:.6f} nats, exempt orbits sane: {ok}')
    for pmax in cutoffs:
        G = math.prod(r['p'] for r in rows if r['gate'] and r['p'] <= pmax)
        holds = G * B.denominator > B.numerator
        print(f'  primes 11..{pmax}: {sum(1 for r in rows if r["gate"] and r["p"] <= pmax)} gates, '
              f'margin {math.log(G) - lb:.6f} nats, exact inequality holds: {holds}')
        if pmax == cutoffs[-1]:
            ok &= holds
    return ok


def compare(k, rows, gates_file):
    pub = {r['p']: r for r in json.load(open(gates_file))['primes']}
    mine = {r['p']: r for r in rows}
    bad = [p for p in pub if p not in mine or mine[p]['gate'] != pub[p]['gate']
           or mine[p]['level1'] != pub[p]['level1_orbits']]
    extra = [p for p in mine if p not in pub]
    print(f'  k={k}: {len(pub)} published primes, verdict/orbit-count mismatches: {bad}, '
          f'unpublished primes: {extra}')
    return not bad and not extra


if __name__ == '__main__':
    r6, r7 = load(sys.argv[1]), load(sys.argv[2])
    ok = True
    for k, rows in ((6, r6), (7, r7)):
        print(f'k={k}: {len(rows)} primes, non-gates {sorted(r["p"] for r in rows if not r["gate"])}')
    ok &= check(6, r6, [997])
    ok &= check(7, r7, [1129, 1153])
    if len(sys.argv) > 3:
        ok &= compare(6, r6, os.path.join(sys.argv[3], 'gates_k6.json'))
        ok &= compare(7, r7, os.path.join(sys.argv[3], 'gates_k7.json'))
    print('REFEREE CHECK', 'OK' if ok else 'FAILED')
    sys.exit(0 if ok else 1)

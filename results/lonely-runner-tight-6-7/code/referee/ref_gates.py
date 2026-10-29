"""Second implementation (referee): gate sieve.  Usage:
    python ref_gates.py K PMIN PMAX OUT.jsonl
A prime p > k+1 is a TIGHT GATE for k speeds iff no level-one orbit outside the exempt
set (orbits of the known tight k-tuples mod p, skipping a tuple with an entry = 0 mod p)
has a closed-improper, non-Lemma-G lift mod (k+1)p.
Resumable: primes already in OUT are skipped.
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_level1 as L1
import ref_lift as LF

KNOWN = {6: [(1, 2, 3, 4, 5, 6)],
         7: [(1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]}


def isprime(n):
    return n > 1 and all(n % q for q in range(2, int(n ** .5) + 1))


def exempt(k, p):
    ex = set()
    for T in KNOWN[k]:
        if any(a % p == 0 for a in T):
            continue
        ex.add(L1.canon([a % p for a in T], p))
    return ex


def gate(k, p):
    t0 = time.time()
    orbs, st = L1.level_one_orbits(k, p)
    ex = exempt(k, p)
    surv = []
    for o in sorted(orbs):
        if o in ex:
            continue
        h = LF.improper_lifts(o, p, k, limit=1)
        if h:
            surv.append([list(o), list(h[0])])
    # sanity: every exempt orbit present at level one and has an improper lift
    ex_ok = all(o in orbs and LF.improper_lifts(o, p, k, limit=1) for o in ex)
    return dict(k=k, p=p, level1=len(orbs), exempt=sorted(ex), exempt_sane=bool(ex_ok),
                nonexempt=len(orbs - ex), survivors=len(surv), surv_examples=surv[:3],
                gate=(len(surv) == 0), nodes=st['nodes'], secs=round(time.time() - t0, 2))


if __name__ == '__main__':
    k, pmin, pmax, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    done = set()
    if os.path.exists(out):
        for line in open(out):
            done.add(json.loads(line)['p'])
    for p in range(max(pmin, k + 3), pmax + 1):
        if not isprime(p) or p in done:
            continue
        r = gate(k, p)
        with open(out, 'a') as f:
            f.write(json.dumps(r) + '\n')
        print(json.dumps({x: r[x] for x in ('k', 'p', 'level1', 'nonexempt', 'survivors',
                                            'gate', 'exempt_sane', 'secs')}), flush=True)

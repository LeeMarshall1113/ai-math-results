"""Referee validation: (a) pruned level-one DFS == unpruned brute force (all multisets
containing class 1) at small primes, k = 6 and 7; (b) meet-in-the-middle lift == naive
per-digit-vector count on every level-one orbit at p = 11, 13, 17, 19 (k = 6) and
p = 11, 13 (k = 7); (c) full brute-force gate verdicts at those primes."""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_level1 as L1, ref_lift as LF, ref_gates as G

ok = True
for k, ps in ((6, [p for p in range(11, 98) if G.isprime(p)]),
              (7, [p for p in range(11, 62) if G.isprime(p)])):
    for p in ps:
        t = time.time()
        a, _ = L1.level_one_orbits(k, p)
        b = L1.brute_level_one_orbits(k, p)
        print(f'level-one k={k} p={p}: dfs {len(a)} brute {len(b)} equal={a == b} ({time.time()-t:.1f}s)', flush=True)
        ok &= a == b
for k, ps in ((6, [11, 13, 17, 19]), (7, [11, 13])):
    for p in ps:
        orbs = L1.brute_level_one_orbits(k, p)
        ex = G.exempt(k, p)
        bad = 0
        surv = 0
        for o in orbs:
            n1 = len(LF.improper_lifts(o, p, k, want_all=True))
            n2 = LF.naive_improper_count(o, p, k)
            bad += n1 != n2
            surv += (o not in ex) and n2 > 0
        print(f'lift k={k} p={p}: {len(orbs)} orbits, lift mismatches {bad}, brute verdict gate={surv == 0}', flush=True)
        ok &= bad == 0
print('REFEREE VALIDATION', 'OK' if ok else 'FAILED')

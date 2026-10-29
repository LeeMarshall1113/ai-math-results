"""Brute-force check of the k = 7 closed lift at small p: for every multiset of level-1
classes containing 1, test all 8^7 lifts mod 8p at ALL times (closed criterion, Lemma G
parity filter, no pruning); compare the set of orbits with a closed-improper lift to the
pipeline (lr.level_one_reps + tight7._lift7)."""
import os, sys, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tight7, lr
import numpy as np
from numba import njit

@njit(cache=True)
def brute_one(cls, p):
    k = cls.shape[0]; N = 8 * p
    pinv = 0
    for a in range(1, 9):
        if (p * a) % 8 == 1:
            pinv = a
    x = np.zeros(k, dtype=np.int64)
    for code in range(8 ** k):
        c = code; odd = 0
        for i in range(k):
            d = c % 8; c //= 8
            if d % 2 == 1:
                odd += 1
            x[i] = cls[i] + p * ((((d - cls[i]) % 8) * pinv) % 8)
        if odd < 2:
            continue
        proper = False
        for t in range(1, N // 2 + 1):
            ok = True
            for i in range(k):
                r = (t * x[i]) % N
                if r > N - r:
                    r = N - r
                if 8 * r <= N:
                    ok = False; break
            if ok:
                proper = True; break
        if not proper:
            return True
    return False

for p in [int(a) for a in sys.argv[1:]]:
    t = time.time(); H = (p - 1) // 2
    brute = set()
    for rest in itertools.combinations_with_replacement(range(1, H + 1), 6):
        tup = (1,) + rest
        if brute_one(np.array(tup, dtype=np.int64), p):
            brute.add(tuple(int(a) for a in lr.canon(np.array([tup]), p)[0]))
    reps, _ = lr.level_one_reps(7, p)
    out = np.zeros((1, 7), dtype=np.int64)
    pipe = {tuple(int(a) for a in r) for r in reps if tight7._lift7(np.array(r, dtype=np.int64), p, 1, out) > 0}
    print(f'k=7 p={p}: brute {len(brute)} orbits, pipeline {len(pipe)}, equal={brute == pipe}, '
          f'known tight orbits in brute: {len(tight7.exempt_reps(p, tight7.KNOWN) & brute)} ({time.time()-t:.0f}s)', flush=True)

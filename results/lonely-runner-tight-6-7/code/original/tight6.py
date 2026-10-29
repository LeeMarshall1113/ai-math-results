"""Tight-gate sieve for six speeds (seven runners).

Goal (Perarnau-Serra 2025, Problem 1, case of six speeds): every tight sextuple
(ML = 1/7) is a dilate of (1,2,3,4,5,6).

Definitions (k = 6 speeds, N = 7p, p prime > 7):
  * a residue tuple x (mod 7p, no x_i = 0 mod p) is CLOSED-IMPROPER if every time
    t in Z/N has some i with ||t x_i / N|| <= 1/7, i.e. dist(t x_i mod N, 0) <= p.
    (A real tight tuple v with p not dividing prod v reduces to a closed-improper x,
    since t/N is a real time.)
  * p is a TIGHT GATE if every closed-improper tuple mod 7p is, mod p, in the orbit
    of (1,..,6) under unit scaling, permutation and per-coordinate sign.

Pipeline:
  level 1 (times j/p): closed and open covers coincide because 7*dist = p is
  impossible for p != 7, so the repo's validated lr.level_one_reps(6, p) gives every
  level-1 orbit. The tight orbit is removed; every other orbit is lifted to mod 7p
  (all 7^6 digit vectors, reduced by the kernel units w = 1 mod p) and tested with
  the CLOSED criterion. Times t = 0 mod 7 are level-1 times, already covered.

Run:  python tight6.py validate
      python tight6.py sieve PMAX [OUT.json]
"""
import os, sys, math, time, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numba import njit
import lr

K = 6


@njit(cache=True)
def _lift_closed(cls, p, pinv7, use_sym, out, maxout):
    """Count closed-improper lifts of level-1 classes cls (mod p) to mod 7p.
    Returns number found (canonical under kernel units if use_sym)."""
    k = cls.shape[0]
    N = 7 * p
    d = np.zeros(k, dtype=np.int64)
    x = np.zeros(k, dtype=np.int64)
    nfound = 0
    total = 7 ** k
    for code in range(total):
        c = code
        for i in range(k):
            d[i] = c % 7
            c //= 7
        if use_sym:
            minimal = True
            for w in range(2, 7):
                for i in range(k):
                    a = (w * d[i]) % 7
                    if a != d[i]:
                        if a < d[i]:
                            minimal = False
                        break
                if not minimal:
                    break
            if not minimal:
                continue
        # Lemma G: if 7 divides all but at most one coordinate, a primitive real tuple
        # with these residues has ML >= 1/6 (LRC(5) + a shift by j/7), so it is not tight.
        nz = 0
        for i in range(k):
            if d[i] != 0:
                nz += 1
        if nz <= 1:
            continue
        for i in range(k):
            x[i] = cls[i] + p * ((((d[i] - cls[i]) % 7) * pinv7) % 7)
        proper = False
        for t in range(1, N // 2 + 1):
            if t % 7 == 0:
                continue
            ok = True
            for i in range(k):
                r = (t * x[i]) % N
                if r > N - r:
                    r = N - r
                if r <= p:          # covered in the closed sense: ||t x/N|| <= 1/7
                    ok = False
                    break
            if ok:
                proper = True
                break
        if not proper:
            if nfound < maxout:
                for i in range(k):
                    out[nfound, i] = x[i]
            nfound += 1
    return nfound


def tight_rep(p):
    return tuple(int(a) for a in lr.canon(np.array([[1, 2, 3, 4, 5, 6]]), p)[0])


def tight_gate(p, use_sym=True):
    t0 = time.perf_counter()
    reps, st = lr.level_one_reps(K, p)
    trep = tight_rep(p)
    reps_t = [tuple(int(a) for a in r) for r in reps]
    nontight = [r for r in reps_t if r != trep]
    pinv7 = pow(p, -1, 7)
    surv = []
    out = np.zeros((64, K), dtype=np.int64)
    nsurv = 0
    for r in nontight:
        n = _lift_closed(np.array(r, dtype=np.int64), p, pinv7, use_sym, out, 64)
        if n:
            nsurv += n
            surv.append((r, [list(map(int, out[j])) for j in range(min(n, 3))]))
    return dict(p=p, level1_orbits=len(reps_t), tight_in_level1=(trep in reps_t),
                nontight=len(nontight), lifted_survivors=nsurv, gate=(nsurv == 0),
                examples=surv[:3], secs=round(time.perf_counter() - t0, 3))


# ------------------------------------------------------------------ validation
def brute_closed_improper_mod7p(p):
    """Brute force: every multiset of 6 level-1 classes (1..H) containing class 1,
    every digit vector mod 7 (no symmetry), closed criterion at ALL times 1..N/2.
    Returns the set of canonical level-1 class tuples having a closed-improper lift."""
    H = (p - 1) // 2
    N = 7 * p
    pinv7 = pow(p, -1, 7)
    found = set()
    out = np.zeros((1, K), dtype=np.int64)
    for rest in itertools.combinations_with_replacement(range(1, H + 1), K - 1):
        tup = (1,) + rest
        n = _lift_closed_alltimes(np.array(tup, dtype=np.int64), p, pinv7)
        if n:
            found.add(tuple(int(a) for a in lr.canon(np.array([tup]), p)[0]))
    return found


@njit(cache=True)
def _lift_closed_alltimes(cls, p, pinv7):
    k = cls.shape[0]
    N = 7 * p
    d = np.zeros(k, dtype=np.int64)
    x = np.zeros(k, dtype=np.int64)
    nfound = 0
    for code in range(7 ** k):
        c = code
        for i in range(k):
            d[i] = c % 7
            c //= 7
        # Lemma G: if 7 divides all but at most one coordinate, a primitive real tuple
        # with these residues has ML >= 1/6 (LRC(5) + a shift by j/7), so it is not tight.
        nz = 0
        for i in range(k):
            if d[i] != 0:
                nz += 1
        if nz <= 1:
            continue
        for i in range(k):
            x[i] = cls[i] + p * ((((d[i] - cls[i]) % 7) * pinv7) % 7)
        proper = False
        for t in range(1, N // 2 + 1):      # all times, including multiples of 7
            ok = True
            for i in range(k):
                r = (t * x[i]) % N
                if r > N - r:
                    r = N - r
                if r <= p:
                    ok = False
                    break
            if ok:
                proper = True
                break
        if not proper:
            nfound += 1
    return nfound


def ml_exact(v):
    """Exact ML of an integer tuple (Lemma 2.2 of Cordella 2609.03444), as a Fraction."""
    from fractions import Fraction
    best = Fraction(0)
    Ls = set()
    for i in range(len(v)):
        for j in range(i, len(v)):
            for L in (v[i] + v[j], abs(v[i] - v[j])):
                if L:
                    Ls.add(L)
    for L in Ls:
        for m in range(L):
            val = min(min((m * a) % L, L - (m * a) % L) for a in v)
            if val * best.denominator > best.numerator * L:
                best = Fraction(val, L)
    return best


def validate():
    ok = True
    # (1) pipeline vs brute force: which level-1 orbits have closed-improper lifts
    for p in (11, 13):
        t = time.time()
        brute = brute_closed_improper_mod7p(p)
        reps, _ = lr.level_one_reps(K, p)
        pinv7 = pow(p, -1, 7)
        out = np.zeros((1, K), dtype=np.int64)
        pipe = set()
        for r in reps:
            if _lift_closed(np.array(r, dtype=np.int64), p, pinv7, True, out, 1):
                pipe.add(tuple(int(a) for a in r))
        same = brute == pipe
        print(f'validate p={p}: brute {len(brute)} orbits with closed-improper lift, '
              f'pipeline {len(pipe)}, equal={same}, tight rep {tight_rep(p)} in brute: '
              f'{tight_rep(p) in brute}  ({time.time()-t:.1f}s)', flush=True)
        ok &= same and (tight_rep(p) in brute)
    # (2) symmetry reduction does not change existence (p = 17, 19, 23)
    for p in (17, 19, 23):
        reps, _ = lr.level_one_reps(K, p)
        pinv7 = pow(p, -1, 7)
        out = np.zeros((1, K), dtype=np.int64)
        a = [bool(_lift_closed(np.array(r, dtype=np.int64), p, pinv7, True, out, 1)) for r in reps]
        b = [bool(_lift_closed(np.array(r, dtype=np.int64), p, pinv7, False, out, 1)) for r in reps]
        print(f'validate sym p={p}: {a == b}', flush=True)
        ok &= a == b
    # (3) real tuples: (1..6) and its dilates are tight; a few near-tight checks
    from fractions import Fraction
    assert ml_exact((1, 2, 3, 4, 5, 6)) == Fraction(1, 7)
    assert ml_exact((2, 4, 6, 8, 10, 12)) == Fraction(1, 7)
    assert ml_exact((1, 2, 3, 4, 5)) == Fraction(1, 6)
    assert ml_exact((1, 3, 4, 5, 9)) == Fraction(1, 6)
    assert ml_exact((1, 2, 3, 4, 5, 7, 12)) == Fraction(1, 8)   # Wills
    assert ml_exact((1, 4, 5, 6, 7, 11, 13)) == Fraction(1, 8)  # Goddyn-Wong
    print('validate real ML values: ok', flush=True)
    print('VALIDATION', 'OK' if ok else 'FAILED')
    return ok


def sieve(pmax, outpath):
    rows = []
    mass = 0.0
    for p in range(11, pmax + 1):
        if not all(p % q for q in range(2, int(p ** 0.5) + 1)):
            continue
        r = tight_gate(p)
        if r['gate']:
            mass += math.log(p)
        r['gate_mass_so_far'] = round(mass, 3)
        rows.append(r)
        print(json.dumps(r), flush=True)
    with open(outpath, 'w') as f:
        json.dump(rows, f, indent=1)
    return rows


def validate2(primes=(17, 19, 23, 29)):
    """Stronger check at primes with non-tight survivors: the set of level-1 orbits
    with a closed-improper lift, pipeline (lr.level_one_reps + lift skipping t = 0 mod 7)
    versus brute force (all multisets containing 1, all times)."""
    ok = True
    for p in primes:
        t = time.time()
        brute = brute_closed_improper_mod7p(p)
        reps, _ = lr.level_one_reps(K, p)
        pinv7 = pow(p, -1, 7)
        out = np.zeros((1, K), dtype=np.int64)
        pipe = set()
        for r in reps:
            if _lift_closed(np.array(r, dtype=np.int64), p, pinv7, True, out, 1):
                pipe.add(tuple(int(a) for a in r))
        print(f'validate2 p={p}: brute {len(brute)}, pipeline {len(pipe)}, equal={brute == pipe}'
              f' ({time.time()-t:.1f}s)', flush=True)
        ok &= brute == pipe
    print('VALIDATION2', 'OK' if ok else 'FAILED')
    return ok


if __name__ == '__main__':
    if sys.argv[1] == 'validate':
        sys.exit(0 if validate() else 1)
    elif sys.argv[1] == 'validate2':
        sys.exit(0 if validate2() else 1)
    elif sys.argv[1] == 'sieve':
        pmax = int(sys.argv[2])
        outp = sys.argv[3] if len(sys.argv) > 3 else os.path.join(
            os.path.dirname(os.path.abspath(__file__)), f'tight6_gates_{pmax}.json')
        sieve(pmax, outp)

"""Tight-gate sieve for seven speeds (eight runners).

Same logic as tight6.py with k = 7, N = 8p, factor chain 2 -> 4 -> 8.
A residue tuple mod Lp (L | 8) is CLOSED-IMPROPER if every time t has some i with
8 * dist(t x_i mod Lp, 0) <= Lp.  For L = 1, 2, 4 equality is impossible (p odd), so
closed = open there and the repo's level-1 orbits apply; the closed criterion bites at L = 8.
Lemma G (k = 7): if 2 divides all but at most one coordinate, a primitive real tuple has
ML >= min(1/7, 1/4) > 1/8 (LRC(6) + shift by j/2), so it is not tight: such lifts are skipped.

Exempt orbits: level-1 orbits of the known tight septuples passed on the command line
(default: (1..7), Wills (1,2,3,4,5,7,12), Goddyn-Wong (1,4,5,6,7,11,13)).
A prime is a tight gate if no non-exempt level-1 orbit has a closed-improper lift at level 8.

Run: python tight7.py PMIN PMAX [OUT.json]
"""
import os, sys, math, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numba import njit
import lr

K = 7
KNOWN = [(1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]


@njit(cache=True)
def _covered_all(x, N, Lp_closed, k, skip_mod):
    """True iff every t in 1..N/2 (t not = 0 mod skip_mod) has a covered coordinate.
    Covered: 8*dist <= N if Lp_closed else 8*dist < N."""
    for t in range(1, N // 2 + 1):
        if t % skip_mod == 0:
            continue
        ok = True
        for i in range(k):
            r = (t * x[i]) % N
            if r > N - r:
                r = N - r
            if Lp_closed:
                if 8 * r <= N:
                    ok = False
                    break
            else:
                if 8 * r < N:
                    ok = False
                    break
        if ok:
            return False
    return True


@njit(cache=True)
def _lift7(cls, p, maxout, out, closed_final=True):
    """DFS over residues mod 2, 4, 8 (bits b0,b1,b2 per coordinate).  Returns the number
    of closed-improper tuples mod 8p lifting cls (no symmetry reduction)."""
    k = cls.shape[0]
    x = np.zeros(k, dtype=np.int64)
    found = 0
    res2 = np.zeros(k, dtype=np.int64)
    res4 = np.zeros(k, dtype=np.int64)
    res8 = np.zeros(k, dtype=np.int64)
    for c2 in range(2 ** k):
        for i in range(k):
            res2[i] = (c2 >> i) & 1
        if not _lemma_g_ok(res2, 2, k):
            continue
        _crt_vec(cls, res2, 2, p, x)
        if not _covered_all(x, 2 * p, False, k, 2):
            continue
        for c4 in range(2 ** k):
            for i in range(k):
                res4[i] = res2[i] + 2 * ((c4 >> i) & 1)
            _crt_vec(cls, res4, 4, p, x)
            if not _covered_all(x, 4 * p, False, k, 2):
                continue
            for c8 in range(2 ** k):
                for i in range(k):
                    res8[i] = res4[i] + 4 * ((c8 >> i) & 1)
                _crt_vec(cls, res8, 8, p, x)
                if _covered_all(x, 8 * p, closed_final, k, 2):
                    if found < maxout:
                        for i in range(k):
                            out[found, i] = x[i]
                    found += 1
    return found


@njit(cache=True)
def _lemma_g_ok(res, L, k):
    """False iff all but at most one coordinate are even (then not tight)."""
    odd = 0
    for i in range(k):
        if res[i] % 2 == 1:
            odd += 1
    return odd >= 2


@njit(cache=True)
def _crt_vec(cls, res, L, p, x):
    # x = cls mod p, x = res mod L, 0 <= x < L p
    pinv = 0
    for a in range(1, L + 1):
        if (p * a) % L == 1 % L:
            pinv = a
            break
    for i in range(cls.shape[0]):
        x[i] = cls[i] + p * ((((res[i] - cls[i]) % L) * pinv) % L)


def exempt_reps(p, known):
    out = set()
    for T in known:
        if any(a % p == 0 for a in T):
            continue
        out.add(tuple(int(a) for a in lr.canon(np.array([T]), p)[0]))
    return out


def tight_gate7(p, known=KNOWN):
    t0 = time.perf_counter()
    reps, st = lr.level_one_reps(K, p)
    reps_t = [tuple(int(a) for a in r) for r in reps]
    ex = exempt_reps(p, known)
    nonex = [r for r in reps_t if r not in ex]
    out = np.zeros((4, K), dtype=np.int64)
    surv = []
    tot = 0
    for r in nonex:
        n = _lift7(np.array(r, dtype=np.int64), p, 4, out)
        if n:
            tot += n
            surv.append((r, [list(map(int, out[j])) for j in range(min(n, 2))]))
    return dict(p=p, level1_orbits=len(reps_t), exempt_present=sorted(ex & set(reps_t)),
                nonexempt=len(nonex), lifted_survivors=tot, gate=(tot == 0),
                examples=surv[:3], secs=round(time.perf_counter() - t0, 2))


if __name__ == '__main__':
    pmin, pmax = int(sys.argv[1]), int(sys.argv[2])
    rows, mass = [], 0.0
    outp = sys.argv[3] if len(sys.argv) > 3 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), f'tight7_gates_{pmin}_{pmax}.json')
    for p in range(max(pmin, 11), pmax + 1):
        if any(p % q == 0 for q in range(2, int(p ** .5) + 1)):
            continue
        r = tight_gate7(p)
        mass += math.log(p) if r['gate'] else 0.0
        r['gate_mass_so_far'] = round(mass, 2)
        rows.append(r)
        print(json.dumps(r), flush=True)
        with open(outp, 'w') as f:
            json.dump(rows, f, indent=1)

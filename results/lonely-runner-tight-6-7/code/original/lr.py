"""Own implementation of the Rosenfeld / Trakulthongchai / Sungkawichai-Trakulthongchai
finite check for the lonely runner conjecture (LRC).

Convention: k = number of moving (nonzero-speed) runners = number of speeds;
"k+1 runners" in the literature's titles.  LRC(k): for distinct positive integers
v_1..v_k there is t with ||t v_i|| >= 1/(k+1) for all i.

Gate test for a prime p (Lemma 5 of Trakulthongchai 2511.22427, Lemma 6 of Rosenfeld):
  if the improper set I(k, k+1, p) is empty, every counterexample has p | v_1...v_k.

Pipeline here:
  level 1 : enumerate unit/sign/permutation orbit representatives of k-tuples of
            nonzero residues mod p whose "cover sets" cover all times 1..(p-1)/2.
            v covers t  iff  (k+1) * dist(t*v mod p, 0) < p.
            Representatives contain the class 1.  Search = exact set-cover enumeration,
            branching on the uncovered time with fewest live candidates.
            A cover of size s < k is a leaf with k-s free "padding" coordinates.
  lifting : level l -> c*l (c prime factor chain of k+1). Coordinate = (class mod p,
            residue mod l). A level-L tuple is improper iff
              (a) no i with gcd(x_j : j != i, L p) > 1, and
              (b) no t in [1, Lp-1] with dist(t x_i mod Lp, 0) >= p for all i.
            (Lemma 7 of 2511.22427: improper at level m implies its projection to
             level l | m is improper.)
  gate closed iff nothing survives at level k+1.
"""
import math
import numpy as np
from numba import njit

# ---------------------------------------------------------------- level one

def cover_masks(k, p):
    """bool matrix M[v, t] for v,t in 1..H (index 0 unused): v covers t."""
    H = (p - 1) // 2
    v = np.arange(H + 1)[:, None]
    t = np.arange(H + 1)[None, :]
    r = (v * t) % p
    d = np.minimum(r, p - r)
    M = (k + 1) * d < p
    M[0, :] = False
    M[:, 0] = False
    return M


def pack(M):
    """pack rows of bool matrix (cols 1..H) into uint64 words."""
    H = M.shape[1] - 1
    W = (H + 63) // 64
    out = np.zeros((M.shape[0], W), dtype=np.uint64)
    for t in range(1, H + 1):
        w, b = divmod(t - 1, 64)
        out[:, w] |= (M[:, t].astype(np.uint64) << np.uint64(b))
    return out


@njit(cache=True)
def _popc(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return int((x * np.uint64(0x0101010101010101)) >> np.uint64(56))


@njit(cache=True)
def level_one(k, H, cov, covT, csize, out, max_out, use_bound=False):
    """Enumerate covering sets S (containing class 1, |S|<=k) of times 1..H.

    cov[v]   : packed cover set of class v  (v = 1..H)
    covT[t]  : list of classes covering t, as packed bitset over classes (bit v-1)
    Returns (n_leaves, n_nodes); leaves written to out[:, :k] (0 = free padding slot).
    Exact enumeration: each leaf set appears once (standard exclusion branching).
    """
    W = cov.shape[1]
    WC = covT.shape[1]
    # stack-based DFS
    depth_max = k + 1
    unc = np.zeros((depth_max, W), dtype=np.uint64)      # uncovered times
    alive = np.zeros((depth_max, WC), dtype=np.uint64)   # classes still allowed
    chosen = np.zeros(depth_max, dtype=np.int64)
    # candidates at each depth
    cand = np.zeros((depth_max, H), dtype=np.int64)
    ncand = np.zeros(depth_max, dtype=np.int64)
    pos = np.zeros(depth_max, dtype=np.int64)
    topr = np.zeros(k + 1, dtype=np.int64)
    n_leaves = 0
    n_nodes = 0
    # root: choose class 1
    for w in range(W):
        unc[1, w] = ~cov[1, w]
    # mask off bits beyond H
    for t in range(H, W * 64):
        unc[1, t // 64] &= ~(np.uint64(1) << np.uint64(t % 64))
    for w in range(WC):
        alive[1, w] = ~np.uint64(0)
    alive[1, 0] &= ~np.uint64(1)  # class 1 used
    for v in range(H, WC * 64):
        alive[1, v // 64] &= ~(np.uint64(1) << np.uint64(v % 64))
    chosen[1] = 1
    d = 1  # number chosen
    ncand[d] = -1
    while d >= 1:
        if ncand[d] == -1:
            n_nodes += 1
            # node expansion
            nunc = 0
            for w in range(W):
                nunc += _popc(unc[d, w])
            if nunc == 0:
                if n_leaves < max_out:
                    for i in range(k):
                        out[n_leaves, i] = chosen[i + 1] if i + 1 <= d else 0
                n_leaves += 1
                d -= 1
                continue
            r = k - d
            if r == 0 or nunc > r * csize:
                d -= 1
                continue
            if use_bound:
                # top-r bound: sum of the r largest |cov(v) & unc| over alive v
                for j in range(r):
                    topr[j] = 0
                for ww in range(WC):
                    x = alive[d, ww]
                    while x != 0:
                        low = x & (~x + np.uint64(1))
                        b = _popc(low - np.uint64(1))
                        v = ww * 64 + b + 1
                        g = 0
                        for w in range(W):
                            g += _popc(cov[v, w] & unc[d, w])
                        # insert into descending topr
                        if g > topr[r - 1]:
                            j = r - 1
                            while j > 0 and topr[j - 1] < g:
                                topr[j] = topr[j - 1]
                                j -= 1
                            topr[j] = g
                        x ^= low
                s_ = 0
                for j in range(r):
                    s_ += topr[j]
                if s_ < nunc:
                    d -= 1
                    continue
            # pick uncovered time with fewest alive candidates
            best_t = -1
            best_c = 1 << 30
            for w in range(W):
                x = unc[d, w]
                while x != 0:
                    low = x & (~x + np.uint64(1))
                    b = _popc(low - np.uint64(1))
                    t = w * 64 + b  # time index (time t+1)
                    c = 0
                    for ww in range(WC):
                        c += _popc(covT[t, ww] & alive[d, ww])
                    if c < best_c:
                        best_c = c
                        best_t = t
                    x ^= low
            if best_c == 0:
                d -= 1
                continue
            n = 0
            for ww in range(WC):
                x = covT[best_t, ww] & alive[d, ww]
                while x != 0:
                    low = x & (~x + np.uint64(1))
                    b = _popc(low - np.uint64(1))
                    cand[d, n] = ww * 64 + b + 1
                    n += 1
                    x ^= low
            ncand[d] = n
            pos[d] = 0
        if pos[d] >= ncand[d]:
            ncand[d] = -1
            d -= 1
            continue
        c = cand[d, pos[d]]
        pos[d] += 1
        # child: choose c; then exclude c from this node's alive set for later siblings
        nd = d + 1
        for w in range(W):
            unc[nd, w] = unc[d, w] & ~cov[c, w]
        for w in range(WC):
            alive[nd, w] = alive[d, w]
        vb = c - 1
        alive[nd, vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))
        alive[d, vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))
        chosen[nd] = c
        ncand[nd] = -1
        d = nd
    return n_leaves, n_nodes


def level_one_run(k, p, max_out=2_000_000, use_bound=False):
    H = (p - 1) // 2
    M = cover_masks(k, p)
    cov = pack(M)                      # rows v: times
    covT = pack(np.ascontiguousarray(M.T))  # rows t: classes
    covT = covT[1:]                    # index by t-1
    csize = int(M[1].sum())
    out = np.zeros((max_out, k), dtype=np.int64)
    nl, nn = level_one(k, H, cov, covT, csize, out, max_out, use_bound)
    return out[:min(nl, max_out)], nl, nn


# ---------------------------------------------------------------- lifting

@njit(cache=True)
def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=True)
def improper(x, L, p, k):
    """x: residues mod L*p (none divisible by p). True iff improper at level L."""
    N = L * p
    if L > 1:
        for i in range(k):
            g = N
            for j in range(k):
                if j != i:
                    g = _gcd(g, x[j])
            if g > 1:
                return False
    for t in range(1, N):
        ok = True
        for i in range(k):
            r = (t * x[i]) % N
            if r > N - r:
                r = N - r
            if r < p:
                ok = False
                break
        if ok:
            return False
    return True


@njit(cache=True)
def lift(classes, resid, L, c, p, k, out_c, out_r, max_out):
    """Lift level-L survivors (classes mod p, residues mod L) to level c*L.
    Returns number of improper lifts (written to out_*)."""
    n = 0
    M = c * L
    x = np.zeros(k, dtype=np.int64)
    digs = np.zeros(k, dtype=np.int64)
    total = c ** k
    for s in range(classes.shape[0]):
        for code in range(total):
            cc = code
            for i in range(k):
                digs[i] = cc % c
                cc //= c
            for i in range(k):
                r = resid[s, i] + digs[i] * L       # residue mod M
                # CRT: x = class mod p, r mod M
                v = classes[s, i]
                # find x in [0, M p) with x%p==v, x%M==r
                # x = v + p*j, j in [0,M): need (v+p*j)%M == r
                for j in range(M):
                    if (v + p * j) % M == r:
                        x[i] = v + p * j
                        break
            if improper(x, M, p, k):
                if n < max_out:
                    for i in range(k):
                        out_c[n, i] = classes[s, i]
                        out_r[n, i] = resid[s, i] + digs[i] * L
                n += 1
    return n


def expand_padding(leaves, k, H):
    """Leaves with 0 (free) slots -> all completions with any classes 1..H (sorted dedupe)."""
    rows = set()
    for row in leaves:
        fixed = [int(v) for v in row if v != 0]
        free = k - len(fixed)
        if free == 0:
            rows.add(tuple(sorted(fixed)))
            continue
        # combinations with replacement of free extra classes
        from itertools import combinations_with_replacement
        for extra in combinations_with_replacement(range(1, H + 1), free):
            rows.add(tuple(sorted(fixed + list(extra))))
    return np.array(sorted(rows), dtype=np.int64).reshape(-1, k)


def factor_chain(n):
    f, d = [], 2
    while n > 1:
        while n % d == 0:
            f.append(d)
            n //= d
        d += 1
    return f


def gate(k, p, verbose=False):
    """Return (closed, stats). closed=True means I(k,k+1,p) empty -> p | product."""
    import time
    stats = {}
    t0 = time.perf_counter()
    leaves, nl, nn = level_one_run(k, p)
    H = (p - 1) // 2
    stats['level1_leaves'] = nl
    stats['level1_nodes'] = nn
    stats['t_level1'] = time.perf_counter() - t0
    t1 = time.perf_counter()
    cls = expand_padding(leaves, k, H)
    stats['level1_tuples'] = len(cls)
    res = np.zeros_like(cls)
    L = 1
    surv = [len(cls)]
    for c in factor_chain(k + 1):
        cap = max(1, len(cls)) * (c ** k)
        cap = min(cap, 2_000_000)
        oc = np.zeros((cap, k), dtype=np.int64)
        orr = np.zeros((cap, k), dtype=np.int64)
        n = lift(cls, res, L, c, p, k, oc, orr, cap)
        if n > cap:
            surv.append(f'>{cap}')
            stats['survivors'] = surv
            stats['t_lift'] = time.perf_counter() - t1
            return False, stats
        cls, res = oc[:n].copy(), orr[:n].copy()
        L *= c
        surv.append(n)
        if n == 0:
            break
    stats['survivors'] = surv
    stats['t_lift'] = time.perf_counter() - t1
    return len(cls) == 0, stats


# ---------------------------------------------------------------- streaming lift (no storage)
from numba import prange


@njit(cache=True)
def _crt(v, r, M, p):
    for j in range(M):
        x = v + p * j
        if x % M == r:
            return x
    return -1


@njit(cache=True)
def _stream_one(cls, chain, p, k, counts, witness_first):
    """DFS over the lift chain for one level-1 representative. counts[l] += survivors at
    chain level l. Returns number of final survivors."""
    nlev = chain.shape[0]
    # levels: L_0=1, L_{j+1} = L_j*chain[j]
    Ls = np.ones(nlev + 1, dtype=np.int64)
    for j in range(nlev):
        Ls[j + 1] = Ls[j] * chain[j]
    res = np.zeros((nlev + 1, k), dtype=np.int64)
    code = np.zeros(nlev + 1, dtype=np.int64)
    x = np.zeros(k, dtype=np.int64)
    final = 0
    lev = 1
    code[1] = 0
    while lev >= 1:
        c = chain[lev - 1]
        tot = c ** k
        if code[lev] >= tot:
            lev -= 1
            continue
        cc = code[lev]
        code[lev] += 1
        Lp = Ls[lev - 1]
        M = Ls[lev]
        for i in range(k):
            d = cc % c
            cc //= c
            res[lev, i] = res[lev - 1, i] + d * Lp
            x[i] = _crt(cls[i], res[lev, i], M, p)
        if improper(x, M, p, k):
            counts[lev] += 1
            if lev == nlev:
                final += 1
            else:
                lev += 1
                code[lev] = 0
    return final


@njit(parallel=True, cache=True)
def stream_lift(classes, chain, p, k):
    n = classes.shape[0]
    nlev = chain.shape[0]
    allc = np.zeros((n, nlev + 1), dtype=np.int64)
    for s in prange(n):
        cnt = np.zeros(nlev + 1, dtype=np.int64)
        _stream_one(classes[s], chain, p, k, cnt, 0)
        for j in range(nlev + 1):
            allc[s, j] = cnt[j]
    return allc.sum(axis=0)


def gate_stream(k, p, chain=None):
    import time
    t0 = time.perf_counter()
    leaves, nl, nn = level_one_run(k, p)
    H = (p - 1) // 2
    t_l1 = time.perf_counter() - t0
    cls = expand_padding(leaves, k, H)
    ch = np.array(chain if chain else factor_chain(k + 1), dtype=np.int64)
    t1 = time.perf_counter()
    counts = stream_lift(cls, ch, p, k)
    counts[0] = len(cls)
    st = dict(k=k, p=p, level1_leaves=int(nl), level1_nodes=int(nn), level1_tuples=len(cls),
              chain=[int(c) for c in ch], survivors=[int(c) for c in counts],
              t_level1=t_l1, t_lift=time.perf_counter() - t1)
    return counts[-1] == 0, st


# ---------------------------------------------------------------- Knuth tree-size estimator
@njit(cache=True)
def _expand(k, H, d, unc, alive, cov, covT, csize, topr, cand):
    """Evaluate a node of the level-one tree. Returns (status, ncand):
    status 1 = leaf (cover found), 0 = pruned, 2 = internal with ncand children in cand."""
    W = cov.shape[1]
    WC = covT.shape[1]
    nunc = 0
    for w in range(W):
        nunc += _popc(unc[w])
    if nunc == 0:
        return 1, 0
    r = k - d
    if r == 0 or nunc > r * csize:
        return 0, 0
    for j in range(r):
        topr[j] = 0
    for ww in range(WC):
        x = alive[ww]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            v = ww * 64 + b + 1
            g = 0
            for w in range(W):
                g += _popc(cov[v, w] & unc[w])
            if g > topr[r - 1]:
                j = r - 1
                while j > 0 and topr[j - 1] < g:
                    topr[j] = topr[j - 1]
                    j -= 1
                topr[j] = g
            x ^= low
    s_ = 0
    for j in range(r):
        s_ += topr[j]
    if s_ < nunc:
        return 0, 0
    best_t = -1
    best_c = 1 << 30
    for w in range(W):
        x = unc[w]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            t = w * 64 + b
            c = 0
            for ww in range(WC):
                c += _popc(covT[t, ww] & alive[ww])
            if c < best_c:
                best_c = c
                best_t = t
            x ^= low
    if best_c == 0:
        return 0, 0
    n = 0
    for ww in range(WC):
        x = covT[best_t, ww] & alive[ww]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            cand[n] = ww * 64 + b + 1
            n += 1
            x ^= low
    return 2, n


@njit(parallel=True, cache=True)
def knuth(k, H, cov, covT, csize, nsamp, seed):
    """Knuth (1975) random-path estimates of (#nodes, #leaves) of the level-one tree.
    Returns arrays of per-sample estimates (unbiased each)."""
    W = cov.shape[1]
    WC = covT.shape[1]
    est_nodes = np.zeros(nsamp)
    est_leaves = np.zeros(nsamp)
    for s in prange(nsamp):
        np.random.seed(seed + s)
        unc = np.zeros(W, dtype=np.uint64)
        alive = np.zeros(WC, dtype=np.uint64)
        topr = np.zeros(k + 1, dtype=np.int64)
        cand = np.zeros(H, dtype=np.int64)
        for w in range(W):
            unc[w] = ~cov[1, w]
        for t in range(H, W * 64):
            unc[t // 64] &= ~(np.uint64(1) << np.uint64(t % 64))
        for w in range(WC):
            alive[w] = ~np.uint64(0)
        alive[0] &= ~np.uint64(1)
        for v in range(H, WC * 64):
            alive[v // 64] &= ~(np.uint64(1) << np.uint64(v % 64))
        d = 1
        weight = 1.0
        tot = 0.0
        while True:
            tot += weight
            st, n = _expand(k, H, d, unc, alive, cov, covT, csize, topr, cand)
            if st == 1:
                est_leaves[s] = weight
                break
            if st == 0:
                break
            i = np.random.randint(0, n)
            for j in range(i + 1):   # siblings before i are excluded, i is chosen
                vb = cand[j] - 1
                alive[vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))
            c = cand[i]
            for w in range(W):
                unc[w] &= ~cov[c, w]
            weight *= n
            d += 1
        est_nodes[s] = tot
    return est_nodes, est_leaves


def knuth_run(k, p, nsamp=20000, seed=1):
    H = (p - 1) // 2
    M = cover_masks(k, p)
    cov = pack(M)
    covT = pack(np.ascontiguousarray(M.T))[1:]
    csize = int(M[1].sum())
    return knuth(k, H, cov, covT, csize, nsamp, seed)


# ================================================================ v2 (2026-10-21)
# State-of-the-art level one, reimplemented from the papers' descriptions:
#  * Allikvere 2609.02604 v2, Sec. 4: two branches + private-time lemma.
#      n = (p-1)/2, m = floor((p-1)/(k+1)) (= size of every cover set),
#      q_s = ceil(max(0, 2n - s m) / s).  A cover by <= s classes has a class with
#      >= q_s private times (covered by no other class of the cover).
#    Branch A: k distinct classes, irredundant (every class keeps a private time);
#              the class with >= q_k private times is moved to 1 by a unit.
#    Branch B: support <= k-1 classes, or k distinct but redundant.  Then there is a
#              covering S' (|S'| <= k-1, S' = support or support minus a removable class)
#              whose class c with >= q_{k-1} private times has multiplicity <= 2 in M;
#              M = M' + x with M' a (k-1)-multiset on S', mult_{M'}(c) = 1.  c -> 1.
#              The irredundant subcover S00 of S' still contains c with >= q_{k-1}
#              private times, so the DFS enumerates irredundant covers with |S| <= k-1
#              and then pads: M = S + E + x, E (k-1-|S| classes != 1) keeping
#              priv(1) >= q_{k-1}, x arbitrary.
#    In the DFS both "every chosen class keeps a private time" and "class 1 keeps >= q
#    private times" are monotone (private sets only shrink), so they prune, and a class
#    whose addition would violate them is removed from the alive set of the subtree.
#  * Sungkawichai-Trakulthongchai 2604.23906 Sec. 5.1: unit-orbit equivalence
#    (permutation, sign, unit scaling); representatives 1 = v1 <= ... <= vk <= (p-1)/2.
#    canon() below picks the lexicographically least such representative, so the
#    output of both branches is deduplicated to exactly one tuple per orbit.
#  * Lifting (S-T Def. 2.1): level l, times t in (1/(lp))Z, threshold ||t v|| >= 1/(k+1)
#    at EVERY level, i.e. (k+1) * dist(t x mod lp) >= l p; gcd taken with l.
#    (The v1 lift above used dist >= p at every level: sound but much weaker below
#    level k+1.)  Symmetry: units w = 1 mod p, w = 1 mod L act on lifts of a level-L
#    tuple; only lifts minimal under that kernel are kept.

def q_bound(k, p, s):
    n = (p - 1) // 2
    m = (p - 1) // (k + 1)
    return -(-max(0, 2 * n - s * m) // s)


@njit(cache=True)
def _eval_node(kmax, d, q, exact, H, unc, alive, once, chosen, cov, covT, csize, topr,
               cand, kbuf):
    """Evaluate a node (d classes chosen: chosen[1..d], chosen[1] == 1).
    alive is modified in place (classes that would break the private-time conditions
    are removed for the whole subtree).
    Returns (status, ncand): 1 = leaf, 0 = pruned, 2 = internal, 3 = cover of size < kmax
    in exact mode (belongs to the other branch)."""
    W = cov.shape[1]
    WC = covT.shape[1]
    # private-time conditions
    p1 = 0
    for w in range(W):
        p1 += _popc(cov[1, w] & once[w])
    if p1 < q:
        return 0, 0
    for i in range(2, d + 1):
        u = chosen[i]
        g = 0
        for w in range(W):
            g += _popc(cov[u, w] & once[w])
        if g == 0:
            return 0, 0
    nunc = 0
    for w in range(W):
        nunc += _popc(unc[w])
    if nunc == 0:
        if exact and d < kmax:
            return 3, 0
        return 1, 0
    r = kmax - d
    if r == 0 or nunc > r * csize:
        return 0, 0
    # remove classes that cover all private times of some chosen u != 1
    for i in range(2, d + 1):
        u = chosen[i]
        first = True
        anyk = np.uint64(0)
        for w in range(W):
            x = cov[u, w] & once[w]
            while x != 0:
                low = x & (~x + np.uint64(1))
                b = _popc(low - np.uint64(1))
                t = w * 64 + b
                anyk = np.uint64(0)
                if first:
                    first = False
                    for ww in range(WC):
                        kbuf[ww] = covT[t, ww] & alive[ww]
                        anyk |= kbuf[ww]
                else:
                    for ww in range(WC):
                        kbuf[ww] &= covT[t, ww]
                        anyk |= kbuf[ww]
                if anyk == 0:
                    break
                x ^= low
            if anyk == 0:
                break
        if anyk != 0:
            for ww in range(WC):
                alive[ww] &= ~kbuf[ww]
    slack1 = p1 - q
    for j in range(r):
        topr[j] = 0
    for ww in range(WC):
        x = alive[ww]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            v = ww * 64 + b + 1
            x ^= low
            if slack1 < csize:
                h = 0
                for w in range(W):
                    h += _popc(cov[v, w] & cov[1, w] & once[w])
                if h > slack1:
                    alive[ww] &= ~low
                    continue
            g = 0
            for w in range(W):
                g += _popc(cov[v, w] & unc[w])
            if g > topr[r - 1]:
                j = r - 1
                while j > 0 and topr[j - 1] < g:
                    topr[j] = topr[j - 1]
                    j -= 1
                topr[j] = g
    s_ = 0
    for j in range(r):
        s_ += topr[j]
    if s_ < nunc:
        return 0, 0
    best_t = -1
    best_c = 1 << 30
    for w in range(W):
        x = unc[w]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            t = w * 64 + b
            c = 0
            for ww in range(WC):
                c += _popc(covT[t, ww] & alive[ww])
            if c < best_c:
                best_c = c
                best_t = t
                if c == 0:
                    break
            x ^= low
        if best_c == 0:
            break
    if best_c == 0:
        return 0, 0
    n = 0
    for ww in range(WC):
        x = covT[best_t, ww] & alive[ww]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            cand[n] = ww * 64 + b + 1
            n += 1
            x ^= low
    return 2, n


@njit(cache=True)
def _root(H, cov, W, WC, unc, alive, once, twice):
    for w in range(W):
        unc[w] = ~cov[1, w]
        once[w] = cov[1, w]
        twice[w] = np.uint64(0)
    for t in range(H, W * 64):
        unc[t // 64] &= ~(np.uint64(1) << np.uint64(t % 64))
    for w in range(WC):
        alive[w] = ~np.uint64(0)
    alive[0] &= ~np.uint64(1)
    for v in range(H, WC * 64):
        alive[v // 64] &= ~(np.uint64(1) << np.uint64(v % 64))


@njit(cache=True)
def _child(c, cov, W, WC, unc, alive, once, twice, unc2, alive2, once2, twice2):
    for w in range(W):
        cw = cov[c, w]
        unc2[w] = unc[w] & ~cw
        tw = twice[w] | (once[w] & cw)
        twice2[w] = tw
        once2[w] = (once[w] | cw) & ~tw
    for w in range(WC):
        alive2[w] = alive[w]
    vb = c - 1
    alive2[vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))


@njit(cache=True)
def level_one2(k, H, cov, covT, csize, kmax, q, exact, out, max_out):
    """Exact DFS over irredundant covers S (1 in S, |S| <= kmax, class 1 with >= q private
    times; exact=True: only |S| == kmax).  Leaves -> out rows (0-padded).
    Returns (n_leaves, n_nodes)."""
    W = cov.shape[1]
    WC = covT.shape[1]
    D = kmax + 2
    unc = np.zeros((D, W), dtype=np.uint64)
    once = np.zeros((D, W), dtype=np.uint64)
    twice = np.zeros((D, W), dtype=np.uint64)
    alive = np.zeros((D, WC), dtype=np.uint64)
    chosen = np.zeros(D, dtype=np.int64)
    cand = np.zeros((D, H + 1), dtype=np.int64)
    kbuf = np.zeros(WC, dtype=np.uint64)
    ncand = np.full(D, -1, dtype=np.int64)
    pos = np.zeros(D, dtype=np.int64)
    topr = np.zeros(k + 2, dtype=np.int64)
    _root(H, cov, W, WC, unc[1], alive[1], once[1], twice[1])
    chosen[1] = 1
    n_leaves = 0
    n_nodes = 0
    d = 1
    while d >= 1:
        if ncand[d] == -1:
            n_nodes += 1
            st, n = _eval_node(kmax, d, q, exact, H, unc[d], alive[d], once[d], chosen,
                               cov, covT, csize, topr, cand[d], kbuf)
            if st == 1:
                if n_leaves < max_out:
                    for i in range(out.shape[1]):
                        out[n_leaves, i] = chosen[i + 1] if i + 1 <= d else 0
                n_leaves += 1
                d -= 1
                continue
            if st != 2:
                d -= 1
                continue
            ncand[d] = n
            pos[d] = 0
        if pos[d] >= ncand[d]:
            ncand[d] = -1
            d -= 1
            continue
        c = cand[d, pos[d]]
        pos[d] += 1
        nd = d + 1
        _child(c, cov, W, WC, unc[d], alive[d], once[d], twice[d],
               unc[nd], alive[nd], once[nd], twice[nd])
        vb = c - 1
        alive[d, vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))
        chosen[nd] = c
        ncand[nd] = -1
        d = nd
    return n_leaves, n_nodes


@njit(parallel=True, cache=True)
def knuth2(k, H, cov, covT, csize, kmax, q, exact, nsamp, seed):
    """Knuth random-path estimates (nodes, leaves) for the level_one2 tree."""
    W = cov.shape[1]
    WC = covT.shape[1]
    est_nodes = np.zeros(nsamp)
    est_leaves = np.zeros(nsamp)
    for s in prange(nsamp):
        np.random.seed(seed + s)
        unc = np.zeros(W, dtype=np.uint64)
        once = np.zeros(W, dtype=np.uint64)
        twice = np.zeros(W, dtype=np.uint64)
        alive = np.zeros(WC, dtype=np.uint64)
        unc2 = np.zeros(W, dtype=np.uint64)
        once2 = np.zeros(W, dtype=np.uint64)
        twice2 = np.zeros(W, dtype=np.uint64)
        alive2 = np.zeros(WC, dtype=np.uint64)
        kbuf = np.zeros(WC, dtype=np.uint64)
        chosen = np.zeros(kmax + 2, dtype=np.int64)
        topr = np.zeros(k + 2, dtype=np.int64)
        cand = np.zeros(H + 1, dtype=np.int64)
        _root(H, cov, W, WC, unc, alive, once, twice)
        chosen[1] = 1
        d = 1
        weight = 1.0
        tot = 0.0
        while True:
            tot += weight
            st, n = _eval_node(kmax, d, q, exact, H, unc, alive, once, chosen,
                               cov, covT, csize, topr, cand, kbuf)
            if st == 1:
                est_leaves[s] = weight
                break
            if st != 2:
                break
            i = np.random.randint(0, n)
            for j in range(i):   # earlier siblings are excluded in the chosen subtree
                vb = cand[j] - 1
                alive[vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))
            c = cand[i]
            _child(c, cov, W, WC, unc, alive, once, twice, unc2, alive2, once2, twice2)
            for w in range(W):
                unc[w] = unc2[w]
                once[w] = once2[w]
                twice[w] = twice2[w]
            for w in range(WC):
                alive[w] = alive2[w]
            d += 1
            chosen[d] = c
            weight *= n
        est_nodes[s] = tot
    return est_nodes, est_leaves


def _tables(k, p):
    H = (p - 1) // 2
    M = cover_masks(k, p)
    cov = pack(M)
    covT = pack(np.ascontiguousarray(M.T))[1:]
    csize = int(M[1].sum())
    return H, M, cov, covT, csize


def knuth2_run(k, p, branch, nsamp=20000, seed=1):
    H, M, cov, covT, csize = _tables(k, p)
    if branch == 'A':
        return knuth2(k, H, cov, covT, csize, k, q_bound(k, p, k), True, nsamp, seed)
    return knuth2(k, H, cov, covT, csize, k - 1, q_bound(k, p, k - 1), False, nsamp, seed)


@njit(cache=True)
def canon_rows(rows, p, inv):
    """Unit-orbit normal form: for each row (k classes in 1..H, sorted multiset), the
    lexicographic minimum over units u = c^{-1} (c in the row) of sorted(fold(u * row)).
    Signs are folded, so u and -u give the same tuple."""
    n, k = rows.shape
    out = np.zeros_like(rows)
    y = np.zeros(k, dtype=np.int64)
    best = np.zeros(k, dtype=np.int64)
    for s in range(n):
        have = False
        for a in range(k):
            if a > 0 and rows[s, a] == rows[s, a - 1]:
                continue
            u = inv[rows[s, a]]
            for i in range(k):
                r = (u * rows[s, i]) % p
                y[i] = r if r <= p - r else p - r
            y.sort()
            better = not have
            if have:
                for i in range(k):
                    if y[i] != best[i]:
                        better = y[i] < best[i]
                        break
            if better:
                for i in range(k):
                    best[i] = y[i]
                have = True
        for i in range(k):
            out[s, i] = best[i]
    return out


def canon(rows, p):
    rows = np.asarray(rows, dtype=np.int64)
    if len(rows) == 0:
        return rows.reshape(0, rows.shape[-1] if rows.ndim == 2 else 0)
    rows = np.sort(rows, axis=1)
    inv = np.zeros(p, dtype=np.int64)
    for c in range(1, p):
        inv[c] = pow(c, -1, p)
    return np.unique(canon_rows(rows, p, inv), axis=0)


def _padB(leaves, k, p, q, M):
    """Branch B padding: M = S + E + x (see header)."""
    from itertools import combinations_with_replacement
    H = (p - 1) // 2
    rows = []
    cov1 = M[1]
    for row in leaves:
        S = [int(v) for v in row if v != 0]
        e = k - 1 - len(S)
        others = [v for v in S if v != 1]
        if e == 0:
            Es = [()]
        else:
            covo = np.zeros(H + 1, dtype=np.int64)
            for v in others:
                covo += M[v]
            Es = []
            for E in combinations_with_replacement(range(2, H + 1), e):
                cc = covo.copy()
                for v in set(E) - set(others):
                    cc += M[v]
                if int((cov1 & (cc == 0)).sum()) >= q:
                    Es.append(E)
        if not Es:
            continue
        B = np.array([S + list(E) for E in Es], dtype=np.int64)
        X = np.repeat(B, H, axis=0)
        xs = np.tile(np.arange(1, H + 1, dtype=np.int64), len(B))[:, None]
        rows.append(np.hstack([X, xs]))
    if not rows:
        return np.zeros((0, k), dtype=np.int64)
    return np.vstack(rows)


def level_one_reps(k, p, max_out=5_000_000):
    """One representative (normal form) per unit orbit of I(k, p, 1).
    Returns (reps, stats)."""
    import time
    t0 = time.perf_counter()
    H, M, cov, covT, csize = _tables(k, p)
    qA, qB = q_bound(k, p, k), q_bound(k, p, k - 1)
    outA = np.zeros((max_out, k), dtype=np.int64)
    nA, nnA = level_one2(k, H, cov, covT, csize, k, qA, True, outA, max_out)
    outB = np.zeros((max_out, k), dtype=np.int64)
    nB, nnB = level_one2(k, H, cov, covT, csize, k - 1, qB, False, outB, max_out)
    if nA > max_out or nB > max_out:
        raise RuntimeError('max_out too small')
    t1 = time.perf_counter()
    padded = _padB(outB[:nB, :k - 1], k, p, qB, M)
    allrows = np.vstack([outA[:nA], padded]) if len(padded) else outA[:nA]
    reps = canon(allrows, p) if len(allrows) else np.zeros((0, k), dtype=np.int64)
    st = dict(k=k, p=p, qA=qA, qB=qB, nodesA=int(nnA), leavesA=int(nA), nodesB=int(nnB),
              leavesB=int(nB), padded=int(len(padded)), reps=int(len(reps)),
              t_search=t1 - t0, t_pad=time.perf_counter() - t1)
    return reps, st


# ---------------------------------------------------------------- v2 lift
@njit(cache=True)
def improper2(x, L, p, k, c_last, wit, nwit):
    """S-T Def. 2.1 at level L: x residues mod L*p, none divisible by p.  True iff improper.
    Times t that are multiples of c_last are skipped: they are times of the parent level
    L/c_last, where the parent is known to be improper.  t and N-t are equivalent.
    wit[0:nwit[0]] is a move-to-front cache of recent witnesses (per level, per thread)."""
    N = L * p
    if L > 1:
        for i in range(k):
            g = L
            for j in range(k):
                if j != i:
                    g = _gcd(g, x[j])
                    if g == 1:
                        break
            if g > 1:
                return False
    K1 = k + 1
    nw = nwit[0]
    for a in range(nw):
        t = wit[a]
        ok = True
        for i in range(k):
            r = (t * x[i]) % N
            if r > N - r:
                r = N - r
            if K1 * r < N:
                ok = False
                break
        if ok:
            for b in range(a, 0, -1):
                wit[b] = wit[b - 1]
            wit[0] = t
            return False
    for t in range(1, N // 2 + 1):
        if t % c_last == 0:
            continue
        ok = True
        for i in range(k):
            r = (t * x[i]) % N
            if r > N - r:
                r = N - r
            if K1 * r < N:
                ok = False
                break
        if ok:
            m = nw if nw < wit.shape[0] else wit.shape[0] - 1
            for b in range(m, 0, -1):
                wit[b] = wit[b - 1]
            wit[0] = t
            if nw < wit.shape[0]:
                nwit[0] = nw + 1
            return False
    return True


@njit(cache=True)
def _kernel_min(res, L, Lp, c, k):
    """True iff res (mod L) is lexicographically minimal among w*res, w = 1 + j*Lp unit mod L."""
    for j in range(1, c):
        w = 1 + j * Lp
        if _gcd(w, L) != 1:
            continue
        for i in range(k):
            a = (w * res[i]) % L
            if a != res[i]:
                if a < res[i]:
                    return False
                break
    return True


@njit(cache=True)
def _stream_one2(cls, chain, p, k, counts, code_lo, code_hi, wit, nwit):
    """DFS of canonical improper lifts of one level-1 representative along chain.
    Level-1 -> level chain[0] digit codes restricted to [code_lo, code_hi).
    counts[l] += canonical improper tuples at chain level l.  Returns final count."""
    nlev = chain.shape[0]
    Ls = np.ones(nlev + 1, dtype=np.int64)
    for j in range(nlev):
        Ls[j + 1] = Ls[j] * chain[j]
    pinv = np.zeros(nlev + 1, dtype=np.int64)
    for j in range(1, nlev + 1):
        L = Ls[j]
        for a in range(1, L + 1):
            if (p * a) % L == 1 % L:
                pinv[j] = a
                break
    res = np.zeros((nlev + 1, k), dtype=np.int64)
    code = np.zeros(nlev + 1, dtype=np.int64)
    hi = np.zeros(nlev + 1, dtype=np.int64)
    x = np.zeros(k, dtype=np.int64)
    final = 0
    lev = 1
    code[1] = code_lo
    hi[1] = code_hi
    while lev >= 1:
        c = chain[lev - 1]
        if code[lev] >= hi[lev]:
            lev -= 1
            continue
        cc = code[lev]
        code[lev] += 1
        Lp = Ls[lev - 1]
        M = Ls[lev]
        for i in range(k):
            dg = cc % c
            cc //= c
            res[lev, i] = res[lev - 1, i] + dg * Lp
        if not _kernel_min(res[lev], M, Lp, c, k):
            continue
        for i in range(k):
            x[i] = cls[i] + p * ((((res[lev, i] - cls[i]) % M) * pinv[lev]) % M)
        if improper2(x, M, p, k, c, wit[lev], nwit[lev]):
            counts[lev] += 1
            if lev == nlev:
                final += 1
            else:
                lev += 1
                code[lev] = 0
                hi[lev] = chain[lev - 1] ** k
    return final


@njit(parallel=True, cache=True)
def stream_lift2(classes, chain, p, k):
    """Per-representative survivor counts at each chain level (canonical lifts)."""
    n = classes.shape[0]
    nlev = chain.shape[0]
    allc = np.zeros((n, nlev + 1), dtype=np.int64)
    for s in prange(n):
        cnt = np.zeros(nlev + 1, dtype=np.int64)
        wit = np.zeros((nlev + 1, 16), dtype=np.int64)
        nwit = np.zeros((nlev + 1, 1), dtype=np.int64)
        _stream_one2(classes[s], chain, p, k, cnt, 0, chain[0] ** k, wit, nwit)
        for j in range(nlev + 1):
            allc[s, j] = cnt[j]
        allc[s, 0] = 1
    return allc


@njit(parallel=True, cache=True)
def stream_lift2_one(cls, chain, p, k, nchunks):
    """One representative, parallel over chunks of the first-level digit codes."""
    nlev = chain.shape[0]
    tot = chain[0] ** k
    allc = np.zeros((nchunks, nlev + 1), dtype=np.int64)
    for s in prange(nchunks):
        lo = tot * s // nchunks
        hi = tot * (s + 1) // nchunks
        cnt = np.zeros(nlev + 1, dtype=np.int64)
        wit = np.zeros((nlev + 1, 16), dtype=np.int64)
        nwit = np.zeros((nlev + 1, 1), dtype=np.int64)
        _stream_one2(cls, chain, p, k, cnt, lo, hi, wit, nwit)
        for j in range(nlev + 1):
            allc[s, j] = cnt[j]
    out = allc.sum(axis=0)
    out[0] = 1
    return out


def gate2(k, p, chain=None):
    """v2 gate: level-one orbit representatives, then canonical streaming lift.
    closed=True  <=>  every representative has no improper lift at the last chain level."""
    import time
    t0 = time.perf_counter()
    reps, st = level_one_reps(k, p)
    t1 = time.perf_counter()
    ch = np.array(chain if chain else factor_chain(k + 1), dtype=np.int64)
    if len(reps):
        per = stream_lift2(reps, ch, p, k)
        counts = per.sum(axis=0)
        surv = reps[per[:, -1] > 0]
    else:
        counts = np.zeros(len(ch) + 1, dtype=np.int64)
        surv = reps
    st.update(chain=[int(c) for c in ch], survivors=[int(c) for c in counts],
              surviving_reps=[list(map(int, r)) for r in surv[:5]],
              t_level1=t1 - t0, t_lift=time.perf_counter() - t1)
    return bool(counts[-1] == 0), st


# ---------------------------------------------------------------- v3 unified level one
# One DFS over covers S (1 in S, |S| <= k), as v1, with:
#  * private-time normalisation (Allikvere Lemma 4.1 applied to the irredundant subcover
#    S00 of supp(M), |S00| = s <= k: some class has >= q_s >= q_k private times; that
#    class is moved to 1).  Prune: priv(1) >= q_k; optional prune: every chosen class
#    keeps >= 1 private time (S00 is irredundant), leaves then are irredundant covers.
#  * children visited in decreasing gain |cov(c) & unc|, with the child-gain cutoff:
#    child c (earlier siblings excluded) needs  gain(c) + top_{r-1}(gains of the other
#    alive classes) >= #uncovered, gains taken at the parent (upper bounds for the child).
#    The bound is monotone along the sorted siblings, so the loop stops at the first
#    failure; such children are never evaluated.
#  * zero-gain classes are dropped from the alive set of the subtree.
#  * leaves with |S| < k are padded with every multiset of k - |S| classes, then reduced
#    to one normal form per unit orbit (canon).
@njit(cache=True)
def _eval3(kmax, d, q, irred, unc, alive, once, chosen, cov, covT, csize,
           hist, cand, gain):
    """Returns (status, ncand, nunc). status 1 leaf, 0 pruned, 2 internal."""
    W = cov.shape[1]
    WC = covT.shape[1]
    if q > 0:
        p1 = 0
        for w in range(W):
            p1 += _popc(cov[1, w] & once[w])
        if p1 < q:
            return 0, 0, 0
    if irred:
        for i in range(2, d + 1):
            u = chosen[i]
            g = 0
            for w in range(W):
                g += _popc(cov[u, w] & once[w])
            if g == 0:
                return 0, 0, 0
    nunc = 0
    for w in range(W):
        nunc += _popc(unc[w])
    if nunc == 0:
        return 1, 0, 0
    r = kmax - d
    if r == 0 or nunc > r * csize:
        return 0, 0, nunc
    # time with fewest alive candidates
    best_t = -1
    best_c = 1 << 30
    for w in range(W):
        x = unc[w]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            t = w * 64 + b
            c = 0
            for ww in range(WC):
                c += _popc(covT[t, ww] & alive[ww])
            if c < best_c:
                best_c = c
                best_t = t
                if c == 0:
                    break
            x ^= low
        if best_c == 0:
            break
    if best_c == 0:
        return 0, 0, nunc
    # gains of all alive classes -> histogram; drop zero-gain classes
    for g in range(csize + 1):
        hist[g] = 0
    for ww in range(WC):
        x = alive[ww]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            v = ww * 64 + b + 1
            x ^= low
            g = 0
            for w in range(W):
                g += _popc(cov[v, w] & unc[w])
            if g == 0:
                alive[ww] &= ~low
            else:
                hist[g] += 1
    # top-r bound
    s_ = 0
    need = r
    for g in range(csize, 0, -1):
        h = hist[g]
        if h >= need:
            s_ += need * g
            need = 0
            break
        s_ += h * g
        need -= h
    if s_ < nunc:
        return 0, 0, nunc
    # candidates sorted by decreasing gain
    n = 0
    for ww in range(WC):
        x = covT[best_t, ww] & alive[ww]
        while x != 0:
            low = x & (~x + np.uint64(1))
            b = _popc(low - np.uint64(1))
            v = ww * 64 + b + 1
            x ^= low
            g = 0
            for w in range(W):
                g += _popc(cov[v, w] & unc[w])
            j = n
            while j > 0 and gain[j - 1] < g:
                gain[j] = gain[j - 1]
                cand[j] = cand[j - 1]
                j -= 1
            gain[j] = g
            cand[j] = v
            n += 1
    return 2, n, nunc


@njit(cache=True)
def _cut_ok(g, r1, hist, csize, nunc):
    """g + top_{r1}(hist with one copy of gain g removed) >= nunc ?"""
    hist[g] -= 1
    s_ = g
    need = r1
    if need > 0:
        for gg in range(csize, 0, -1):
            h = hist[gg]
            if h >= need:
                s_ += need * gg
                need = 0
                break
            s_ += h * gg
            need -= h
    hist[g] += 1
    return s_ >= nunc


@njit(cache=True)
def level_one3(k, H, cov, covT, csize, kmax, q, irred, out, max_out, max_nodes=-1):
    """Leaves = covers S (1 in S, |S| <= kmax) passing the prunes; rows 0-padded.
    Returns (n_leaves, n_nodes_evaluated)."""
    W = cov.shape[1]
    WC = covT.shape[1]
    D = kmax + 2
    unc = np.zeros((D, W), dtype=np.uint64)
    once = np.zeros((D, W), dtype=np.uint64)
    twice = np.zeros((D, W), dtype=np.uint64)
    alive = np.zeros((D, WC), dtype=np.uint64)
    chosen = np.zeros(D, dtype=np.int64)
    cand = np.zeros((D, H + 1), dtype=np.int64)
    gain = np.zeros((D, H + 1), dtype=np.int64)
    hist = np.zeros((D, csize + 1), dtype=np.int64)
    ncand = np.full(D, -1, dtype=np.int64)
    nuncs = np.zeros(D, dtype=np.int64)
    pos = np.zeros(D, dtype=np.int64)
    _root(H, cov, W, WC, unc[1], alive[1], once[1], twice[1])
    chosen[1] = 1
    n_leaves = 0
    n_nodes = 0
    d = 1
    while d >= 1:
        if ncand[d] == -1:
            n_nodes += 1
            if max_nodes > 0 and n_nodes > max_nodes:
                break
            st, n, nu = _eval3(kmax, d, q, irred, unc[d], alive[d], once[d], chosen,
                               cov, covT, csize, hist[d], cand[d], gain[d])
            if st == 1:
                if n_leaves < max_out:
                    for i in range(out.shape[1]):
                        out[n_leaves, i] = chosen[i + 1] if i + 1 <= d else 0
                n_leaves += 1
                d -= 1
                continue
            if st != 2:
                d -= 1
                continue
            ncand[d] = n
            nuncs[d] = nu
            pos[d] = 0
        if pos[d] >= ncand[d]:
            ncand[d] = -1
            d -= 1
            continue
        j = pos[d]
        c = cand[d, j]
        g = gain[d, j]
        if not _cut_ok(g, kmax - d - 1, hist[d], csize, nuncs[d]):
            ncand[d] = -1
            d -= 1
            continue
        pos[d] += 1
        nd = d + 1
        _child(c, cov, W, WC, unc[d], alive[d], once[d], twice[d],
               unc[nd], alive[nd], once[nd], twice[nd])
        vb = c - 1
        alive[d, vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))
        hist[d, g] -= 1           # c is excluded for the later siblings
        chosen[nd] = c
        ncand[nd] = -1
        d = nd
    return n_leaves, n_nodes


@njit(parallel=True, cache=True)
def knuth3(k, H, cov, covT, csize, kmax, q, irred, nsamp, seed):
    """Knuth estimates of (evaluated nodes, leaves) of the level_one3 tree."""
    W = cov.shape[1]
    WC = covT.shape[1]
    est_nodes = np.zeros(nsamp)
    est_leaves = np.zeros(nsamp)
    for s in prange(nsamp):
        np.random.seed(seed + s)
        unc = np.zeros(W, dtype=np.uint64)
        once = np.zeros(W, dtype=np.uint64)
        twice = np.zeros(W, dtype=np.uint64)
        alive = np.zeros(WC, dtype=np.uint64)
        unc2 = np.zeros(W, dtype=np.uint64)
        once2 = np.zeros(W, dtype=np.uint64)
        twice2 = np.zeros(W, dtype=np.uint64)
        alive2 = np.zeros(WC, dtype=np.uint64)
        chosen = np.zeros(kmax + 2, dtype=np.int64)
        cand = np.zeros(H + 1, dtype=np.int64)
        gain = np.zeros(H + 1, dtype=np.int64)
        hist = np.zeros(csize + 1, dtype=np.int64)
        _root(H, cov, W, WC, unc, alive, once, twice)
        chosen[1] = 1
        d = 1
        weight = 1.0
        tot = 0.0
        while True:
            tot += weight
            st, n, nu = _eval3(kmax, d, q, irred, unc, alive, once, chosen,
                               cov, covT, csize, hist, cand, gain)
            if st == 1:
                est_leaves[s] = weight
                break
            if st != 2:
                break
            # number of children passing the cutoff (sequential exclusion)
            nv = 0
            for j in range(n):
                if not _cut_ok(gain[j], kmax - d - 1, hist, csize, nu):
                    break
                hist[gain[j]] -= 1
                nv += 1
            if nv == 0:
                break
            i = np.random.randint(0, nv)
            for j in range(i):
                vb = cand[j] - 1
                alive[vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))
            c = cand[i]
            _child(c, cov, W, WC, unc, alive, once, twice, unc2, alive2, once2, twice2)
            for w in range(W):
                unc[w] = unc2[w]
                once[w] = once2[w]
                twice[w] = twice2[w]
            for w in range(WC):
                alive[w] = alive2[w]
            d += 1
            chosen[d] = c
            weight *= nv
        est_nodes[s] = tot
    return est_nodes, est_leaves


def knuth3_run(k, p, nsamp=20000, seed=1, q=True, irred=True):
    H, M, cov, covT, csize = _tables(k, p)
    qq = q_bound(k, p, k) if q else 0
    return knuth3(k, H, cov, covT, csize, k, qq, irred, nsamp, seed)


def _pad_all(leaves, k, H):
    """Leaves with s < k classes -> every multiset completion (numpy, chunked)."""
    from itertools import combinations_with_replacement
    rows = []
    full = leaves[(leaves != 0).sum(axis=1) == k]
    if len(full):
        rows.append(full)
    part = leaves[(leaves != 0).sum(axis=1) < k]
    for row in part:
        S = [int(v) for v in row if v != 0]
        e = k - len(S)
        E = np.array(list(combinations_with_replacement(range(1, H + 1), e)), dtype=np.int64)
        rows.append(np.hstack([np.tile(np.array(S, dtype=np.int64), (len(E), 1)), E]))
    if not rows:
        return np.zeros((0, k), dtype=np.int64)
    return np.vstack(rows)


def level_one_reps3(k, p, max_out=5_000_000, q=True, irred=True):
    """v3: one normal-form representative per unit orbit of I(k, p, 1)."""
    import time
    t0 = time.perf_counter()
    H, M, cov, covT, csize = _tables(k, p)
    qq = q_bound(k, p, k) if q else 0
    out = np.zeros((max_out, k), dtype=np.int64)
    nl, nn = level_one3(k, H, cov, covT, csize, k, qq, irred, out, max_out)
    if nl > max_out:
        raise RuntimeError('max_out too small')
    t1 = time.perf_counter()
    rows = _pad_all(out[:nl], k, H)
    reps = canon(rows, p) if len(rows) else np.zeros((0, k), dtype=np.int64)
    st = dict(k=k, p=p, q=qq, irred=irred, nodes=int(nn), leaves=int(nl),
              padded=int(len(rows)), reps=int(len(reps)),
              t_search=t1 - t0, t_pad=time.perf_counter() - t1)
    return reps, st


def gate3(k, p, chain=None):
    """v3 level one (unified search + normal form) followed by the v2 canonical lift."""
    import time
    t0 = time.perf_counter()
    reps, st = level_one_reps3(k, p)
    t1 = time.perf_counter()
    ch = np.array(chain if chain else factor_chain(k + 1), dtype=np.int64)
    if len(reps):
        per = stream_lift2(reps, ch, p, k)
        counts = per.sum(axis=0)
        surv = reps[per[:, -1] > 0]
    else:
        counts = np.zeros(len(ch) + 1, dtype=np.int64)
        surv = reps
    st.update(chain=[int(c) for c in ch], survivors=[int(c) for c in counts],
              n_surviving_reps=int(len(surv)),
              surviving_reps=[list(map(int, r)) for r in surv[:5]],
              t_level1=t1 - t0, t_lift=time.perf_counter() - t1)
    return bool(counts[-1] == 0), st


def tight_rep(k, p):
    return np.array(sorted(min(i % p, p - i % p) for i in range(1, k + 1)), dtype=np.int64)


# ---------------------------------------------------------------- direct level-L cover search
# For one level-1 representative cls (classes mod p) and a level L (gcd(L, p) = 1), an
# improper level-L lift is a choice of residue r_i mod L per coordinate such that the
# items (i, r_i) cover every time t in 1..Lp/2 with t != 0 mod L (multiples of L are
# level-1 times, covered for every choice), where (i, r) covers t iff
# (k+1) dist(t x mod Lp) < Lp with x = CRT(cls_i mod p, r mod L); plus the gcd clause
# (for every prime d | L, at least two coordinates not divisible by d).  Searched like level one
# (fewest-candidate time, exclusion branching, per-coordinate top-r gain bound), with
# no intermediate levels.
def _item_tables(cls, L, p, k):
    N = L * p
    T = np.array([t for t in range(1, N // 2 + 1) if t % L], dtype=np.int64)
    nT = len(T)
    xs = np.zeros((k, L), dtype=np.int64)
    for i in range(k):
        for r in range(L):
            xs[i, r] = next(v for v in range(int(cls[i]), N, p) if v % L == r)
    X = xs.reshape(-1)
    R = (T[None, :] * X[:, None]) % N
    D = np.minimum(R, N - R)
    C = (k + 1) * D < N                       # items x times
    Cp = np.zeros((k * L, nT + 1), dtype=bool)
    Cp[:, 1:] = C
    icov = pack(Cp)
    Ct = np.zeros((nT, k * L + 1), dtype=bool)
    Ct[:, 1:] = C.T
    itT = pack(Ct)
    return T, xs, C, icov, itT


@njit(cache=True)
def _cover_task(k, L, icov, itT, unc0, alive0, assign0, max_sol, sols, nsol, stats, lprimes):
    """Sequential DFS from a given state. stats[0] nodes, stats[1] leaves (partial
    assignments covering all times), stats[2] leaves passing the gcd clause."""
    WT = icov.shape[1]
    WI = itT.shape[1]
    D = k + 2
    unc = np.zeros((D, WT), dtype=np.uint64)
    alive = np.zeros((D, WI), dtype=np.uint64)
    assign = np.zeros((D, k), dtype=np.int64)
    cand = np.zeros((D, k * L), dtype=np.int64)
    gain = np.zeros((D, k * L), dtype=np.int64)
    ncand = np.full(D, -1, dtype=np.int64)
    pos = np.zeros(D, dtype=np.int64)
    cmax = np.zeros(k, dtype=np.int64)
    for w in range(WT):
        unc[0, w] = unc0[w]
    for w in range(WI):
        alive[0, w] = alive0[w]
    for i in range(k):
        assign[0, i] = assign0[i]
    d = 0
    while d >= 0:
        if ncand[d] == -1:
            stats[0] += 1
            nunc = 0
            for w in range(WT):
                nunc += _popc(unc[d, w])
            nfree = 0
            for i in range(k):
                if assign[d, i] < 0:
                    nfree += 1
            if nunc == 0:
                stats[1] += 1
                # gcd clause: for every prime d | L, at least two coordinates not
                # divisible by d (free coordinates are set to residue 1)
                gok = True
                for a in range(lprimes.shape[0]):
                    dd = lprimes[a]
                    nd_ = nfree
                    for i in range(k):
                        if assign[d, i] >= 0 and assign[d, i] % dd != 0:
                            nd_ += 1
                    if nd_ < 2:
                        gok = False
                if gok:
                    if nsol[0] < max_sol:
                        for i in range(k):
                            sols[nsol[0], i] = assign[d, i]
                    nsol[0] += 1
                    stats[2] += 1
                d -= 1
                continue
            if nfree == 0:
                d -= 1
                continue
            for i in range(k):
                cmax[i] = 0
            for ww in range(WI):
                x = alive[d, ww]
                while x != 0:
                    low = x & (~x + np.uint64(1))
                    b = _popc(low - np.uint64(1))
                    it = ww * 64 + b
                    x ^= low
                    g = 0
                    for w in range(WT):
                        g += _popc(icov[it, w] & unc[d, w])
                    i = it // L
                    if g > cmax[i]:
                        cmax[i] = g
            cmax.sort()
            s_ = 0
            for j in range(nfree):
                s_ += cmax[k - 1 - j]
            if s_ < nunc:
                d -= 1
                continue
            best_t = -1
            best_c = 1 << 30
            for w in range(WT):
                x = unc[d, w]
                while x != 0:
                    low = x & (~x + np.uint64(1))
                    b = _popc(low - np.uint64(1))
                    t = w * 64 + b
                    x ^= low
                    c = 0
                    for ww in range(WI):
                        c += _popc(itT[t, ww] & alive[d, ww])
                    if c < best_c:
                        best_c = c
                        best_t = t
                        if c == 0:
                            break
                if best_c == 0:
                    break
            if best_c == 0:
                d -= 1
                continue
            n = 0
            for ww in range(WI):
                x = itT[best_t, ww] & alive[d, ww]
                while x != 0:
                    low = x & (~x + np.uint64(1))
                    b = _popc(low - np.uint64(1))
                    it = ww * 64 + b
                    x ^= low
                    g = 0
                    for w in range(WT):
                        g += _popc(icov[it, w] & unc[d, w])
                    j = n
                    while j > 0 and gain[d, j - 1] < g:
                        gain[d, j] = gain[d, j - 1]
                        cand[d, j] = cand[d, j - 1]
                        j -= 1
                    gain[d, j] = g
                    cand[d, j] = it
                    n += 1
            ncand[d] = n
            pos[d] = 0
        if pos[d] >= ncand[d]:
            ncand[d] = -1
            d -= 1
            continue
        it = cand[d, pos[d]]
        pos[d] += 1
        nd = d + 1
        i = it // L
        for w in range(WT):
            unc[nd, w] = unc[d, w] & ~icov[it, w]
        for w in range(WI):
            alive[nd, w] = alive[d, w]
        for r in range(L):
            jt = i * L + r
            alive[nd, jt // 64] &= ~(np.uint64(1) << np.uint64(jt % 64))
        for j in range(k):
            assign[nd, j] = assign[d, j]
        assign[nd, i] = it % L
        alive[d, it // 64] &= ~(np.uint64(1) << np.uint64(it % 64))
        ncand[nd] = -1
        d = nd
    return 0


@njit(parallel=True, cache=True)
def _cover_par(k, L, icov, itT, tasks_unc, tasks_alive, tasks_assign, max_sol, lprimes):
    nt = tasks_unc.shape[0]
    allstats = np.zeros((nt, 3), dtype=np.int64)
    sols = np.zeros((nt, max_sol, k), dtype=np.int64)
    nsols = np.zeros((nt, 1), dtype=np.int64)
    for j in prange(nt):
        _cover_task(k, L, icov, itT, tasks_unc[j], tasks_alive[j], tasks_assign[j],
                    max_sol, sols[j], nsols[j], allstats[j], lprimes)
    return allstats, sols, nsols


def _bits_set(words, n):
    out = np.zeros(n, dtype=bool)
    for j in range(n):
        out[j] = (int(words[j // 64]) >> (j % 64)) & 1
    return out


def level_cover(cls, L, p, k, max_sol=4, split_depth=2):
    """Exhaustive search for improper level-L lifts of representative cls.
    Returns (stats, list of improper full tuples mod L*p; free coordinates set to an odd
    residue)."""
    T, xs, C, icov, itT = _item_tables(cls, L, p, k)
    nT = len(T)
    WT, WI = icov.shape[1], itT.shape[1]
    nI = k * L
    unc = np.zeros(WT, dtype=np.uint64)
    for t in range(nT):
        unc[t // 64] |= np.uint64(1) << np.uint64(t % 64)
    alive = np.zeros(WI, dtype=np.uint64)
    for it in range(nI):
        alive[it // 64] |= np.uint64(1) << np.uint64(it % 64)
    states = [(unc, alive, -np.ones(k, dtype=np.int64))]
    for _ in range(split_depth):   # same branching rule, unrolled to make parallel tasks
        new = []
        for (u, a, asg) in states:
            ub = _bits_set(u, nT)
            ab = _bits_set(a, nI)
            if not ub.any() or (asg >= 0).all():
                new.append((u, a, asg))
                continue
            cnt = C[ab][:, ub].sum(axis=0)
            ti = np.flatnonzero(ub)[np.argmin(cnt)]
            cands = [it for it in range(nI) if ab[it] and C[it, ti]]
            a2 = a.copy()
            for it in cands:
                i = it // L
                a3 = a2.copy()
                for r in range(L):
                    jt = i * L + r
                    a3[jt // 64] &= ~(np.uint64(1) << np.uint64(jt % 64))
                g = asg.copy()
                g[i] = it % L
                new.append((u & ~icov[it], a3, g))
                a2[it // 64] &= ~(np.uint64(1) << np.uint64(it % 64))
        states = new
    if not states:
        return dict(tasks=0, nodes=0, cover_leaves=0, improper_leaves=0), []
    perm = np.random.RandomState(0).permutation(len(states))   # balance prange chunks
    states = [states[j] for j in perm]
    TU = np.array([s[0] for s in states])
    TA = np.array([s[1] for s in states])
    TG = np.array([s[2] for s in states])
    lprimes = np.array(sorted(set(factor_chain(L))), dtype=np.int64)
    st, sols, ns = _cover_par(k, L, icov, itT, TU, TA, TG, max_sol, lprimes)
    out = []
    for j in range(len(states)):
        for s in range(min(int(ns[j, 0]), max_sol)):
            asg = sols[j, s]
            r = [int(asg[i]) if asg[i] >= 0 else 1 for i in range(k)]
            out.append([int(xs[i, r[i]]) for i in range(k)])
    return dict(tasks=len(states), nodes=int(st[:, 0].sum()), cover_leaves=int(st[:, 1].sum()),
                improper_leaves=int(st[:, 2].sum())), out


@njit(parallel=True, cache=True)
def _cover_knuth(k, L, icov, itT, nT, nsamp, seed):
    """Knuth random-path estimate of the node count of the level_cover tree (from root)."""
    WT = icov.shape[1]
    WI = itT.shape[1]
    est = np.zeros(nsamp)
    leaves = np.zeros(nsamp)
    for s in prange(nsamp):
        np.random.seed(seed + s)
        unc = np.zeros(WT, dtype=np.uint64)
        alive = np.zeros(WI, dtype=np.uint64)
        assign = -np.ones(k, dtype=np.int64)
        cand = np.zeros(k * L, dtype=np.int64)
        cmax = np.zeros(k, dtype=np.int64)
        for t in range(nT):
            unc[t // 64] |= np.uint64(1) << np.uint64(t % 64)
        for it in range(k * L):
            alive[it // 64] |= np.uint64(1) << np.uint64(it % 64)
        weight = 1.0
        tot = 0.0
        while True:
            tot += weight
            nunc = 0
            for w in range(WT):
                nunc += _popc(unc[w])
            nfree = 0
            for i in range(k):
                if assign[i] < 0:
                    nfree += 1
            if nunc == 0:
                leaves[s] = weight
                break
            if nfree == 0:
                break
            for i in range(k):
                cmax[i] = 0
            for ww in range(WI):
                x = alive[ww]
                while x != 0:
                    low = x & (~x + np.uint64(1))
                    b = _popc(low - np.uint64(1))
                    it = ww * 64 + b
                    x ^= low
                    g = 0
                    for w in range(WT):
                        g += _popc(icov[it, w] & unc[w])
                    if g > cmax[it // L]:
                        cmax[it // L] = g
            cmax.sort()
            s_ = 0
            for j in range(nfree):
                s_ += cmax[k - 1 - j]
            if s_ < nunc:
                break
            best_t = -1
            best_c = 1 << 30
            for w in range(WT):
                x = unc[w]
                while x != 0:
                    low = x & (~x + np.uint64(1))
                    b = _popc(low - np.uint64(1))
                    t = w * 64 + b
                    x ^= low
                    c = 0
                    for ww in range(WI):
                        c += _popc(itT[t, ww] & alive[ww])
                    if c < best_c:
                        best_c = c
                        best_t = t
            if best_c == 0:
                break
            n = 0
            for ww in range(WI):
                x = itT[best_t, ww] & alive[ww]
                while x != 0:
                    low = x & (~x + np.uint64(1))
                    b = _popc(low - np.uint64(1))
                    cand[n] = ww * 64 + b
                    n += 1
                    x ^= low
            # the exact search orders candidates by gain; order only matters through the
            # exclusions, so sort by gain here too
            gains = np.zeros(n, dtype=np.int64)
            for j in range(n):
                g = 0
                for w in range(WT):
                    g += _popc(icov[cand[j], w] & unc[w])
                gains[j] = -g
            order = np.argsort(gains, kind='mergesort')
            i0 = np.random.randint(0, n)
            for j in range(i0):
                jt = cand[order[j]]
                alive[jt // 64] &= ~(np.uint64(1) << np.uint64(jt % 64))
            it = cand[order[i0]]
            i = it // L
            for w in range(WT):
                unc[w] &= ~icov[it, w]
            for r in range(L):
                jt = i * L + r
                alive[jt // 64] &= ~(np.uint64(1) << np.uint64(jt % 64))
            assign[i] = it % L
            weight *= n
        est[s] = tot
    return est, leaves


def cover_knuth(cls, L, p, k, nsamp=20000, seed=1):
    T, xs, C, icov, itT = _item_tables(cls, L, p, k)
    return _cover_knuth(k, L, icov, itT, len(T), nsamp, seed)


# ---------------------------------------------------------------- v3, parallel
@njit(cache=True)
def _l1_dfs(k, H, cov, covT, csize, kmax, q, irred, d0, unc0, alive0, once0, twice0,
            chosen0, frontier, F_unc, F_alive, F_once, F_twice, F_chosen, nF,
            out, max_out, nl_io):
    """The level_one3 DFS started from a state at depth d0.  If frontier > d0, nodes at
    depth == frontier are not evaluated but stored in F_* (nF[0] counts them).
    Returns n_nodes evaluated; nl_io[0] += leaves (rows appended to out)."""
    W = cov.shape[1]
    WC = covT.shape[1]
    D = kmax + 2
    unc = np.zeros((D, W), dtype=np.uint64)
    once = np.zeros((D, W), dtype=np.uint64)
    twice = np.zeros((D, W), dtype=np.uint64)
    alive = np.zeros((D, WC), dtype=np.uint64)
    chosen = np.zeros(D, dtype=np.int64)
    cand = np.zeros((D, H + 1), dtype=np.int64)
    gain = np.zeros((D, H + 1), dtype=np.int64)
    hist = np.zeros((D, csize + 1), dtype=np.int64)
    ncand = np.full(D, -1, dtype=np.int64)
    nuncs = np.zeros(D, dtype=np.int64)
    pos = np.zeros(D, dtype=np.int64)
    for w in range(W):
        unc[d0, w] = unc0[w]
        once[d0, w] = once0[w]
        twice[d0, w] = twice0[w]
    for w in range(WC):
        alive[d0, w] = alive0[w]
    for i in range(d0 + 1):
        chosen[i] = chosen0[i]
    n_nodes = 0
    d = d0
    while d >= d0:
        if ncand[d] == -1:
            if frontier > d0 and d == frontier:
                j = nF[0]
                if j < F_unc.shape[0]:
                    for w in range(W):
                        F_unc[j, w] = unc[d, w]
                        F_once[j, w] = once[d, w]
                        F_twice[j, w] = twice[d, w]
                    for w in range(WC):
                        F_alive[j, w] = alive[d, w]
                    for i in range(d + 1):
                        F_chosen[j, i] = chosen[i]
                nF[0] += 1
                d -= 1
                continue
            n_nodes += 1
            st, n, nu = _eval3(kmax, d, q, irred, unc[d], alive[d], once[d], chosen,
                               cov, covT, csize, hist[d], cand[d], gain[d])
            if st == 1:
                j = nl_io[0]
                if j < max_out:
                    for i in range(out.shape[1]):
                        out[j, i] = chosen[i + 1] if i + 1 <= d else 0
                nl_io[0] += 1
                d -= 1
                continue
            if st != 2:
                d -= 1
                continue
            ncand[d] = n
            nuncs[d] = nu
            pos[d] = 0
        if pos[d] >= ncand[d]:
            ncand[d] = -1
            d -= 1
            continue
        j = pos[d]
        c = cand[d, j]
        g = gain[d, j]
        if not _cut_ok(g, kmax - d - 1, hist[d], csize, nuncs[d]):
            ncand[d] = -1
            d -= 1
            continue
        pos[d] += 1
        nd = d + 1
        _child(c, cov, W, WC, unc[d], alive[d], once[d], twice[d],
               unc[nd], alive[nd], once[nd], twice[nd])
        vb = c - 1
        alive[d, vb // 64] &= ~(np.uint64(1) << np.uint64(vb % 64))
        hist[d, g] -= 1
        chosen[nd] = c
        ncand[nd] = -1
        d = nd
    return n_nodes


@njit(parallel=True, cache=True)
def _l1_tasks(k, H, cov, covT, csize, kmax, q, irred, fd, F_unc, F_alive, F_once, F_twice,
              F_chosen, nT, max_out_each):
    nodes = np.zeros(nT, dtype=np.int64)
    nls = np.zeros((nT, 1), dtype=np.int64)
    outs = np.zeros((nT, max_out_each, k), dtype=np.int64)
    dummyU = np.zeros((1, cov.shape[1]), dtype=np.uint64)
    dummyA = np.zeros((1, covT.shape[1]), dtype=np.uint64)
    dummyC = np.zeros((1, kmax + 2), dtype=np.int64)
    for j in prange(nT):
        nF = np.zeros(1, dtype=np.int64)
        nodes[j] = _l1_dfs(k, H, cov, covT, csize, kmax, q, irred, fd, F_unc[j], F_alive[j],
                           F_once[j], F_twice[j], F_chosen[j], -1, dummyU, dummyA, dummyU,
                           dummyU, dummyC, nF, outs[j], max_out_each, nls[j])
    return nodes, nls, outs


def level_one3_par(k, p, frontier=3, max_out_each=256, q=True, irred=True, seed=0):
    """Parallel v3 level one: serial DFS to the frontier depth, then prange over the
    frontier states (shuffled for balance).  Returns (leaves array, n_leaves, n_nodes)."""
    H, M, cov, covT, csize = _tables(k, p)
    qq = q_bound(k, p, k) if q else 0
    W, WC = cov.shape[1], covT.shape[1]
    unc = np.zeros(W, dtype=np.uint64); once = np.zeros(W, dtype=np.uint64)
    twice = np.zeros(W, dtype=np.uint64); alive = np.zeros(WC, dtype=np.uint64)
    _root(H, cov, W, WC, unc, alive, once, twice)
    chosen = np.zeros(k + 2, dtype=np.int64); chosen[1] = 1
    cap = 300_000
    F_unc = np.zeros((cap, W), dtype=np.uint64); F_once = np.zeros((cap, W), dtype=np.uint64)
    F_twice = np.zeros((cap, W), dtype=np.uint64); F_alive = np.zeros((cap, WC), dtype=np.uint64)
    F_chosen = np.zeros((cap, k + 2), dtype=np.int64)
    nF = np.zeros(1, dtype=np.int64)
    out0 = np.zeros((max_out_each, k), dtype=np.int64)
    nl0 = np.zeros(1, dtype=np.int64)
    n0 = _l1_dfs(k, H, cov, covT, csize, k, qq, irred, 1, unc, alive, once, twice, chosen,
                 frontier, F_unc, F_alive, F_once, F_twice, F_chosen, nF, out0, max_out_each, nl0)
    nT = int(nF[0])
    if nT > cap:
        raise RuntimeError('frontier too large')
    perm = np.random.RandomState(seed).permutation(nT)
    nodes, nls, outs = _l1_tasks(k, H, cov, covT, csize, k, qq, irred, frontier,
                                 F_unc[perm], F_alive[perm], F_once[perm], F_twice[perm],
                                 F_chosen[perm], nT, max_out_each)
    if (nls[:, 0] > max_out_each).any() or nl0[0] > max_out_each:
        raise RuntimeError('max_out_each too small')
    rows = [out0[:nl0[0]]] + [outs[j, :nls[j, 0]] for j in range(nT) if nls[j, 0]]
    leaves = np.vstack(rows) if rows else np.zeros((0, k), dtype=np.int64)
    return leaves, int(nl0[0] + nls.sum()), int(n0 + nodes.sum()), nT


def level_one_reps3_par(k, p, frontier=3, **kw):
    import time
    t0 = time.perf_counter()
    leaves, nl, nn, nT = level_one3_par(k, p, frontier, **kw)
    t1 = time.perf_counter()
    rows = _pad_all(leaves, k, (p - 1) // 2)
    reps = canon(rows, p) if len(rows) else np.zeros((0, k), dtype=np.int64)
    return reps, dict(k=k, p=p, nodes=nn, leaves=nl, tasks=nT, padded=int(len(rows)),
                      reps=int(len(reps)), t_search=t1 - t0, t_pad=time.perf_counter() - t1)

"""F1 repair (2026-10-31): the corrected fundamental domain D* and the extra 1/10 start boxes it needs.

D* = { p in [0,1]^6 : p1, p2, p3 <= 1/2,  |s_z - 1| <= |s_y - 1| <= |s_x - 1| },
     s_z = p2 + p4, s_y = p1 + p6, s_x = p3 + p5.

Claim: every orbit of the 48 cube symmetries on configurations p meets D*.
Proof.  p encodes the six contact points (0,p1,p2),(1,p1,p2),(p3,0,p4),(p3,1,p4),(p5,p6,0),(p5,p6,1).  For axis i,
s_i is the sum of the i-th coordinates of the contact points on the two face pairs not normal to e_i.
 (a) The reflection x_i -> 1 - x_i maps the contact set to a contact set of the same form; it replaces the two
     i-coordinates (u, v) by (1-u, 1-v), so s_i -> 2 - s_i and |s_i - 1| is unchanged; the other s_j are untouched.
     Concretely: reflecting y flips (p1, p6), reflecting z flips (p2, p4), reflecting x flips (p3, p5).
 (b) A coordinate permutation pi maps the face pair normal to e_i to the pair normal to e_pi(i) and carries
     i-coordinates to pi(i)-coordinates, so s'_{pi(i)} = s_i: the triple (s_x, s_y, s_z) is permuted by pi, and every
     one of the 6 permutations of the triple is realised.
 Given p, first apply the permutation that sorts the triple (|s_x-1|, |s_y-1|, |s_z-1|) into
 |s_z-1| <= |s_y-1| <= |s_x-1|.  Then apply, independently, the y-, z- and x-reflections to those of p1, p2, p3 that
 exceed 1/2 (p1 is moved only by the y-reflection, p2 only by z, p3 only by x).  By (a) the reflections keep the
 sorted order of |s_i - 1|, and afterwards p1, p2, p3 <= 1/2.  The image lies in D*.  QED
 (D* is not claimed minimal: orbits on the boundary meet it more than once, which is harmless for a covering.)

Start boxes.  A 1/10 box B = prod [c_i/10, (c_i+1)/10] with c in [0,5)^3 x [0,10)^3 meets the closed set D* iff
max(m_y, m_z) <= min(M_y, M_x), where [m_i, M_i] is the range of |s_i - 1| over B (s_x, s_y, s_z depend on disjoint
coordinates, so they vary independently over B and each |s_i - 1| takes every value in [m_i, M_i]).  Then every
point of D* lies in a closed start box that meets D*, so  D* is covered by {boxes meeting D*}.
EXT = {boxes meeting D*} minus NEED (NEED = the existing 39,380 boxes, as in certcheck.py / certify2f.py).

Usage: python dstar_domain.py [--write dstar_ext_boxes.json] [--write-min dstar_ext_min.json]   (exact integer/rational arithmetic, 1 thread)"""
import itertools, random, sys, json, os
from fractions import Fraction as Fr

def pts(p, one=1):
    """contact points with their face labels (axis, side)."""
    p1, p2, p3, p4, p5, p6 = p
    return [((0, p1, p2), 0, 0), ((one, p1, p2), 0, 1), ((p3, 0, p4), 1, 0), ((p3, one, p4), 1, 1),
            ((p5, p6, 0), 2, 0), ((p5, p6, one), 2, 1)]

G = [(perm, fl) for perm in itertools.permutations(range(3)) for fl in itertools.product((0, 1), repeat=3)]

def image(g, p, one=1):
    """g = (perm, fl): x'_i = x_{perm[i]}, then x'_i -> one - x'_i where fl[i].  The face labels are carried along
    (a point on face (j, s) goes to face (i, s xor fl[i]) with perm[i] = j), so ties/zeros are never ambiguous."""
    perm, fl = g; inv = {perm[i]: i for i in range(3)}; d = {}
    for q, ax, sd in pts(p, one):
        r = [q[perm[i]] for i in range(3)]
        r = tuple(one - r[i] if fl[i] else r[i] for i in range(3))
        i = inv[ax]; s2 = sd ^ fl[i]
        assert r[i] == (one if s2 else 0)
        if s2 == 0: d[i] = r
    return (d[0][1], d[0][2], d[1][0], d[1][2], d[2][0], d[2][1])

def inDstar(x, one=1):   # one = the value of "1" (scaled integer grids use one = N)
    h2 = one  # compare 2*p_i <= one
    return all(2 * y <= h2 for y in x[:3]) and abs(x[1] + x[3] - one) <= abs(x[0] + x[5] - one) <= abs(x[2] + x[4] - one)

def canon(p, one=1):
    """the constructive group element of the proof; returns (g, image)."""
    s = {0: p[2] + p[4], 1: p[0] + p[5], 2: p[1] + p[3]}            # s_x, s_y, s_z
    # want new axes: x gets largest |s-1|, z smallest.  new coordinate i reads old coordinate perm[i].
    order = sorted(range(3), key=lambda i: abs(s[i] - one))          # ascending
    perm = (order[2], order[1], order[0])
    q = image((perm, (0, 0, 0)), p, one)
    fl = [0, 0, 0]
    if 2 * q[0] > one: fl[1] = 1     # p1 is a y-coordinate
    if 2 * q[1] > one: fl[2] = 1     # p2 is a z-coordinate
    if 2 * q[2] > one: fl[0] = 1     # p3 is an x-coordinate
    g2 = ((0, 1, 2), tuple(fl))
    return image(g2, q, one)

def in_UD(k):
    k1, k2, k3, k4, k5, k6 = k
    return k2 + k4 <= k1 + k6 + 1 <= k3 + k5 + 2

def need_set():
    return set(c for c in itertools.product(range(5), range(5), range(5), range(10), range(10), range(10))
               if any(in_UD([2 * c[i] + e[i] for i in range(6)]) for e in itertools.product((0, 1), repeat=6)))

def absrange(lo, hi):  # range of |s - 10| for s in [lo, hi] (units of 1/10)
    a, b = lo - 10, hi - 10
    m = 0 if a <= 0 <= b else min(abs(a), abs(b))
    return m, max(abs(a), abs(b))

def meets_dstar(c):
    sz = (c[1] + c[3], c[1] + c[3] + 2); sy = (c[0] + c[5], c[0] + c[5] + 2); sx = (c[2] + c[4], c[2] + c[4] + 2)
    mz, Mz = absrange(*sz); my, My = absrange(*sy); mx, Mx = absrange(*sx)
    return max(my, mz) <= min(My, Mx)

def meets_dstar_open(c):
    """the open box meets the interior of D* (strict version; used to pick orbit representatives)."""
    sz = (c[1] + c[3], c[1] + c[3] + 2); sy = (c[0] + c[5], c[0] + c[5] + 2); sx = (c[2] + c[4], c[2] + c[4] + 2)
    mz, Mz = absrange(*sz); my, My = absrange(*sy); mx, Mx = absrange(*sx)
    return max(my, mz) < min(My, Mx)

def meets_dstar_brute(c, k=2):
    """independent check: some point of the grid (c + j/k)/10, j = 0..k, in each coordinate, lies in D*."""
    N = 10 * k
    axes = [[c[i] * k + j for j in range(k + 1)] for i in range(6)]
    return any(inDstar(x, N) for x in itertools.product(*axes))

def gbox(g, c):
    """image of the closed 1/10 box with lower corner c/10 under g (g acts affinely coordinatewise, so it is again
    a 1/10 grid box); computed on the box centre (2c+1)/20."""
    im = image(g, tuple(2 * x + 1 for x in c), 20)
    return tuple((x - 1) // 2 for x in im)

def orbit_reps(NEED, EXT):
    """Box-orbit reduction.  Every p in [0,1]^6 lies in a closed 1/10 box B of [0,1]^6.  If some g has g(B) in NEED,
    g(p) is already covered.  Otherwise p's orbit is covered once ONE box of B's orbit is certified.  So it suffices
    to certify one representative per box orbit that contains no NEED box.  (Conversely such an orbit is genuinely
    missing: interior points of B have all their images in interiors of non-NEED boxes.)  The representative is the
    lexicographically smallest box of the orbit whose interior meets int D*, so the result is a subset of EXT."""
    seen = set(); reps = []
    for c in itertools.product(range(10), repeat=6):
        if c in seen: continue
        orb = set(gbox(g, c) for g in G); seen |= orb
        if orb & NEED: continue
        r = sorted(x for x in orb if x in EXT and meets_dstar_open(x))
        assert r, 'orbit without a D*-meeting box: contradicts the D* proof'
        reps.append(r[0])
    return sorted(reps)

def main():
    rng = random.Random(20261031)
    # 1. referee counterexample
    p = [Fr(459, 500), Fr(1481, 5000), Fr(1689, 10000), Fr(9529, 10000), Fr(9359, 10000), Fr(3079, 10000)]
    ims = [image(g, p) for g in G]
    print('referee example: images in D*', sum(inDstar(x) for x in ims), '; canon in D*', inDstar(canon(p)))
    # 2. random rationals: brute force over all 48 images AND the constructive element
    bad = 0; N = 20000
    for _ in range(N):
        den = rng.choice([7, 10, 12, 97, 10000])
        q = [Fr(rng.randint(0, den), den) for _ in range(6)]
        ok48 = any(inDstar(image(g, q)) for g in G); okc = inDstar(canon(q))
        bad += (not ok48) + (not okc)
    print('random rationals N=%d: failures %d' % (N, bad))
    # 3. exhaustive integer grids (all points k/n, k = 0..n, including every boundary/tie case)
    for n, mode in ((6, 'all48'), (10, 'canon')):
        cnt = 0; fail = 0
        for x in itertools.product(range(n + 1), repeat=6):
            cnt += 1
            if mode == 'all48': ok = any(inDstar(image(g, x, n), n) for g in G)
            else: ok = inDstar(canon(x, n), n)
            fail += not ok
        print('grid 1/%d (%s): %d points, failures %d' % (n, mode, cnt, fail))
    # 4. start boxes
    NEED = need_set()
    MEET = set(c for c in itertools.product(range(5), range(5), range(5), range(10), range(10), range(10))
               if meets_dstar(c))
    EXT = sorted(MEET - NEED)
    print('need (existing)', len(NEED), 'meet D*', len(MEET), 'meet D* and in need', len(MEET & NEED),
          'EXT = meet D* not in need', len(EXT), 'need not meeting D*', len(NEED - MEET))
    # independent brute-force cross-check of meets_dstar on all candidate boxes (grid 1/20 inside each box)
    mism = 0
    for c in itertools.product(range(5), range(5), range(5), range(10), range(10), range(10)):
        if (c in MEET) != meets_dstar_brute(c): mism += 1
    print('meets_dstar vs brute-force (1/20 grid) mismatches', mism)
    # coverage spot check: random points of D* lie in a closed box of NEED u EXT
    U = NEED | set(EXT); miss = 0; hit = 0
    while hit < 20000:
        q = [Fr(rng.randint(0, 1000), 1000) for _ in range(6)]
        q = canon(q)
        hit += 1
        cand = [[c for c in range((4 if i < 3 else 9) + 1) if Fr(c, 10) <= x <= Fr(c + 1, 10)] for i, x in enumerate(q)]
        miss += not any(c in U for c in itertools.product(*cand))
    print('random D* points not in NEED u EXT:', miss, 'of', hit)
    rngs = [(min(c[i] for c in EXT), max(c[i] for c in EXT)) for i in range(6)]
    print('EXT index ranges (box c covers [c/10,(c+1)/10]):', rngs)
    REPS = orbit_reps(NEED, set(EXT))
    print('box orbits with no NEED box (= start boxes actually needed):', len(REPS))
    print('REPS index ranges:', [(min(c[i] for c in REPS), max(c[i] for c in REPS)) for i in range(6)])
    if '--write-min' in sys.argv:
        fn = sys.argv[sys.argv.index('--write-min') + 1]
        json.dump({'domain': 'one representative per 1/10 box orbit (48 cube symmetries) containing no U^D_20 need box',
                   'grid': '1/10', 'rule': 'lexicographically smallest box of the orbit whose interior meets int D*', 'n': len(REPS),
                   'boxes': [list(c) for c in REPS]}, open(fn, 'w'))
        print('wrote', fn)
    if '--write' in sys.argv:
        fn = sys.argv[sys.argv.index('--write') + 1]
        json.dump({'domain': 'D* = p1,p2,p3<=1/2, |s_z-1|<=|s_y-1|<=|s_x-1|', 'grid': '1/10',
                   'rule': 'closed 1/10 box meets D* and is not in the U^D_20 need set', 'n': len(EXT),
                   'boxes': [list(c) for c in EXT]}, open(fn, 'w'))
        print('wrote', fn)

if __name__ == '__main__':
    main()

"""Independent replay checker for the H_3 <= 13 certificate (newmethod/fullcert).

Written from the spec (NEWMETHOD.md, hadwiger.tex, Prymak arXiv 2112.10698) and the
data format only; no production code was read before this file produced results.

What a leaf asserts (my reading of the spec)
--------------------------------------------
Leaf box P = prod [a_i, a_i+d_i] in configuration space [0,1]^6.
Q_P = intersection of O_v over the vertices v of P, O_v = conv(columns of A_v).
R_P = union of the six rectangles {A_p e_c : p in P}.
  plain leaf ("k":"p"): 13 translates T_1..T_13 + Q_P cover E u R_P (E = cube edges).
  chord leaf ("k":"c", axis i, w_s = e_j + s e_k, {j<k} = other axes): for EACH s in {+1,-1}
     12 translates of Q_P and one translate of (1-lam) Q_P + lam t [0, w_s] cover E u R_P,
     lam = 7/10, t = 569/1000.
Data: "q" = integers over 1e9; I use ONLY q[0:39] as the 13 translation vectors
(and, for chord leaves, q[56:64] merely as *candidate* mu values, see below). The
breakpoints / edge assignment / structure "st" are NOT used: instead each translate's
exact intersection interval with each edge is computed and the union of the 13 intervals
must cover the whole edge [0,1]. Each rectangle must lie inside a single translate
(all 4 corners in it; the translate is convex). This is a structure-free and strictly
more lenient check than Prymak's L(P,tau), but it is sound: it directly proves
E u R_P is covered by the 13 closed bodies.

Chord body interval on an edge line: for fixed mu the body contains
T + mu*lam*t*w + (1-lam)Q_P, whose edge interval is exact; the projection of a convex set
is convex, so the hull of the intervals over any finite set M of mu in [0,1] is contained
in the true interval (inner approximation => sound). M = {0, 1, k/32, stored mu's}.
Rectangle corners in the chord body: exact 1-D feasibility in mu.
Which translate is the chord body: slot v0 (vertex translates in vertex order) is tried
first; any other single slot is also accepted (one chord body + 12 Q_P bodies is all the
reduction needs) and reported.

Q_P: by the sign-pattern lemma (hadwiger.tex Lemma 16; re-derived: h_{O_v}(u) is a max of
terms in disjoint parameter groups, so min over the box is attained at the 8 vertices
w(sigma)), Q_P = intersection of 8 O_w. Option use64=True uses all 64 vertices instead
(cross-check). Facets of O_w: all planes through 3 of its 6 points that support all 6
(O_w is full-dimensional: the three column pairs differ by e_1, e_2, e_3).
Redundant facets are harmless (intersection).

Arithmetic
----------
Exact path: Python Fractions / ints throughout (exact_leaf).
Fast path (fast_plain_batch, numpy): all geometry is scaled to integers:
  G = 640 (every leaf corner is a multiple of 1/640 -- asserted), D = 1e9 (q denominator).
  scaled points X = G*x are ints in [0,640]; facet normals n = (B-A)x(C-A) have
  |n_i| <= 2*640^2 < 2^20; h = n.A has |h| <= 3*2^20*640 < 2^31.
  Point x=v+s e_k (edge) in T+Q, T=q/D:  G*D*n_k*s <= R,  R := D*(h - G n.v) + G n.q.
  |D*(h - G n.v)| <= 1e9*(2^31+3*2^20*640) < 4.1e18, |G n.q| <= 640*3*2^20*|q|_max;
  we assert |q| < 1.5e9 so |G n.q| < 3.1e18; |R| < 7.2e18 < 2^63 = 9.22e18: exact int64.
  Rectangle corner x (G*x int in [0,640]) in T+Q:  D*n.(Gx) - G n.q <= D h, same bounds.
  These comparisons are exact integer comparisons.
  Edge endpoints s_f = R/(G*D*n_k) are formed in float64: R -> float (rel err <= 2^-53),
  G*D*n_k -> float (rel err <= 2^-53), division (rel err <= 2^-53), so
  |s^_f - s_f| <= eps|s_f| with eps = 3.4e-16 (and s^_f, s_f have the same sign).
  hi = min_f s_f over n_k>0. If s <= hi^ - tau with tau = 1e-12*max(1,|hi^|) then for each f
  s_f >= s^_f - 2eps|s^_f| >= hi^ - 2eps|hi^| (x - 2eps|x| is increasing) > hi^ - tau >= s,
  so s satisfies every upper constraint; symmetric for lo. The addition hi^ - tau is itself
  rounded (error <= 2^-53(|hi^|+tau)), which tau = 1e-12*max(1,|hi^|) absorbs with room.
  Hence [lo^+tau, hi^-tau] is a rigorous INNER interval (computed as floats, which are
  exact reals); union coverage of [0,1] by float intervals is then decided by exact float
  comparisons. Constraints with n_k = 0 are decided exactly (R >= 0 in int64).
  If the fast test is not conclusive the leaf goes to the exact path.
  STRICT mode: the inner interval [lo^+tau, hi^-tau] already lies in the OPEN interval
  (s <= hi^-tau < s_f for every upper facet f, shown above), so the float sweep is unchanged;
  n_k = 0 constraints need R > 0 and rectangle corners need strict inequalities.
"""
from fractions import Fraction as Fr
import itertools
import numpy as np

G = 640
DEN = 10**9
LAM = Fr(7, 10)
T_CH = Fr(569, 1000)
EDGES = [(v, k) for k in range(3) for v in itertools.product((0, 1), repeat=3) if v[k] == 0]
TRIPLES = list(itertools.combinations(range(6), 3))


def columns(p):
    p1, p2, p3, p4, p5, p6 = p
    return [(0, p1, p2), (1, p1, p2), (p3, 0, p4), (p3, 1, p4), (p5, p6, 0), (p5, p6, 1)]


def rect_corners(lo, hi):
    """corners of the 6 rectangles R_P cap F (lists of 3-tuples), from box lo/hi (6-vectors)."""
    out = []
    pairs = [(0, 1), (2, 3), (4, 5)]  # parameter pairs for columns (1,2),(3,4),(5,6)
    for c in range(6):
        i, j = pairs[c // 2]
        pts = []
        for x in (lo[i], hi[i]):
            for y in (lo[j], hi[j]):
                p = [None] * 6
                p[i], p[j] = x, y
                for m in range(6):
                    if p[m] is None:
                        p[m] = lo[m]
                pts.append(columns(p)[c])
        out.append(pts)
    return out


def body_vertices(lo, hi, use64=False):
    """parameter vectors w whose O_w intersect to Q_P."""
    if use64:
        return [tuple(hi[m] if b[m] else lo[m] for m in range(6)) for b in itertools.product((0, 1), repeat=6)]
    out = []
    for sx, sy, sz in itertools.product((0, 1), repeat=3):
        # sigma_x controls p3,p5 (x-coords); sigma_y controls p1,p6; sigma_z controls p2,p4
        sel = [sy, sz, sx, sz, sx, sy]
        out.append(tuple(hi[m] if sel[m] else lo[m] for m in range(6)))
    return out


def facets_of(points):
    """H-rep {n.x <= h} of conv(points) (6 points, full-dim), from supporting triples."""
    F = set()
    for a, b, c in TRIPLES:
        A, B, C = points[a], points[b], points[c]
        u = [B[i] - A[i] for i in range(3)]
        v = [C[i] - A[i] for i in range(3)]
        n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        if n == (0, 0, 0):
            continue
        h = n[0] * A[0] + n[1] * A[1] + n[2] * A[2]
        vals = [n[0] * P[0] + n[1] * P[1] + n[2] * P[2] - h for P in points]
        if all(x <= 0 for x in vals):
            F.add((n, h))
        elif all(x >= 0 for x in vals):
            F.add(((-n[0], -n[1], -n[2]), -h))
    return F


def QP_facets_exact(lo, hi, use64=False):
    F = set()
    for w in body_vertices(lo, hi, use64):
        F |= facets_of(columns(w))
    return list(F)


# ---------------- exact path ----------------
# STRICT mode (set hcheck.STRICT = True): every covered point must lie in the INTERIOR of its body.
# Why these and only these become strict: the conclusion needed is
#   E u V_p  subset of  union_i (T_i + int O_p)   [one body: T + int B_lam(p)],
# using only Lemma 2.4(i) Q_P subset O_p (so int Q_P subset int O_p), not Lemma 2.4(ii).
#  * {x : n.x < h for every listed facet} is contained in int Q_P (listed rows are valid
#    inequalities of Q_P), so strict facet inequalities = interior membership (sufficient).
#  * Rectangle R_P cap F in T + int Q_P  <=  its 4 corners strictly inside (int Q_P is convex).
#  * Edge: {s : v+s e_k in T + int Q_P} is an OPEN interval (lo,hi) (constraints with n_k = 0
#    need r > 0); the union of these open intervals must cover the CLOSED [0,1]: first lo < 0,
#    each next lo < current reach, final reach > 1 (covers_unit_exact below).
#  * Chord body: T + (1-lam) int Q_P + lam t [0,w] is open and lies in T + (1-lam) int O_p +
#    lam t [0,w], which is inside int K up to translation (Lemma inner). A point x is in it iff
#    some mu in the CLOSED [0,1] makes every facet inequality strict; on an edge, the hull of the
#    open intervals over finitely many mu lies in the (convex) projection, so it is a sound
#    inner approximation as before. mu itself may be 0 or 1 (the segment is closed).
#  * Breakpoints / vertices need nothing separate: they are points of edges.
STRICT = False


def interval_exact(Fs, alpha, T, v, k):
    """exact {s : v + s e_k in T + alpha*Q} (closed; open in STRICT mode); (lo, hi) or None.
    None as lo/hi means unbounded."""
    lo, hi = None, None
    for n, h in Fs:
        r = alpha * h - sum(n[i] * (v[i] - T[i]) for i in range(3))
        if n[k] == 0:
            if r < 0 or (STRICT and r == 0):
                return None
            continue
        s = r / n[k]
        if n[k] > 0:
            hi = s if hi is None or s < hi else hi
        else:
            lo = s if lo is None or s > lo else lo
    if lo is not None and hi is not None and (lo > hi or (STRICT and lo == hi)):
        return None
    return (lo, hi)


def point_in_exact(Fs, alpha, T, x):
    if STRICT:
        return all(sum(n[i] * (x[i] - T[i]) for i in range(3)) < alpha * h for n, h in Fs)
    return all(sum(n[i] * (x[i] - T[i]) for i in range(3)) <= alpha * h for n, h in Fs)


def point_in_chord_exact(Fs, T, cw, x):
    """exists mu in [0,1]: x - T - mu*cw in (1-lam) Q  (in (1-lam) int Q in STRICT mode)."""
    alpha = 1 - LAM
    L, U = None, None  # bounds on mu from facets (strict in STRICT mode)
    for n, h in Fs:
        r = alpha * h - sum(n[i] * (x[i] - T[i]) for i in range(3))
        b = -sum(n[i] * cw[i] for i in range(3))  # constraint: b*mu <= r  (< r if STRICT)
        if b == 0:
            if r < 0 or (STRICT and r == 0):
                return False
        elif b > 0:
            U = r / b if U is None or r / b < U else U
        else:
            L = r / b if L is None or r / b > L else L
    if not STRICT:
        lo = Fr(0) if L is None else max(Fr(0), L)
        hi = Fr(1) if U is None else min(Fr(1), U)
        return lo <= hi
    # exists mu in [0,1] with L < mu < U
    if L is not None and U is not None and not L < U:
        return False
    if L is not None and not L < 1:
        return False
    if U is not None and not U > 0:
        return False
    return True


def covers_unit_exact(intervals):
    BIG = Fr(10**9)
    iv = [((-BIG if lo is None else lo), (BIG if hi is None else hi)) for lo, hi in intervals]
    if not STRICT:
        iv.sort()
        reach = Fr(0)
        covered_zero = False
        for lo, hi in iv:
            if lo > reach:
                break
            if hi >= reach:
                reach = hi
                covered_zero = True
        return covered_zero and reach >= 1
    # open intervals: point cur must lie strictly inside some interval; advance to its hi
    cur = Fr(0)
    while True:
        best = max((hi for lo, hi in iv if lo < cur), default=None)
        if best is None or best <= cur:
            return False
        cur = best
        if cur > 1:
            return True


def leaf_box(leaf):
    a = [Fr(x) for x in leaf['a']]
    d = [Fr(x) for x in leaf['d']]
    return a, [a[i] + d[i] for i in range(6)]


def check_cover_exact(Fs, lo, hi, translates, chord=None):
    """translates: list of 13 Fraction 3-vectors. chord: None or (slot, cw, mus).
    Returns (ok, reason)."""
    alpha1 = Fr(1)
    # edges
    for v, k in EDGES:
        ivs = []
        for idx, T in enumerate(translates):
            if chord is not None and idx == chord[0]:
                slot, cw, mus = chord
                los, his = [], []
                for mu in mus:
                    Tp = tuple(T[i] + mu * cw[i] for i in range(3))
                    I = interval_exact(Fs, 1 - LAM, Tp, v, k)
                    if I is not None:
                        los.append(I[0]); his.append(I[1])
                if los:
                    lo_ = None if any(x is None for x in los) else min(los)
                    hi_ = None if any(x is None for x in his) else max(his)
                    ivs.append((lo_, hi_))
            else:
                I = interval_exact(Fs, alpha1, T, v, k)
                if I is not None:
                    ivs.append(I)
        if not covers_unit_exact(ivs):
            return False, 'edge %s+s*e%d not covered' % (v, k)
    # rectangles
    for c, corners in enumerate(rect_corners(lo, hi)):
        ok = False
        for idx, T in enumerate(translates):
            if chord is not None and idx == chord[0]:
                if all(point_in_chord_exact(Fs, T, chord[1], x) for x in corners):
                    ok = True
                    break
            elif all(point_in_exact(Fs, alpha1, T, x) for x in corners):
                ok = True
                break
        if not ok:
            return False, 'rectangle of column %d not in any single translate' % (c + 1)
    return True, ''


def translates_of(q):
    return [tuple(Fr(q[3 * m + i], DEN) for i in range(3)) for m in range(13)]


VERT_ORDERS = {
    # candidate vertex-index -> coordinates conventions (only used to pick which slot to try first)
    'x-major': lambda i: ((i >> 2) & 1, (i >> 1) & 1, i & 1),
}


def exact_leaf(leaf, use64=False):
    """full exact verification of one leaf. returns (ok, info)."""
    lo, hi = leaf_box(leaf)
    Fs = QP_facets_exact(lo, hi, use64)
    if leaf['k'] == 'p':
        q = leaf['c']['q']
        assert len(q) >= 39
        ok, why = check_cover_exact(Fs, lo, hi, translates_of(q))
        return ok, why
    assert leaf['k'] == 'c'
    ax = leaf['ax']
    j, k = [m for m in range(3) if m != ax]
    info = {}
    for sname, s in (('plus', 1), ('minus', -1)):
        cert = leaf[sname]
        q = cert['q']
        w = [0, 0, 0]; w[j] = 1; w[k] = s
        cw = tuple(LAM * T_CH * w[i] for i in range(3))
        mus = sorted(set([Fr(0), Fr(1)] + [Fr(m, 32) for m in range(33)] +
                         [Fr(x, DEN) for x in q[56:] if 0 <= x <= DEN]))
        tr = translates_of(q)
        v0 = cert['st'][1]
        order = [v0] + [m for m in range(13) if m != v0]
        okany = False
        for slot in order:
            ok, why = check_cover_exact(Fs, lo, hi, tr, chord=(slot, cw, mus))
            if ok:
                okany = True
                info[sname] = slot
                break
        if not okany:
            return False, 'sign %s: no slot works (%s)' % (sname, why)
    return True, info


# ---------------- fast path (plain leaves, numpy int64) ----------------

_TRI = np.array(TRIPLES, dtype=np.int64)  # (20,3)
_EV = np.array([e[0] for e in EDGES], dtype=np.int64)  # (12,3)
_EK = np.array([e[1] for e in EDGES], dtype=np.int64)
_SIG = np.array([[sy, sz, sx, sz, sx, sy] for sx, sy, sz in itertools.product((0, 1), repeat=3)], dtype=np.int64)  # (8,6)
EPS_TAU = 1e-12


def _cols_np(P):
    """P: (...,6) int (scaled by G). returns (...,6,3) scaled column points."""
    z = np.zeros(P.shape[:-1], dtype=np.int64)
    g = np.full(P.shape[:-1], G, dtype=np.int64)
    p1, p2, p3, p4, p5, p6 = [P[..., m] for m in range(6)]
    cols = [(z, p1, p2), (g, p1, p2), (p3, z, p4), (p3, g, p4), (p5, p6, z), (p5, p6, g)]
    return np.stack([np.stack(c, axis=-1) for c in cols], axis=-2)


def facets_np(LO, HI):
    """LO,HI: (L,6) int64 scaled. returns n (L,160,3), h (L,160), valid (L,160) bool."""
    L = LO.shape[0]
    W = np.where(_SIG[None, :, :] == 1, HI[:, None, :], LO[:, None, :])  # (L,8,6)
    pts = _cols_np(W)  # (L,8,6,3)
    A = pts[:, :, _TRI[:, 0], :]; B = pts[:, :, _TRI[:, 1], :]; C = pts[:, :, _TRI[:, 2], :]  # (L,8,20,3)
    u = B - A; v = C - A
    n = np.stack([u[..., 1] * v[..., 2] - u[..., 2] * v[..., 1],
                  u[..., 2] * v[..., 0] - u[..., 0] * v[..., 2],
                  u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]], axis=-1)
    h = (n * A).sum(-1)  # (L,8,20)
    vals = np.einsum('lbtk,lbpk->lbtp', n, pts) - h[..., None]  # (L,8,20,6)
    allle = (vals <= 0).all(-1)
    allge = (vals >= 0).all(-1)
    nz = (n != 0).any(-1)
    valid = nz & (allle | allge)
    flip = nz & allge & ~allle
    n = np.where(flip[..., None], -n, n)
    h = np.where(flip, -h, h)
    return n.reshape(L, 160, 3), h.reshape(L, 160), valid.reshape(L, 160)


def fast_plain_batch(LO, HI, Q):
    """LO,HI: (L,6) int64 (scaled by G); Q: (L,13,3) int64 numerators over DEN.
    returns bool array (L,): True = proven covered (rigorously); False = inconclusive."""
    assert np.abs(Q).max() < 1.5e9
    n, h, valid = facets_np(LO, HI)
    assert np.abs(n).max() < 2**20 and np.abs(h).max() < 2**31
    L = LO.shape[0]
    ok = np.ones(L, dtype=bool)
    nq = np.einsum('lfk,ltk->lft', n, Q)  # (L,160,13): |.| < 3*2^20*1.5e9 ~ 4.7e15
    Gnq = G * nq
    BIG = np.int64(0)
    # ---- edges
    for e in range(12):
        v = _EV[e]; k = _EK[e]
        hv = h - G * (n @ v)  # (L,160)
        R = DEN * hv[:, :, None] + Gnq  # (L,160,13) exact int64 (see bounds)
        nk = n[:, :, k][:, :, None]
        validb = np.broadcast_to(valid[:, :, None], R.shape)
        # n_k == 0 constraints: need R >= 0
        zero_bad = (validb & (nk == 0) & ((R <= 0) if STRICT else (R < 0))).any(1)  # (L,13) -> body misses line
        den = (float(G) * float(DEN)) * nk.astype(np.float64)
        with np.errstate(divide='ignore', invalid='ignore'):
            s = R.astype(np.float64) / den
        up = np.where(validb & (nk > 0), s, np.inf).min(1)  # (L,13)
        dn = np.where(validb & (nk < 0), s, -np.inf).max(1)
        tau_u = EPS_TAU * np.maximum(1.0, np.abs(up))
        tau_d = EPS_TAU * np.maximum(1.0, np.abs(dn))
        hi_in = np.where(np.isfinite(up), up - tau_u, np.inf)
        lo_in = np.where(np.isfinite(dn), dn + tau_d, -np.inf)
        good = (~zero_bad) & (lo_in <= hi_in)
        lo_in = np.where(good, lo_in, np.inf)
        hi_in = np.where(good, hi_in, -np.inf)
        # sweep: sort by lo
        order = np.argsort(lo_in, axis=1)
        los = np.take_along_axis(lo_in, order, 1)
        his = np.take_along_axis(hi_in, order, 1)
        reach = np.zeros(L)
        started = np.zeros(L, dtype=bool)
        alive = np.ones(L, dtype=bool)
        for m in range(13):
            lm = los[:, m]; hm = his[:, m]
            gap = alive & (lm > reach) & (reach < 1.0)
            alive &= ~gap
            ext = alive & (hm >= reach)
            reach = np.where(ext, np.maximum(reach, hm), reach)
            started |= ext
        ok &= alive & started & (reach >= 1.0)
        # note: alive=False after a gap means intervals starting later can't help (sorted).
    # ---- rectangles: exact int64
    rect_ok = np.zeros((L, 6), dtype=bool)
    for c, cornersel in enumerate(_rect_corner_np(LO, HI)):
        # cornersel: (L,4,3) scaled ints in [0,640]
        nx = np.einsum('lfk,lck->lfc', n, cornersel)  # (L,160,4)
        lhs = DEN * nx[:, :, :, None] - Gnq[:, :, None, :]  # (L,160,4,13)
        rhs = (DEN * h)[:, :, None, None]
        sat = ((lhs < rhs) if STRICT else (lhs <= rhs)) | ~valid[:, :, None, None]
        inside = sat.all(1).all(1)  # (L,13)
        rect_ok[:, c] = inside.any(1)
    ok &= rect_ok.all(1)
    return ok


def _rect_corner_np(LO, HI):
    L = LO.shape[0]
    pairs = [(0, 1), (2, 3), (4, 5)]
    out = []
    for c in range(6):
        i, j = pairs[c // 2]
        pts = []
        for xi in (LO[:, i], HI[:, i]):
            for yj in (LO[:, j], HI[:, j]):
                P = LO.copy()
                P[:, i] = xi; P[:, j] = yj
                pts.append(_cols_np(P)[:, c, :])
        out.append(np.stack(pts, axis=1))
    return out


# ---------------- tiling of a start box ----------------

def tiling_ok(sb_lo, sb_hi, LO, HI):
    """sb_lo/hi: int 6-tuples (scaled); LO,HI: (L,6) int arrays. Proves that the leaves tile the
    start box: recursively find an axis-parallel hyperplane (any axis, the midpoint of the current
    box) that no leaf crosses; recurse; each final cell must equal exactly one leaf.
    (Midpoint-of-any-axis is more general than the spec's longest-side rule.)
    Returns (ok, msg)."""
    LO = np.asarray(LO); HI = np.asarray(HI)
    if (LO < np.asarray(sb_lo)).any() or (HI > np.asarray(sb_hi)).any() or (HI <= LO).any():
        return False, 'leaf outside start box or degenerate'
    stack = [(np.array(sb_lo), np.array(sb_hi), np.arange(LO.shape[0]))]
    while stack:
        lo, hi, idx = stack.pop()
        if idx.size == 0:
            return False, 'gap: empty cell %s-%s' % (lo.tolist(), hi.tolist())
        if idx.size == 1:
            i = idx[0]
            if (LO[i] == lo).all() and (HI[i] == hi).all():
                continue
        # leaf equal to the cell?
        eq = ((LO[idx] == lo) & (HI[idx] == hi)).all(1)
        if eq.any():
            return False, 'cell %s equals a leaf but has %d leaves (overlap/dup)' % (lo.tolist(), idx.size)
        done = False
        for ax in np.argsort(-(hi - lo), kind='stable'):
            if (hi[ax] - lo[ax]) % 2:
                continue
            mid = (lo[ax] + hi[ax]) // 2
            left = HI[idx, ax] <= mid
            right = LO[idx, ax] >= mid
            if (left | right).all():
                l2 = hi.copy(); l2[ax] = mid
                r2 = lo.copy(); r2[ax] = mid
                stack.append((lo, l2, idx[left]))
                stack.append((r2, hi, idx[right]))
                done = True
                break
        if not done:
            return False, 'cell %s-%s cannot be split cleanly (%d leaves)' % (lo.tolist(), hi.tolist(), idx.size)
    return True, ''

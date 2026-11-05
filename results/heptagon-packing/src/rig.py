"""Rigorous vectorised interval arithmetic in numpy (no libm trust).

Every elementary operation is IEEE round-to-nearest; enclosures are widened outward by one
ulp per operation (nextafter), which is rigorous for +, -, *, /, sqrt. Trigonometric and
arctangent enclosures do NOT call np.cos/np.sin/np.arctan for the bounds: they use argument
reduction with a double-double enclosure of pi/2 and Taylor polynomials with explicit
Lagrange/alternating remainders, all evaluated in this interval arithmetic.
(np.cos etc. are used only to choose branch points/monotonic pieces, never for bounds.)
"""
import numpy as np
INF = np.inf
_EPS = 2.0**-51          # >= 2 ulps relative for normal numbers
_ETA = 2.0**-1000        # absolute floor (covers zero/subnormal results)
def dn(a):
    # a - (|a| 2^-51 + 2^-1000): exact result <= a - 2 ulp(a) (normal a), so round-to-nearest
    # cannot return a value above a - ulp(a)/2 ... hence strictly below the exact a. Rigorous.
    return a - (np.abs(a)*_EPS + _ETA)
def up(a):
    return a + (np.abs(a)*_EPS + _ETA)

class R:
    """interval array [lo, hi]"""
    __slots__ = ('lo', 'hi')
    def __init__(s, lo, hi=None):
        if hi is None: hi = lo
        s.lo = np.asarray(lo, float); s.hi = np.asarray(hi, float)
    def __add__(s, o):
        if isinstance(o, R): return R(dn(s.lo+o.lo), up(s.hi+o.hi))
        return R(dn(s.lo+o), up(s.hi+o))
    __radd__ = __add__
    def __sub__(s, o):
        if isinstance(o, R): return R(dn(s.lo-o.hi), up(s.hi-o.lo))
        return R(dn(s.lo-o), up(s.hi-o))
    def __rsub__(s, o): return R(dn(o-s.hi), up(o-s.lo))
    def __neg__(s): return R(-s.hi, -s.lo)
    def __mul__(s, o):
        if not isinstance(o, R):
            o = np.asarray(o, float)
            a, b = s.lo*o, s.hi*o
            return R(dn(np.minimum(a, b)), up(np.maximum(a, b)))
        p1, p2, p3, p4 = s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi
        return R(dn(np.minimum(np.minimum(p1, p2), np.minimum(p3, p4))),
                 up(np.maximum(np.maximum(p1, p2), np.maximum(p3, p4))))
    __rmul__ = __mul__
    def sq(s):
        a, b = s.lo*s.lo, s.hi*s.hi
        lo = np.where((s.lo <= 0) & (s.hi >= 0), 0.0, dn(np.minimum(a, b)))
        return R(lo, up(np.maximum(a, b)))
    def inv(s):
        pos = s.lo > 0; neg = s.hi < 0
        lo = np.where(pos | neg, dn(1.0/np.where(pos | neg, s.hi, 1.0)), -INF)
        hi = np.where(pos | neg, up(1.0/np.where(pos | neg, s.lo, 1.0)), INF)
        return R(lo, hi)
    def sqrt(s):
        return R(dn(np.sqrt(np.maximum(s.lo, 0.0))), up(np.sqrt(np.maximum(s.hi, 0.0))))
    def mag(s): return np.maximum(np.abs(s.lo), np.abs(s.hi))
    def hull(s, o): return R(np.minimum(s.lo, o.lo), np.maximum(s.hi, o.hi))
    def __getitem__(s, k): return R(s.lo[k], s.hi[k])

# double-double enclosure of pi/2 and pi
def _enclose_pio2():
    # pi/2 = 1.57079632679489661923132169163975...
    lo = np.float64(1.5707963267948966)    # < pi/2
    hi = np.nextafter(lo, INF)              # > pi/2
    return R(lo, hi)
PIO2 = _enclose_pio2()
PI = PIO2*2.0

_FACT = [1.0]
for _k in range(1, 40): _FACT.append(_FACT[-1]*_k)

def _sin_taylor(r, N=12):
    """r: R with |r| <= 0.8. sin r = sum_{k<N} (-1)^k r^{2k+1}/(2k+1)!, remainder |r|^{2N+1}/(2N+1)!"""
    r2 = r.sq(); acc = R(np.zeros_like(r.lo))
    for k in range(N-1, -1, -1):
        c = (-1)**k / _FACT[2*k+1]
        acc = acc*r2 + R(dn(c), up(c))
    acc = acc*r
    m = r.mag(); rem = up(m**(2*N+1)/_FACT[2*N+1]*(1+1e-12))
    return R(dn(acc.lo - rem), up(acc.hi + rem))

def _cos_taylor(r, N=12):
    r2 = r.sq(); acc = R(np.zeros_like(r.lo))
    for k in range(N-1, -1, -1):
        c = (-1)**k / _FACT[2*k]
        acc = acc*r2 + R(dn(c), up(c))
    m = r.mag(); rem = up(m**(2*N)/_FACT[2*N]*(1+1e-12))
    return R(dn(acc.lo - rem), up(acc.hi + rem))

def _point_sincos(x):
    """rigorous enclosures of sin(x), cos(x) for float array x (|x| < 1e3)."""
    k = np.round(x/1.5707963267948966)
    r = R(x) - PIO2*k          # interval reduction
    s, c = _sin_taylor(r), _cos_taylor(r)
    q = np.mod(k, 4).astype(int)
    S = R(np.where(q == 0, s.lo, np.where(q == 1, c.lo, np.where(q == 2, -s.hi, -c.hi))),
          np.where(q == 0, s.hi, np.where(q == 1, c.hi, np.where(q == 2, -s.lo, -c.lo))))
    C = R(np.where(q == 0, c.lo, np.where(q == 1, -s.hi, np.where(q == 2, -c.hi, s.lo))),
          np.where(q == 0, c.hi, np.where(q == 1, -s.lo, np.where(q == 2, -c.lo, s.hi))))
    return S, C

def rcos(a):
    """cos over interval a (any width; width >= 2pi gives [-1,1])."""
    sl, cl = _point_sincos(a.lo); sh, ch = _point_sincos(a.hi)
    lo = np.minimum(cl.lo, ch.lo); hi = np.maximum(cl.hi, ch.hi)
    # interior extrema: max at 2k pi, min at (2k+1) pi  (branch detection with a safe margin)
    tp = 2*np.pi
    has_max = np.floor(a.hi/tp + 1e-12) >= np.ceil(a.lo/tp - 1e-12)
    has_min = np.floor((a.hi - np.pi)/tp + 1e-12) >= np.ceil((a.lo - np.pi)/tp - 1e-12)
    hi = np.where(has_max, 1.0, hi); lo = np.where(has_min, -1.0, lo)
    wide = (a.hi - a.lo) >= 6.2
    return R(np.where(wide, -1.0, np.maximum(lo, -1.0)), np.where(wide, 1.0, np.minimum(hi, 1.0)))

def rsin(a):
    return rcos(a - PIO2)

def _atan_point(t):
    """rigorous atan(t) enclosure for float array t >= 0."""
    t = np.asarray(t, float)
    big = t > 1.0
    u = R(np.where(big, t, 1.0)).inv()
    v = R(np.where(big, u.lo, t), np.where(big, u.hi, t))     # v in [0,1]
    # two halvings: atan(v) = 2 atan(v / (1 + sqrt(1 + v^2)))
    for _ in range(2):
        v = v * ((v.sq() + 1.0).sqrt() + 1.0).inv()
    # series on |v| <= 0.2: sum_{k<N} (-1)^k v^{2k+1}/(2k+1), alternating remainder
    N = 14; v2 = v.sq(); acc = R(np.zeros_like(v.lo))
    for k in range(N-1, -1, -1):
        c = (-1)**k/(2*k+1)
        acc = acc*v2 + R(dn(c), up(c))
    acc = acc*v
    m = v.mag(); rem = up(m**(2*N+1)/(2*N+1)*(1+1e-12))
    at = R(dn(acc.lo - rem), up(acc.hi + rem))*4.0
    res = R(np.where(big, (PIO2 - at).lo, at.lo), np.where(big, (PIO2 - at).hi, at.hi))
    return res

def ratan(a):
    """atan over interval a with a.lo >= 0 (all our uses), monotone."""
    lo = np.maximum(a.lo, 0.0)
    return R(_atan_point(lo).lo, _atan_point(np.maximum(a.hi, 0.0)).hi)

if __name__ == '__main__':
    import mpmath as mp
    mp.mp.dps = 40
    xs = np.linspace(-12, 12, 20001)
    S, C = _point_sincos(xs)
    worst = 0
    for x, sl, sh, cl, ch in zip(xs, S.lo, S.hi, C.lo, C.hi):
        s, c = mp.sin(mp.mpf(x)), mp.cos(mp.mpf(x))
        assert sl <= s <= sh and cl <= c <= ch, x
        worst = max(worst, sh - sl, ch - cl)
    print('sin/cos enclosures OK on 20001 points; max width', worst)
    ts = np.concatenate([np.linspace(0, 5, 10001), [1e-12, 1.0, 1e6]])
    A = _atan_point(ts)
    for t, lo, hi in zip(ts, A.lo, A.hi):
        v = mp.atan(mp.mpf(t)); assert lo <= v <= hi, t
    print('atan enclosures OK; max width', (A.hi - A.lo).max())
    # interval cos over random intervals vs dense sampling
    rng = np.random.default_rng(0); lo = rng.uniform(-10, 10, 3000); hi = lo + rng.uniform(0, 1.5, 3000)
    E = rcos(R(lo, hi))
    for l, h, el, eh in zip(lo, hi, E.lo, E.hi):
        g = np.linspace(l, h, 50); vals = [float(mp.cos(mp.mpf(v))) for v in g]
        assert el <= min(vals) and max(vals) <= eh
    print('interval cos OK')

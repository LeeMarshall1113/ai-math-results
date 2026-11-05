"""Pilot interval branch-and-bound (Step 2).

Lemma to certify (pilot, analogue of Hales-Kusner Lemma 35):
  every nonobtuse triangle whose vertices are the centres of three pairwise
  non-overlapping regular heptagons (circumradius 1) has area >= TARGET.

Coordinates: c0 = 0, c1 = (d, 0), c2 = (x, y); orientations t0, t1, t2 (mod 2pi/7).
WLOG (relabel) c0c1 is a longest edge, and (mirror) x <= d/2. Then x in [0, d/2],
|c2| <= d, |c2 - c1| <= d, every edge >= 2u (else the inscribed discs overlap).

A box is discarded when one of these holds on the whole box (interval arithmetic):
  (A) area lower bound >= TARGET;
  (B) a centre-geometry constraint fails (some edge < 2u, |c2|>d, obtuse at c2);
  (C) some vertex of one heptagon lies strictly inside another heptagon (overlap).
Otherwise it is bisected along its widest coordinate.

Interval arithmetic: numpy float64 with outward rounding by nextafter after every
operation; cos/sin padded by 4e-16 absolute. This is NOT a formally rigorous
implementation (numpy's libm cos/sin carry no ulp guarantee); it is a cost pilot.
For a certificate, swap in python-flint arb or C with directed rounding.
"""
import numpy as np, sys, time, json

U = np.cos(np.pi/7); P7 = 2*np.pi/7; TWO_U = 2*U
INF = np.inf
def dn(a): return np.nextafter(a, -INF)
def up(a): return np.nextafter(a, INF)

class I:
    __slots__ = ('lo', 'hi')
    def __init__(s, lo, hi): s.lo, s.hi = lo, hi
    def __add__(s, o):
        if isinstance(o, I): return I(dn(s.lo+o.lo), up(s.hi+o.hi))
        return I(dn(s.lo+o), up(s.hi+o))
    __radd__ = __add__
    def __sub__(s, o):
        if isinstance(o, I): return I(dn(s.lo-o.hi), up(s.hi-o.lo))
        return I(dn(s.lo-o), up(s.hi-o))
    def __rsub__(s, o): return I(dn(o-s.hi), up(o-s.lo))
    def __neg__(s): return I(-s.hi, -s.lo)
    def __mul__(s, o):
        if not isinstance(o, I):
            o = float(o)
            return I(dn(s.lo*o), up(s.hi*o)) if o >= 0 else I(dn(s.hi*o), up(s.lo*o))
        p = np.stack([s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi])
        return I(dn(p.min(0)), up(p.max(0)))
    __rmul__ = __mul__
    def sq(s):
        lo2, hi2 = s.lo*s.lo, s.hi*s.hi
        lo = np.where((s.lo <= 0) & (s.hi >= 0), 0.0, np.minimum(lo2, hi2))
        return I(dn(lo), up(np.maximum(lo2, hi2)))

def icos(a):
    """cos over interval a (width < pi assumed)."""
    lo, hi = a.lo, a.hi
    cl, ch = np.cos(lo), np.cos(hi)
    mn = np.minimum(cl, ch); mx = np.maximum(cl, ch)
    # maxima at 2k pi, minima at (2k+1) pi
    k = np.ceil(lo/(2*np.pi)); has_max = 2*np.pi*k <= hi
    k2 = np.ceil((lo-np.pi)/(2*np.pi)); has_min = 2*np.pi*k2+np.pi <= hi
    mx = np.where(has_max, 1.0, mx); mn = np.where(has_min, -1.0, mn)
    return I(np.maximum(mn-4e-16, -1.0), np.minimum(mx+4e-16, 1.0))

def isin(a):
    return icos(a - np.pi/2)

def hept_verts(cx, cy, t):
    """intervals of the 7 vertices"""
    V = []
    for k in range(7):
        a = t + (2*np.pi*k/7)
        V.append((cx + icos(a), cy + isin(a)))
    return V

def strictly_inside(V, cx, cy, t, px, py):
    """bool array: point (px,py) strictly inside heptagon (cx,cy,t) for the whole box.
    Use the half-plane form: n_k . (p - c) < u for all k, n_k = (cos(t+(2k+1)pi/7), sin(...))."""
    ok = None
    for k in range(7):
        a = t + ((2*k+1)*np.pi/7)
        val = icos(a)*(px - cx) + isin(a)*(py - cy)
        c = val.hi < U*(1-1e-15)
        ok = c if ok is None else (ok & c)
    return ok

def evaluate(B, target):
    """B: array (n, 6, 2) of boxes [d, x, y, t0, t1, t2]. Returns discard mask."""
    d = I(B[:, 0, 0], B[:, 0, 1]); x = I(B[:, 1, 0], B[:, 1, 1]); y = I(B[:, 2, 0], B[:, 2, 1])
    t = [I(B[:, 3+i, 0], B[:, 3+i, 1]) for i in range(3)]
    # (A) area
    area_lo = dn(dn(B[:, 0, 0]*B[:, 2, 0])*0.5)
    disc = area_lo >= target
    # (B) centre geometry
    r02 = x.sq() + y.sq(); r12 = (x - d).sq() + y.sq(); dd = d.sq()
    disc |= r02.hi < TWO_U**2*(1-1e-15)          # edge c0c2 too short
    disc |= r12.hi < TWO_U**2*(1-1e-15)
    disc |= r02.lo > dd.hi                       # c0c1 not longest
    disc |= r12.lo > dd.hi
    dot2 = (-x)*(d - x) + y*y                    # (c0-c2).(c1-c2): obtuse at c2 if < 0
    disc |= dot2.hi < 0
    disc |= dn(B[:, 1, 0]*2.0) > B[:, 0, 1]      # x > d/2 (mirror reduction)
    # (C) overlaps (vertex of one strictly inside another)
    zero = I(np.zeros(len(B)), np.zeros(len(B)))
    C = [(zero, zero), (d, zero), (x, y)]
    V = [hept_verts(C[i][0], C[i][1], t[i]) for i in range(3)]
    for i in range(3):
        for j in range(3):
            if i == j: continue
            for (px, py) in V[j]:
                disc |= strictly_inside(V[i], C[i][0], C[i][1], t[i], px, py)
    return disc

def run(target, region, max_boxes=10**9, tlimit=1e9, minw=1e-7, log_every=10):
    lo = np.array(region['lo']); hi = np.array(region['hi'])
    stack = np.stack([lo, hi], -1)[None]
    processed = 0; t0 = time.time(); depth_hist = {}; stuck = []
    last = t0
    while len(stack):
        batch = stack[-20000:]; stack = stack[:-20000]
        disc = evaluate(batch, target)
        processed += len(batch)
        keep = batch[~disc]
        if len(keep):
            w = keep[:, :, 1] - keep[:, :, 0]
            j = np.argmax(w, 1)
            small = w[np.arange(len(keep)), j] < minw
            if small.any():
                stuck.append(keep[small][:5]); keep = keep[~small]
                j = j[~small]
            mid = 0.5*(keep[np.arange(len(keep)), j, 0] + keep[np.arange(len(keep)), j, 1])
            L = keep.copy(); R = keep.copy()
            L[np.arange(len(keep)), j, 1] = mid; R[np.arange(len(keep)), j, 0] = mid
            stack = np.concatenate([stack, L, R]) if len(stack) else np.concatenate([L, R])
        if time.time() - last > log_every:
            print(f'  processed {processed:,} open {len(stack):,} [{time.time()-t0:.0f}s]', flush=True); last = time.time()
        if processed > max_boxes or time.time()-t0 > tlimit:
            return dict(done=False, processed=processed, open=len(stack), secs=time.time()-t0, stuck=len(stuck))
    return dict(done=True, processed=processed, open=0, secs=time.time()-t0, stuck=sum(len(s) for s in stuck))

if __name__ == '__main__':
    target = float(sys.argv[1]); frac = float(sys.argv[2]); tl = float(sys.argv[3])
    # full domain; t0 restricted to the first `frac` of its period (pilot sub-region)
    region = dict(lo=[TWO_U, 0.0, 1.2, 0.0, 0.0, 0.0],
                  hi=[2.25, 1.125, 1.75, P7*frac, P7, P7])
    r = run(target, region, tlimit=tl)
    r.update(target=target, frac=frac)
    print(json.dumps(r))

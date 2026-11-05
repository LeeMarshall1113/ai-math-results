"""Interval branch-and-bound PILOT for the main nonobtuse inequality

    b(T) = area(T) + f(e01) + f(e12) + f(e20) >= ACRIT + margin_target

over a sub-region of canonical nonobtuse triangles (c0 = 0, c1 = (d, 0), c2 = (x, y),
orientations t0, t1, t2; c0c1 a longest edge; x <= d/2 by mirror symmetry of b).

Bounding: max(natural interval extension, mean-value form with interval forward-mode
gradients). Pruning: non-overlap (vertex strictly inside), edge < 2u, obtuse, c0c1 not
longest, x > d/2, circumradius > 2.

Rounding: numpy float64 with nextafter outward rounding per op, trig padded by 4e-16.
COST PILOT ONLY (numpy libm has no ulp guarantee). Measures boxes/s and box counts.
"""
import numpy as np, sys, time, json
sys.path.insert(0, 'src')
from iv_pilot import I, icos, isin, dn, up, hept_verts, strictly_inside, U, P7, TWO_U

PI = np.pi
RMIN = 2*U; RMAX = 4.0
ACRIT = 1.5326754446782214

# ---------- interval forward-mode AD (value + 6 partials) ----------
NV = 6
def izero(n): return I(np.zeros(n), np.zeros(n))
def iconst(c, n): return I(np.full(n, float(c)), np.full(n, float(c)))

class D:
    __slots__ = ('v', 'g')
    def __init__(s, v, g): s.v, s.g = v, g
    @staticmethod
    def var(iv, k, n):
        g = [izero(n) for _ in range(NV)]; g[k] = iconst(1.0, n); return D(iv, g)
    @staticmethod
    def const(c, n):
        return D(iconst(c, n), [izero(n) for _ in range(NV)])
    def __add__(s, o):
        if isinstance(o, D): return D(s.v + o.v, [a + b for a, b in zip(s.g, o.g)])
        return D(s.v + o, s.g)
    __radd__ = __add__
    def __sub__(s, o):
        if isinstance(o, D): return D(s.v - o.v, [a - b for a, b in zip(s.g, o.g)])
        return D(s.v - o, s.g)
    def __rsub__(s, o): return D(o - s.v, [-a for a in s.g])
    def __neg__(s): return D(-s.v, [-a for a in s.g])
    def __mul__(s, o):
        if isinstance(o, D):
            return D(s.v*o.v, [s.v*b + o.v*a for a, b in zip(s.g, o.g)])
        return D(s.v*o, [a*o for a in s.g])
    __rmul__ = __mul__

def iinv(a):
    assert (a.lo > 0).all() or (a.hi < 0).all() or True
    lo = np.where(a.lo > 0, dn(1/a.hi), np.where(a.hi < 0, dn(1/a.hi), -np.inf))
    hi = np.where(a.lo > 0, up(1/a.lo), np.where(a.hi < 0, up(1/a.lo), np.inf))
    return I(lo, hi)

def isqrt(a):
    return I(dn(np.sqrt(np.maximum(a.lo, 0))), up(np.sqrt(np.maximum(a.hi, 0))))

def iatan(a):
    return I(np.arctan(a.lo) - 4e-16, np.arctan(a.hi) + 4e-16)

def Dinv(s):
    iv = iinv(s.v); d = -(iv*iv)
    return D(iv, [d*a for a in s.g])

def Dsqrt(s):
    r = isqrt(s.v); d = iinv(r*2.0)
    return D(r, [d*a for a in s.g])

def Datan(s):
    v = iatan(s.v); d = iinv(s.v.sq() + 1.0)
    return D(v, [d*a for a in s.g])

def Dcos(s):
    c = icos(s.v); sn = isin(s.v)
    return D(c, [(-sn)*a for a in s.g])

def Dsin(s):
    c = icos(s.v); sn = isin(s.v)
    return D(sn, [c*a for a in s.g])

# ---------- transfer f ----------
class FIv:
    def __init__(self, path, DR=3, DF=2):
        self.c = np.load(path).reshape(DR+1, 2*DF+1, 2*DF+1); self.DR, self.DF = DR, DF
    def g(self, r, a, b, n):
        """sum_{k,ja,jb} c * T_k(x(r)) * Fa_ja(7a) * Fb_jb(7b), factorised."""
        x = (r*2.0 - (RMIN+RMAX))*(1.0/(RMAX-RMIN))
        Tr = [D.const(1.0, n), x]
        for k in range(2, self.DR+1): Tr.append(x*Tr[-1]*2.0 - Tr[-2])
        A7, B7 = a*7.0, b*7.0
        Fa = [D.const(1.0, n)]; Fb = [D.const(1.0, n)]
        for k in range(1, self.DF+1):
            Fa += [Dcos(A7*float(k)), Dsin(A7*float(k))]; Fb += [Dcos(B7*float(k)), Dsin(B7*float(k))]
        tot = None
        for k in range(self.DR+1):
            inner_a = None
            for ja in range(2*self.DF+1):
                lin = None
                for jb in range(2*self.DF+1):
                    cc = float(self.c[k, ja, jb])
                    if cc == 0.0: continue
                    t = Fb[jb]*cc
                    lin = t if lin is None else lin + t
                if lin is None: continue
                t = Fa[ja]*lin
                inner_a = t if inner_a is None else inner_a + t
            if inner_a is None: continue
            t = Tr[k]*inner_a
            tot = t if tot is None else tot + t
        return tot
    def f(self, r, a, b, n):
        p7 = PI/7
        return (self.g(r, a, b, n) - self.g(r, b + p7, a + p7, n)
                + self.g(r, (-b) + p7, (-a) + p7, n) - self.g(r, -a, -b, n))

def b_eval(F, B, use_D=True):
    """B (n,6,2) boxes [d,x,y,t0,t1,t2]. returns interval of b (natural, with gradient)"""
    n = len(B)
    X = [D.var(I(B[:, k, 0].copy(), B[:, k, 1].copy()), k, n) for k in range(6)]
    d, x, y, t0, t1, t2 = X
    area = d*y*0.5
    # edge 0->1: r = d, phi = 0
    e01 = (d, t0, t1)
    # edge 1->2: v = (x-d, y), phi = pi/2 + atan((d-x)/y)
    iy = Dinv(y)
    r12 = Dsqrt((x - d)*(x - d) + y*y)
    phi12 = Datan((d - x)*iy) + PI/2
    e12 = (r12, t1 - phi12, t2 - phi12)
    # edge 2->0: v = (-x, -y), phi = -pi/2 - atan(x/y)
    r20 = Dsqrt(x*x + y*y)
    phi20 = -(Datan(x*iy)) - PI/2
    e20 = (r20, t2 - phi20, t0 - phi20)
    tot = area
    for e in (e01, e12, e20):
        tot = tot + F.f(e[0], e[1], e[2], n)
    return tot

def lower_bound(F, B):
    n = len(B)
    full = b_eval(F, B)
    nat = full.v.lo
    # mean-value: b(c) + sum grad(B) * [-h, h]
    C = B.mean(2); Cb = np.stack([C, C], 2)
    bc = b_eval(F, Cb).v
    h = (B[:, :, 1] - B[:, :, 0])/2
    mv = bc.lo.copy()
    for k in range(6):
        gk = full.g[k]
        mv = mv - up(np.maximum(np.abs(gk.lo), np.abs(gk.hi))*h[:, k])
    return np.maximum(nat, mv)

def infeasible(B):
    """geometric pruning; returns mask of boxes certainly not valid configurations"""
    d = I(B[:, 0, 0], B[:, 0, 1]); x = I(B[:, 1, 0], B[:, 1, 1]); y = I(B[:, 2, 0], B[:, 2, 1])
    t = [I(B[:, 3+i, 0], B[:, 3+i, 1]) for i in range(3)]
    r02 = x.sq() + y.sq(); r12 = (x - d).sq() + y.sq(); dd = d.sq()
    m = r02.hi < TWO_U**2*(1-1e-15)
    m |= r12.hi < TWO_U**2*(1-1e-15)
    m |= r02.lo > dd.hi
    m |= r12.lo > dd.hi
    m |= ((-x)*(d - x) + y*y).hi < 0
    m |= dn(B[:, 1, 0]*2.0) > B[:, 0, 1]
    # circumradius > 2  <=>  d^2 r02 r12 > 64 area^2 = 16 d^2 y^2  <=> r02 r12 > 16 y^2
    m |= (r02*r12).lo > (y.sq()*16.0).hi
    zero = I(np.zeros(len(B)), np.zeros(len(B)))
    Cc = [(zero, zero), (d, zero), (x, y)]
    V = [hept_verts(Cc[i][0], Cc[i][1], t[i]) for i in range(3)]
    if SAT:
        for i, j in ((0, 1), (1, 2), (0, 2)):
            gap_max = None
            for owner in (i, j):
                for k in range(7):
                    a = t[owner] + ((2*k+1)*np.pi/7)
                    nx, ny = icos(a), isin(a)
                    proj_i = [nx*px + ny*py for (px, py) in V[i]]
                    proj_j = [nx*px + ny*py for (px, py) in V[j]]
                    # gap of j beyond i along n, and of i beyond j along -n
                    g1 = np.min([q.hi for q in proj_j], 0) - np.max([p.lo for p in proj_i], 0)
                    g2 = np.min([-p.lo for p in proj_i], 0) - np.max([-q.hi for q in proj_j], 0)
                    g = np.maximum(g1, g2)
                    gap_max = g if gap_max is None else np.maximum(gap_max, g)
            m |= gap_max < 0
    else:
        for i in range(3):
            for j in range(3):
                if i == j: continue
                for (px, py) in V[j]:
                    m |= strictly_inside(V[i], Cc[i][0], Cc[i][1], t[i], px, py)
    return m

DUMP = None
import os
SAT = os.environ.get('HEPT_SAT') == '1'
DUAL = os.environ.get('HEPT_DUAL') == '1'
def run(F, region, target, tlimit, minw=1e-6, batch=4000, log_every=30):
    lo = np.array(region['lo'], float); hi = np.array(region['hi'], float)
    stack = np.stack([lo, hi], -1)[None]
    if os.environ.get('HEPT_RESUME'):
        stack = np.load(os.environ['HEPT_RESUME']); print('resumed', len(stack), 'open boxes', flush=True)
    processed = 0; t0 = time.time(); last = t0; stuck = 0; minlb_open = []
    while len(stack):
        bt = stack[-batch:]; stack = stack[:-batch]
        inf = infeasible(bt)
        keep = bt[~inf]
        if len(keep):
            lb = lower_bound(F, keep)
            keep = keep[lb < target]
        if len(keep) and DUAL:
            lb2 = dual_lower_bound(F, keep)
            keep = keep[lb2 < target]
        processed += len(bt)
        if len(keep):
            w = keep[:, :, 1] - keep[:, :, 0]
            # scale angles and lengths alike
            j = np.argmax(w, 1)
            small = w[np.arange(len(keep)), j] < minw
            stuck += int(small.sum()); keep = keep[~small]; j = j[~small]
            mid = 0.5*(keep[np.arange(len(keep)), j, 0] + keep[np.arange(len(keep)), j, 1])
            L = keep.copy(); R = keep.copy()
            L[np.arange(len(keep)), j, 1] = mid; R[np.arange(len(keep)), j, 0] = mid
            stack = np.concatenate([stack, L, R]) if len(stack) else np.concatenate([L, R])
        if time.time() - last > log_every:
            el = time.time()-t0
            print(f'  processed {processed:,} open {len(stack):,} rate {processed/el:,.0f}/s [{el:.0f}s]', flush=True); last = time.time()
        if time.time()-t0 > tlimit:
            if DUMP: np.save(DUMP, stack)
            return dict(done=False, processed=processed, open=len(stack), secs=time.time()-t0, stuck=stuck)
    return dict(done=True, processed=processed, open=0, secs=time.time()-t0, stuck=stuck)

# ---------------- dual (Lagrangian) bounding near contact surfaces ----------------
PAIRS = ((0, 1), (1, 2), (0, 2))

def forced_halfaxes(B):
    """For each pair: index of the unique half-axis whose gap upper bound is >= 0 on the box
    (else -1). Half-axis code = owner_flag*14 + k*2 + side (owner_flag 0: owner = i, 1: owner = j;
    side 0: other on +n side, 1: other on -n side)."""
    n = len(B)
    d = I(B[:, 0, 0], B[:, 0, 1]); x = I(B[:, 1, 0], B[:, 1, 1]); y = I(B[:, 2, 0], B[:, 2, 1])
    t = [I(B[:, 3+q, 0], B[:, 3+q, 1]) for q in range(3)]
    zero = I(np.zeros(n), np.zeros(n))
    Cc = [(zero, zero), (d, zero), (x, y)]
    V = [hept_verts(Cc[q][0], Cc[q][1], t[q]) for q in range(3)]
    out = []
    for i, j in PAIRS:
        cnt = np.zeros(n, int); code = np.full(n, -1)
        for of, (o, oth) in enumerate(((i, j), (j, i))):
            for k in range(7):
                a = t[o] + ((2*k+1)*np.pi/7)
                nx, ny = icos(a), isin(a)
                po = [nx*px + ny*py for (px, py) in V[o]]
                pt = [nx*px + ny*py for (px, py) in V[oth]]
                g1 = np.min([q.hi for q in pt], 0) - np.max([p.lo for p in po], 0)
                ok = g1 >= 0     # (the reversed side is the same quantity; outward normals suffice for SAT)
                cnt += ok; code = np.where(ok, of*14 + k*2, code)
        out.append(np.where(cnt == 1, code, -1))
    return out

def h_D(X, pair, code, mvert):
    """D-valued valid inequality h >= 0 (per box arrays code, mvert)."""
    d, x, y, t0, t1, t2 = X; n = len(code)
    zero = D.const(0.0, n)
    C = [(zero, zero), (d, zero), (x, y)]; T = [t0, t1, t2]
    i, j = pair
    res = None
    for of, (o, oth) in enumerate(((i, j), (j, i))):
        k = (code % 14)//2; side = code % 2
        off_n = (2*k + 1)*np.pi/7; off_v = 2*np.pi*mvert/7
        a = T[o] + off_n; th = T[oth] + off_v
        qx = C[oth][0] + Dcos(th) - C[o][0]; qy = C[oth][1] + Dsin(th) - C[o][1]
        proj = Dcos(a)*qx + Dsin(a)*qy
        sgn = np.where(side == 0, 1.0, -1.0); rho = np.where(side == 0, U, 1.0)
        h = proj*I(sgn, sgn) - rho
        sel = (code // 14) == of
        res = h if res is None else Dwhere(sel, h, res)
    return res

def Dwhere(m, A, Bd):
    def iw(a, b): return I(np.where(m, a.lo, b.lo), np.where(m, a.hi, b.hi))
    return D(iw(A.v, Bd.v), [iw(a, b) for a, b in zip(A.g, Bd.g)])

def mid_vertex(Bc, pair, code):
    """vertex of the 'other' heptagon minimising the signed projection at the box centre."""
    n = len(Bc); i, j = pair
    d, x, y = Bc[:, 0], Bc[:, 1], Bc[:, 2]; T = [Bc[:, 3], Bc[:, 4], Bc[:, 5]]
    C = [(np.zeros(n), np.zeros(n)), (d, np.zeros(n)), (x, y)]
    of = code // 14; k = (code % 14)//2; side = code % 2
    o = np.where(of == 0, i, j); oth = np.where(of == 0, j, i)
    best = np.zeros(n, int); bv = np.full(n, np.inf)
    for m in range(7):
        val = np.zeros(n)
        for q in range(3):
            for r in range(3):
                sel = (o == q) & (oth == r)
                if not sel.any(): continue
                a = T[q][sel] + (2*k[sel]+1)*np.pi/7; th = T[r][sel] + 2*np.pi*m/7
                qx = C[r][0][sel] + np.cos(th) - C[q][0][sel]; qy = C[r][1][sel] + np.sin(th) - C[q][1][sel]
                pr = np.cos(a)*qx + np.sin(a)*qy
                val[sel] = np.where(side[sel] == 0, pr, -pr)
        upd = val < bv; best = np.where(upd, m, best); bv = np.minimum(bv, val)
    return best

def dual_lower_bound(F, B):
    n = len(B)
    codes = forced_halfaxes(B)
    Xb = [D.var(I(B[:, k, 0].copy(), B[:, k, 1].copy()), k, n) for k in range(6)]
    C = B.mean(2); Cb = np.stack([C, C], 2)
    Xc = [D.var(I(Cb[:, k, 0].copy(), Cb[:, k, 1].copy()), k, n) for k in range(6)]
    bB = b_eval(F, B); bC = b_eval(F, Cb)
    gb = np.stack([0.5*(g.lo+g.hi) for g in bC.g], 1)
    hs = []
    for p, pair in enumerate(PAIRS):
        code = codes[p]
        if (code >= 0).sum() == 0: hs.append(None); continue
        cc = np.where(code >= 0, code, 0)
        mv = mid_vertex(C, pair, cc)
        hB = h_D(Xb, pair, cc, mv); hC = h_D(Xc, pair, cc, mv)
        hs.append((code >= 0, hB, hC))
    # multipliers mu_p >= 0 by projecting grad b onto grad h (sequential nonnegative fit)
    phiB = bB; phiC_lo = bC.v.lo.copy(); resid = gb.copy()
    mus = []
    for item in hs:
        if item is None: continue
        act, hB, hC = item
        gh = np.stack([0.5*(g.lo+g.hi) for g in hC.g], 1)
        mu = np.where(act, np.maximum(0.0, (resid*gh).sum(1)/np.maximum((gh*gh).sum(1), 1e-30)), 0.0)
        resid = resid - mu[:, None]*gh
        phiB = phiB - hB*I(mu, mu)
        phiC_lo = phiC_lo - up(np.maximum(mu*hC.v.hi, mu*hC.v.lo))   # subtract upper bound of mu*h(c)
    h = (B[:, :, 1] - B[:, :, 0])/2
    mvb = phiC_lo.copy()
    for k in range(6):
        gk = phiB.g[k]
        mvb = mvb - up(np.maximum(np.abs(gk.lo), np.abs(gk.hi))*h[:, k])
    return np.maximum(phiB.v.lo, mvb)


if __name__ == '__main__':
    cpath = sys.argv[1]; target = float(sys.argv[2]); tl = float(sys.argv[3])
    reg = json.loads(sys.argv[4])
    F = FIv(cpath)
    if len(sys.argv) > 5: DUMP = sys.argv[5]
    r = run(F, reg, target, tl)
    r.update(target=target, region=reg, rate=r['processed']/max(r['secs'], 1e-9))
    print(json.dumps(r))



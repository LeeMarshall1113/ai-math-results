"""Rigorous interval branch-and-bound for the main nonobtuse inequality

    b(T) = area(T) + f(e01) + f(e12) + f(e20) >= TARGET

over canonical nonobtuse triangles (c0 = 0, c1 = (d,0), c2 = (x,y), orientations t0..t2;
c0c1 a longest edge; x <= d/2 by the mirror symmetry of b).

All bounds use rig.R (outward-rounded numpy intervals with self-made rigorous trig).
Bounding ladder for each box (cheap to expensive):
  1. geometric pruning (edges < 2u, obtuse, c0c1 not longest, x > d/2, circumradius > 2,
     interval separating-axis overlap test);
  2. natural interval extension and first-order mean-value form (forward AD, order 1);
  3. second-order Taylor form: b(c) + grad b(c).delta + 1/2 delta^T H(B) delta (AD order 2);
  4. Lagrangian dual: as 3. for b - sum mu_p h_p, where h_p >= 0 is the vertex-on-edge
     inequality forced by feasibility when a pair can be separated by only one SAT axis on B.
"""
import numpy as np, sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rig import R, rcos, rsin, ratan, PIO2, PI, dn, up

U_ = R(np.float64(0.9009688679024191), np.nextafter(np.float64(0.9009688679024191), np.inf))  # cos(pi/7)
_ul = 0.9009688679024191
U_ = R(np.nextafter(_ul, -np.inf), np.nextafter(_ul, np.inf))
TWO_U_HI = float(np.nextafter(2*_ul, np.inf)) * (1 + 1e-15)
PI7 = PI*(1.0/7.0)          # interval enclosure of pi/7
P7I = PI*(2.0/7.0)          # 2pi/7
NV = 6

# ------------------------- AD (order 1 or 2) over R -------------------------
def _z(n): return R(np.zeros(n))
class A:
    __slots__ = ('v', 'g', 'H', 'n', 'o')
    def __init__(s, v, g, H, n, o): s.v, s.g, s.H, s.n, s.o = v, g, H, n, o
    @staticmethod
    def var(iv, k, n, o):
        g = [_z(n) for _ in range(NV)]; g[k] = R(np.ones(n))
        H = {(i, j): _z(n) for i in range(NV) for j in range(i, NV)} if o == 2 else None
        return A(iv, g, H, n, o)
    @staticmethod
    def const(c, n, o):
        cv = c if isinstance(c, R) else R(np.full(n, float(c)))
        H = {(i, j): _z(n) for i in range(NV) for j in range(i, NV)} if o == 2 else None
        return A(cv, [_z(n) for _ in range(NV)], H, n, o)
    def _lift(s, o): return o if isinstance(o, A) else A.const(o, s.n, s.o)
    def __add__(s, o):
        if not isinstance(o, A): return A(s.v + o, s.g, s.H, s.n, s.o)
        H = {k: s.H[k] + o.H[k] for k in s.H} if s.o == 2 else None
        return A(s.v + o.v, [a + b for a, b in zip(s.g, o.g)], H, s.n, s.o)
    __radd__ = __add__
    def __neg__(s):
        H = {k: -h for k, h in s.H.items()} if s.o == 2 else None
        return A(-s.v, [-a for a in s.g], H, s.n, s.o)
    def __sub__(s, o):
        if not isinstance(o, A): return A(s.v - o, s.g, s.H, s.n, s.o)
        return s + (-o)
    def __rsub__(s, o): return (-s) + o
    def scale(s, c):
        """c: float, array or R (constant)"""
        H = {k: h*c for k, h in s.H.items()} if s.o == 2 else None
        return A(s.v*c, [a*c for a in s.g], H, s.n, s.o)
    def __mul__(s, o):
        if not isinstance(o, A): return s.scale(o)
        v = s.v*o.v
        g = [s.v*b + o.v*a for a, b in zip(s.g, o.g)]
        H = None
        if s.o == 2:
            H = {}
            for (i, j) in s.H:
                H[(i, j)] = s.H[(i, j)]*o.v + o.H[(i, j)]*s.v + s.g[i]*o.g[j] + s.g[j]*o.g[i]
        return A(v, g, H, s.n, s.o)
    __rmul__ = __mul__
    def unary(s, f0, f1, f2):
        g = [f1*a for a in s.g]; H = None
        if s.o == 2:
            H = {(i, j): f2*s.g[i]*s.g[j] + f1*s.H[(i, j)] for (i, j) in s.H}
        return A(f0, g, H, s.n, s.o)

def Acos(u):
    c, sn = rcos(u.v), rsin(u.v)
    return u.unary(c, -sn, -c)
def Asin(u):
    c, sn = rcos(u.v), rsin(u.v)
    return u.unary(sn, c, -sn)
def Ainv(u):
    iv = u.v.inv(); return u.unary(iv, -(iv.sq()), iv.sq()*iv*2.0)
def Asqrt(u):
    r = u.v.sqrt(); ir = r.inv()
    return u.unary(r, ir*0.5, -(ir*ir*ir)*0.25)
def Aatan(u):
    d = (u.v.sq() + 1.0).inv()
    return u.unary(ratan(u.v), d, -(u.v*d.sq())*2.0)

# ------------------------- transfer f -------------------------
RMIN_ = U_*2.0; RMAX = 4.0
class FTrans:
    def __init__(s, path, DR, DF, mu=None, mu_idx=9):
        s.c = np.load(path).astype(float).reshape(DR+1, 2*DF+1, 2*DF+1); s.DR, s.DF = DR, DF
        s.cR = {}
        if mu is not None:     # exactification: coefficient mu_idx becomes the interval c + [mu_lo, mu_hi]
            k, ja, jb = np.unravel_index(mu_idx, s.c.shape)
            base = s.c[k, ja, jb]
            s.cR[(k, ja, jb)] = R(base) + R(mu[0], mu[1])
    def g(s, r, a, b):
        n, o = r.n, r.o
        # x = (2r - RMIN - RMAX)/(RMAX - RMIN), with RMIN an interval
        den = (R(np.full(n, RMAX)) - RMIN_).inv()
        x = (r.scale(2.0) - (RMIN_ + RMAX)).scale(den)
        Tr = [A.const(1.0, n, o), x]
        for k in range(2, s.DR+1): Tr.append((x*Tr[-1]).scale(2.0) - Tr[-2])
        def mults(t):
            c1, s1 = Acos(t), Asin(t)
            C, S = [A.const(1.0, n, o), c1], [None, s1]
            for k in range(2, s.DF+1):
                C.append(C[-1]*c1 - S[-1]*s1); S.append(S[-1]*c1 + C[-2]*s1)
            out = [C[0]]
            for k in range(1, s.DF+1): out += [C[k], S[k]]
            return out
        Fa = mults(a.scale(7.0)); Fb = mults(b.scale(7.0))
        tot = None
        for k in range(s.DR+1):
            ia = None
            for ja in range(2*s.DF+1):
                lin = None
                for jb in range(2*s.DF+1):
                    cc = s.cR.get((k, ja, jb), s.c[k, ja, jb])
                    if not isinstance(cc, R) and cc == 0.0: continue
                    t = Fb[jb].scale(cc); lin = t if lin is None else lin + t
                if lin is None: continue
                t = Fa[ja]*lin; ia = t if ia is None else ia + t
            if ia is None: continue
            t = Tr[k]*ia; tot = t if tot is None else tot + t
        return tot
    def f(s, r, a, b):
        return (s.g(r, a, b) - s.g(r, b + PI7, a + PI7) + s.g(r, (-b) + PI7, (-a) + PI7) - s.g(r, -a, -b))

def b_A(F, Bx, o):
    """Bx: list of 6 R (box or point). returns A of b with AD order o."""
    n = len(Bx[0].lo)
    X = [A.var(Bx[k], k, n, o) for k in range(6)]
    d, x, y, t0, t1, t2 = X
    area = (d*y).scale(0.5)
    iy = Ainv(y)
    r12 = Asqrt((x - d)*(x - d) + y*y)
    phi12 = Aatan((d - x)*iy) + PIO2
    r20 = Asqrt(x*x + y*y)
    phi20 = -(Aatan(x*iy)) - PIO2
    tot = area + F.f(d, t0, t1) + F.f(r12, t1 - phi12, t2 - phi12) + F.f(r20, t2 - phi20, t0 - phi20)
    return tot

# ------------------------- geometry -------------------------
def hverts(cx, cy, t):
    V = []
    for k in range(7):
        a = t + P7I*float(k)
        V.append((cx + rcos(a), cy + rsin(a)))
    return V

def sat_axis_gaps(B):
    """per pair: list over the 14 outward edge normals of upper bounds of the gap."""
    n = len(B)
    d = R(B[:, 0, 0], B[:, 0, 1]); x = R(B[:, 1, 0], B[:, 1, 1]); y = R(B[:, 2, 0], B[:, 2, 1])
    t = [R(B[:, 3+q, 0], B[:, 3+q, 1]) for q in range(3)]
    z = R(np.zeros(n)); Cc = [(z, z), (d, z), (x, y)]
    V = [hverts(Cc[q][0], Cc[q][1], t[q]) for q in range(3)]
    out = {}
    for i, j in ((0, 1), (1, 2), (0, 2)):
        gl = []
        for (o, oth) in ((i, j), (j, i)):
            for k in range(7):
                a = t[o] + PI7*float(2*k+1)
                nx, ny = rcos(a), rsin(a)
                po = [nx*px + ny*py for (px, py) in V[o]]
                pt = [nx*px + ny*py for (px, py) in V[oth]]
                gl.append(np.min([q.hi for q in pt], 0) - np.max([p.lo for p in po], 0))
        out[(i, j)] = gl
    return out

def infeasible(B, gaps=None):
    d = R(B[:, 0, 0], B[:, 0, 1]); x = R(B[:, 1, 0], B[:, 1, 1]); y = R(B[:, 2, 0], B[:, 2, 1])
    r02 = x.sq() + y.sq(); r12 = (x - d).sq() + y.sq(); dd = d.sq()
    E2 = dn(dn(U_.lo*2)**2)
    m = r02.hi < E2                    # edge c0c2 < 2u for sure
    m |= r12.hi < E2
    m |= r02.lo > dd.hi
    m |= r12.lo > dd.hi
    if OBT is None:
        m |= ((-x)*(d - x) + y*y).hi < 0
    m |= dn(B[:, 1, 0]*2.0) > B[:, 0, 1]
    m |= (r02*r12).lo > (y.sq()*16.0).hi
    if gaps is None: gaps = sat_axis_gaps(B)
    for key, gl in gaps.items():
        m |= np.max(gl, 0) < 0
    return m, gaps

# ------------------------- bounds -------------------------
def box_R(B): return [R(B[:, k, 0], B[:, k, 1]) for k in range(6)]
def ctr_R(B):
    c = B.mean(2); return [R(c[:, k]) for k in range(6)], c

def half_w(B):
    """outward-rounded half-widths about the float centre c = B.mean(2):
    h_k >= max(hi_k - c_k, c_k - lo_k) (audit FLAW-1)."""
    c = B.mean(2)
    return up(np.maximum(up(B[:, :, 1] - c), up(c - B[:, :, 0])))

def lb_first(F, B):
    full = b_A(F, box_R(B), 1)
    Cr, _ = ctr_R(B); bc = b_A(F, Cr, 1).v
    h = half_w(B)
    mv = bc.lo.copy()
    for k in range(6): mv = dn(mv - up(full.g[k].mag()*h[:, k]))
    return np.maximum(full.v.lo, mv)

def taylor2_lb(valc, gradc, Hbox, h):
    """lower bound of v(c) + g.delta + 1/2 delta^T H delta over |delta_i| <= h_i."""
    lb = valc.lo.copy()
    for k in range(6): lb = dn(lb - up(gradc[k].mag()*h[:, k]))
    for (i, j), Hij in Hbox.items():
        if i == j:
            lb = dn(lb + dn(0.5*np.minimum(Hij.lo, 0.0)*up(h[:, i]*h[:, i])))
        else:
            lb = dn(lb - up(Hij.mag()*up(h[:, i]*h[:, j])))
    return lb

def lb_second(F, B):
    Cr, _ = ctr_R(B)
    c1 = b_A(F, Cr, 1); H2 = b_A(F, box_R(B), 2)
    h = half_w(B)
    return taylor2_lb(c1.v, c1.g, H2.H, h), c1, H2

DUAL2 = os.environ.get('HEPT_DUAL2') == '1'
STUCK_DUMP = os.environ.get('HEPT_STUCK'); STUCK_LIST = []
STOP_REQ = [False]       # graceful stop: SIGTERM -> dump open boxes and exit as TIMEOUT
def _on_term(signum, frame): STOP_REQ[0] = True
try:
    import signal as _sg
    _sg.signal(_sg.SIGTERM, _on_term)
except Exception:
    pass
OBT = None   # obtuse-transfer table (set from HEPT_OBT); see obtuse_extra()

class ObtTable:
    def __init__(self, Dpath, binpath):
        b = json.load(open(binpath)); self.rb = np.array(b['rb']); self.cb = np.array(b['cb'])
        self.D = np.load(Dpath); self.RTH = self.rb[0]
        self.rev = np.maximum.accumulate(self.D[:, ::-1], axis=1)[:, ::-1]   # rev[i,j] = max_{j'>=j} D[i,j']
        self.nr, self.nc = self.D.shape
    def rbins(self, r):
        i0 = np.clip(np.searchsorted(self.rb, r.lo, side='right')-1, 0, self.nr-1)
        i1 = np.clip(np.searchsorted(self.rb, r.hi, side='right')-1, 0, self.nr-1)
        return i0, i1
    def dhat_ub(self, r, ca):
        # upper bound of max D over obtuse neighbours compatible with (edge length r, receiver cos ca)
        n = len(r.lo); out = np.zeros(n)
        act = (r.hi >= self.RTH) & (ca.hi > 0)
        if not act.any(): return out
        i0, i1 = self.rbins(r)
        jj = np.clip(np.searchsorted(self.cb[1:], -ca.hi, side='left'), 0, self.nc-1)
        for i in range(self.nr):
            sel = act & (i0 <= i) & (i1 >= i)
            if sel.any(): out = np.where(sel, np.maximum(out, self.rev[i, jj]), out)
        return out
    def d_lb(self, r, cg):
        # lower bound of D(bin(r), bin(cos gamma)) over the box
        n = len(r.lo); out = np.full(n, np.inf)
        i0, i1 = self.rbins(r)
        j0 = np.clip(np.searchsorted(self.cb, cg.lo, side='right')-1, 0, self.nc-1)
        j1 = np.clip(np.searchsorted(self.cb, cg.hi, side='right')-1, 0, self.nc-1)
        for i in range(self.nr):
            for j in range(self.nc):
                sel = (i0 <= i) & (i1 >= i) & (j0 <= j) & (j1 >= j)
                if sel.any(): out = np.where(sel, np.minimum(out, self.D[i, j]), out)
        return np.where(np.isfinite(out), out, 0.0)


class ObtTable4:
    """Continuous (bilinear-on-nodes) obtuse-transfer table D4.
    Nodes Dn[i, j] at (rb[i], cb[j]); in cell (i, j) D is the bilinear interpolant of its 4 nodes.
    Bilinear functions attain extremes over a sub-rectangle at its corners, so corner evaluation
    (padded by 1e-12 for float rounding) gives rigorous bounds."""
    def __init__(self, Dpath, binpath):
        b = json.load(open(binpath)); self.rb = np.array(b['rb']); self.cb = np.array(b['cb'])
        self.Dn = np.load(Dpath); self.RTH = self.rb[0]
        self.nr, self.nc = len(self.rb)-1, len(self.cb)-1
        self.revn = np.maximum.accumulate(self.Dn[:, ::-1], axis=1)[:, ::-1]   # revn[i,j] = max_{j'>=j} Dn[i,j']
    def _cell(self, grid, v, n):
        return np.clip(np.searchsorted(grid, v, side='right')-1, 0, n-1)
    def _bil(self, i, j, t, s_):
        D00 = self.Dn[i, j]; D10 = self.Dn[i+1, j]; D01 = self.Dn[i, j+1]; D11 = self.Dn[i+1, j+1]
        return (1-t)*(1-s_)*D00 + t*(1-s_)*D10 + (1-t)*s_*D01 + t*s_*D11
    def _rect(self, rlo, rhi, clo, chi, want_max):
        """exact bound of the bilinear table over rectangles [rlo,rhi]x[clo,chi] (arrays), cell by cell."""
        n = len(rlo); out = np.full(n, -np.inf if want_max else np.inf)
        valid = (rhi >= rlo) & (chi >= clo)
        for i in range(self.nr):
            ra, rb_ = self.rb[i], self.rb[i+1]
            mi = valid & (rhi >= ra) & (rlo <= rb_)
            if not mi.any(): continue
            tl = np.clip((np.maximum(rlo, ra) - ra)/(rb_ - ra), 0, 1); th = np.clip((np.minimum(rhi, rb_) - ra)/(rb_ - ra), 0, 1)
            for j in range(self.nc):
                ca_, cb_ = self.cb[j], self.cb[j+1]
                m = mi & (chi >= ca_) & (clo <= cb_)
                if not m.any(): continue
                if not (self.Dn[i:i+2, j:j+2] != 0).any():
                    vals = np.zeros(n)
                else:
                    sl = np.clip((np.maximum(clo, ca_) - ca_)/(cb_ - ca_), 0, 1); sh = np.clip((np.minimum(chi, cb_) - ca_)/(cb_ - ca_), 0, 1)
                    v4 = [self._bil(i, j, tl, sl), self._bil(i, j, tl, sh), self._bil(i, j, th, sl), self._bil(i, j, th, sh)]
                    vals = np.maximum.reduce(v4) if want_max else np.minimum.reduce(v4)
                out = np.where(m, np.maximum(out, vals) if want_max else np.minimum(out, vals), out)
        return out
    def dhat_ub(self, r, ca):
        """upper bound of the charge across an edge of length r whose opposite angle has cosine ca:
        max of D over r' in r and c in [-ca.hi, 0] (Delaunay-compatible obtuse neighbours)."""
        n = len(r.lo)
        act = (r.hi >= self.RTH) & (ca.hi > 0)
        if not act.any(): return np.zeros(n)
        rlo = np.maximum(r.lo, self.RTH); rhi = np.minimum(r.hi, self.rb[-1])
        clo = np.maximum(-ca.hi, self.cb[0]); chi = np.full(n, self.cb[-1])
        v = self._rect(rlo, rhi, clo, chi, True)
        return np.where(act & np.isfinite(v), np.maximum(v, 0.0) + 1e-12, 0.0)
    def d_lb(self, r, cg):
        """lower bound of D(r, cos gamma) over the box (obtuse part, cg <= 0)."""
        n = len(r.lo)
        rlo = np.maximum(r.lo, self.rb[0]); rhi = np.minimum(r.hi, self.rb[-1])
        clo = np.maximum(cg.lo, self.cb[0]); chi = np.minimum(cg.hi, self.cb[-1])
        v = self._rect(rlo, rhi, clo, chi, False)
        return np.where(np.isfinite(v), np.maximum(v - 1e-12, 0.0), 0.0)

def _cos_at(P, Q1, Q2):
    ux, uy = Q1[0] - P[0], Q1[1] - P[1]; vx, vy = Q2[0] - P[0], Q2[1] - P[1]
    num = ux*vx + uy*vy
    den = ((ux.sq() + uy.sq())*(vx.sq() + vy.sq())).sqrt()
    return num*den.inv()

def obtuse_extra(B):
    """per-box upper bound of the obtuse-transfer correction to subtract from b - acrit."""
    n = len(B)
    d = R(B[:, 0, 0], B[:, 0, 1]); x = R(B[:, 1, 0], B[:, 1, 1]); y = R(B[:, 2, 0], B[:, 2, 1])
    z = R(np.zeros(n)); c0 = (z, z); c1 = (d, z); c2 = (x, y)
    r01 = d; r12 = ((x - d).sq() + y.sq()).sqrt(); r20 = (x.sq() + y.sq()).sqrt()
    cg = _cos_at(c2, c0, c1)          # angle at c2, opposite edge 01 (the longest edge)
    ca12 = _cos_at(c0, c1, c2)        # angle at c0, opposite edge 12
    ca20 = _cos_at(c1, c2, c0)        # angle at c1, opposite edge 20
    inc_short = OBT.dhat_ub(r12, ca12) + OBT.dhat_ub(r20, ca20)
    certain_obt = cg.hi < 0
    certain_non = cg.lo >= 0
    cg_obt = R(cg.lo, np.minimum(cg.hi, 0.0)); cg_non = R(np.maximum(cg.lo, 0.0), cg.hi)
    extra_obt = inc_short - OBT.d_lb(r01, cg_obt)
    extra_non = inc_short + OBT.dhat_ub(r01, cg_non)
    return np.where(certain_obt, extra_obt, np.where(certain_non, extra_non, np.maximum(extra_obt, extra_non)))
CHARTS = None
PAIRS = ((0, 1), (1, 2), (0, 2))
def h_A(X, pair, oflag, k, m):
    """vertex m of 'other' on outward-normal side of edge k of 'owner' (per-box arrays)."""
    d, x, y, t0, t1, t2 = X; n = d.n; o_ = d.o
    zero = A.const(0.0, n, o_); C = [(zero, zero), (d, zero), (x, y)]; T = [t0, t1, t2]
    i, j = pair; res = None
    for of, (own, oth) in enumerate(((i, j), (j, i))):
        offn = PI7*((2*k+1).astype(float)); offv = P7I*(m.astype(float))
        a = T[own] + offn; th = T[oth] + offv
        qx = C[oth][0] + Acos(th) - C[own][0]; qy = C[oth][1] + Asin(th) - C[own][1]
        hval = Acos(a)*qx + Asin(a)*qy - U_
        sel = oflag == of
        res = hval if res is None else Awhere(sel, hval, res)
    return res

def Rwhere(m, a, b): return R(np.where(m, a.lo, b.lo), np.where(m, a.hi, b.hi))
def Awhere(m, P, Q):
    H = {k: Rwhere(m, P.H[k], Q.H[k]) for k in P.H} if P.o == 2 else None
    return A(Rwhere(m, P.v, Q.v), [Rwhere(m, a, b) for a, b in zip(P.g, Q.g)], H, P.n, P.o)

def lb_dual(F, B, gaps, c1, H2):
    n = len(B); h = half_w(B)
    Cr, cpt = ctr_R(B)
    gb = np.stack([0.5*(g.lo + g.hi) for g in c1.g], 1)
    Xc = [A.var(Cr[k], k, n, 1) for k in range(6)]
    Xb = [A.var(R(B[:, k, 0], B[:, k, 1]), k, n, 2) for k in range(6)]
    phic, phiB = c1, H2; resid = gb.copy(); used = False
    for pair in PAIRS:
        gl = np.array(gaps[pair])                       # (14, n)
        ok = gl >= 0; cnt = ok.sum(0); forced = cnt == 1
        if not forced.any(): continue
        ax = np.argmax(ok, 0); oflag = ax//7; k = ax % 7
        # vertex of 'other' minimising the projection at the centre (float choice only)
        m = _min_vertex(cpt, pair, oflag, k)
        hc = h_A(Xc, pair, oflag, k, m); hB = h_A(Xb, pair, oflag, k, m)
        gh = np.stack([0.5*(g.lo + g.hi) for g in hc.g], 1)
        mu = np.where(forced, np.maximum(0.0, (resid*gh).sum(1)/np.maximum((gh*gh).sum(1), 1e-30)), 0.0)
        if not (mu > 0).any(): continue
        resid = resid - mu[:, None]*gh; used = True
        muR = R(mu)
        phic = phic - hc.scale(muR); phiB = phiB - hB.scale(muR)
    if not used: return np.full(n, -np.inf)
    return taylor2_lb(phic.v, phic.g, phiB.H, h)

def _min_vertex(c, pair, oflag, k):
    n = len(c); i, j = pair
    C = [(np.zeros(n), np.zeros(n)), (c[:, 0], np.zeros(n)), (c[:, 1], c[:, 2])]; T = [c[:, 3], c[:, 4], c[:, 5]]
    own = np.where(oflag == 0, i, j); oth = np.where(oflag == 0, j, i)
    best = np.zeros(n, int); bv = np.full(n, np.inf)
    for mm in range(7):
        val = np.full(n, np.inf)
        for q in range(3):
            for r in range(3):
                sel = (own == q) & (oth == r)
                if not sel.any(): continue
                a = T[q][sel] + (2*k[sel]+1)*np.pi/7; th = T[r][sel] + 2*np.pi*mm/7
                qx = C[r][0][sel] + np.cos(th) - C[q][0][sel]; qy = C[r][1][sel] + np.sin(th) - C[q][1][sel]
                val[sel] = np.cos(a)*qx + np.sin(a)*qy
        upd = val < bv; best = np.where(upd, mm, best); bv = np.minimum(bv, val)
    return best

# ------------------------- driver -------------------------
def run(F, region, target, tlimit, minw=1e-7, batch=3000, log_every=60, dump=None, use2=True):
    lo = np.array(region['lo'], float); hi = np.array(region['hi'], float)
    stack = np.stack([lo, hi], -1)[None]
    if os.environ.get('HEPT_RESUME'): stack = np.load(os.environ['HEPT_RESUME'])
    st = dict(processed=0, infeasible=0, first=0, second=0, dual=0, split=0, stuck=0)
    t0 = time.time(); last = t0
    ckpt = os.environ.get('HEPT_CKPT'); ck_every = float(os.environ.get('HEPT_CKPT_EVERY', '600')); last_ck = t0
    while len(stack):
        if ckpt and time.time() - last_ck > ck_every:
            np.save(ckpt + '.tmp.npy', stack); os.replace(ckpt + '.tmp.npy', ckpt); last_ck = time.time()
            print(f'  checkpoint {len(stack):,} open boxes -> {ckpt}', flush=True)
        bt = stack[-batch:]; stack = stack[:-batch]; st['processed'] += len(bt)
        inf, gaps = infeasible(bt); st['infeasible'] += int(inf.sum())
        tg = np.full(len(bt), target) if OBT is None else up(target + up(obtuse_extra(bt)))
        if CHARTS is not None:
            ex = CHARTS.inside(bt) & ~inf; st['chart'] = st.get('chart', 0) + int(ex.sum()); inf = inf | ex
        keep = bt[~inf]; tg = tg[~inf]; gaps = {k: [g[~inf] for g in v] for k, v in gaps.items()}
        if len(keep):
            lb = lb_first(F, keep); ok = lb >= tg; st['first'] += int(ok.sum())
            keep = keep[~ok]; tg = tg[~ok]; gaps = {k: [g[~ok] for g in v] for k, v in gaps.items()}
        if len(keep) and use2:
            lb2, c1, H2 = lb_second(F, keep); ok = lb2 >= tg; st['second'] += int(ok.sum())
            rest = ~ok
            if rest.any():
                sub = keep[rest]; tgs = tg[rest]; g2 = {k: [g[rest] for g in v] for k, v in gaps.items()}
                c1s = _subA(c1, rest); H2s = _subA(H2, rest)
                lb3 = lb_dual(F, sub, g2, c1s, H2s)
                if DUAL2:
                    lb3 = np.maximum(lb3, lb_dual2(F, sub, c1s, H2s))
                ok3 = lb3 >= tgs; st['dual'] += int(ok3.sum())
                keep = sub[~ok3]
            else:
                keep = keep[:0]
        if len(keep):
            w = keep[:, :, 1] - keep[:, :, 0]; jj = np.argmax(w, 1)
            small = w[np.arange(len(keep)), jj] < minw
            st['stuck'] += int(small.sum())
            if small.any() and STUCK_DUMP is not None: STUCK_LIST.append(keep[small][:200])
            keep = keep[~small]; jj = jj[~small]
            mid = 0.5*(keep[np.arange(len(keep)), jj, 0] + keep[np.arange(len(keep)), jj, 1])
            L = keep.copy(); Rr = keep.copy()
            L[np.arange(len(keep)), jj, 1] = mid; Rr[np.arange(len(keep)), jj, 0] = mid
            st['split'] += len(keep)
            stack = np.concatenate([stack, L, Rr]) if len(stack) else np.concatenate([L, Rr])
        if time.time() - last > log_every:
            el = time.time() - t0
            print(f'  {st} open {len(stack):,} rate {st["processed"]/el:,.0f}/s [{el:.0f}s]', flush=True); last = time.time()
        if time.time() - t0 > tlimit or STOP_REQ[0]:
            if dump: np.save(dump, stack)
            if STUCK_DUMP and STUCK_LIST: np.save(STUCK_DUMP, np.concatenate(STUCK_LIST))
            st.update(done=False, open=len(stack), secs=time.time()-t0); return st
    if STUCK_DUMP and STUCK_LIST: np.save(STUCK_DUMP, np.concatenate(STUCK_LIST))
    st.update(done=True, open=0, secs=time.time()-t0); return st

def _subA(a, m):
    H = {k: a.H[k][m] for k in a.H} if a.o == 2 else None
    return A(a.v[m], [g[m] for g in a.g], H, int(m.sum()), a.o)

# ---- dual bounding v2: vertex-in-slab inequalities (valid whenever the vertex is strictly
# inside the other six half-planes of the owner on the whole box: then h < 0 => overlap) ----
def _cand_float(cpt, pair):
    """float centre values: for owner/other orientation, each vertex m of other: nearest edge k
    of owner and h value."""
    n = len(cpt); i, j = pair
    Cx = [np.zeros(n), cpt[:, 0], cpt[:, 1]]; Cy = [np.zeros(n), np.zeros(n), cpt[:, 2]]; T = [cpt[:, 3], cpt[:, 4], cpt[:, 5]]
    out = []
    for of, (own, oth) in enumerate(((i, j), (j, i))):
        for m in range(7):
            th = T[oth] + 2*np.pi*m/7
            qx = Cx[oth] + np.cos(th) - Cx[own]; qy = Cy[oth] + np.sin(th) - Cy[own]
            best_k = np.zeros(n, int); best_v = np.full(n, -np.inf)
            for k in range(7):
                a = T[own] + (2*k+1)*np.pi/7
                v = np.cos(a)*qx + np.sin(a)*qy - _ul
                upd = v > best_v; best_k = np.where(upd, k, best_k); best_v = np.maximum(best_v, v)
            out.append((of, m, best_k, best_v))
    return out

def _slab_ok(B, pair, of, k, m):
    """interval check on box: vertex m of other strictly inside the 6 half-planes l != k of owner."""
    n = len(B); i, j = pair
    d = R(B[:, 0, 0], B[:, 0, 1]); x = R(B[:, 1, 0], B[:, 1, 1]); y = R(B[:, 2, 0], B[:, 2, 1])
    t = [R(B[:, 3+q, 0], B[:, 3+q, 1]) for q in range(3)]
    z = R(np.zeros(n)); Cc = [(z, z), (d, z), (x, y)]
    own = np.where(of == 0, i, j) if isinstance(of, np.ndarray) else (i if of == 0 else j)
    oth = (j if of == 0 else i)
    th = t[oth] + P7I*(m.astype(float))
    qx = Cc[oth][0] + rcos(th) - Cc[own][0]; qy = Cc[oth][1] + rsin(th) - Cc[own][1]
    ok = np.ones(n, bool)
    for l in range(1, 7):
        kk = (k + l) % 7
        a = t[own] + PI7*((2*kk+1).astype(float))
        v = rcos(a)*qx + rsin(a)*qy - U_
        ok &= v.hi < 0
    return ok

def lb_dual2(F, B, c1, H2, hmax=0.03):
    n = len(B); h = half_w(B)
    Cr, cpt = ctr_R(B)
    gb = np.stack([0.5*(g.lo + g.hi) for g in c1.g], 1)
    Xc = [A.var(Cr[q], q, n, 1) for q in range(6)]
    Xb = [A.var(R(B[:, q, 0], B[:, q, 1]), q, n, 2) for q in range(6)]
    phic, phiB = c1, H2; resid = gb.copy(); used = False
    for pair in PAIRS:
        for (of, m, k, hv) in _cand_float(cpt, pair):
            act = np.abs(hv) < hmax
            if not act.any(): continue
            mm = np.full(n, m)
            act &= _slab_ok(B, pair, of, k, mm)
            if not act.any(): continue
            ofa = np.full(n, of)
            hc = h_A(Xc, pair, ofa, k, mm); hB = h_A(Xb, pair, ofa, k, mm)
            gh = np.stack([0.5*(g.lo + g.hi) for g in hc.g], 1)
            mu = np.where(act, np.maximum(0.0, (resid*gh).sum(1)/np.maximum((gh*gh).sum(1), 1e-30)), 0.0)
            if not (mu > 0).any(): continue
            resid = resid - mu[:, None]*gh; used = True
            muR = R(mu); phic = phic - hc.scale(muR); phiB = phiB - hB.scale(muR)
    if not used: return np.full(n, -np.inf)
    return taylor2_lb(phic.v, phic.g, phiB.H, h)


def _provenance(paths):
    import hashlib, platform
    print('python', sys.version.replace(chr(10), ' '), '| platform', platform.platform(), flush=True)
    print('numpy', np.__version__, np.__file__, flush=True)
    for p_ in paths:
        if p_ and os.path.exists(p_):
            print('sha256', p_, hashlib.sha256(open(p_, 'rb').read()).hexdigest(), flush=True)
    for k_ in sorted(os.environ):
        if k_.startswith('HEPT_'): print('env', k_, '=', os.environ[k_], flush=True)

if __name__ == '__main__':
    cpath, DR, DF = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    _here = os.path.dirname(os.path.abspath(__file__))
    _files = [os.path.join(_here, 'rig.py'), os.path.join(_here, 'bnb2.py'), os.path.join(_here, 'chart_excl.py'), cpath,
              os.path.join(os.path.dirname(cpath), 'dl_ref.json'), os.environ.get('HEPT_RESUME')]
    for _k in ('HEPT_OBT', 'HEPT_OBT4'):
        if os.environ.get(_k): _files += os.environ[_k].split(',')
    if os.environ.get('HEPT_CHARTS'): _files.append(os.environ['HEPT_CHARTS'])
    _provenance(_files)
    print('argv', json.dumps(sys.argv[1:]), flush=True)
    target = float(sys.argv[4]); tl = float(sys.argv[5]); reg = json.loads(sys.argv[6])
    dump = sys.argv[7] if len(sys.argv) > 7 else None
    mu = json.loads(os.environ['HEPT_MU']) if os.environ.get('HEPT_MU') else None
    F = FTrans(cpath, DR, DF, mu=mu)
    if os.environ.get('HEPT_OBT'):
        OBT = ObtTable(*os.environ['HEPT_OBT'].split(','))
    if os.environ.get('HEPT_OBT4'):
        OBT = ObtTable4(*os.environ['HEPT_OBT4'].split(','))
    if os.environ.get('HEPT_CHARTS'):
        from chart_excl import Charts
        CHARTS = Charts(os.environ['HEPT_CHARTS'], os.path.join(os.path.dirname(cpath), 'dl_ref.json'))
    r = run(F, reg, target, tl, dump=dump)
    r.update(target=target, region=reg, rate=r['processed']/max(r['secs'], 1e-9))
    print(json.dumps(r))



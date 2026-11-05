"""Extend the explicit local neighbourhood along the flex curve (tube of charts), ball arithmetic.

Chart at a certified point z_c on the contact curve C = {g = 0}:
  zeta = J_c (z - z_c), J_c = [DG(z_c); v], Phi(z) = (g(z), v.z).
  D_c = [0, eta]^5 x [s_c - eta, s_c + eta], eta = (1 - q) rho.
Certified on the zeta-box of radius rho (subdivided n^6):
  (A) touching vertices stay inside the other 6 half-planes (so g_k < 0 forces overlap);
  (B) ||DPhi J_c^{-1} - I||_inf <= q < 1 (Banach: Phi^{-1}(D_c) exists, unique, inside chart);
  (C) w = grad b . DPhi^{-1} has w_1..w_5 > 0 (b increases as contacts open);
  (D) beta'' >= B2 on the chart, and beta(s_c), beta'(s_c) enclosed at z_c.
Then for every non-overlapping z in the chart with Phi(z) in D_c:
  b(z) >= min_{|d| <= eta} beta(s_c) + beta'(s_c) d + B2 d^2 / 2  =: m_c.
For the central chart (s_c = 0) beta'(0) = 0 exactly and m_c = acrit.
"""
import sys, itertools, time, json
sys.path.insert(0, 'src')
import numpy as np
from flint import arb, arb_series, arb_mat, ctx
import local_radius_arb as LR
from local_radius_arb import b_series, g_series, halfplanes_other, CONTACTS, coef, ORD
ctx.prec = 128

P, LAM, V, J0, _ = LR.setup(arb(13)/100)
import os
FPATH = os.environ.get('HEPT_F')
if FPATH:
    from f_arb import FArb, exactify
    _F0 = FArb(FPATH)
    _F, MU, MU_IDX, BASE = exactify(_F0, V)
    print('exactify: grad b . v before =', BASE.str(8), ' mu =', MU.str(8), ' on basis', MU_IDX, flush=True)
    def BFUN(z, p=None): return _F.b(z)
else:
    def BFUN(z, p=None): return b_series(z, p)
Vf = np.array([float(x.mid()) for x in V])

def zser_c(zc, Minv, cen, rad, dz=None):
    zeta = [arb(cen[j], rad[j]) for j in range(6)]
    out = []
    for q in range(6):
        c0 = zc[q] + sum((Minv[q, j]*zeta[j] for j in range(6)), arb(0))
        c1 = arb(0) if dz is None else sum((Minv[q, j]*dz[j] for j in range(6)), arb(0))
        out.append(arb_series([c0, c1], prec=ORD))
    return out

def jac_at(zc_float):
    z = [arb_series([arb(zc_float[q]), arb(0)], prec=ORD) for q in range(6)]
    J = arb_mat(6, 6)
    for j in range(6):
        e = [arb(0)]*6; e[j] = arb(1)
        zz = [arb_series([arb(zc_float[q]), e[q]], prec=ORD) for q in range(6)]
        for kk, c in enumerate(CONTACTS): J[kk, j] = coef(g_series(zz, c), 1)
    for q in range(6): J[5, q] = V[q]
    return J

def curve_point(s, z_init):
    """float Newton for g(z) = 0, v.z = s"""
    z = z_init.copy()
    for _ in range(30):
        zz = [arb_series([arb(z[q]), arb(0)], prec=ORD) for q in range(6)]
        F = np.array([float(coef(g_series(zz, c), 0).mid()) for c in CONTACTS] + [Vf @ z - s])
        J = np.array([[float(x.mid()) for x in row] for row in jac_at(z).tolist()])
        dz = np.linalg.solve(J, -F); z = z + dz
        if np.abs(dz).max() < 1e-15: break
    return z

def chart(zc_float, s_c, rho, n, central=False):
    zc = [arb(x) for x in zc_float]
    Jc = jac_at(zc_float); Minv = Jc.inv()
    # tight certified curve point near zc: Banach on tiny box, target (0, s_c) -> shift
    # compute Phi(zc) exactly (balls); the true curve point z* = zc + Minv*zeta*, |zeta*| <= tiny
    zz = [arb_series([zc[q], arb(0)], prec=ORD) for q in range(6)]
    Phi_c = [coef(g_series(zz, c), 0) for c in CONTACTS] + [sum((V[q]*zc[q] for q in range(6)), arb(0)) - s_c]
    off = max(float(abs(x).upper()) for x in Phi_c)
    tiny = max(4*off, 1e-30)
    h = 2*rho/n; centers = [-rho + h*(i+0.5) for i in range(n)]
    HP = -1e9; Q = 0.0; W = [1e9]*5; B2 = 1e9
    for cell in itertools.product(range(n), repeat=6):
        cen = [centers[i] for i in cell]; rad = [h/2*(1+1e-9)]*6     # inflated: sub-boxes overlap, no seams (audit FLAW-3)
        z0 = zser_c(zc, Minv, cen, rad)
        for c in CONTACTS:
            for val in halfplanes_other(z0, c): HP = max(HP, float(coef(val, 0).upper()))
        gb = []; DP = arb_mat(6, 6)
        for j in range(6):
            e = [arb(0)]*6; e[j] = arb(1)
            z = zser_c(zc, Minv, cen, rad, e)
            gb.append(coef(BFUN(z, P), 1))
            for kk, c in enumerate(CONTACTS): DP[kk, j] = coef(g_series(z, c), 1)
            DP[5, j] = sum((V[q]*Minv[q, j] for q in range(6)), arb(0))
        for i in range(6):
            srow = arb(0)
            for j in range(6): srow += abs(DP[i, j] - (1 if i == j else 0))
            Q = max(Q, float(srow.upper()))
        rhs = arb_mat(6, 1)
        for j in range(6): rhs[j, 0] = gb[j]
        w = DP.transpose().solve(rhs)
        W = [min(W[k], float(w[k, 0].lower())) for k in range(5)]
        e6 = arb_mat(6, 1); e6[5, 0] = arb(1)
        zp = DP.solve(e6); zpl = [zp[j, 0] for j in range(6)]
        d2 = arb_mat(6, 1)
        for kk, c in enumerate(CONTACTS):
            d2[kk, 0] = -2*coef(g_series(zser_c(zc, Minv, cen, rad, zpl), c), 2)
        zpp = DP.solve(d2)
        bpp = 2*coef(BFUN(zser_c(zc, Minv, cen, rad, zpl), P), 2) + sum((gb[j]*zpp[j, 0] for j in range(6)), arb(0))
        B2 = min(B2, float(bpp.lower()))
    # beta(s_c), beta'(s_c) on the tiny certified box around the curve point
    cen0 = [0.0]*6; rad0 = [tiny]*6
    zt = zser_c(zc, Minv, cen0, rad0)
    beta_c = coef(BFUN(zt, P), 0)
    gb = []; DP = arb_mat(6, 6)
    for j in range(6):
        e = [arb(0)]*6; e[j] = arb(1)
        z = zser_c(zc, Minv, cen0, rad0, e)
        gb.append(coef(BFUN(z, P), 1))
        for kk, c in enumerate(CONTACTS): DP[kk, j] = coef(g_series(z, c), 1)
        DP[5, j] = sum((V[q]*Minv[q, j] for q in range(6)), arb(0))
    rhs = arb_mat(6, 1)
    for j in range(6): rhs[j, 0] = gb[j]
    w = DP.transpose().solve(rhs); betap = w[5, 0]
    # audit FLAW-3: Banach needs |t - Phi(z_c)| <= (1-q) rho for every target t in D_c. With the chart
    # window centred at Phi(z_c) (g in [0, eta], s within eta of v.z_c) and |Phi(z_c) - (0, s_c)| <= off,
    # eta = (1-q) rho - 2 off (rounded down) is sufficient.
    eta = float(np.nextafter((1-Q)*rho - 2*off - 1e-30, -np.inf)) if Q < 1 else 0.0
    ok = Q < 1 and HP < 0 and min(W) > 0 and B2 > 0
    acrit = coef(BFUN([arb_series([arb(0), arb(0)], prec=ORD)]*6, P), 0)
    if central:
        m = float(acrit.lower())  # exact: beta'(0) = 0, beta'' > 0  => min = acrit
        excess = 0.0
    else:
        bl = float((beta_c - acrit).lower()); bp = float(abs(betap).upper())
        # min over |d|<=eta of bl - bp|d| + B2 d^2/2  (lower bound)
        cands = [0.0, eta] + ([bp/B2] if B2 > 0 and bp/B2 < eta else [])
        excess = min(bl - bp*d + 0.5*B2*d*d for d in cands)
        m = None
    Jlo = [[float(Jc[i, j].lower()) for j in range(6)] for i in range(6)]
    Jhi = [[float(Jc[i, j].upper()) for j in range(6)] for i in range(6)]
    return dict(off=off, zc=[float(x) for x in zc_float], Jlo=Jlo, Jhi=Jhi, s_c=s_c, rho=rho, q=Q, eta=eta, hp=HP, w_min=min(W), B2=B2, ok=ok,
                beta_minus_acrit_c=float((beta_c-acrit).mid()), betap_c=float(betap.mid()), excess_lower=excess,
                tiny=tiny)

if __name__ == '__main__':
    rho = float(sys.argv[1]); n = int(sys.argv[2]); smax = float(sys.argv[3])
    t0 = time.time(); out = []
    c0 = chart(np.zeros(6), 0.0, rho, n, central=True); out.append(c0)
    print(json.dumps(c0), flush=True)
    for sign in (+1, -1):
        s = 0.0; z = np.zeros(6); last_eta = c0['eta']
        while abs(s) < smax:
            step = 0.9*last_eta   # overlap with previous chart's s-range
            s_new = s + sign*step
            z = curve_point(s_new, z)
            r = chart(z, s_new, rho, n)
            print(json.dumps(r), flush=True); out.append(r)
            if not r['ok']: break
            s = s_new; last_eta = min(r['eta'], last_eta) if r['eta'] > 0 else last_eta
            # each chart covers [s - eta, s + eta]; next center spaced by 0.9*eta keeps overlap
    json.dump(out, open(f'runs/tube_rho{rho}_n{n}' + ('_f' if FPATH else '') + '.json', 'w'))
    print('secs', round(time.time()-t0, 1))

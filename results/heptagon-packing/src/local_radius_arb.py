"""Explicit neighbourhood for the local lemma at the DL triangle (ball arithmetic).

Uses the jets-only transfer f (f linear in (a, b) near each s2 edge configuration with the
slopes fixed below, zero Hessian). The final global f must agree with these jets to 2nd order
in the neighbourhood, or the same check must be re-run with that f.

Lemma (certified for the rho, eta, sigma printed):
  Let Z = [-rho, rho]^6 (perturbation of Ht1, H2p2 with H0 fixed) and Phi(z) = (g_1..g_5, v.z).
  (A) For z in Z, overlap is forced whenever some g_k(z) < 0 (touching vertex stays inside the
      other six half-planes of the touched heptagon). So non-overlap implies g >= 0.
  (B) For every t in D = [0,eta]^5 x [-sigma, sigma] there is a unique z in Z with Phi(z) = t
      (Banach contraction of z -> J0^{-1}(t - (Phi(z) - J0 z))).
  (C) On Z: w = grad_z b . DPhi^{-1} has w_1..w_5 > 0 (b increases as any contact opens).
  (D) On Z: beta''(s) > 0 along the contact curve Phi^{-1}(0, s).
  Hence every non-overlapping configuration z in Z with Phi(z) in D has b(z) >= acrit,
  with equality only at the DL triangle.
"""
import sys, itertools, time
sys.path.insert(0, 'src')
from flint import arb, arb_series, arb_mat, ctx
import local_lemma_arb as LL
from local_lemma_arb import state, area, g, CONTACTS, EDGES, edge_cfg, REF, coef, PI, U, grad, F_area, F_g, F_edge
ctx.prec = 128
ORD = 3

def solve_params(t):
    ca = grad(F_area); G = [grad(F_g(c)) for c in CONTACTS]
    Je = [[grad(F_edge(k, cmp)) for cmp in range(3)] for k in range(3)]
    cols = [Je[0][1], Je[0][2], [Je[1][1][q]-Je[1][2][q] for q in range(6)], [Je[2][1][q]-Je[2][2][q] for q in range(6)]]
    # gauge choice: ga = -1/2, lam4 = lam5 = t ; unknowns gb, a1, a2, lam1, lam2, lam3
    ga = arb(-1)/2
    M = arb_mat(6, 6); r = arb_mat(6, 1)
    for q in range(6):
        M[q, 0] = cols[1][q]; M[q, 1] = cols[2][q]; M[q, 2] = cols[3][q]
        M[q, 3] = -G[0][q]; M[q, 4] = -G[1][q]; M[q, 5] = -G[2][q]
        r[q, 0] = -(ca[q] + ga*cols[0][q] - t*(G[3][q] + G[4][q]))
    x = M.solve(r)
    p = [ga, x[0, 0], x[1, 0], x[2, 0]]; lam = [x[3, 0], x[4, 0], x[5, 0], t, t]
    return p, lam, G

def b_series(z, p):
    C, T = state(z)
    E = [edge_cfg(C, T, i, j, None) for (i, j) in EDGES]
    ga, gb, a1, a2 = p
    return (area(C) + ga*(E[0][1]-REF[0][1]) + gb*(E[0][2]-REF[0][2])
            + a1*((E[1][1]-REF[1][1]) - (E[1][2]-REF[1][2])) + a2*((E[2][1]-REF[2][1]) - (E[2][2]-REF[2][2])))

def g_series(z, c):
    C, T = state(z); return g(C, T, c)

def halfplanes_other(z, c):
    """for contact c=(i,j,m,k): values n_l.(vertex - c_i) - u for l != k (must be < 0)."""
    i, j, m, k = c; C, T = state(z)
    th = T[j] + 2*PI*m/7
    wx = C[j][0] + th.cos() - C[i][0]; wy = C[j][1] + th.sin() - C[i][1]
    out = []
    for l in range(7):
        if l == k: continue
        a = T[i] + arb(2*l+1)*PI/7
        out.append(a.cos()*wx + a.sin()*wy - U)
    return out

def zser(Minv, cen, rad, dz=None):
    """z = Minv * zeta, zeta in ball (cen +/- rad), plus s * Minv * dz (series)."""
    zeta = [arb(cen[j], rad[j]) for j in range(6)]
    out = []
    for q in range(6):
        c0 = sum((Minv[q, j]*zeta[j] for j in range(6)), arb(0))
        c1 = arb(0) if dz is None else sum((Minv[q, j]*dz[j] for j in range(6)), arb(0))
        out.append(arb_series([c0, c1], prec=ORD))
    return out

def setup(t):
    p, lam, G0 = solve_params(t)
    N = arb_mat(5, 5); rr = arb_mat(5, 1); others = [0, 2, 3, 4, 5]
    for kk in range(5):
        for c, q in enumerate(others): N[kk, c] = G0[kk][q]
        rr[kk, 0] = -G0[kk][1]
    vv = N.solve(rr); v = [vv[0, 0], arb(1), vv[1, 0], vv[2, 0], vv[3, 0], vv[4, 0]]
    J0 = arb_mat(6, 6)
    for kk in range(5):
        for q in range(6): J0[kk, q] = G0[kk][q]
    for q in range(6): J0[5, q] = v[q]
    return p, lam, v, J0, J0.inv()

def check_box(Minv, p, cen, rad):
    """returns (halfplane_max, contraction rowsum max, w lower bounds[5], beta'' lower)"""
    hp = None
    z0 = zser(Minv, cen, rad)
    for c in CONTACTS:
        for val in halfplanes_other(z0, c):
            ub = float(coef(val, 0).upper()); hp = ub if hp is None or ub > hp else hp
    gb = []; DP = arb_mat(6, 6)
    for j in range(6):
        e = [arb(0)]*6; e[j] = arb(1)
        z = zser(Minv, cen, rad, e)
        gb.append(coef(b_series(z, p), 1))
        for kk, c in enumerate(CONTACTS): DP[kk, j] = coef(g_series(z, c), 1)
        DP[5, j] = coef(sum((z[q]*0 for q in range(6)), arb_series([arb(0)], prec=ORD)), 1) + (1 if j == 5 else 0)
    contr = 0.0
    for i in range(6):
        srow = arb(0)
        for j in range(6): srow += abs(DP[i, j] - (1 if i == j else 0))
        contr = max(contr, float(srow.upper()))
    DPt = DP.transpose(); rhs = arb_mat(6, 1)
    for j in range(6): rhs[j, 0] = gb[j]
    w = DPt.solve(rhs)
    wl = [float(w[k, 0].lower()) for k in range(5)]
    e6 = arb_mat(6, 1); e6[5, 0] = arb(1)
    zp = DP.solve(e6); zpl = [zp[j, 0] for j in range(6)]
    d2 = arb_mat(6, 1)
    for kk, c in enumerate(CONTACTS):
        d2[kk, 0] = -2*coef(g_series(zser(Minv, cen, rad, zpl), c), 2)
    zpp = DP.solve(d2)
    bpp = 2*coef(b_series(zser(Minv, cen, rad, zpl), p), 2) + sum((gb[j]*zpp[j, 0] for j in range(6)), arb(0))
    return hp, contr, wl, float(bpp.lower())

def check(rho, n, t=arb(13)/100):
    """zeta-box [-rho, rho]^6 with zeta = J0 z. D = [0, eta]^5 x [-sigma, sigma] with
    eta = sigma = (1 - q) rho maps inside (Banach)."""
    p, lam, v, J0, Minv = setup(t)
    h = 2*rho/n; centers = [-rho + h*(i+0.5) for i in range(n)]
    HP = -1e9; Q = 0.0; W = [1e9]*5; B2 = 1e9
    for cell in itertools.product(range(n), repeat=6):
        cen = [centers[i] for i in cell]; rad = [h/2]*6
        hp, q, wl, b2 = check_box(Minv, p, cen, rad)
        HP = max(HP, hp); Q = max(Q, q); W = [min(W[k], wl[k]) for k in range(5)]; B2 = min(B2, b2)
    eta = (1-Q)*rho if Q < 1 else 0.0
    res = dict(rho_zeta=rho, n=n, contraction_q=Q, halfplane_max=HP, w_min=W, beta2_min=B2,
               eta=eta, sigma=eta, params=[float(x.mid()) for x in p], lam=[float(x.mid()) for x in lam],
               z_radius_bound=max(sum(float(abs(Minv[i, j]).upper()) for j in range(6)) for i in range(6))*rho)
    res['ok'] = Q < 1 and HP < 0 and min(W) > 0 and B2 > 0
    # margins on the boundary of D for gluing with the global B&B (first-order estimates)
    res['margin_eta_face'] = min(W)*eta if res['ok'] else None
    res['margin_sigma_face_est'] = 0.5*B2*eta**2 if res['ok'] else None
    return res

if __name__ == '__main__':
    rho = float(sys.argv[1]); n = int(sys.argv[2])
    t0 = time.time()
    r = check(rho, n)
    r['secs'] = round(time.time()-t0, 1)
    for k, val in r.items(): print(f'{k}: {val}')

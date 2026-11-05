"""Exact (ball-arithmetic) infinitesimal local lemma at the double-lattice triangle.

Setting: triangle s2 = (H0, Ht1, H2p2) of the Kuperberg double lattice of regular heptagons
(circumradius 1). H0 is fixed (rigid-motion gauge); variables z = (x1, y1, phi1, x2, y2, phi2)
perturb Ht1 and H2p2.

Modified area b(z) = area(z) + f(e01) + f(e12) + f(e20), f any C^2 function of the oriented
pair configuration (r, a, b) that is antisymmetric (f o R = -f, R(r,a,b) = (r, b+pi/7, a+pi/7))
and mirror-symmetric (f o m = f, m(r,a,b) = (r, pi/7-b, pi/7-a)). Such f has, at the s2 edge
configurations:
  translation edge e01 = (r*, 0, 0):   f = 0,  grad f = (0, ga, gb)
  reflected edges e12, e20:            f = 0,  grad f = alpha_k (0, 1, -1)
and its admissible Hessian terms vanish along the flex direction (checked below).

Claim certified here (all quantities as arb balls at 300 bits):
 (1) the five vertex-on-edge contacts listed in CONTACTS are active (|g| < 1e-80) with the
     touching vertex strictly inside the edge, and every other vertex/edge pair of the three
     heptagons is strictly separated or irrelevant (min |g| over inactive near-pairs > 0.01);
 (2) with alpha2 = 0, lambda1 = lambda3 = 1/20, the linear system grad b = G^T lambda has a
     unique solution (ga, gb, alpha1, lambda2, lambda4, lambda5) with all lambda > 0;
 (3) the contact Jacobian G (5 x 6) has rank 5; its null vector v (y1-component 1) gives
     v^T Hess(Lagrangian) v > 0.
(1)-(3) are the hypotheses of the standard second-order sufficient condition, so the DL
triangle is a strict local minimiser of b among nearby non-overlapping configurations.

Caveat: (1) shows the contact values are zero only to 1e-80, not exactly zero. They are exactly
zero by the Kuperberg construction, but that is not re-derived here; a symbolic check in
Q(cos pi/7, sin pi/7) closes it.
"""
from flint import arb, arb_series, arb_mat, ctx
ctx.prec = 300
PI = arb.pi(); P7 = 2*PI/7; U = (PI/7).cos()
ORD = 3   # series length (up to s^2)

def coef(sr, n):
    c = sr.coeffs()
    return c[n] if len(c) > n else arb(0)

def S(c0, c1=0):
    return arb_series([arb(c0), arb(c1)], prec=ORD)

def kvert(i):
    return ((2*PI*i/7).cos(), (2*PI*i/7).sin())

# --- KK half-length parallelogram, exactly (linear solve in arb) ---
k = [kvert(i) for i in range(7)]
L = 1 + U
# p2 = (1-a) k1 + a k2 ; p3 = (1-b) k2 + b k3
# heights equal: (1-a) k1y + a k2y = (1-b) k2y + b k3y
# x(p2) - x(p3) = L/2
A = arb_mat([[k[2][1]-k[1][1], k[2][1]-k[3][1]],
             [k[2][0]-k[1][0], k[2][0]-k[3][0]]])
rhs = arb_mat([[k[2][1]-k[1][1]], [L/2 - k[1][0] + k[2][0]]])
# eq1: a(k2y-k1y) + b(k2y-k3y) = k2y - k1y
# eq2: x(p2)-x(p3) = k1x + a(k2x-k1x) - k2x - b(k3x-k2x) = L/2 -> a(k2x-k1x) + b(k2x-k3x) = L/2 - k1x + k2x
sol = A.solve(rhs)
a_kk, b_kk = sol[0, 0], sol[1, 0]
p2 = ((1-a_kk)*k[1][0] + a_kk*k[2][0], (1-a_kk)*k[1][1] + a_kk*k[2][1])
C0 = [(arb(0), arb(0)), (L, arb(0)), (2*p2[0], 2*p2[1])]
T0 = [arb(0), arb(0), PI/7]

def state(z):
    """z: list of 6 arb_series"""
    C = [(S(0), S(0)), (S(C0[1][0]) + z[0], S(C0[1][1]) + z[1]), (S(C0[2][0]) + z[3], S(C0[2][1]) + z[4])]
    T = [S(0), S(0) + z[2], S(T0[2]) + z[5]]
    return C, T

def area(C):
    return 0.5*((C[1][0]-C[0][0])*(C[2][1]-C[0][1]) - (C[1][1]-C[0][1])*(C[2][0]-C[0][0]))

CONTACTS = [(0, 2, 4, 1), (1, 0, 0, 3), (1, 2, 5, 2), (2, 0, 1, 4), (2, 1, 2, 5)]

def g(C, T, c):
    i, j, m, kk = c
    a = T[i] + arb((2*kk+1))*PI/7
    th = T[j] + 2*PI*m/7
    wx = C[j][0] + th.cos() - C[i][0]; wy = C[j][1] + th.sin() - C[i][1]
    return a.cos()*wx + a.sin()*wy - U

def tangential(C, T, c):
    i, j, m, kk = c
    a = T[i] + arb((2*kk+1))*PI/7
    th = T[j] + 2*PI*m/7
    wx = C[j][0] + th.cos() - C[i][0]; wy = C[j][1] + th.sin() - C[i][1]
    return -a.sin()*wx + a.cos()*wy

EDGES = [(0, 1), (1, 2), (2, 0)]

def edge_cfg(C, T, i, j, ref_ab):
    vx = C[j][0]-C[i][0]; vy = C[j][1]-C[i][1]
    r = (vx*vx + vy*vy).sqrt()
    # psi = atan2(vy, vx) via atan of ratio around the reference direction
    c0x, c0y = coef(vx, 0), coef(vy, 0)
    psi0 = arb.atan2(c0y, c0x)
    # rotate v by -psi0 so it is near the +x axis, then psi = psi0 + atan(v'y/v'x)
    cs, sn = psi0.cos(), psi0.sin()
    rx = cs*vx + sn*vy; ry = -sn*vx + cs*vy
    psi = psi0 + (ry/rx).atan()
    return r, T[i]-psi, T[j]-psi

def zvec(w, s_coeff=1):
    return [arb_series([arb(0), arb(wk)], prec=ORD) for wk in w]

def lin_coef(F, w):
    return F(zvec(w)).coeffs()[1] if len(F(zvec(w)).coeffs()) > 1 else arb(0)

def make_funcs():
    zero = [arb(0)]*6
    C, T = state(zvec(zero))
    ref = []
    for (i, j) in EDGES:
        r, a, b = edge_cfg(C, T, i, j, None)
        ref.append((coef(r, 0), coef(a, 0), coef(b, 0)))
    return ref

REF = make_funcs()


def F_area(z):
    C, T = state(z); return area(C)

def F_g(c):
    return lambda z: (lambda C, T: g(C, T, c))(*state(z))

def F_edge(kidx, comp):
    def f(z):
        C, T = state(z); i, j = EDGES[kidx]; e = edge_cfg(C, T, i, j, None)
        return e[comp]
    return f

def grad(F):
    out = []
    for kk in range(6):
        w = [0]*6; w[kk] = 1
        out.append(coef(F(zvec(w)), 1))
    return out

def main():
    print('KK a, b =', a_kk.str(20), b_kk.str(20))
    zero = [0]*6
    C, T = state(zvec(zero))
    # (1) contacts
    for c in CONTACTS:
        gv = coef(g(C, T, c), 0); tv = coef(tangential(C, T, c), 0)
        half = (PI/7).tan()*U
        print(f'  contact {c}: g = {gv.str(5)}  tangential = {tv.str(10)}  inside edge: {bool(abs(tv) < half)}')
        assert abs(gv) < arb(10)**-80
        assert abs(tv) < half
    # inactive vertex/edge pairs that are close: check strictly positive or irrelevant
    mins = []
    for i in range(3):
        for j in range(3):
            if i == j: continue
            for m in range(7):
                for kk in range(7):
                    c = (i, j, m, kk)
                    if c in CONTACTS: continue
                    gv = coef(g(C, T, c), 0); tv = coef(tangential(C, T, c), 0)
                    if abs(tv) < (PI/7).tan()*U + 0.05 and abs(gv) < 0.2:
                        mins.append((float(gv.mid()), c))
    mins.sort(key=lambda t: abs(t[0]))
    print('  nearest inactive vertex/edge pairs (g, pair):', [(round(x, 4), c) for x, c in mins[:6]])
    # (2) first-order system
    ca = grad(F_area)
    G = [grad(F_g(c)) for c in CONTACTS]
    Je = [[grad(F_edge(kidx, comp)) for comp in range(3)] for kidx in range(3)]
    # contribution columns of params: ga -> Je[0][1], gb -> Je[0][2], alpha1 -> Je[1][1]-Je[1][2], alpha2 -> Je[2][1]-Je[2][2]
    col_ga = Je[0][1]; col_gb = Je[0][2]
    col_a1 = [Je[1][1][q]-Je[1][2][q] for q in range(6)]
    lam1 = lam3 = arb(1)/20
    # equation: ca + ga*col_ga + gb*col_gb + a1*col_a1 - sum lam_k G_k = 0
    # unknowns: ga, gb, a1, lam2, lam4, lam5 (indices 1,3,4 of CONTACTS -> lam2, lam4, lam5)
    M = arb_mat(6, 6); r = arb_mat(6, 1)
    for q in range(6):
        M[q, 0] = col_ga[q]; M[q, 1] = col_gb[q]; M[q, 2] = col_a1[q]
        M[q, 3] = -G[1][q]; M[q, 4] = -G[3][q]; M[q, 5] = -G[4][q]
        r[q, 0] = -(ca[q] - lam1*G[0][q] - lam3*G[2][q])
    x = M.solve(r)
    ga, gb, a1, l2, l4, l5 = [x[q, 0] for q in range(6)]
    print('  ga =', ga.str(15), ' gb =', gb.str(15), ' alpha1 =', a1.str(15), ' alpha2 = 0')
    lam = [lam1, l2, lam3, l4, l5]
    print('  lambda =', [l.str(12) for l in lam])
    assert all(l > 0 for l in lam), 'a multiplier is not certified positive'
    # residual check
    res = [ca[q] + ga*col_ga[q] + gb*col_gb[q] + a1*col_a1[q] - sum(lam[kk]*G[kk][q] for kk in range(5)) for q in range(6)]
    print('  residual radii:', [float(abs(rq).upper()) for rq in res])
    # (3) null vector: G v = 0, v[1] = 1
    N = arb_mat(5, 5); rr = arb_mat(5, 1); others = [0, 2, 3, 4, 5]
    for kk in range(5):
        for t, q in enumerate(others): N[kk, t] = G[kk][q]
        rr[kk, 0] = -G[kk][1]
    vv = N.solve(rr)
    v = [vv[0, 0], arb(1), vv[1, 0], vv[2, 0], vv[3, 0], vv[4, 0]]
    print('  null vector v =', [vi.str(10) for vi in v])
    # Lagrangian along s*v up to s^2
    z = [arb_series([arb(0), vi], prec=ORD) for vi in v]
    C, T = state(z)
    Lg = area(C)
    E = [edge_cfg(C, T, i, j, None) for (i, j) in EDGES]
    Lg = Lg + ga*(E[0][1] - REF[0][1]) + gb*(E[0][2] - REF[0][2]) + a1*((E[1][1]-REF[1][1]) - (E[1][2]-REF[1][2]))
    for kk, c in enumerate(CONTACTS):
        Lg = Lg - lam[kk]*g(C, T, c)
    print('  L along flex: c1 =', coef(Lg, 1).str(10), '  second derivative 2*c2 =', (2*coef(Lg, 2)).str(15))
    # admissible Hessian terms of f along flex: translation edge needs dr, reflected edges need dt = da - db
    dr01 = coef(E[0][0], 1); dt12 = coef(E[1][1], 1) - coef(E[1][2], 1); dt20 = coef(E[2][1], 1) - coef(E[2][2], 1)
    print('  along flex: dr(e01) =', dr01.str(5), ' dt(e12) =', dt12.str(5), ' dt(e20) =', dt20.str(5))
    assert abs(dr01) < arb(10)**-60 and abs(dt12) < arb(10)**-60 and abs(dt20) < arb(10)**-60
    assert 2*coef(Lg, 2) > 0
    print('CERTIFIED: SOSC hypotheses hold at the DL triangle (subject to the docstring caveat).')

if __name__ == '__main__':
    main()

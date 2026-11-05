"""Kill test: is there an antisymmetric pair-transfer f with
   area(T) + f(e01) + f(e12) + f(e20) >= acrit  for all nonobtuse Delaunay triangles T?
Cutting-plane LP over a finite basis. Non-rigorous (float, sampled)."""
import numpy as np, sys, time, json
sys.path.insert(0, 'src')
from vec import *
from scipy.optimize import linprog, minimize
ACRIT = 1.5326754446782214; RMIN = 2*U; RMAX = 4.0
DR, DF = int(sys.argv[1]), int(sys.argv[2]); ROUNDS = int(sys.argv[3]); TL = float(sys.argv[4])
NONOB = (sys.argv[5] == 'nonobtuse') if len(sys.argv) > 5 else True
rng = np.random.default_rng(11)
BOX = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
NW, NR = 50, 25
ETA = float(sys.argv[7]) if len(sys.argv) > 7 else 1e-3

def gbasis(E):
    """E (n,3) -> (n,nb) basis values of g(r,a,b)"""
    r, a, b = E[:, 0], E[:, 1]*7, E[:, 2]*7
    x = (2*r-RMIN-RMAX)/(RMAX-RMIN)
    Tr = [np.ones_like(x), x]
    for k in range(2, DR+1): Tr.append(2*x*Tr[-1]-Tr[-2])
    Fa = [np.ones_like(a)] + [f(k*a) for k in range(1, DF+1) for f in (np.cos, np.sin)]
    Fb = [np.ones_like(b)] + [f(k*b) for k in range(1, DF+1) for f in (np.cos, np.sin)]
    return np.stack([t*fa*fb for t in Tr for fa in Fa for fb in Fb], 1)

def rev(E):
    return np.stack([E[:, 0], np.mod(E[:, 2]+np.pi/7, P7), np.mod(E[:, 1]+np.pi/7, P7)], 1)

def fbasis(E):
    mir = np.stack([E[:, 0], np.pi/7 - E[:, 2], np.pi/7 - E[:, 1]], 1)
    Sm = np.stack([E[:, 0], -E[:, 1], -E[:, 2]], 1)
    return gbasis(E) - gbasis(rev(E)) + gbasis(mir) - gbasis(Sm)

def rows(C, T):
    a, R, dmin, s, E = tri_feats(C, T)
    M = sum(fbasis(e) for e in E)
    return a, M

def feasible(C, T):
    a, R, dmin, s, E = tri_feats(C, T)
    ok = (a > 0) & (R <= 2) & (s >= -1e-9)
    if NONOB: ok &= dmin >= 0
    return ok

def sample(n):
    out_C, out_T = [], []
    while sum(len(c) for c in out_C) < n:
        m = 20000
        C = np.zeros((m, 3, 2)); C[:, 1, 0] = rng.uniform(RMIN, 3.6, m)
        ang = rng.uniform(np.pi/8, np.pi/1.8, m); l = rng.uniform(RMIN, 3.6, m)
        C[:, 2] = np.stack([l*np.cos(ang), l*np.sin(ang)], 1)
        T = rng.uniform(0, P7, (m, 3))
        # shrink centres toward centroid until first contact (bisection), to sample near-contacts
        cen = C.mean(1, keepdims=True); lo = np.full(m, 0.3); hi = np.full(m, 1.0)
        ok = feasible(C, T) | True
        for _ in range(18):
            mid = (lo+hi)/2; Cm = cen + (C-cen)*mid[:, None, None]
            a, R, dmin, s, E = tri_feats(Cm, T)
            good = s >= 0
            hi = np.where(good, mid, hi); lo = np.where(good, lo, mid)
        Cs = cen + (C-cen)*hi[:, None, None]
        # random extra expansion to cover non-contact configs too
        Cs = cen + (Cs-cen)*np.where(rng.random(m) < 0.5, 1.0, rng.uniform(1, 1.6, m))[:, None, None]
        f = feasible(Cs, T)
        out_C.append(Cs[f]); out_T.append(T[f])
    return np.concatenate(out_C)[:n], np.concatenate(out_T)[:n]

def x2ct(x):
    C = np.array([[[0, 0], [x[0], 0], [x[1], x[2]]]]); T = np.array([x[3:6]]); return C, T

def bval(x, c):
    C, T = x2ct(x); a, M = rows(C, T); return a[0] + M[0] @ c

def bcons(x):
    C, T = x2ct(x); a, R, dmin, s, E = tri_feats(C, T)
    g = [s[0], 2-R[0], a[0]]
    if NONOB:
        for k in range(3):
            P = C[0]; u, v = P[(k+1) % 3]-P[k], P[(k+2) % 3]-P[k]; g.append(u @ v)
    return np.array(g)

def solve(A, Mrows, box):
    c, t = solve0(A, Mrows, box)
    # phase 2: keep t >= t* - ETA, minimise L1 norm of c (smallest transfer on the optimal face)
    nb = Mrows.shape[1]; tt = t - ETA
    # variables c = p - q, p,q >= 0 ; constraint: -M (p - q) <= a - tt
    Aub = np.hstack([-Mrows, Mrows])
    res = linprog(np.ones(2*nb), A_ub=Aub, b_ub=A - tt, bounds=[(0, box)]*(2*nb), method='highs')
    if res.status != 0: return c, t
    return res.x[:nb] - res.x[nb:], t

def solve0(A, Mrows, box):
    nb = Mrows.shape[1]
    # variables: c (nb), t ; maximise t : t - M c <= a
    cobj = np.zeros(nb+1); cobj[-1] = -1
    Aub = np.hstack([-Mrows, np.ones((len(A), 1))])
    res = linprog(cobj, A_ub=Aub, b_ub=A, bounds=[(-box, box)]*nb+[(None, None)], method='highs')
    return res.x[:nb], res.x[-1]

if __name__ == '__main__':
    t0 = time.time()
    C, T = sample(40000)
    from hept import kk_double_lattice
    d = kk_double_lattice(); p = [np.zeros(2), d['t1'], 2*d['p2']]
    Cd = np.array([[p[0], p[1], p[2]], [p[0], 2*d['p6'], p[1]]]); Td = np.array([[0, 0, np.pi/7], [0, np.pi/7, 0]])
    Cd = Cd - Cd[:, :1]
    for q in range(2):
        ph = np.arctan2(Cd[q, 1, 1], Cd[q, 1, 0]); Rm = np.array([[np.cos(-ph), -np.sin(-ph)], [np.sin(-ph), np.cos(-ph)]])
        Cd[q] = Cd[q] @ Rm.T; Td[q] = Td[q] - ph
    # rotate so first edge on x-axis (the LP rows are invariant anyway)
    C = np.concatenate([C, Cd]); T = np.concatenate([T, Td])
    A, M = rows(C, T)
    print('DL triangle areas', A[-2:])
    print(f'nb={M.shape[1]} samples={len(A)} min area sample={A.min():.5f}', flush=True)
    # also seed with DL triangle and its rotations/relabels
    log = []
    for rd in range(ROUNDS):
        c, t = solve(A, M, BOX)
        bs = A + M @ c; worst = np.concatenate([np.argsort(bs)[:NW], rng.choice(len(bs), NR, replace=False)])
        # local minimisation of b with f fixed
        newx = []
        for w in worst:
            x0 = np.array([C[w, 1, 0], C[w, 2, 0], C[w, 2, 1], *T[w]])
            r = minimize(lambda x: bval(x, c), x0, method='SLSQP', constraints=[{'type': 'ineq', 'fun': bcons}],
                         options={'maxiter': 200, 'ftol': 1e-12})
            if bcons(r.x).min() > -1e-8: newx.append(r.x)
            if time.time()-t0 > TL: break
        newx = np.array(newx)
        Cn = np.zeros((len(newx), 3, 2)); Cn[:, 1, 0] = newx[:, 0]; Cn[:, 2] = newx[:, 1:3]; Tn = newx[:, 3:6]
        An, Mn = rows(Cn, Tn); bn = An + Mn @ c
        print(f'round {rd}: LP t={t:.6f} (t/acrit={t/ACRIT:.5f})  true min b after local search={bn.min():.6f} ({bn.min()/ACRIT:.5f})  [{time.time()-t0:.0f}s]', flush=True)
        log.append(dict(round=rd, lp_t=float(t), local_min_b=float(bn.min())))
        C = np.concatenate([C, Cn]); T = np.concatenate([T, Tn]); A = np.concatenate([A, An]); M = np.concatenate([M, Mn])
        if time.time()-t0 > TL: break
    c, t = solve(A, M, BOX)
    print(f'final LP t={t:.6f} t/acrit={t/ACRIT:.6f}')
    np.save(f'runs/lpsym_c_DR{DR}_DF{DF}_{"nonob" if NONOB else "any"}.npy', c)
    json.dump(log, open(f'runs/lpsym_log_DR{DR}_DF{DF}_{"nonob" if NONOB else "any"}.json', 'w'))

"""Second-order check along the contact-preserving flex, for minimal-norm first-order transfer params."""
import numpy as np, sys
sys.path.insert(0, 'src')
from local_dl import *
z0 = np.zeros(6); CT = find_contacts(); nl = len(CT)
G = np.array([grad(lambda z: gfun(z, c), z0) for c in CT]); ca = grad(area, z0); J = jac_edges(z0)
Fm = np.array([J[0][1], J[0][2], J[1][1]-J[1][2], J[2][1]-J[2][2]]).T
u, s, vt = np.linalg.svd(G); v = vt[5]

def solve(lmin):
    # min sum|p| s.t. ca + Fm p = G^T lam, lam >= lmin ; p = pp - pm
    npar = 4
    Aeq = np.hstack([Fm, -Fm, -G.T]); beq = -ca
    cost = np.concatenate([np.ones(2*npar), np.zeros(nl)])
    res = linprog(cost, A_eq=Aeq, b_eq=beq, bounds=[(0, None)]*(2*npar)+[(lmin, None)]*nl, method='highs')
    return (res.x[:npar]-res.x[npar:2*npar], res.x[2*npar:]) if res.status == 0 else (None, None)

def b_of(z, p, Hq):
    E = edges(z); dE = [E[k]-EREF[k] for k in range(3)]
    ga, gb, a1, a2 = p
    val = area(z)
    val += ga*dE[0][1] + gb*dE[0][2]
    val += a1*(dE[1][1]-dE[1][2]) + a2*(dE[2][1]-dE[2][2])
    # quadratic terms: e01: q1 dr*da + q2 dr*db ; e12: q3 dt*dr + q4 dt*ds ; e20: q5 dt*dr + q6 dt*ds
    q = Hq
    val += q[0]*dE[0][0]*dE[0][1] + q[1]*dE[0][0]*dE[0][2]
    for k, (qa, qb) in ((1, (q[2], q[3])), (2, (q[4], q[5]))):
        dt = dE[k][1]-dE[k][2]; ds = dE[k][1]+dE[k][2]
        val += dt*(qa*dE[k][0] + qb*ds)
    return val

def lagr(z, p, lam, Hq):
    return b_of(z, p, Hq) - sum(l*gfun(z, c) for l, c in zip(lam, CT))

def d2(fn, v, h=1e-4):
    return (fn(h*v) - 2*fn(0*v) + fn(-h*v))/h**2

for lmin in (0.01, 0.05, 0.1):
    p, lam = solve(lmin)
    if p is None: print('lmin', lmin, 'infeasible'); continue
    gb_ = grad(lambda z: b_of(z, p, np.zeros(6)), z0)
    print(f'lmin={lmin}: p={np.round(p, 5)} lam={np.round(lam, 5)}  |grad b - G^T lam|={np.linalg.norm(gb_-G.T@lam):.2e}')
    print('   second derivative of Lagrangian along flex (no quadratic f terms):', d2(lambda z: lagr(z, p, lam, np.zeros(6)), v))
    # sensitivity to each quadratic param
    sens = [d2(lambda z: lagr(z, p, lam, np.eye(6)[k]), v) - d2(lambda z: lagr(z, p, lam, np.zeros(6)), v) for k in range(6)]
    print('   d(second deriv)/d q_k:', np.round(sens, 5))
    print('   edge perturbations along flex (dr,da,db) per edge:', [np.round(J[k] @ v, 5) for k in range(3)])

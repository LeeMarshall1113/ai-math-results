"""Local lemma at the double lattice, single-triangle level (float exploration).
Triangle s2 = (H0, Ht1, H2p2). Variables z: perturbations (x,y,phi) of Ht1 and H2p2 (H0 fixed).
b(z) = area + f(e01) + f(e12) + f(e20), f antisymmetric (f o R = -f) and mirror-symmetric (f o m = f).
Jets allowed at the s2 edge configs:
  e01 translation edge (r,0,0):  f = 0, grad = (0, ga, gb)        [from f o S = -f, S(r,a,b)=(r,-a,-b)]
  e12, e20 point-reflected edges: f = 0, grad = alpha_k (0, 1, -1)  [from f o R = -f]
"""
import numpy as np, sys
sys.path.insert(0, 'src')
from hept import *
from scipy.optimize import linprog
P7 = 2*np.pi/7
d = kk_double_lattice()
C0 = [np.zeros(2), d['t1'], 2*d['p2']]; T0 = [0.0, 0.0, np.pi/7]

def state(z):
    C = [C0[0], C0[1]+z[0:2], C0[2]+z[3:5]]; T = [T0[0], T0[1]+z[2], T0[2]+z[5]]
    return C, T

def edge_cfg(C, T, i, j, ref):
    v = C[j]-C[i]; phi = np.arctan2(v[1], v[0])
    a = T[i]-phi; b = T[j]-phi
    # unwrap near the reference values
    a = ref[1] + ((a-ref[1]+P7/2) % P7 - P7/2); b = ref[2] + ((b-ref[2]+P7/2) % P7 - P7/2)
    return np.array([np.linalg.norm(v), a, b])

C, T = state(np.zeros(6))
EDGES = [(0, 1), (1, 2), (2, 0)]
EREF = []
for i, j in EDGES:
    v = C[j]-C[i]; phi = np.arctan2(v[1], v[0])
    EREF.append(np.array([np.linalg.norm(v), np.mod(T[i]-phi, P7), np.mod(T[j]-phi, P7)]))

def edges(z):
    C, T = state(z); return [edge_cfg(C, T, i, j, EREF[k]) for k, (i, j) in enumerate(EDGES)]

def area(z):
    C, _ = state(z); return tri_area(*C)

# contacts: (i, j, vertex index of H_j, edge index of H_i): g = n_k(t_i).(v - c_i) - u
def find_contacts(tol=1e-9):
    C, T = state(np.zeros(6)); out = []
    for i in range(3):
        for j in range(3):
            if i == j: continue
            Vj = verts(C[j], T[j])
            for vi, v in enumerate(Vj):
                for k in range(7):
                    a = T[i]+(2*k+1)*np.pi/7; n = np.array([np.cos(a), np.sin(a)])
                    g = n @ (v-C[i]) - U
                    if abs(g) < tol:
                        # must also lie within the edge (tangential coordinate)
                        tng = np.array([-n[1], n[0]]) @ (v-C[i])
                        if abs(tng) <= np.tan(np.pi/7)*U + 1e-9:
                            out.append((i, j, vi, k, tng))
    return out

def gfun(z, c):
    i, j, vi, k, _ = c; C, T = state(z)
    v = verts(C[j], T[j])[vi]; a = T[i]+(2*k+1)*np.pi/7
    return np.array([np.cos(a), np.sin(a)]) @ (v-C[i]) - U

def grad(fn, z, h=1e-6):
    g = np.zeros(6)
    for k in range(6):
        e = np.zeros(6); e[k] = h; g[k] = (fn(z+e)-fn(z-e))/(2*h)
    return g

def jac_edges(z, h=1e-6):
    J = np.zeros((3, 3, 6))
    for k in range(6):
        e = np.zeros(6); e[k] = h
        J[:, :, k] = (np.array(edges(z+e))-np.array(edges(z-e)))/(2*h)
    return J

if __name__ == '__main__':
    z0 = np.zeros(6)
    CT = find_contacts()
    print('edge configs', np.round(EREF, 6))
    print('contacts (i, j, vertex of j, edge of i, tangential pos):')
    for c in CT: print('  ', c[:4], round(c[4], 5))
    G = np.array([grad(lambda z: gfun(z, c), z0) for c in CT])
    print('contact Jacobian rank', np.linalg.matrix_rank(G, 1e-8), 'of', G.shape)
    ca = grad(area, z0)
    J = jac_edges(z0)
    # f gradient terms: e01: (0,ga,gb); e12: alpha1 (0,1,-1); e20: alpha2 (0,1,-1)
    cols = [J[0][1], J[0][2], J[1][1]-J[1][2], J[2][1]-J[2][2]]   # d/dz of each param's contribution
    Fm = np.array(cols).T   # 6 x 4
    # feasibility LP: ca + Fm p = G^T lam, lam >= eps ; maximise eps
    nl = len(CT); npar = 4
    # vars: p (4, free), lam (nl), eps ; equalities: Fm p - G^T lam = -ca ; lam_k - eps >= 0
    Aeq = np.hstack([Fm, -G.T, np.zeros((6, 1))]); beq = -ca
    Aub = np.hstack([np.zeros((nl, npar)), -np.eye(nl), np.ones((nl, 1))]); bub = np.zeros(nl)
    cost = np.zeros(npar+nl+1); cost[-1] = -1
    res = linprog(cost, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq,
                  bounds=[(-10, 10)]*npar+[(0, 10)]*nl+[(None, 1)], method='highs')
    print('LP status', res.status, res.message)
    if res.status == 0:
        p = res.x[:npar]; lam = res.x[npar:npar+nl]
        print('max min-multiplier eps =', res.x[-1])
        print('f params (ga, gb, alpha1, alpha2) =', np.round(p, 6))
        print('multipliers', np.round(lam, 6))
    # null space of contacts
    u, s, vt = np.linalg.svd(G); r = (s > 1e-8).sum(); N = vt[r:]
    print('null space dim', len(N))
    for v in N:
        print('  null dir', np.round(v, 5), ' darea.v', round(ca @ v, 8))
    # without transfer: is plain area first-order stable?
    res0 = linprog(np.zeros(nl), A_eq=-G.T, b_eq=-ca, bounds=[(0, None)]*nl, method='highs')
    print('plain area: grad in contact cone?', res0.status == 0)

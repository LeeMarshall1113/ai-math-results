"""Float exploration: augmented-Lagrangian local certificate F = L - acrit + c*sum g_k^2 >= 0."""
import numpy as np, sys
sys.path.insert(0, 'src')
from local_dl import find_contacts, gfun, state, area, edges, EREF
from scipy.optimize import minimize
CT = find_contacts()
p = np.array([-0.5, -0.07971591511176641, 0.038369436798363965, -0.3090046002347956])
lam = np.array([0.6279860426230446, 0.6375923879211624, 0.14052208653213275, 0.13, 0.13])
def b(z):
    E = edges(z); d = [E[k]-EREF[k] for k in range(3)]
    return area(z) + p[0]*d[0][1] + p[1]*d[0][2] + p[2]*(d[1][1]-d[1][2]) + p[3]*(d[2][1]-d[2][2])
A0 = b(np.zeros(6))
def gs(z): return np.array([gfun(z, c) for c in CT])
def L(z): return b(z) - lam @ gs(z)
def hess(f, z, h=1e-4):
    H = np.zeros((6, 6)); I = np.eye(6)*h
    for i in range(6):
        for j in range(6):
            H[i, j] = (f(z+I[i]+I[j]) - f(z+I[i]-I[j]) - f(z-I[i]+I[j]) + f(z-I[i]-I[j]))/(4*h*h)
    return H
def jac(z, h=1e-6):
    return np.array([[(gfun(z+h*e, c)-gfun(z-h*e, c))/(2*h) for e in np.eye(6)] for c in CT])
HL = hess(L, np.zeros(6)); G = jac(np.zeros(6))
print('eig H_L', np.round(np.linalg.eigvalsh(HL), 4))
print('sing G', np.round(np.linalg.svd(G)[1], 4))
for c in (1, 3, 10, 30, 100):
    print('c', c, 'min eig H_L + 2c G^T G =', round(np.linalg.eigvalsh(HL + 2*c*G.T@G).min(), 5))
# how far does F stay >= 0 (min over sphere of radius r, sampled + local min)?
def F(z, c): return L(z) - A0 + c*(gs(z)**2).sum()
rng = np.random.default_rng(1)
for c in (3, 10, 30):
    for r in (0.005, 0.01, 0.02, 0.05):
        worst = 1e9
        for _ in range(300):
            u = rng.normal(size=6); u /= np.abs(u).max(); z = r*u
            worst = min(worst, F(z, c))
        mumin = min(lam.min() - c*np.abs(gs(r*rng.uniform(-1, 1, 6))).max() for _ in range(300))
        print(f'  c={c} r_inf={r}: min F on sampled box boundary {worst:.3e}, min mu sampled {mumin:.3f}')

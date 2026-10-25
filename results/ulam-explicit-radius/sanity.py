"""Floating-point sanity checks of PROOF.md (NOT part of the proof).

1. Lemma 3 construction on a concrete body: lattice is admissible and
   det(I+Q') <= 1 + A lam(U) + kappa(L) G lam(U).
2. Lemma 6 duality constant: ||A lam||_{L1(SO3)} / ||lam||_1 vs 1/N_eff.
"""
import numpy as np
from scipy.special import sph_harm_y
from scipy.spatial.transform import Rotation

rng = np.random.default_rng(1)
S6 = np.array([(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]) / np.sqrt(2)
S12 = np.vstack([S6, -S6])


def dual(k):
    rho = np.zeros(6); rho[k] = 1
    s12, s13, s23 = rho[0] + rho[1], rho[2] + rho[3], rho[4] + rho[5]
    Q = np.zeros((3, 3))
    Q[0, 0] = (s12 + s13 - s23) / 2; Q[1, 1] = (s12 + s23 - s13) / 2; Q[2, 2] = (s13 + s23 - s12) / 2
    Q[0, 1] = Q[1, 0] = (rho[0] - rho[1]) / 2
    Q[0, 2] = Q[2, 0] = (rho[2] - rho[3]) / 2
    Q[1, 2] = Q[2, 1] = (rho[4] - rho[5]) / 2
    return Q


DUAL = np.array([dual(k) for k in range(6)])


def realY(l, coef, xyz):
    """real even function: sum_m coef[m] * Re/Im Y_l^m."""
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    th = np.arccos(np.clip(z, -1, 1)); ph = np.arctan2(y, x)
    out = np.zeros(x.shape)
    for j, m in enumerate(range(-l, l + 1)):
        Y = sph_harm_y(l, abs(m), th, ph)
        out += coef[j] * (Y.real if m >= 0 else Y.imag)
    return out


def sphere_samples(n):
    v = rng.normal(size=(n, 3)); return v / np.linalg.norm(v, axis=1, keepdims=True)


# ------------------------------------------------------------ test 1
l = 10
coef = rng.normal(size=2 * l + 1)
P = sphere_samples(200000)
Yv = realY(l, coef, P)
t = 2e-5 / np.abs(Yv).max()          # sup |lam| = 2e-5 (well inside convexity)
lam = lambda xyz: t * realY(l, coef, xyz)
# Lipschitz constant of lam (numerical, finite differences on samples)
h = 1e-4
grads = []
for _ in range(3):
    d = rng.normal(size=P.shape); d -= (d * P).sum(1, keepdims=True) * P
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    Pd = P + h * d; Pd /= np.linalg.norm(Pd, axis=1, keepdims=True)
    grads.append(np.abs(lam(Pd) - lam(P)) / h)
L = 1.2 * max(g.max() for g in grads)
eps = np.abs(lam(P)).max() * 1.01
print(f"test 1: degree-{l} body, sup|lam|={eps:.2e}, Lip~{L:.2e}")
worst_margin, worst_det = np.inf, -np.inf
for trial in range(300):
    U = Rotation.random(random_state=rng).as_matrix()
    X = S6 @ U.T                       # U x_i
    lam_i = lam(X); r_i = np.exp(lam_i); rho = r_i - 1
    tau = np.zeros(6)
    for _ in range(100):              # fixed point iteration of Psi
        Qp = np.tensordot(rho + tau, DUAL, 1)
        perp = np.array([np.linalg.norm(Qp @ x - (x @ Qp @ x) * x) for x in S6])
        tau = np.exp(eps) * (np.exp(L * np.exp(eps) * perp) - 1)
    Qp = np.tensordot(rho + tau, DUAL, 1)
    pts = (S6 @ (np.eye(3) + Qp).T) @ U.T   # U (I+Q') x_i
    norms = np.linalg.norm(pts, axis=1)
    rdir = np.exp(lam(pts / norms[:, None]))
    worst_margin = min(worst_margin, (norms - rdir).min())
    detv = np.linalg.det(np.eye(3) + Qp)
    G = np.abs(lam_i).sum()
    cbar = np.exp(2 * eps) * L * np.sqrt(6)
    kappa = eps * np.exp(eps) / 4 + np.exp(eps) * cbar / (2 * (1 - cbar)) + 72 * eps * np.exp(2 * eps)
    bound = 1 + 0.5 * lam_i.sum() + kappa * G
    worst_det = max(worst_det, detv - bound)
print(f"   min_i (|p_i| - r(dir p_i)) over trials = {worst_margin:.3e}  (must be >= 0)")
print(f"   max (det(I+Q') - bound) over trials   = {worst_det:.3e}  (must be <= 0)")

# ------------------------------------------------------------ test 2
print("test 2: ||A lam||_{L1(SO3)} / ||lam||_1  (proof uses >= 1/N_eff = 1/180.7 = 0.00553)")
Us = Rotation.random(40000, random_state=rng).as_matrix()
for l in (4, 10, 12, 22):
    for rep in range(2):
        coef = rng.normal(size=2 * l + 1)
        lam1 = np.mean(np.abs(realY(l, coef, P)))
        pts = np.einsum('nab,kb->nka', Us, S6)
        A = 0.5 * realY(l, coef, pts).sum(1)
        print(f"   l={l:2d}: ratio = {np.mean(np.abs(A)) / lam1:.4f}")

"""Exclusion of boxes lying inside the certified local charts around the DL triangle.

A canonical box (d, x, y, t0, t1, t2) is mapped (rigorously, rig.R) to the chart frame of the
s2 triangle (H0 fixed at the origin with orientation 0; Ht1, H2p2 perturbed by z).
Canonical DL labelling: P0 = H2p2, P1 = H0, P2 = Ht1.
The box is discarded if, for some certified chart c, every z in the image satisfies
  |J_c (z - z_c)|_inf <= rho_c,  g_k(z) <= eta_c (k = 1..5),  |(J_c (z - z_c))_6| <= eta_c.
Every non-overlapping configuration there has b >= acrit. This is NOT a per-chart statement
(a single non-central chart's own bound m_c can be negative). It follows from the tube theorem
(PROOFS-DRAFT 6.2, Thm 6.1): (A)-(D) per chart, (E) eta = (1-q) rho - 2 off, (G) gluing and (Ov)
window overlap (certified by tube_glue_arb.py), so beta is one convex function with beta'(0) = 0;
then monotonicity along w > 0 inside each chart. g_k < 0 forces overlap inside a chart by (A).
"""
import json, numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rig import R, rcos, rsin, PI, PIO2
U_lo, U_hi = np.nextafter(0.9009688679024191, -np.inf), np.nextafter(0.9009688679024191, np.inf)
UR = R(U_lo, U_hi)
PI7 = PI*(1.0/7.0); P7I = PI*(2.0/7.0)
W = 1e-15
KROT = -2      # ph = -t1 - 2*(2pi/7) for the canonical DL labelling

class Charts:
    def __init__(self, tube_json, ref_json):
        D = json.load(open(tube_json)); ref = json.load(open(ref_json))
        self.ch = [c for c in D if c['ok']]
        C0 = ref['C0']
        self.t1ref = (R(C0[1][0][0]-W, C0[1][0][1]+W), R(C0[1][1][0]-W, C0[1][1][1]+W))
        self.p2ref = (R(C0[2][0][0]-W, C0[2][0][1]+W), R(C0[2][1][0]-W, C0[2][1][1]+W))
        self.contacts = [tuple(c) for c in ref['contacts']]
        self.zc = np.array([c['zc'] for c in self.ch])
        self.Jlo = np.array([c['Jlo'] for c in self.ch]); self.Jhi = np.array([c['Jhi'] for c in self.ch])
        self.rho = np.array([c['rho'] for c in self.ch]); self.eta = np.array([c['eta'] for c in self.ch])
        # canonical DL (for a cheap pre-filter)
        self.dl = np.array([1.9149801, 0.8896, 1.6007, 0.7941123, 0.3453134, 0.3453134])

    def to_z(self, B):
        n = len(B)
        d = R(B[:, 0, 0], B[:, 0, 1]); x = R(B[:, 1, 0], B[:, 1, 1]); y = R(B[:, 2, 0], B[:, 2, 1])
        t0 = R(B[:, 3, 0], B[:, 3, 1]); t1 = R(B[:, 4, 0], B[:, 4, 1]); t2 = R(B[:, 5, 0], B[:, 5, 1])
        al = -t1 + P7I*float(KROT)                # true rotation back to the chart frame (k fixed for the DL)
        ca, sa = rcos(al), rsin(al)
        def rot(px, py):                          # R(al) (p - P1), P1 = (d, 0)
            qx, qy = px - d, py
            return ca*qx - sa*qy, sa*qx + ca*qy
        zero = R(np.zeros(n))
        Q2 = rot(x, y); Q0 = rot(zero, zero)
        z = [Q2[0] - self.t1ref[0], Q2[1] - self.t1ref[1], t2 - t1,
             Q0[0] - self.p2ref[0], Q0[1] - self.p2ref[1], t0 - t1 - PI7]
        return z

    def g(self, z, c):
        i, j, m, k = c
        C = [(R(0.0), R(0.0)), (self.t1ref[0] + z[0], self.t1ref[1] + z[1]), (self.p2ref[0] + z[3], self.p2ref[1] + z[4])]
        T = [R(0.0), z[2], PI7 + z[5]]
        a = T[i] + PI7*float(2*k+1); th = T[j] + P7I*float(m)
        wx = C[j][0] + rcos(th) - C[i][0]; wy = C[j][1] + rsin(th) - C[i][1]
        return rcos(a)*wx + rsin(a)*wy - UR

    def inside(self, B):
        n = len(B); res = np.zeros(n, bool)
        near = np.all(np.abs(B.mean(2) - self.dl) < 0.12, 1)
        if not near.any(): return res
        idx = np.where(near)[0]; Bn = B[idx]
        z = self.to_z(Bn)
        gs = [self.g(z, c) for c in self.contacts]
        gmax = np.max([g.hi for g in gs], 0)
        zmid = np.stack([0.5*(q.lo + q.hi) for q in z], 1)                     # (m,6)
        Jmid = 0.5*(self.Jlo + self.Jhi)                                        # (C,6,6)
        zeta6 = np.einsum('cj,mcj->mc', Jmid[:, 5, :], zmid[:, None, :] - self.zc[None, :, :])
        cand = np.abs(zeta6) < 2*self.eta[None, :] + 1e-3                      # coarse float pre-selection only
        for ci in np.where(cand.any(0))[0]:
            ok = (gmax <= self.eta[ci]) & cand[:, ci]
            if not ok.any(): continue
            dz = [z[q] - self.zc[ci, q] for q in range(6)]
            for i in range(6):
                acc = R(np.zeros(len(idx)))
                for j in range(6):
                    acc = acc + dz[j]*R(np.full(len(idx), self.Jlo[ci, i, j]), np.full(len(idx), self.Jhi[ci, i, j]))
                lim = self.rho[ci] if i < 5 else self.eta[ci]
                ok &= (acc.hi <= lim) & (acc.lo >= -lim)
            res[idx[ok]] = True
        return res

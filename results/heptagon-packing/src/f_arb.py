"""Evaluate the fitted global transfer f (lp_transfer2 basis) on arb_series inputs, and the
modified area b for the s2 chart variables. Basis must match lp_transfer2.gbasis/fbasis."""
import sys
import numpy as np
sys.path.insert(0, 'src')
from flint import arb, arb_series
from local_lemma_arb import state, area, edge_cfg, EDGES, PI, U

RMIN = 2*U; RMAX = arb(4)

class FArb:
    def __init__(self, coef_path, DR=3, DF=2):
        self.c = [arb(float(x)) for x in np.load(coef_path)]
        self.DR, self.DF = DR, DF

    def g(self, r, a, b):
        x = (2*r - RMIN - RMAX)/(RMAX - RMIN)
        Tr = [x*0 + 1, x]
        for k in range(2, self.DR+1): Tr.append(2*x*Tr[-1] - Tr[-2])
        A7, B7 = 7*a, 7*b
        Fa = [A7*0 + 1]; Fb = [B7*0 + 1]
        for k in range(1, self.DF+1):
            Fa += [(k*A7).cos(), (k*A7).sin()]; Fb += [(k*B7).cos(), (k*B7).sin()]
        vals = [t*fa*fb for t in Tr for fa in Fa for fb in Fb]
        return vals

    def f(self, r, a, b):
        p7 = PI/7
        terms = [self.g(r, a, b), self.g(r, b + p7, a + p7), self.g(r, p7 - b, p7 - a), self.g(r, -a, -b)]
        out = None
        for i, ci in enumerate(self.c):
            bi = terms[0][i] - terms[1][i] + terms[2][i] - terms[3][i]
            out = ci*bi if out is None else out + ci*bi
        return out

    def b(self, z):
        C, T = state(z)
        val = area(C)
        for (i, j) in EDGES:
            r, a, bb = edge_cfg(C, T, i, j, None)
            val = val + self.f(r, a, bb)
        return val


def exactify(F, v, nb=None):
    """Return (F', mu, idx): F' = F + mu * basis_idx, with mu (a ball enclosing one exact real)
    chosen so that grad b_{F'}(0) . v = 0 exactly (stationarity along the flex)."""
    from local_lemma_arb import coef
    from flint import arb_series
    ORD = 3
    def dir_deriv(bfun):
        z = [arb_series([arb(0), v[q]], prec=ORD) for q in range(6)]
        return coef(bfun(z), 1)
    base = dir_deriv(F.b)
    nb = len(F.c)
    best = None
    for i in range(nb):
        Fi = FArb.__new__(FArb); Fi.DR, Fi.DF = F.DR, F.DF
        Fi.c = [arb(0)]*nb; Fi.c[i] = arb(1)
        # derivative of the f-part only (area part cancels in difference)
        di = dir_deriv(lambda z: Fi.b(z)) - dir_deriv(lambda z: FArb_zero(F).b(z))
        if best is None or abs(di.mid()) > abs(best[1].mid()):
            best = (i, di)
    i, di = best
    mu = -base/di
    F2 = FArb.__new__(FArb); F2.DR, F2.DF = F.DR, F.DF
    F2.c = list(F.c); F2.c[i] = F2.c[i] + mu
    return F2, mu, i, base

def FArb_zero(F):
    Z = FArb.__new__(FArb); Z.DR, Z.DF = F.DR, F.DF; Z.c = [arb(0)]*len(F.c); return Z

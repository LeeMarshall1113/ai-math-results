"""lpfast with the float LP sent straight to HiGHS (hlp.solve_ub): the same dense matrix, same options, so the
solution vector is bit-identical to scipy linprog's; the item matrices are filled in place instead of
np.vstack-ing ~60 freshly allocated blocks.  items_for / round_solution / exact_ok are lpfast's own."""
import os, math
import numpy as np
from lpfast import *            # noqa: F401,F403  (CUBE, EDGES, items_for, exact_ok, LAM, T, DEN, ...)
from lpfast import LAM, T, DEN
import hlp

def float_solve(items, nv, nmu, Nf, cf, s_vec):
    m = len(cf); nt = nv + nmu + 1; k_items = len(items)
    A = np.zeros((m * k_items, nt)); b = np.empty(m * k_items)
    mu = nv
    sv = np.array([float(x) for x in s_vec]) if s_vec is not None else None
    lamf = float(LAM); cfl = cf * (1 - lamf)
    for j, (base, lin, Ti, mg) in enumerate(items):
        Aj = A[j * m:(j + 1) * m]
        bb = np.array([float(x) for x in base])
        b[j * m:(j + 1) * m] = (cfl if mg else cf) + Nf @ bb
        for k, u in lin: Aj[:, k] -= Nf @ np.array(u, float)
        Aj[:, Ti:Ti + 3] += Nf
        if mg:
            Aj[:, mu] += Nf @ (float(LAM * T) * sv); mu += 1
        Aj[:, -1] = 1.0
    c = np.zeros(nt); c[-1] = -1
    lb = np.full(nt, -hlp.INF); ub = np.full(nt, hlp.INF)
    lb[nv:nv + nmu] = 0.0; ub[nv:nv + nmu] = 1.0; ub[-1] = 1.0
    ok, x = hlp.solve_ub(c, A, b, lb, ub)
    if not ok: return None, None
    return x[-1], x[:-1]

_exact_ok_slow = exact_ok   # lpfast's all-python-int version (kept for the uncertain entries / testing)
FILTER_REL = 1e-12          # float filter margin, relative to the sum of |terms|; true rounding error is < 1e-15

STRICT = os.environ.get('CERT_STRICT') == '1'   # 2026-10-31 (Lean): demand every row value > 0 (open translates)

def exact_ok(items, nv, nmu, q, Hint, s_vec):
    """Same answer as lpfast.exact_ok, decided by a float filter with a rigorous error margin where possible:
    every value v = c0 D + n.z (or 3 c0 D + 10 n.z) is a sum of 4 products of integers; in float64 its error is below
    ~6e-16 * (|c0| D + |n|.|z|).  Entries with |v_float| > FILTER_REL * that magnitude have the sign of v_float;
    the rest (normally none) are recomputed in exact python integers, exactly as lpfast does."""
    Zs = []; Ds = []; M = []
    mu = nv
    for base, lin, Ti, mg in items:
        L = 1
        for x in base: L = L * x.denominator // math.gcd(L, x.denominator)
        D = DEN * L
        z = []
        for i in range(3):
            v = base[i].numerator * (D // base[i].denominator)
            for k, u in lin: v += q[k] * u[i] * L
            v -= q[Ti + i] * L
            z.append(v)
        if mg:
            lt = LAM * T
            z = [zi * lt.denominator - q[mu] * L * lt.numerator * s_vec[i] for i, zi in enumerate(z)]
            D = D * lt.denominator; mu += 1
        Zs.append(z); Ds.append(D); M.append(mg)
    H = Hint if isinstance(Hint, np.ndarray) else np.array(Hint, dtype=np.int64)
    c0 = H[:, 0].astype(float); N = H[:, 1:].astype(float)
    Zf = np.array(Zs, dtype=float).T; Df = np.array(Ds, dtype=float); Mv = np.array(M)
    cw = np.where(Mv, 3.0, 1.0); nw = np.where(Mv, 10.0, 1.0)
    v = c0[:, None] * (Df * cw)[None, :] + (N @ Zf) * nw[None, :]
    mag = np.abs(c0)[:, None] * (Df * cw)[None, :] + (np.abs(N) @ np.abs(Zf)) * nw[None, :]
    tol = FILTER_REL * mag
    if (v < -tol).any(): return False
    unsure = np.abs(v) <= tol
    if not unsure.any(): return True
    rows, cols = np.nonzero(unsure)
    for r, j in zip(rows.tolist(), cols.tolist()):
        h = [int(x) for x in H[r]]; z = Zs[j]; D = Ds[j]
        nz = h[1] * z[0] + h[2] * z[1] + h[3] * z[2]
        val = 3 * h[0] * D + 10 * nz if M[j] else h[0] * D + nz
        if val < 0 or (STRICT and val == 0): return False
    return True

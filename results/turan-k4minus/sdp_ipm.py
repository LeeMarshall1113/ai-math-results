"""Primal-dual interior-point solver (HKM direction, Mehrotra predictor-corrector) for the
flag-algebra SDP of prep.py.

Certificate side (P):  min lambda
   s.t. for every H:  sum_b <A_bH, Q_b> + z_H - lambda + sum_r mu_r S_rH = -e_H,
        Q_b PSD, z >= 0, lambda, mu free.
i.e.  e_H + sum_b <A_bH,Q_b> + sum_r mu_r S_rH <= lambda  for all H.
Moment side (D): max sum_H e_H phi_H, phi >= 0, sum phi = 1, sum phi_H A_bH PSD, S phi = 0
(with y = -phi).

usage: python sdp_ipm.py data7.npz [stat] [out.npz]
"""
import sys, time, os
sys.path.append(os.path.expandvars(r"C:/Users/hackathon/AppData/Roaming/Python/Python314/site-packages"))
import numpy as np, scipy.sparse as sp, scipy.linalg as sla

def log(*a):
    print(*a, flush=True)

LOADINFO = {}
PRODUCT = False
NODEG = False

def load(path, use_stat, interval=None):
    D = np.load(path)
    N = int(D["N"]); e = D["edges"].astype(float) / (N * (N - 1) * (N - 2) / 6)
    m = len(e)
    blocks = []
    for b in range(int(D["nblocks"])):
        T = D[f"b{b}_trip"]; den = float(D[f"b{b}_den"]); n = len(D[f"b{b}_flags"])
        h, i, j, c = T[:, 0], T[:, 1], T[:, 2], T[:, 3].astype(float) / den
        At = sp.csr_matrix((c, (h, i * n + j)), shape=(m, n * n))
        At.sum_duplicates()
        blocks.append(dict(n=n, At=At, AtT=At.T.tocsr(), s=int(D[f"b{b}_s"])))
    Bcols = [-np.ones(m)]
    if use_stat:
        S = D["stat420"].astype(float) / 420.0
        # drop linearly dependent rows
        q, r, piv = sla.qr(S.T, mode="economic", pivoting=True)
        rk = int((np.abs(np.diag(r)) > 1e-9 * abs(r[0, 0])).sum())
        sel = np.sort(piv[:rk]); S = S[sel]; LOADINFO["statsel"] = sel
        log("stationarity rows used", S.shape[0])
        Bcols += list(S)
    B = np.array(Bcols).T
    # free columns F (lambda, linear-stat mu) handled in an augmented Schur system;
    # LP columns C (w >= 0) with costs cw: interval multipliers p_g, q_g
    C = []; cw = []
    if interval is not None:
        t1, t2 = interval
        Ds = np.load(path.replace("data", "stat"))
        den = float(Ds["den"]); D_ = Ds["cntD"] / den; G_ = Ds["cntG"] / den
        for Dg, Gg in ([] if NODEG else zip(D_, G_)):
            C += [Dg - t1 * Gg, -(Dg - t2 * Gg)]; cw += [0.0, 0.0]
        if PRODUCT:
            pd = float(Ds["pden"])
            for PE, P in zip(Ds["cntPE"] / pd, Ds["cntP"] / pd):
                C += [PE - t1 * P, -(PE - t2 * P)]; cw += [0.0, 0.0]
            log("product rows", len(Ds["cntP"]))
        log("interval rows", len(D_), "t in", interval)
    C = np.array(C).T if C else np.zeros((m, 0))
    return N, e, blocks, B, C, np.array(cw)

def sym(A):
    return 0.5 * (A + A.T)

def opA(blk, X):
    return blk["At"] @ X.reshape(-1)

def opAT(blk, y):
    n = blk["n"]
    return sym((blk["AtT"] @ y).reshape(n, n))

def schur_block(blk, X, Si, chunk=96):
    n = blk["n"]; At = blk["At"]; m = At.shape[0]
    M = np.zeros((m, m))
    for c0 in range(0, m, chunk):
        idx = range(c0, min(m, c0 + chunk))
        W = np.zeros((n * n, len(idx)))
        ip, ind, dat = At.indptr, At.indices, At.data
        for t, k in enumerate(idx):
            lo, hi = ip[k], ip[k + 1]
            if lo == hi:
                continue
            ii, jj = np.divmod(ind[lo:hi], n)
            R, inv = np.unique(np.concatenate([ii, jj]), return_inverse=True)
            L = hi - lo; Am = np.zeros((len(R), len(R)))
            np.add.at(Am, (inv[:L], inv[L:]), dat[lo:hi])
            Wk = (X[:, R] @ Am) @ Si[R, :]
            W[:, t] = Wk.reshape(-1)
        M[:, c0:c0 + len(idx)] = At @ W
    return sym(M)

def maxstep(X, dX):
    try:
        L = np.linalg.cholesky(X)
    except np.linalg.LinAlgError:
        return 0.0
    Li = sla.solve_triangular(L, np.eye(len(X)), lower=True)
    ev = np.linalg.eigvalsh(sym(Li @ dX @ Li.T))
    lm = ev[0]
    return 1e30 if lm >= 0 else -1.0 / lm

def solve(N, e, blocks, F, C, cw, maxit=120, tol=1e-9, floor=None):
    """free variables u (columns F, cost cu = (1,0,...)), LP variables w >= 0 (columns C, cost cw)."""
    m = len(e); b = -e; p = C.shape[1]; nf = F.shape[1]
    if floor is not None:
        # lambda = floor + s0, s0 >= 0 (last LP column): bounded even when [t1,t2] is excluded
        b = -(e - floor)
    cu = np.zeros(nf)
    if nf: cu[0] = 1.0
    X = [np.eye(k["n"]) for k in blocks]; S = [np.eye(k["n"]) for k in blocks]
    z = np.ones(m); s = np.ones(m); y = np.zeros(m); u = np.zeros(nf)
    w = np.ones(p); sw = np.ones(p)
    ntot = sum(k["n"] for k in blocks) + m + p
    t0 = time.time()
    for it in range(maxit):
        rp = b - sum(opA(k, Xb) for k, Xb in zip(blocks, X)) - z - C @ w - F @ u
        Rd = [-opAT(k, y) - Sb for k, Sb in zip(blocks, S)]
        rdz = -y - s
        rdw = cw - C.T @ y - sw
        rdu = cu - F.T @ y
        mu = (sum(np.sum(Xb * Sb) for Xb, Sb in zip(X, S)) + z @ s + w @ sw) / ntot
        pobj = cw @ w + cu @ u; dobj = b @ y
        pinf = np.linalg.norm(rp) / (1 + np.linalg.norm(b))
        dinf = (np.sqrt(sum(np.sum(R * R) for R in Rd) + rdz @ rdz) + np.linalg.norm(rdw) + np.linalg.norm(rdu))
        gap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj))
        log(f"it {it:2d} pobj {pobj:.10f} dobj {dobj:.10f} gap {gap:.2e} pinf {pinf:.2e} dinf {dinf:.2e} mu {mu:.2e} t {time.time()-t0:.0f}s")
        merit = max(gap, pinf, dinf)
        if it == 0 or merit < best[0]:
            best = (merit, [x.copy() for x in X], z.copy(), w.copy(), u.copy(), y.copy())
        if merit < tol or (it > 5 and merit > 100 * best[0]) or mu < 1e-14:
            break
        Si = [sym(np.linalg.inv(Sb)) for Sb in S]
        M = np.diag(z / s)
        if p:
            M += (C * (w / sw)) @ C.T
        for k, Xb, Sib in zip(blocks, X, Si):
            tb = time.time(); M += schur_block(k, Xb, Sib)
            if it == 0: log(f"   schur block n={k['n']} {time.time()-tb:.1f}s")
        try:
            cf = sla.cho_factor(M, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            try:
                M[np.diag_indices(m)] += 1e-11 * np.abs(np.diag(M)).max()
                cf = sla.cho_factor(M, lower=True, check_finite=False)
            except np.linalg.LinAlgError:
                log("Schur complement not PD; stopping at current iterate")
                break
        MiF = sla.cho_solve(cf, F, check_finite=False) if nf else np.zeros((m, 0))
        FMF = F.T @ MiF

        def direction(sig, corrX=None, corrz=None, corrw=None):
            K = []
            for idx, (Xb, Sib, Rb) in enumerate(zip(X, Si, Rd)):
                T = sig * mu * Sib - Xb - Xb @ Rb @ Sib
                if corrX is not None:
                    T = T - corrX[idx]
                K.append(sym(T))
            kz = sig * mu / s - z - (z / s) * rdz
            kw = sig * mu / sw - w - (w / sw) * rdw
            if corrz is not None:
                kz = kz - corrz; kw = kw - corrw
            r1 = rp - sum(opA(k, Kb) for k, Kb in zip(blocks, K)) - kz - C @ kw
            Mir = sla.cho_solve(cf, r1, check_finite=False)
            du = np.linalg.solve(FMF, F.T @ Mir - rdu) if nf else np.zeros(0)
            dy = Mir - MiF @ du
            dS = [Rb - opAT(k, dy) for k, Rb in zip(blocks, Rd)]
            dX = [sym(Kb + Xb @ opAT(k, dy) @ Sib) for k, Kb, Xb, Sib in zip(blocks, K, X, Si)]
            ds = rdz - dy
            dz = kz + (z / s) * dy
            dsw = rdw - C.T @ dy
            dw = kw + (w / sw) * (C.T @ dy)
            return dX, dz, dw, du, dy, dS, ds, dsw

        def lpstep(v, dv):
            neg = dv < 0
            return np.min(-v[neg] / dv[neg]) if neg.any() else 1e30

        def steps(dX, dz, dw, dS, ds, dsw):
            ap = min([maxstep(Xb, d) for Xb, d in zip(X, dX)] + [lpstep(z, dz), lpstep(w, dw)])
            ad = min([maxstep(Sb, d) for Sb, d in zip(S, dS)] + [lpstep(s, ds), lpstep(sw, dsw)])
            return min(1.0, ap), min(1.0, ad)

        dX, dz, dw, du, dy, dS, ds, dsw = direction(0.0)
        ap, ad = steps(dX, dz, dw, dS, ds, dsw)
        mua = (sum(np.sum((Xb + ap * a) * (Sb + ad * c)) for Xb, a, Sb, c in zip(X, dX, S, dS))
               + (z + ap * dz) @ (s + ad * ds) + (w + ap * dw) @ (sw + ad * dsw)) / ntot
        sig = min(1.0, (mua / mu) ** 3)
        corrX = [a @ c @ Sib for a, c, Sib in zip(dX, dS, Si)]
        dX, dz, dw, du, dy, dS, ds, dsw = direction(sig, corrX, dz * ds / s, dw * dsw / sw)
        ap, ad = steps(dX, dz, dw, dS, ds, dsw)
        ap = min(1.0, 0.95 * ap); ad = min(1.0, 0.95 * ad)
        log(f"   sig {sig:.2e} ap {ap:.3f} ad {ad:.3f}")
        X = [Xb + ap * d for Xb, d in zip(X, dX)]; z = z + ap * dz; w = w + ap * dw; u = u + ap * du
        S = [Sb + ad * d for Sb, d in zip(S, dS)]; s = s + ad * ds; y = y + ad * dy; sw = sw + ad * dsw
    merit, X, z, w, u, y = best
    log("best merit", merit)
    return X, z, w, u, y

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("data"); ap.add_argument("--stat", action="store_true")
    ap.add_argument("--interval", nargs=2, type=float); ap.add_argument("--out")
    ap.add_argument("--product", action="store_true")
    ap.add_argument("--no-degree", action="store_true")
    ap.add_argument("--floor", type=float, help="lambda >= floor (lambda becomes floor + LP slack)")
    a = ap.parse_args()
    PRODUCT = a.product; NODEG = a.no_degree
    t = time.time()
    N, e, blocks, F, C, cw = load(a.data, a.stat, a.interval)
    log("loaded", N, "m", len(e), "blocks", [k["n"] for k in blocks], "free", F.shape[1], "lp cols", C.shape[1], f"{time.time()-t:.0f}s")
    if a.floor is not None:
        assert F.shape[1] == 1, "--floor only without linear stat rows"
        F = F[:, 1:]; C = np.column_stack([C, -np.ones(len(e))]); cw = np.append(cw, 1.0)
    X, z, w, u, y = solve(N, e, blocks, F, C, cw, floor=a.floor)
    if a.floor is not None:
        lam = a.floor + w[-1]; w = w[:-1]; C = C[:, :-1]; cw = cw[:-1]; u = np.array([lam])
    val = e + sum(opA(k, Xb) for k, Xb in zip(blocks, X)) + C @ w + F[:, 1:] @ u[1:]
    log("lambda", u[0], "max_H (float)", val.max(), "dual", -e @ y)
    if a.out:
        np.savez(a.out, **{f"Q{b}": Xb for b, Xb in enumerate(X)}, w=w, u=u, phi=-y, cw=cw,
                 statsel=LOADINFO.get("statsel", np.zeros(0, int)), stat=int(a.stat),
                 interval=np.array(a.interval if a.interval else []), product=int(a.product), nodeg=int(a.no_degree))

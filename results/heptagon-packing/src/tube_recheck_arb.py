"""Audit-2 FLAW-4 recheck of a tube file with OUTWARD-widened exported data (arb, 200 bits).

The exported floats (J entries, q, off, hp, w_min, B2) come from float(arb), which rounds to
nearest (error <= 1/2 ulp). Here every exported number is widened by 2 ulps in the safe direction
and all chart conditions are rechecked in arb:

 (A)  hp+ < 0,  (B) q+ < 1,  (C) w- > 0,  (D) B2- > 0                  [per chart]
 (K') inflated-cube lemma: chart sub-boxes have radius (h/2)(1 + 1e-9) with h = 2 rho / n (n = 2),
      so their union contains the cube {|zeta| <= rho (1 + 5e-10)} = K'. (A)-(D) hold on K'.
 (E') for every z the B&B excludes with the exported J' (|J'(z - z_c)| <= rho, g <= eta,
      |v'.(z - z_c)| <= eta): with D = max entry |J - J'| <= 2 ulp(|J'|) and
      d = 6 D ||J'^{-1}|| rho (a bound on |(J - J')(z - z_c)|), require
        rho + d <= rho (1 + 5e-10)           (z lies in K'), and
        eta + max(off+, d) <= (1 - q+) rho (1 + 5e-10)   (Banach targets within reach in K').
 (G)  gluing with widened J (as in tube_glue_arb.py).
 (Ov) consecutive windows [v.z_c - eta_c, v.z_c + eta_c] intersect (widened v).
 (W)  window check: the next chart's centre parameter s* = v.z_{c'} lies in I_c, i.e.
      |v.(z_{c'} - z_c)| <= eta_c.
Prints one line per failure and a summary; exit code 0 iff everything passes.
"""
import json, sys, math
import numpy as np
from flint import arb, arb_mat, ctx
ctx.prec = 200
INF = math.inf
def up2(x): return float(np.nextafter(np.nextafter(x, INF), INF))
def dn2(x): return float(np.nextafter(np.nextafter(x, -INF), -INF))
def ulp(x): return float(np.spacing(abs(x))) if x != 0 else 5e-324

D = json.load(open(sys.argv[1])); ok = [r for r in D if r['ok']]
cen = ok[0]; assert cen['s_c'] == 0.0
plus = [r for r in ok[1:] if r['s_c'] > 0]; minus = [r for r in ok[1:] if r['s_c'] < 0]
INFL = arb(1) + arb(5)/10**10

def Jball(r):
    M = arb_mat(6, 6); Dmax = 0.0
    for i in range(6):
        for j in range(6):
            lo, hi = dn2(r['Jlo'][i][j]), up2(r['Jhi'][i][j])
            M[i, j] = (arb(lo) + arb(hi))/2 + arb(0, (hi - lo)/2*(1 + 1e-12) + 1e-320)
            Dmax = max(Dmax, 2*max(ulp(r['Jlo'][i][j]), ulp(r['Jhi'][i][j])))
    return M, Dmax
def Jmid(r):
    M = arb_mat(6, 6)
    for i in range(6):
        for j in range(6): M[i, j] = (arb(r['Jlo'][i][j]) + arb(r['Jhi'][i][j]))/2
    return M
def rowmax(M):
    return max((sum((abs(M[i, j]) for j in range(6)), arb(0)) for i in range(6)), key=lambda x: float(x.upper()))

fails = 0; worstE = 0.0; worstG = 0.0; worstOv = -INF; worstW = 0.0
cache = {}
def prep(r):
    k = id(r)
    if k not in cache:
        Jb, Dm = Jball(r); Jinv = Jmid(r).inv()
        cache[k] = (Jb, Dm, Jinv)
    return cache[k]
for r in ok:
    q, off, rho, eta = up2(r['q']), up2(r['off']), arb(r['rho']), arb(r['eta'])
    cond = [('A', up2(r['hp']) < 0), ('B', q < 1), ('C', dn2(r['w_min']) > 0), ('D', dn2(r['B2']) > 0)]
    Jb, Dm, Jinv = prep(r)
    d = arb(6)*arb(Dm)*rowmax(Jinv)*rho
    e1 = (rho + d) <= rho*INFL
    reach = eta + (arb(off) if off > float(d.upper()) else d)
    e2 = reach <= (1 - arb(q))*rho*INFL
    worstE = max(worstE, float((reach/((1 - arb(q))*rho*INFL)).upper()))
    cond += [("K'", bool(e1)), ("E'", bool(e2))]
    for name, v in cond:
        if not v:
            fails += 1; print('FAIL', name, 's_c', r['s_c'])
def sval(r, J): return sum((J[5, j]*arb(r['zc'][j]) for j in range(6)), arb(0))
npairs = 0
for walk in (plus, minus):
    prev = cen
    for nxt in walk:
        Jc, _, _ = prep(prev); Jn, _, Jn_inv_mid = prep(nxt)
        Jn_inv = Jn.inv()
        dz = [arb(nxt['zc'][j]) - arb(prev['zc'][j]) for j in range(6)]
        Jdz = [sum((Jc[i, j]*dz[j] for j in range(6)), arb(0)) for i in range(6)]
        lhs = max(abs(x) for x in Jdz) + rowmax(Jc*Jn_inv)*arb(2*up2(nxt['off']))
        okG = bool(lhs < arb(prev['rho'])); worstG = max(worstG, float(lhs.upper())/prev['rho'])
        sp_, sn_ = sval(prev, Jc), sval(nxt, Jn)
        if nxt['s_c'] > prev['s_c']: gap = (sn_ - arb(nxt['eta'])) - (sp_ + arb(prev['eta']))
        else: gap = (sp_ - arb(prev['eta'])) - (sn_ + arb(nxt['eta']))
        okO = bool(gap < 0); worstOv = max(worstOv, float(gap.upper()))
        sstar = abs(sum((Jc[5, j]*dz[j] for j in range(6)), arb(0)))
        okW = bool(sstar <= arb(prev['eta'])); worstW = max(worstW, float(sstar.upper())/prev['eta'])
        npairs += 1
        if not (okG and okO and okW):
            fails += 1; print('FAIL pair', prev['s_c'], '->', nxt['s_c'], 'G', okG, 'Ov', okO, 'W', okW)
        prev = nxt
print(f"charts {len(ok)}  pairs {npairs}  failures {fails}")
print(f"worst E' ratio (reach / ((1-q+) rho (1+5e-10))) = {worstE:.12f}")
print(f"worst G ratio = {worstG:.4f}   worst Ov gap = {worstOv:.3e}   worst W ratio |v.(z_c'-z_c)|/eta_c = {worstW:.4f}")
print('ALL CONDITIONS CERTIFIED (A)-(D), K\', E\', G, Ov, W' if fails == 0 else 'NOT CERTIFIED')
sys.exit(0 if fails == 0 else 1)

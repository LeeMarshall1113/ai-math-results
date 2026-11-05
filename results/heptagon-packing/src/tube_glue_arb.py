"""Certify (G) gluing and (Ov) window overlap for a tube file, in arb ball arithmetic.

Walk structure of local_tube_arb.py's output: [central, +s charts in order, -s charts in order];
each walk starts at the central chart. For each consecutive pair c -> c' in a walk:
 (G)  |J_c (z_{c'} - z_c)|_inf + ||J_c J_{c'}^{-1}||_inf * r_{c'} <= rho_c,
      where r_{c'} = 2 off_{c'} >= off_{c'}/(1-q_{c'}) bounds, in c' chart coordinates, the distance from
      z_{c'} to the exact curve point with Phi = (0, v.z_{c'}) (Banach, q < 1/2). So that curve point lies
      in chart c's cube K_c.
 (Ov) the s-windows I_c = [v.z_c - eta_c, v.z_c + eta_c] of consecutive charts intersect.
J entries are taken as balls spanning the exported [Jlo, Jhi]; v is row 6 of J (identical for all charts).
"""
import json, sys
from flint import arb, arb_mat, ctx
ctx.prec = 200
D = json.load(open(sys.argv[1]))
ok = [r for r in D if r['ok']]
cen = ok[0]; assert cen['s_c'] == 0.0
plus = [r for r in ok[1:] if r['s_c'] > 0]; minus = [r for r in ok[1:] if r['s_c'] < 0]
assert len(plus) + len(minus) == len(ok) - 1

def ball(lo, hi):
    lo, hi = arb(lo), arb(hi)
    return (lo + hi)/2 + arb(0, float(((hi - lo)/2).upper())*(1+1e-12) + 1e-300)

def Jm(r):
    M = arb_mat(6, 6)
    for i in range(6):
        for j in range(6): M[i, j] = ball(r['Jlo'][i][j], r['Jhi'][i][j])
    return M

def infnorm_vec(vs): return max(abs(x) for x in vs)
def infnorm_mat(M):
    best = arb(0)
    for i in range(6):
        srow = sum((abs(M[i, j]) for j in range(6)), arb(0))
        best = srow if srow > best or best == 0 else best
    return max((sum((abs(M[i, j]) for j in range(6)), arb(0)) for i in range(6)), key=lambda x: float(x.upper()))

def sval(r, J):
    return sum((J[5, j]*arb(r['zc'][j]) for j in range(6)), arb(0))

worst_G = 0.0; worst_Ov = None; nfail = 0; npairs = 0
for walk in (plus, minus):
    prev = cen
    for nxt in walk:
        Jc = Jm(prev); Jn = Jm(nxt); Jn_inv = Jn.inv()
        dz = [arb(nxt['zc'][j]) - arb(prev['zc'][j]) for j in range(6)]
        Jdz = [sum((Jc[i, j]*dz[j] for j in range(6)), arb(0)) for i in range(6)]
        K = Jc*Jn_inv
        rad_n = arb(2*nxt['off']) if nxt['off'] > 0 else arb(0)
        lhs = infnorm_vec(Jdz) + infnorm_mat(K)*rad_n
        okG = bool(lhs < arb(prev['rho']))
        worst_G = max(worst_G, float(lhs.upper())/prev['rho'])
        # (Ov): intervals [s - eta, s + eta]
        sp_, sn_ = sval(prev, Jc), sval(nxt, Jn)
        lo_hi = (sp_ + prev['eta']) if nxt['s_c'] > prev['s_c'] else (sn_ + nxt['eta'])
        hi_lo = (sn_ - nxt['eta']) if nxt['s_c'] > prev['s_c'] else (sp_ - prev['eta'])
        gap = hi_lo - lo_hi                      # must be < 0 (overlap)
        okO = bool(gap < 0)
        g = float(gap.upper()); worst_Ov = g if worst_Ov is None or g > worst_Ov else worst_Ov
        npairs += 1
        if not (okG and okO):
            nfail += 1
            print('FAIL pair', prev['s_c'], '->', nxt['s_c'], 'G lhs/rho', float(lhs.upper())/prev['rho'], 'gap', g)
        prev = nxt
# (Ov) also the central window contains 0 and the walks' coverage endpoints
print(f'pairs {npairs}  failures {nfail}  worst (G) lhs/rho = {worst_G:.4f}  worst (Ov) gap upper = {worst_Ov:.3e} (must be < 0)')
print('covered s-interval approx', min(r['s_c'] - r['eta'] for r in ok), max(r['s_c'] + r['eta'] for r in ok))
print('CERTIFIED (G) and (Ov)' if nfail == 0 else 'NOT CERTIFIED')

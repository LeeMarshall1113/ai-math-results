"""Rigorous certificate for the nonzero-denominator / strict-inequality side of src/dl_exact.py
(referee UNCLEAR, 2026-11-05).

dl_exact.py proves identities in Q(c, s) by reducing numerators to 0 modulo
    REL = [s^2 + c^2 - 1, 8 c^3 - 4 c^2 - 4 c + 1]
which is legitimate only if every denominator is nonzero at c = cos(pi/7), s = sin(pi/7), and the
contact-interiority claims need strict inequalities there. This script runs dl_exact.py (same
objects, nothing re-typed) and certifies, with exact ball arithmetic (python-flint arb, 400 bits;
every arb result is a rigorous enclosure):

 (0) isolating intervals. c is the unique root of p(c) = 8c^3 - 4c^2 - 4c + 1 in [4/5, 1]
     (p(4/5) < 0 < p(1), p' > 0 on [4/5, 1], checked in exact Fractions) and cos(pi/7) lies in
     that interval, so c = cos(pi/7) is the root. A rational interval [cl, ch] of width 1e-40 is
     obtained by exact bisection (sign tests in Fractions). s > 0 with s^2 = 1 - c^2, and
     sin(pi/7) > 0 satisfies this, so s = sin(pi/7); [sl, sh] is a rational interval from integer
     square roots. Both enclosures are rigorous; c and s are then arb balls containing them.
 (1) every positive is_zero check: the stored numerator reduction is re-verified by exact
     polynomial identity  num - sum(q_i REL_i) - rem == 0  (sympy expand, independent of the
     reducer), rem == 0; the combined denominator is nonzero; and EVERY denominator base (each
     Pow with negative exponent in the input expression) is nonzero.
 (2) the (a, b) linear system of the KK construction has nonzero determinant (arb_mat det).
 (3) det(Asub) != 0 (Asub = Gs[:, others], the matrix dl_exact.py hands to LUsolve) as an arb_mat
     determinant, so rank DG(0) = 5 and ker DG(0) is one-dimensional, spanned by the vector v
     with y1 = 1 that a rigorous arb solve encloses. (The symbolic LUsolve section of dl_exact.py
     does not finish -- expression swell -- and is not used here or by the proof.)
 (4) every contact-interiority comparison holds STRICTLY: s - tangential > 0 and
     s + tangential > 0 (|tangential| < s = half edge length).
Negative controls (entries tagged 'negctl') are excluded from (1) except that their reductions
are re-verified and must have rem != 0.

Exit code 0 iff everything is certified. Single thread, a few minutes, < 1 GB.
"""
import sys, os, io, time, contextlib, runpy
from fractions import Fraction
from math import isqrt
import sympy as sp
from flint import arb, arb_mat, ctx, fmpq

ctx.prec = 400
T0 = time.time()
fails = []
def check(ok, msg):
    print(('CERTIFIED ' if ok else 'NOT CERTIFIED ') + msg, flush=True)
    if not ok: fails.append(msg)
    return ok

# ---------------- (0) isolating intervals ----------------
def p_f(x): return 8*x**3 - 4*x**2 - 4*x + 1
def dp_f(x): return 24*x**2 - 8*x - 4
lo0, hi0 = Fraction(4, 5), Fraction(1)
check(p_f(lo0) < 0 < p_f(hi0), f'p(4/5) = {p_f(lo0)} < 0 < p(1) = {p_f(hi0)}')
check(dp_f(lo0) > 0 and (dp_f(hi0) > 0), f"p'(4/5) = {dp_f(lo0)} > 0, p'(1) = {dp_f(hi0)} > 0, p'' = 48x - 8 > 0 on [4/5,1]: p strictly increasing, unique root in [4/5,1]")
# cos(pi/7) >= 1 - x^2/2 with x = pi/7 < 3.1416/7 < 0.4489 gives cos(pi/7) > 0.8993 > 4/5, and <= 1
xs = Fraction(31416, 10000)/7
check(1 - xs**2/2 > lo0, f'cos(pi/7) >= 1 - x^2/2 >= {float(1 - xs**2/2):.5f} > 4/5 (x <= 3.1416/7)')
cl, ch = lo0, hi0
W = Fraction(1, 10**40)
while ch - cl > W:
    m = (cl + ch)/2
    if p_f(m) < 0: cl = m
    else: ch = m
check(p_f(cl) < 0 < p_f(ch) or p_f(ch) == 0, f'isolating interval  c in [{cl.numerator}/{cl.denominator}, {ch.numerator}/{ch.denominator}]  (= [{float(cl):.17f}, {float(ch):.17f}], width <= 1e-40)')
N = 10**45
s2lo, s2hi = 1 - ch*ch, 1 - cl*cl
sl = Fraction(isqrt(int(s2lo*N*N)), N)
sh = Fraction(isqrt(int(s2hi*N*N)) + 2, N)
check(sl > 0 and sl*sl <= s2lo and sh*sh >= s2hi, f'isolating interval  s in [{sl.numerator}/{sl.denominator}, {sh.numerator}/{sh.denominator}]  (s>0, s^2 = 1 - c^2; = [{float(sl):.17f}, {float(sh):.17f}])')
def ball(a, b): return arb(fmpq(a.numerator, a.denominator)).union(arb(fmpq(b.numerator, b.denominator)))
CB, SB = ball(cl, ch), ball(sl, sh)
print('c ball', CB.str(30, radius=True), ' s ball', SB.str(30, radius=True), flush=True)

# ---------------- evaluator ----------------
c_sym, s_sym = sp.symbols('c s')
_memo = {}
def ev(e):
    if e in _memo: return _memo[e]
    if e == c_sym: r = CB
    elif e == s_sym: r = SB
    elif e.is_Integer: r = arb(int(e))
    elif e.is_Rational: r = arb(fmpq(int(e.p), int(e.q)))
    elif e.is_Add:
        r = arb(0)
        for a in e.args: r = r + ev(a)
    elif e.is_Mul:
        r = arb(1)
        for a in e.args: r = r*ev(a)
    elif e.is_Pow:
        ex = e.exp
        assert ex.is_Integer, ('non-integer power', e)
        b = ev(e.base); k = int(ex)
        r = b**k if k >= 0 else arb(1)/(b**(-k))
    else:
        raise ValueError(('unsupported node', e))
    _memo[e] = r
    return r
def nonzero(b): return (b > 0) or (b < 0)      # python-flint: certain comparison on the ball
def den_bases(expr):
    out = set()
    for pw in expr.atoms(sp.Pow):
        if pw.exp.is_Integer and int(pw.exp) < 0: out.add(pw.base)
    return out
def mag(b): return b.str(8, radius=False)

# ---------------- run dl_exact.py (same objects) ----------------
print('running src/dl_exact.py (stdout captured) ...', flush=True)
here = os.path.dirname(os.path.abspath(__file__))
# dl_exact.py's final section (symbolic LUsolve null vector) suffers expression swell and did not finish
# in a 50 minute run (NOTES, runs/dl_exact.log: rows 0, 1 only). It is NOT needed: this script
# executes the source up to that section (everything through the construction of Gs) and replaces
# the symbolic null vector by a rigorous arb solve below.
src = open(os.path.join(here, 'dl_exact.py')).read()
cut = src.index('# exact null vector with y1-component 1')
G = {'__name__': 'dl_exact_cert_run'}
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(compile(src[:cut], 'dl_exact.py[:null-vector]', 'exec'), G)
txt = buf.getvalue().splitlines()
check('ALL PASS' in txt, 'dl_exact.py reports ALL PASS (exact identities, 1-3, 6)  [elapsed %.0fs]' % (time.time() - T0))
ZC, INT = G['ZERO_CHECKS'], G['INTERIOR']
REL = G['REL']
print(f'{len(ZC)} is_zero calls recorded, {len(INT)} interiority comparisons recorded', flush=True)

# ---------------- (1) identities and denominators ----------------
nmain = nneg = 0; allbases = set()
for idx, (expr, num, den, qs, rem, tag) in enumerate(ZC):
    resid = sp.expand(num - sum(q*r for q, r in zip(qs, REL)) - rem)
    if not check(resid == 0, f'is_zero #{idx} [{tag}]: num == sum q_i REL_i + rem (exact polynomial identity)'): continue
    if tag == 'negctl':
        nneg += 1
        check(sp.simplify(rem) != 0, f'is_zero #{idx} negative control has rem != 0 (as required)')
        continue
    nmain += 1
    check(sp.simplify(rem) == 0, f'is_zero #{idx} rem == 0, so expression == 0 once denominators are nonzero')
    dv = ev(den)
    check(nonzero(dv), f'is_zero #{idx} combined denominator nonzero: {mag(dv)}  (enclosure excludes 0)')
    allbases |= den_bases(expr)
print(f'{nmain} positive identities, {nneg} negative controls', flush=True)

# every denominator base appearing anywhere in the checked expressions, interiority and null vector
Gs = G['Gs']
others = [0, 2, 3, 4, 5]
Asub = Gs[:, others]; rhs = -Gs[:, 1]
for (_, tang, _) in INT: allbases |= den_bases(tang)
allbases |= den_bases(G['A_']) | den_bases(G['B_'])
for e in G['C'][2]: allbases |= den_bases(e)
bad = 0; worst = None
for b in sorted(allbases, key=sp.default_sort_key):
    bv = ev(sp.expand(b)) if False else ev(b)
    if not nonzero(bv): bad += 1; print('NOT CERTIFIED denominator base', b, mag(bv))
    elif worst is None or abs(bv.mid()) < abs(worst.mid()): worst = bv
check(bad == 0, f'all {len(allbases)} distinct denominator bases (every negative power in the input expressions, A_, B_, tangential coordinates, DL vertices) are nonzero; smallest |value| ~ {mag(abs(worst.mid())) if worst is not None else "n/a"}')

# ---------------- (2) the (a, b) system of the KK construction ----------------
a_, b_ = G['a'], G['b']
p2, p3, L = G['p2'], G['p3'], G['L']
eqs = [p2[1] - p3[1], p2[0] - p3[0] - L/2]
Mab = sp.Matrix([[sp.diff(e, v) for v in (a_, b_)] for e in eqs])
mat = arb_mat(2, 2, [ev(sp.expand(Mab[i, j])) for i in range(2) for j in range(2)])
dab = mat.det()
check(nonzero(dab), f'KK (a, b) linear system determinant nonzero: {mag(dab)}')

# ---------------- (3) Asub invertibility ----------------
print('Asub entries are polynomials in c, s?', all(sp.expand(x).is_polynomial(c_sym, s_sym) for x in Asub), flush=True)
mat5 = arb_mat(5, 5, [ev(sp.expand(Asub[i, j])) for i in range(5) for j in range(5)])
d5 = mat5.det()
check(nonzero(d5), f'det(Asub) (5x5 LUsolve matrix) nonzero: {mag(d5)}; hence rank DG(0) = 5 and ker DG(0) = R v')
Gb = arb_mat(5, 6, [ev(sp.expand(Gs[i, j])) for i in range(5) for j in range(6)])
rb = arb_mat(5, 1, [ev(sp.expand(rhs[i])) for i in range(5)])
sol = mat5.solve(rb)      # rigorous: succeeds only if the interval matrix is nonsingular
vb = [sol[0, 0], arb(1), sol[1, 0], sol[2, 0], sol[3, 0], sol[4, 0]]
res = [sum((Gb[i, j]*vb[j] for j in range(6)), arb(0)) for i in range(5)]
check(all(r.contains(0) for r in res), 'arb null vector v (y1 = 1): G v encloses 0 in all five rows (consistency; v is the unique kernel vector since det(Asub) != 0)')
print('v enclosure:', [x.str(15, radius=False) for x in vb], ' max radius', max(float(x.rad()) for x in vb), flush=True)
try:
    import json
    Vj = json.load(open(os.path.join(here, '..', 'runs', 'charts_w030c.json')))[0]['Jlo'][5]
    dev = max(abs(float(x.mid()) - y) for x, y in zip(vb, Vj))
    print('max |v_arb - v_chart(Jlo row 6)| =', dev, flush=True)
except Exception as ex:
    print('chart comparison skipped:', ex)

# ---------------- (4) contact interiority, strict ----------------
for (label, tang, half) in INT:
    tv = ev(tang); hv = ev(half)
    lo_gap, hi_gap = hv - tv, hv + tv
    check((lo_gap > 0) and (hi_gap > 0), f'interiority {label}: tangential in {mag(tv)}, half-edge s = {mag(hv)}; s - t > 0 ({mag(lo_gap)}), s + t > 0 ({mag(hi_gap)})')

print('elapsed %.0fs' % (time.time() - T0))
if fails:
    print('SOME CHECKS NOT CERTIFIED:'); [print('  ', f) for f in fails]; sys.exit(1)
print('ALL CERTIFIED: denominators nonzero, Asub invertible, contact interiority strict')

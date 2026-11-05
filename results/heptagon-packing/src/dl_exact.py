"""Exact (symbolic) verification of the Kuperberg double-lattice (DL) triangle s2 = (H0, Ht1, H2p2),
in the field Q(c, s), c = cos(pi/7), s = sin(pi/7), with 8c^3 - 4c^2 - 4c + 1 = 0, s^2 = 1 - c^2.

Checks (each is an identity in Q(c, s), verified by clearing denominators and reducing the
numerator to 0 modulo the two defining relations; denominators are checked nonzero):
  (1) the KK half-length parallelogram: p2, p3 at equal height, |p2 p3| = |p1 p4| / 2;
  (2) the five vertex-on-edge contacts used by the chart lemma are exactly 0:
      (0,2,4,1), (1,0,0,3), (1,2,5,2), (2,0,1,4), (2,1,2,5);
  (3) area(s2) = A / (2 delta_DL), with A = (7/2) sin(2pi/7) and
      delta_DL = (2/97)(-111 + 492c - 356c^2) (Kallus' closed form);
  (4) the touching vertices lie strictly inside their edges (numerical, 50 digits; strict
      inequalities with margin > 0.1, so float evaluation of the exact numbers suffices).
"""
import sympy as sp

c, s = sp.symbols('c s')
REL = [s**2 + c**2 - 1, 8*c**3 - 4*c**2 - 4*c + 1]
CNUM = sp.cos(sp.pi/7); SNUM = sp.sin(sp.pi/7)

def T(n):  # cos(n x) as polynomial in cos x
    return sp.chebyshevt(n, c)
def Uc(n):  # sin((n+1) x)/sin x
    return sp.chebyshevu(n, c)
def cos_k(n):   # cos(n pi/7), n integer >= 0
    return sp.expand(T(n))
def sin_k(n):   # sin(n pi/7)
    return sp.expand(s*Uc(n-1)) if n >= 1 else sp.Integer(0)

def is_zero(expr):
    num, den = sp.fraction(sp.together(sp.expand(expr)))
    num = sp.expand(num); den = sp.expand(den)
    _, rem = sp.reduced(num, REL, s, c)
    dval = sp.N(den.subs({c: CNUM, s: SNUM}), 50)
    return sp.simplify(rem) == 0, dval

# heptagon K (circumradius 1): vertex i at angle 2 pi i / 7 = (2i) pi/7
def vert(theta_units, i):
    """vertex i of a heptagon with orientation theta = theta_units * pi/7: angle (theta_units + 2i) pi/7"""
    n = theta_units + 2*i
    return sp.Matrix([cos_k(n), sin_k(n)])
def normal(theta_units, k):
    """outward normal of edge k: angle (theta_units + 2k + 1) pi/7"""
    n = theta_units + 2*k + 1
    return sp.Matrix([cos_k(n), sin_k(n)])

kv = [vert(0, i) for i in range(7)]
u = c
# KK parallelogram: p1 = k0, p4 = (k3 + k4)/2 ; p2 = (1-a) k1 + a k2 ; p3 = (1-b) k2 + b k3
a, b = sp.symbols('a b')
p2 = (1-a)*kv[1] + a*kv[2]; p3 = (1-b)*kv[2] + b*kv[3]
L = 1 + u
sol = sp.solve([sp.Eq(p2[1], p3[1]), sp.Eq(p2[0] - p3[0], L/2)], [a, b], dict=True)[0]
A_, B_ = sp.simplify(sol[a]), sp.simplify(sol[b])
P2 = p2.subs({a: A_, b: B_}); P3 = p3.subs({a: A_, b: B_})
p4 = (kv[3] + kv[4])/2
results = []
results.append(('(1a) p4 = (-c, 0): x', is_zero(p4[0] + c)))
results.append(('(1b) p4 y', is_zero(p4[1])))
results.append(('(1c) p2, p3 same height', is_zero(P2[1] - P3[1])))
results.append(('(1d) |p2 p3| = |p1 p4|/2', is_zero((P2[0] - P3[0]) - L/2)))

# DL triangle s2: H0 at 0 (theta 0), Ht1 at (1 + c, 0) (theta 0), H2p2 at 2 p2 (theta pi/7)
C = [sp.Matrix([0, 0]), sp.Matrix([L, 0]), 2*P2]
TH = [0, 0, 1]     # orientations in units of pi/7
CONTACTS = [(0, 2, 4, 1), (1, 0, 0, 3), (1, 2, 5, 2), (2, 0, 1, 4), (2, 1, 2, 5)]
for (i, j, m, k) in CONTACTS:
    w = C[j] + vert(TH[j], m) - C[i]
    g = (normal(TH[i], k).T * w)[0] - u
    results.append((f'(2) contact {(i, j, m, k)} g = 0', is_zero(g)))
    # (4) tangential position inside the edge: |tangential| < half edge length = tan(pi/7) u = s
    nvec = normal(TH[i], k); tvec = sp.Matrix([-nvec[1], nvec[0]])
    tang = (tvec.T * w)[0]
    tv = sp.N(tang.subs({c: CNUM, s: SNUM}), 50); half = sp.N(SNUM, 50)
    results.append((f'(4) contact {(i, j, m, k)} |tangential| = {sp.N(abs(tv), 12)} < {sp.N(half, 12)} (margin {sp.N(half - abs(tv), 6)})',
                    (bool(abs(tv) < half - sp.Rational(1, 10)), None)))

area = sp.Rational(1, 2)*(C[1][0]*C[2][1] - C[1][1]*C[2][0])
Ahept = sp.Rational(7, 2)*sin_k(2)
delta = sp.Rational(2, 97)*(-111 + 492*c - 356*c**2)
results.append(('(3) area(s2) * 2 delta_DL - A = 0', is_zero(area*2*delta - Ahept)))
print('area(s2) =', sp.N(area.subs({c: CNUM, s: SNUM}), 30), ' delta_DL =', sp.N(delta.subs({c: CNUM}), 30))
ok_all = True
for name, (ok, d) in results:
    print(('PASS ' if ok else 'FAIL ') + name + ('' if d is None else f'   [denominator {sp.N(d, 6)}]'))
    ok_all &= bool(ok)
print('ALL PASS' if ok_all else 'SOME FAILED')

# negative controls: the reducer must reject false identities
for name, e in (('c - 9/10', c - sp.Rational(9, 10)), ('area*2*delta - A + 1/1000', area*2*delta - Ahept + sp.Rational(1, 1000)),
                ('contact (0,2,4,1) g + s/10^6', (normal(TH[0], 1).T*(C[2] + vert(TH[2], 4) - C[0]))[0] - u + s/10**6)):
    print('negative control', name, '-> zero?', is_zero(e)[0], '(must be False)')

# ---- H-DL extras: exact contact Jacobian, exact null vector, translation edge ----
# variables z = (x1, y1, phi1, x2, y2, phi2) perturb Ht1 (index 1) and H2p2 (index 2); H0 fixed.
def rot90(v): return sp.Matrix([-v[1], v[0]])
Grows = []
for (i, j, m, k) in CONTACTS:
    n = normal(TH[i], k); w = C[j] + vert(TH[j], m) - C[i]; vm = vert(TH[j], m)
    row = [sp.Integer(0)]*6
    def off(idx): return {1: 0, 2: 3}[idx]
    # d/dC_j : + n ; d/dC_i : - n ; d/dphi_i : rot90(n).w ; d/dphi_j : n.rot90(vm)
    if j != 0:
        row[off(j)] += n[0]; row[off(j)+1] += n[1]; row[off(j)+2] += (n.T*rot90(vm))[0]
    if i != 0:
        row[off(i)] -= n[0]; row[off(i)+1] -= n[1]; row[off(i)+2] += (rot90(n).T*w)[0]
    Grows.append([sp.expand(e) for e in row])
Gs = sp.Matrix(Grows)
# exact null vector with y1-component 1: solve the 5x5 system on the other components
others = [0, 2, 3, 4, 5]
Asub = Gs[:, others]; rhs = -Gs[:, 1]
vsol = Asub.LUsolve(rhs)
vex = [vsol[0], sp.Integer(1), vsol[1], vsol[2], vsol[3], vsol[4]]
for r_ in range(5):
    e = sum(Gs[r_, q]*vex[q] for q in range(6))
    ok, d = is_zero(e)
    print(('PASS ' if ok else 'FAIL ') + f'(5) G v = 0, row {r_}')
num_v = [float(sp.N(x.subs({c: CNUM, s: SNUM}), 30)) for x in vex]
print('exact null vector v (y1 = 1):', [round(x, 12) for x in num_v])
import json as _json
try:
    V_json = _json.load(open('runs/charts_w030c.json'))[0]['Jlo'][5]
    print('max |v_exact - v_chart(lo)| =', max(abs(a - b) for a, b in zip(num_v, V_json)))
except Exception as ex:
    print('chart file comparison skipped:', ex)
# translation edge 0 -> 1: r = 1 + c, both orientations 0, edge direction phi = 0  ->  (1 + c, 0, 0)
print(('PASS ' if is_zero(C[1][0] - (1 + c))[0] and C[1][1] == 0 else 'FAIL ') + '(6) translation edge configuration is exactly (1 + c, 0, 0)')

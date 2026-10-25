"""Rigorous certificate for the explicit Kallus neighbourhood (see PROOF.md).

Every inequality used in the proof that involves a number is checked here
with ball arithmetic (python-flint / Arb: every real is an interval with
outward-rounded radius).  Exact rationals (fmpq) are used where possible.

Usage:  python certify.py [delta]      (default delta = DELTA_STAR below)

Exit code 0 iff every check passes.
"""
import sys
from fractions import Fraction
from flint import arb, fmpq, ctx

ctx.prec = 256
LMAX = 20000           # degrees l <= LMAX summed exactly (in Arb), tail bounded
DELTA_STAR = fmpq(6, 10**8)

fails = []


def check(cond, msg):
    ok = bool(cond)
    print(("PASS  " if ok else "FAIL  ") + msg)
    if not ok:
        fails.append(msg)


def ub(x):   # rigorous upper end of an arb
    return x.mid() + x.rad()


# ---------------------------------------------------------------- Lemma 2
# D3 minimal vectors (scaled to the unit sphere), one per antipodal pair.
X = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
# all vectors carry a factor 1/sqrt2; quadratic forms pick up 1/2.


def qform(Q, x):   # <x,Qx> for x = v/sqrt2, exact
    return Fraction(sum(Q[a][b] * x[a] * x[b] for a in range(3) for b in range(3)), 2)


def solve_dual(k):
    """Symmetric Q with <x_i,Q x_i> = delta_ik (formulas of Lemma 2)."""
    rho = [Fraction(int(i == k)) for i in range(6)]
    s12, s13, s23 = rho[0] + rho[1], rho[2] + rho[3], rho[4] + rho[5]
    Q = [[None] * 3 for _ in range(3)]
    Q[0][0] = (s12 + s13 - s23) / 2
    Q[1][1] = (s12 + s23 - s13) / 2
    Q[2][2] = (s13 + s23 - s12) / 2
    Q[0][1] = Q[1][0] = (rho[0] - rho[1]) / 2
    Q[0][2] = Q[2][0] = (rho[2] - rho[3]) / 2
    Q[1][2] = Q[2][1] = (rho[4] - rho[5]) / 2
    return Q


duals = [solve_dual(k) for k in range(6)]
ok_dual = all(qform(duals[k], X[i]) == (1 if i == k else 0)
              for k in range(6) for i in range(6))
check(ok_dual, "Lemma 2: dual basis Q^(k) satisfies <x_i,Q^(k)x_i> = delta_ik (exact)")
check(all(sum(duals[k][a][a] for a in range(3)) == Fraction(1, 2) for k in range(6)),
      "Lemma 2: tr Q^(k) = 1/2 for every k  (eutaxy coefficients 1/4) (exact)")


def perp_sq(Q, x):  # |Qx - <x,Qx>x|^2 for x = v/sqrt2, exact rational
    Qv = [Fraction(sum(Q[a][b] * x[b] for b in range(3))) for a in range(3)]  # Q v
    # Qx = Qv/sqrt2 ; <x,Qx> x = qform * v/sqrt2
    q = qform(Q, x)
    w = [Qv[a] - q * x[a] for a in range(3)]   # sqrt2 * (Qx - <x,Qx>x)
    return sum(c * c for c in w) / 2


g1sq = max(perp_sq(duals[k], X[i]) for k in range(6) for i in range(6))
check(g1sq == Fraction(3, 8), "Lemma 2: g1^2 = max_{k,i} |Q^(k)_perp x_i|^2 = 3/8 (exact)")
G2 = max(sum(arb(fmpq(perp_sq(duals[k], X[i]).numerator,
                     perp_sq(duals[k], X[i]).denominator)).sqrt() for i in range(6))
         for k in range(6))
print("      G2 = max_k sum_i |Q^(k)_perp x_i| =", G2, "(= sqrt 6)")
check(abs(G2 - arb(6).sqrt()) < arb(10) ** -60, "Lemma 2: G2 = sqrt(6)")
# operator norm of Q^(k): characteristic polynomial roots; every Q^(k) is a
# signed permutation conjugate of diag-block [[1/2,1/2],[1/2,1/2]] (+) [-1/2]
# -> eigenvalues {1, 0, -1/2}.  Check: Q^3 - (tr)Q^2 + e2 Q - det = 0 with
# tr=1/2, e2=-1/2, det=0.
def matmul(A, B):
    return [[sum(A[i][t] * B[t][j] for t in range(3)) for j in range(3)] for i in range(3)]


ok_eig = True
for Q in duals:
    Q2 = matmul(Q, Q)
    Q3 = matmul(Q2, Q)
    tr = sum(Q[a][a] for a in range(3))
    e2 = (tr * tr - sum(Q2[a][a] for a in range(3))) / 2
    det = (Q[0][0] * (Q[1][1] * Q[2][2] - Q[1][2] * Q[2][1])
           - Q[0][1] * (Q[1][0] * Q[2][2] - Q[1][2] * Q[2][0])
           + Q[0][2] * (Q[1][0] * Q[2][1] - Q[1][1] * Q[2][0]))
    ok_eig &= (tr, e2, det) == (Fraction(1, 2), Fraction(-1, 2), 0)
check(ok_eig, "Lemma 2: char. poly of every Q^(k) is t(t-1)(t+1/2) => ||Q^(k)||_op = 1 (exact)")
qop = arb(1)

# ---------------------------------------------------------------- Lemma 6
# c_l = 1 + 4 P_l(1/2) + P_l(0).  Exact rationals for l <= 60, Arb beyond.
Ph_q = [fmpq(1), fmpq(1, 2)]
P0_q = [fmpq(1), fmpq(0)]
for l in range(1, 60):
    Ph_q.append((fmpq(2 * l + 1) * fmpq(1, 2) * Ph_q[l] - l * Ph_q[l - 1]) / (l + 1))
    P0_q.append((-l * P0_q[l - 1]) / (l + 1))
c_q = {l: 1 + 4 * Ph_q[l] + P0_q[l] for l in range(0, 61, 2)}
check(c_q[2] == 0 and c_q[4] == fmpq(7, 32) and c_q[10] == fmpq(65, 65536),
      "Lemma 6: c_2 = 0, c_4 = 7/32, c_10 = 65/65536 (exact)")

# Exact rational recurrence for all l <= LMAX (Arb three-term recurrences
# lose precision exponentially, so rationals are used and converted once).
Phx, P0x = Ph_q[:2], P0_q[:2]
for l in range(1, LMAX):
    Phx.append((fmpq(2 * l + 1, 2) * Phx[l] - l * Phx[l - 1]) / (l + 1))
    P0x.append((-l * P0x[l - 1]) / (l + 1))
c = {l: arb(1 + 4 * Phx[l] + P0x[l]) for l in range(4, LMAX + 1, 2)}
# consistency with exact values
check(all(abs(c[l] - arb(c_q[l])) < arb(10) ** -50 for l in range(4, 61, 2)),
      "Lemma 6: long recurrence agrees with the l <= 60 table")
cmin_l = min(range(4, LMAX + 1, 2), key=lambda l: c[l].mid())
check(all(c[l] > 0 for l in c),
      f"Lemma 6: c_l > 0 for all even 4 <= l <= {LMAX} (min at l={cmin_l}: {c[cmin_l].str(8)})")
print("      worst other c_l: ",
      sorted(((float(c[l].mid()), l) for l in c if l != 10))[:3])

# S = sum_{l>=4 even} (2l+1) (1-c_l)^6 / (1.5 c_l)
S_fin = arb(0)
for l in range(4, LMAX + 1, 2):
    S_fin += (2 * l + 1) * (1 - c[l]) ** 6 / (arb(3) / 2 * c[l])
# tail: Szego Thm 7.3.3 (Bernstein): sqrt(sin th)|P_n(cos th)| < sqrt(2/(pi n))
pi = arb.pi()
Bconst = (2 / pi).sqrt() * (4 * (2 / arb(3).sqrt()).sqrt() + 1)
cmin_tail = 1 - Bconst / arb(LMAX).sqrt()
check(cmin_tail > 0, f"Lemma 6: tail c_l >= 1 - B/sqrt(l) > 0 for l > {LMAX} (B={Bconst.str(6)})")
S_tail = (2 + arb(1) / LMAX) * Bconst ** 6 / (arb(3) / 2 * cmin_tail) / (2 * arb(LMAX))
S_tot = S_fin + S_tail
print("      S_finite =", S_fin.str(10), " S_tail <=", S_tail.str(6))
TVp = arb(31)    # TV(delta - nu + nu*nu) <= 1 + 5 + 25
Neff = 2 * TVp + S_tot.sqrt()
Neff_ub = ub(Neff)
print("      N_eff = 2*31 + sqrt(S) <=", arb(Neff_ub).str(10))
print("      (fiber-only bound of Kallus's Lemma 6 would be >= 2/c_10 =",
      (2 / c[10]).str(8), ")")

# ---------------------------------------------------------------- Lemma 4 (normalisation)
delta = fmpq(Fraction(sys.argv[1]).numerator, Fraction(sys.argv[1]).denominator) \
    if len(sys.argv) > 1 else DELTA_STAR
d = arb(delta)
eta = (1 + d).log() / 2
L0 = ((1 + d) ** 2 - 1).sqrt()
P2 = arb(10) / (3 * arb(3).sqrt())   # ||pi_2||_{inf->inf} = 5 * int_{-1}^1 |P_2| dt/2


def E1(a):
    Ra = 1 / (1 - a) ** 2 - 1 - 2 * a
    z = 2 * a + Ra
    return Ra / 2 + z ** 2 / (4 * (1 - z))


def theta(a):
    return (a / (1 - a)).asin()


a = P2 * eta
for _ in range(50):
    a = P2 * (eta + E1(a) + L0 * theta(a))
a = arb(ub(a)) * (1 + arb(10) ** -6)       # inflate slightly, then verify
rhs = P2 * (eta + E1(a) + L0 * theta(a))
check(rhs <= a and a < arb(1) / 100,
      f"Lemma 4: Brouwer ball radius a = {a.str(6)} satisfies P2(eta+E1(a)+L0*theta(a)) <= a")
eps = (1 + d).log() + 2 * (a + E1(a))
L = L0 * (1 + a) / (1 - a) + ((1 + a) ** 2 / (1 - a) ** 2 - 1) / 2
print(f"      delta = {d.str(6)}:  eps <= {eps.str(6)},  L <= {L.str(8)}  (L0 = {L0.str(8)})")

# ---------------------------------------------------------------- Lemma 3 (construction)
g1 = arb(fmpq(3, 8)).sqrt()
e = eps.exp()
Ybar = 2 * g1 * 6 * (e - 1)
cbar = e ** 2 * L * G2 * (L * e * Ybar).exp()
check(cbar <= arb(1) / 2, f"Lemma 3: cbar = {cbar.str(8)} <= 1/2")
check((1 - 12 * eps * e) * arb(2).sqrt() >= e,
      "Lemma 3: non-minimal lattice vectors stay outside e^eps B")
kappa = eps * e / 4 + e * cbar / (2 * (1 - cbar)) + 72 * eps * e ** 2 + 288 * eps ** 2 * e ** 3
print("      kappa <=", kappa.str(10), "   1/(12 N_eff) >=", (1 / (12 * arb(Neff_ub))).str(10))

# ---------------------------------------------------------------- Theorem
check(kappa * 12 * arb(Neff_ub) < 1,
      f"THEOREM: 12 * kappa * N_eff = {(kappa * 12 * arb(Neff_ub)).str(8)} < 1 for delta = {delta}")

print()
print("ALL CHECKS PASSED" if not fails else f"{len(fails)} CHECK(S) FAILED")
sys.exit(1 if fails else 0)

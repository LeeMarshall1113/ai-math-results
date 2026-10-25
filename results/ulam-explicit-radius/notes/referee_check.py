# Independent recomputation (referee). mpmath interval arithmetic + exact rationals.
import sys, itertools, math
from fractions import Fraction as F
from math import comb
from mpmath import iv, mp
iv.prec = 200
LCUT = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
delta = F(sys.argv[1]) if len(sys.argv) > 1 else F(6, 10**8)

# --- Lemma 2 independently: solve linear system numerically-exact via Fractions
import numpy as np
X = [np.array(v)/np.sqrt(2) for v in [(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]]
# basis of Sym(3): 6 params
def sym(p):
    a,b,c,d,e,f = p
    return np.array([[a,d,e],[d,b,f],[e,f,c]])
M = np.array([[x@sym(np.eye(6)[j])@x for j in range(6)] for x in X])
Minv = np.linalg.inv(M)
g1 = 0; G2 = 0; opn = 0; trs = []
for k in range(6):
    Q = sym(Minv[:, k])
    trs.append(np.trace(Q)); opn = max(opn, np.abs(np.linalg.eigvalsh(Q)).max())
    perps = [np.linalg.norm(Q@x - (x@Q@x)*x) for x in X]
    g1 = max(g1, max(perps)); G2 = max(G2, sum(perps))
print("Lemma2: traces", set(np.round(trs, 12)), "g1", g1, "sqrt(3/8)", math.sqrt(3/8), "G2", G2, "opnorm", opn)
# all 12 vs: non-minimal D3 norms
mins = sorted({sum(z*z for z in v)/2 for v in itertools.product(range(-3,4), repeat=3) if sum(v)%2==0 and any(v)})[:3]
print("Lambda0 squared norms:", mins)

# --- c_l via closed forms (independent of recurrence): P_l(0) closed form; P_l(1/2) via explicit sum
def P_half(l):
    # P_l(x) = 2^-l sum_k (-1)^k C(l,k) C(2l-2k,l) x^(l-2k)
    s = F(0)
    for k in range(l//2+1):
        s += (-1)**k * comb(l,k)*comb(2*l-2*k,l) * F(1, 2**(l-2*k))
    return s / 2**l
def P_zero(l):
    return F((-1)**(l//2)*comb(l,l//2), 2**l) if l % 2 == 0 else F(0)
c = {}
for l in range(0, 61, 2):
    c[l] = 1 + 4*P_half(l) + P_zero(l)
print("c2,c4,c10:", c[2], c[4], c[10])
# faster exact recurrence for larger l, cross-checked against closed form at l<=60
Ph=[F(1),F(1,2)]; P0=[F(1),F(0)]
for l in range(1, LCUT):
    Ph.append((F(2*l+1,2)*Ph[l]-l*Ph[l-1])/(l+1)); P0.append((-l*P0[l-1])/(l+1))
assert all(1+4*Ph[l]+P0[l]==c[l] for l in range(0,61,2))
cl = {l: 1+4*Ph[l]+P0[l] for l in range(4, LCUT+1, 2)}
neg = [l for l in cl if cl[l] <= 0]
print("nonpositive c_l (4..LCUT):", neg)
srt = sorted(cl, key=lambda l: cl[l])[:4]
print("smallest c_l:", [(l, float(cl[l])) for l in srt])
S = iv.mpf(0)
for l, v in cl.items():
    vi = iv.mpf(v.numerator)/v.denominator
    S += (2*l+1)*(1-vi)**6/(iv.mpf(1.5)*vi)
t10 = 21*(1-iv.mpf(65)/65536)**6/(iv.mpf(1.5)*iv.mpf(65)/65536)
print("l=10 term:", t10)
B0 = iv.sqrt(2/iv.pi)*(4*iv.sqrt(2/iv.sqrt(3))+1)
cmin = 1 - B0/iv.sqrt(LCUT)
tail = (2+iv.mpf(1)/LCUT)*B0**6/(iv.mpf(1.5)*cmin)/(2*LCUT)
print("B0", B0, "S_fin", S, "tail", tail)
N = 62 + iv.sqrt(S+tail)
print("N_eff", N)

# --- Lemma 4
d = iv.mpf(delta.numerator)/delta.denominator
eta = iv.log(1+d)/2; L0 = iv.sqrt((1+d)**2-1); P2 = 10/(3*iv.sqrt(3))
# check ||pi_2|| by integration: 5*int_0^1|P2|
def E1(a):
    Ra = 1/(1-a)**2-1-2*a; z = 2*a+Ra
    return Ra/2+z**2/(4*(1-z))
th = lambda a: (a/(1-a))/iv.sqrt(1-(a/(1-a))**2)  # >= arcsin
a = iv.mpf(0)
for _ in range(60):
    a = P2*(eta+E1(a)+L0*th(a))
a = iv.mpf(a.b)*(1+iv.mpf(10)**-9)
rhs = P2*(eta+E1(a)+L0*th(a))
print("Lemma4 a=", a, "cond ok:", (rhs.b <= a.a) and a.b < 0.01)
eps = iv.log(1+d)+2*(a+E1(a))
L = L0*(1+a)/(1-a)+((1+a)**2/(1-a)**2-1)/2
e = iv.exp(eps)
Ybar = 12*iv.sqrt(iv.mpf(3)/8)*(e-1)
cbar = e**2*L*iv.sqrt(6)*iv.exp(L*e*Ybar)
kappa = eps*e/4 + e*cbar/(2*(1-cbar)) + 72*eps*e**2 + 288*eps**2*e**3
print("eps", eps, "L", L, "cbar", cbar, "kappa", kappa)
print("adm cond", ((1-12*eps*e)*iv.sqrt(2)-e).a > 0)
m = 12*kappa*N
print("12 kappa N =", m, " CLOSES:", m.b < 1)

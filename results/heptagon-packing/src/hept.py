"""Geometry of congruent regular heptagons (circumradius 1). Float only (exploration)."""
import numpy as np
N = 7
U = np.cos(np.pi/7); V = np.sin(np.pi/7)
ANG = 2*np.pi*np.arange(N)/N
K = np.stack([np.cos(ANG), np.sin(ANG)], 1)          # reference heptagon vertices k_i
AREA = 0.5*N*np.sin(2*np.pi/N)                       # area of K
INR = U                                              # inradius

def rot(t):
    c, s = np.cos(t), np.sin(t); return np.array([[c, -s], [s, c]])

def verts(c, t, refl=False):
    """vertices of heptagon centred at c, rotated by t (refl: point-reflected first)."""
    k = -K if refl else K
    return k @ rot(t).T + c

def sep(P, Q):
    """signed separation of convex polygons (SAT over edge normals of both).
    >=0 iff interiors disjoint. Value = max over axes of gap."""
    best = -np.inf
    for A, B in ((P, Q), (Q, P)):
        E = np.roll(A, -1, 0) - A
        Nn = np.stack([E[:, 1], -E[:, 0]], 1)            # outward normals (ccw polygon)
        Nn /= np.linalg.norm(Nn, axis=1)[:, None]
        pa = (A @ Nn.T).max(0); qb = (B @ Nn.T).min(0)
        best = max(best, (qb - pa).max())
    return best

def tri_area(a, b, c):
    return 0.5*((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))

def circumradius(a, b, c):
    la, lb, lc = np.linalg.norm(b-c), np.linalg.norm(c-a), np.linalg.norm(a-b)
    ar = abs(tri_area(a, b, c))
    return la*lb*lc/(4*ar) if ar > 0 else np.inf

def kk_double_lattice():
    """Kuperberg double lattice for the heptagon, following Kallus-Kusner sec 3.2."""
    from scipy.optimize import brentq
    k = K
    p1 = k[0]; p4 = 0.5*(k[3]+k[4])
    L = np.linalg.norm(p1-p4)
    def p2(a): return (1-a)*k[1]+a*k[2]
    def p3(b): return (1-b)*k[2]+b*k[3]
    # same height, horizontal distance L/2
    def b_of_a(a):
        y = p2(a)[1]
        return (k[2][1]-y)/(k[2][1]-k[3][1])
    f = lambda a: (p2(a)[0]-p3(b_of_a(a))[0]) - L/2
    a = brentq(f, 0, 1); b = b_of_a(a)
    P2, P3 = p2(a), p3(b)
    P6 = P2*np.array([1, -1]); P5 = P3*np.array([1, -1])
    par_area = abs(2*tri_area(P2, P3, P6))
    return dict(a=a, b=b, p1=p1, p4=p4, p2=P2, p3=P3, p5=P5, p6=P6, par_area=par_area,
                t1=p1-p4, t2=2*(P2-P6))

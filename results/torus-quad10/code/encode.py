"""CNF encoding of NECESSARY conditions for a realisation of the map M = {4,4}_(3,1)
(faces F_i = (i, i+1, i+8, i+7) mod 10) as a polyhedral torus in R^3.

Variables: for each 4-subset T (sorted) of the 10 points, chi(T) in {+,0,-} with
chi(p1..p4) = sign det[1 p1; 1 p2; 1 p3; 1 p4] (affine chirotope, rank 4).
Every clause below is implied by any genuine realisation (see PROOF-SKELETON.md, Lemma 4),
so UNSAT => no realisation.  Clause groups (toggle with flags for ablation):
  GP   3-term Grassmann-Pluecker sign relations (all sigma of size 2, all 4 others)
  PL   faces coplanar
  CV   strict convexity: chi(a,b,c,e)=chi(b,c,d,e)=chi(c,d,a,e)=chi(d,a,b,e) for all e;
       and abc spans a plane (some e off it) -- implied by AD below
  AD   edge-adjacent faces (a,b,c,d),(b,a,x,y): chi(a,b,c,x) = chi(a,b,c,y) != 0
  INT  forbidden intersection patterns for vertex-adjacent and vertex-disjoint face pairs
  SYM  mirror symmetry breaking (chi -> -chi is a symmetry of every clause above)
"""
import itertools, sys

N = 10
F = [tuple((i + d) % N for d in (0, 1, 8, 7)) for i in range(N)]

class Enc:
    def __init__(self):
        self.nv = 0; self.cls = []; self.var = {}; self.prod = []
        for T in itertools.combinations(range(N), 4):
            p = self.new(); n = self.new(); z = self.new()  # pos, neg, nonzero
            self.var[T] = (p, n, z)
            self.cls += [[-p, -n], [-z, p, n], [-p, z], [-n, z]]
    def new(self):
        self.nv += 1; return self.nv
    def lit(self, T, s):
        """literal for 'chi(T) == s' (s=+1/-1) or 'chi(T) != 0' (s=0 means nonzero literal)"""
        T = list(T); assert len(set(T)) == 4
        par = 0
        for i in range(4):
            for j in range(i + 1, 4):
                if T[i] > T[j]: par ^= 1
        p, n, z = self.var[tuple(sorted(T))]
        if s == 0: return z
        if par: s = -s
        return p if s > 0 else n
    # condition atoms, each returns a list of literals whose conjunction is the condition
    def eq(self, T, s):       # chi(T) == s, s in {+1,-1,0}
        return [-self.lit(T, 0)] if s == 0 else [self.lit(T, s)]
    def nz(self, T): return [self.lit(T, 0)]
    def notv(self, T, s):     # chi(T) != s  (s=+1/-1)
        return [-self.lit(T, s)]
    def forbid(self, *conds):
        c = []
        for cond in conds: c += [-l for l in cond]
        self.cls.append(sorted(set(c)))
    def add(self, c): self.cls.append(c)

def gp(E):
    # chi(s,y1,y2)chi(s,y3,y4) - chi(s,y1,y3)chi(s,y2,y4) + chi(s,y1,y4)chi(s,y2,y3) = 0 (real Pluecker)
    prod = {}
    def P(A, B):   # aux literals (pp, pn): product sign + / -
        key = (A, B)
        if key in prod: return prod[key]
        pp, pn = E.new(), E.new()
        ap, an, bp, bn = E.lit(A, 1), E.lit(A, -1), E.lit(B, 1), E.lit(B, -1)
        E.cls += [[-ap, -bp, pp], [-an, -bn, pp], [-ap, -bn, pn], [-an, -bp, pn],
                  [-pp, ap, an], [-pp, bp, bn], [-pp, -ap, bp], [-pp, -an, bn],
                  [-pn, ap, an], [-pn, bp, bn], [-pn, -ap, bn], [-pn, -an, bp]]
        prod[key] = (pp, pn); E.prod.append((pp, pn, A, B)); return pp, pn
    cnt = 0
    for s in itertools.combinations(range(N), 2):
        rest = [x for x in range(N) if x not in s]
        for y in itertools.combinations(rest, 4):
            y1, y2, y3, y4 = y
            t1 = P(s + (y1, y2), s + (y3, y4))
            a = P(s + (y1, y3), s + (y2, y4)); t2 = (a[1], a[0])  # negated term
            t3 = P(s + (y1, y4), s + (y2, y3))
            for (tp, tn) in (t1, t2, t3):
                E.add([-tp, t1[1], t2[1], t3[1]])
                E.add([-tn, t1[0], t2[0], t3[0]])
            cnt += 1
    return cnt

def edges(f): return [(f[i], f[(i + 1) % 4]) for i in range(4)]
def tris(f): a, b, c, d = f; return [(a, b, c), (a, c, d), (a, b, d), (b, c, d)]

def seg_meets_tri(E, p, q, a, b, c):
    """forbid: closed segment pq meets closed triangle abc (abc assumed non-degenerate), with
    the line pq not contained in plane abc.  Pattern: s1=chi(abcp), s2=chi(abcq), s1*s2<=0,
    not both 0; t1..t3 = chi(pqab),chi(pqbc),chi(pqca) all in {tau,0}."""
    for (s1, s2) in [(1, -1), (-1, 1), (0, 1), (0, -1), (1, 0), (-1, 0)]:
        for tau in (1, -1):
            E.forbid(E.eq((a, b, c, p), s1), E.eq((a, b, c, q), s2),
                     E.notv((p, q, a, b), -tau), E.notv((p, q, b, c), -tau), E.notv((p, q, c, a), -tau))

def pt_in_poly_coplanar(E, x, f, others):
    """forbid: x in plane(f) and x in closed convex polygon f.  For a reference e off plane(f):
    chi(a,b,x,e),... each in {sigma,0} with sigma = chi(a,b,c,e)."""
    a, b, c, d = f
    for e in others:
        if e == x: continue
        for sg in (1, -1):
            E.forbid(E.eq((a, b, c, x), 0), E.eq((a, b, c, e), sg),
                     E.notv((a, b, x, e), -sg), E.notv((b, c, x, e), -sg),
                     E.notv((c, d, x, e), -sg), E.notv((d, a, x, e), -sg))

def seg_cross_coplanar(E, p, q, a, b, others):
    """forbid proper crossing of coplanar segments pq, ab (4 distinct points)."""
    for e in others:
        for s in (1, -1):
            for t in (1, -1):
                E.forbid(E.eq((p, q, a, b), 0), E.eq((p, q, a, e), s), E.eq((p, q, b, e), -s),
                         E.eq((a, b, p, e), t), E.eq((a, b, q, e), -t))

def ray_in_cone(E, v, a, x, y, others):
    """forbid: a in plane(v,x,y) and ray v->a in the closed cone of v->x, v->y (so segment va
    meets triangle vxy beyond v)."""
    for e in others:
        for sg in (1, -1):
            E.forbid(E.eq((v, x, y, a), 0), E.eq((v, x, y, e), sg),
                     E.notv((v, x, a, e), -sg), E.notv((v, a, y, e), -sg))

def build(groups=("GP", "PL", "CV", "AD", "INT", "SYM"), faces=None):
    global F, N
    if faces is not None:
        F = [tuple(f) for f in faces]; N = len({u for f in F for u in f})
    E = Enc(); info = {}; E.sym = None
    if "GP" in groups: info["gp_rel"] = gp(E)
    if "PL" in groups:
        for f in F: E.forbid(E.nz(f))
    if "CV" in groups:
        for f in F:
            a, b, c, d = f
            for e in range(N):
                if e in f: continue
                cyc = [(a, b, c, e), (b, c, d, e), (c, d, a, e), (d, a, b, e)]
                for s in (1, -1, 0):
                    for T1, T2 in zip(cyc, cyc[1:] + cyc[:1]):
                        E.forbid(E.eq(T1, s), E.notv(T2, s) if s else E.nz(T2))
    if "AD" in groups:
        for i, f in enumerate(F):
            for j, g in enumerate(F):
                if i == j: continue
                s = set(f) & set(g)
                if len(s) != 2: continue
                a, b = [u for u in f if u in s]
                c = [u for u in f if u not in s][0]
                x, y = [u for u in g if u not in s]
                E.forbid(E.eq((a, b, c, x), 0)); E.forbid(E.eq((a, b, c, y), 0))
                for sg in (1, -1): E.forbid(E.eq((a, b, c, x), sg), E.notv((a, b, c, y), sg))
                if E.sym is None: E.sym = (a, b, c, x)
    if "INT" in groups:
        for i, j in itertools.permutations(range(N), 2):
            P, Q = F[i], F[j]
            s = set(P) & set(Q)
            if len(s) == 2: continue
            others_P = [e for e in range(N) if e not in P]
            if not s:
                for (p, q) in edges(P):
                    for T in tris(Q): seg_meets_tri(E, p, q, *T)
                for x in Q: pt_in_poly_coplanar(E, x, P, others_P)
                for (p, q) in edges(P):
                    for (a, b) in edges(Q):
                        seg_cross_coplanar(E, p, q, a, b, [e for e in range(N) if e not in (p, q, a, b)])
            else:
                v = s.pop()
                for (p, q) in edges(P):
                    if v in (p, q):
                        a = q if p == v else p
                        k = Q.index(v); x, w, y = Q[(k + 1) % 4], Q[(k + 2) % 4], Q[(k + 3) % 4]
                        for (u1, u2) in [(x, w), (w, y), (x, y)]:
                            ray_in_cone(E, v, a, u1, u2, [e for e in range(N) if e not in (v, a, u1, u2)])
                    else:
                        for T in tris(Q): seg_meets_tri(E, p, q, *T)
                        for (a, b) in edges(Q):
                            if v not in (a, b):
                                seg_cross_coplanar(E, p, q, a, b, [e for e in range(N) if e not in (p, q, a, b)])
                for x in Q:
                    if x != v: pt_in_poly_coplanar(E, x, P, others_P)
    if "SYM" in groups:
        # E.sym = (a,b,c,x) from the first AD pair, so chi(E.sym) != 0 in any realisation;
        # the mirror image flips every sign, so WLOG it is +.
        assert E.sym is not None, "SYM needs AD"
        E.add([E.lit(E.sym, 1)])
    return E, info

def write(E, path):
    with open(path, "w") as fh:
        fh.write(f"p cnf {E.nv} {len(E.cls)}\n")
        for c in E.cls: fh.write(" ".join(map(str, c)) + " 0\n")

if __name__ == "__main__":
    groups = sys.argv[2].split(",") if len(sys.argv) > 2 else ("GP", "PL", "CV", "AD", "INT", "SYM")
    E, info = build(groups)
    write(E, sys.argv[1])
    print(sys.argv[1], "vars", E.nv, "clauses", len(E.cls), info)

def chirotope(X):
    """exact signs from rational/integer coordinates X (list of 3-tuples of Fractions/ints)"""
    from fractions import Fraction
    def det3(u, v, w):
        return (u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0]))
    chi = {}
    for T in itertools.combinations(range(len(X)), 4):
        p = [tuple(Fraction(c) for c in X[t]) for t in T]
        d = det3(*[tuple(p[k][m]-p[0][m] for m in range(3)) for k in (1, 2, 3)])
        chi[T] = (d > 0) - (d < 0)
    return chi

def assignment(E, chi):
    val = {}
    for T, (p, n, z) in E.var.items():
        val[p] = chi[T] > 0; val[n] = chi[T] < 0; val[z] = chi[T] != 0
    for (pp, pn, A, B) in E.prod:
        def sg(T):
            par = sum(1 for i in range(4) for j in range(i+1, 4) if T[i] > T[j]) % 2
            return chi[tuple(sorted(T))] * (-1 if par else 1)
        s = sg(A) * sg(B); val[pp] = s > 0; val[pn] = s < 0
    return val

def violated(E, val, maxshow=5):
    bad = [c for c in E.cls if not any((val[abs(l)] if l > 0 else not val[abs(l)]) for l in c)]
    return bad

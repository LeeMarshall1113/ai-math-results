"""CNF for "K_n has a straight-line drawing whose edges split into t crossing-free colour classes".

Signotope (pseudo-configuration) encoding; ingredients from Heule-Scheucher (TACAS 2024).

Points 1..n sorted by x. X(a,b,c), a<b<c, True = positive (counter-clockwise) orientation.
Every point set in general position (needed: K_n drawn straight-line => no 3 collinear) gives a
signotope; so UNSAT over signotopes => no point set works.

Symmetry breaking (all optional, all sound for point sets):
  p1   : X(1,a,b)=True for all 1<a<b  [HS Lemma: p1 a hull vertex, others sorted around it;
         reached by a projective map that keeps every point on one side of the line sent to
         infinity, hence maps segments to segments and preserves all orientations/crossings].
  lex  : HS reflection predicate on consecutive triples (HS eq. (20)).
  col  : colour value-precedence over a fixed edge order (colours are interchangeable).

Crossing: for a<b<c<d only ac|bd and ad|bc can cross (ab, cd have disjoint x-ranges).
  seg pq crosses rs  <=>  o(p,q,r) != o(p,q,s)  and  o(r,s,p) != o(r,s,q).
  Y(pq,rs) <-> crossing (full definition); Y -> not(C(pq,k) and C(rs,k)) for every colour k.

Euler (optional, --euler): every colour class is a crossing-free geometric graph, so it has at
  most 3n-3-h edges (h = hull size). Encoded as: H_v <-> v is a hull vertex (v not inside any
  triangle; for p1-normal form p1 is a hull vertex), and for each colour k,
  #edges(k) + #hull <= 3n-3, via a sequential-counter / totalizer at-most constraint.
"""
import itertools
import sys


class Enc:
    def __init__(self, n, t, p1=True, lex=True, col=True, euler=False, ydef="full", fix=None, tri=False):
        self.tri = tri
        self.n, self.t = n, t
        self.nv = 0
        self.cl = []
        self.X = {tr: self.new() for tr in itertools.combinations(range(1, n + 1), 3)}
        self.signotope()
        if p1:
            for a in range(2, n + 1):
                for b in range(a + 1, n + 1):
                    self.cl.append([self.X[(1, a, b)]])
        if lex:
            self.lex_reflection()
        self.E = list(itertools.combinations(range(1, n + 1), 2))
        self.C = {(e, k): self.new() for e in self.E for k in range(t)}
        for e in self.E:
            self.cl.append([self.C[(e, k)] for k in range(t)])
            if not tri:
                for k1, k2 in itertools.combinations(range(t), 2):
                    self.cl.append([-self.C[(e, k1)], -self.C[(e, k2)]])
        self.Y = {}
        for a, b, c, d in itertools.combinations(range(1, n + 1), 4):
            for (p, q), (r, s) in (((a, c), (b, d)), ((a, d), (b, c))):
                y = self.new()
                self.Y[((p, q), (r, s))] = y
                self.cross_def(y, p, q, r, s, ydef)
                for k in range(t):
                    self.cl.append([-y, -self.C[((p, q), k)], -self.C[((r, s), k)]])
        if tri:
            self.maximality()
        if col and not tri:
            self.precedence()
        self.euler_mode = int(euler)
        if euler:
            self.euler()
        if fix:
            for lit in fix:
                self.cl.append([lit])

    def new(self):
        self.nv += 1
        return self.nv

    def o(self, p, q, r):
        """literal for orient(p,q,r) > 0 (distinct p,q,r in any order)."""
        s = sorted((p, q, r))
        perm = [s.index(v) for v in (p, q, r)]
        inv = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        v = self.X[tuple(s)]
        return v if inv % 2 == 0 else -v

    def signotope(self):
        X = self.X
        for a, b, c, d in itertools.combinations(range(1, self.n + 1), 4):
            abc, abd, acd, bcd = X[(a, b, c)], X[(a, b, d)], X[(a, c, d)], X[(b, c, d)]
            self.cl += [[-abc, -acd, abd], [abc, acd, -abd],
                        [-abc, -bcd, acd], [abc, bcd, -acd]]

    def cross_def(self, y, p, q, r, s, ydef):
        o1, o2, o3, o4 = self.o(p, q, r), self.o(p, q, s), self.o(r, s, p), self.o(r, s, q)
        # (o1 xor o2) and (o3 xor o4) -> y
        for s1 in (1, -1):
            for s3 in (1, -1):
                self.cl.append([-s1 * o1, s1 * o2, -s3 * o3, s3 * o4, y])
        if ydef == "full":  # y -> o1 xor o2 ; y -> o3 xor o4
            self.cl += [[-y, o1, o2], [-y, -o1, -o2], [-y, o3, o4], [-y, -o3, -o4]]

    def maximality(self):
        """tri mode: C(e,k) means 'e in triangulation T_k' (each layer extended to a triangulation,
        so an edge may lie in several T_k). Maximality: e not in T_k -> some f in T_k crosses e."""
        inv = {e: [] for e in self.E}
        for (e1, e2), y in self.Y.items():
            inv[e1].append((e2, y)); inv[e2].append((e1, y))
        for e in self.E:
            for k in range(self.t):
                cl = [self.C[(e, k)]]
                for f, y in inv[e]:
                    z = self.new()
                    self.cl += [[-z, y], [-z, self.C[(f, k)]]]
                    cl.append(z)
                self.cl.append(cl)

    def lex_reflection(self):
        n, X = self.n, self.X
        m, f = (n + 1) // 2, n // 2
        A = [X[(a, a + 1, a + 2)] for a in range(m - 1, 1, -1)]
        B = [X[(a, a + 1, a + 2)] for a in range(f + 1, n - 1)]
        assert len(A) == len(B)
        eq_prev = None
        for i, (a, b) in enumerate(zip(A, B)):
            pre = [] if eq_prev is None else [-eq_prev]
            self.cl.append(pre + [-a, b])
            if i == len(A) - 1:
                break
            e = self.new()
            self.cl.append(pre + [a, b, e])
            self.cl.append(pre + [-a, -b, e])
            eq_prev = e

    def edge_order(self):
        # edges at p1 first (p1 is a hull vertex; its star is the natural anchor), then the rest.
        return sorted(self.E, key=lambda e: (e[0] != 1, e))

    def precedence(self):
        """Value precedence: colour k+1 may be used on edge e_i only if colour k is used on an
        earlier edge. Any colouring can be renamed (first-occurrence order) to satisfy it."""
        order = self.edge_order()
        t, C = self.t, self.C
        # U[i][k]: colour k used on some edge among order[0..i]
        U = []
        for i, e in enumerate(order):
            row = []
            for k in range(t):
                u = self.new()
                row.append(u)
                self.cl.append([-C[(e, k)], u])
                if i:
                    self.cl.append([-U[i - 1][k], u])
                    self.cl.append([-u, U[i - 1][k], C[(e, k)]])
                else:
                    self.cl.append([-u, C[(e, k)]])
            U.append(row)
            for k in range(1, t):
                if i == 0:
                    self.cl.append([-C[(e, k)]])
                else:
                    self.cl.append([-C[(e, k)], U[i - 1][k - 1]])

    def euler(self):
        n, X = self.n, self.X
        H = {}
        self.Zin = {}
        for v in range(1, n + 1):
            # v interior <=> exists triangle (a,b,c) containing v; H_v = not interior
            inner = []
            others = [w for w in range(1, n + 1) if w != v]
            for a, b, c in itertools.combinations(others, 3):
                # v inside abc <=> o(a,b,v)=o(b,c,v)=o(c,a,v)=o(a,b,c)
                z = self.new()
                lits = [self.o(a, b, v), self.o(b, c, v), self.o(c, a, v)]
                ref = self.o(a, b, c)
                # z <-> (all lits equal ref).  z -> lit==ref ; (all equal) -> z
                for l in lits:
                    self.cl.append([-z, -l, ref])
                    self.cl.append([-z, l, -ref])
                self.cl.append([z] + [-l for l in lits] + [-ref])
                self.cl.append([z] + [l for l in lits] + [ref])
                inner.append(z)
                self.Zin[(v, (a, b, c))] = z
            h = self.new()
            H[v] = h
            for z in inner:
                self.cl.append([-z, -h])
            self.cl.append([h] + inner)
        self.H = H
        bound = 3 * n - 3
        if self.euler_mode == 3:
            # Summed Euler bound ("excess"): sum_k |T_k| <= t(3n-3) - t*h, and every edge is in >= 1
            # layer, so  sum_e (m_e - 1) + t*h <= t(3n-3) - C(n,2),  m_e = #layers containing e.
            # tri mode: m_e - 1 >= [m_e>=2] + [m_e>=3] + ... (literals forced by memberships).
            # col mode: an uncrossed edge can be added to every layer, so it counts t-1 times.
            excess = []
            if self.tri:
                for e in self.E:
                    for j in range(2, self.t + 1):
                        d = self.new()
                        for ks in itertools.combinations(range(self.t), j):
                            self.cl.append([d] + [-self.C[(e, k)] for k in ks])
                        excess.append(d)
            else:
                inv = {e: [] for e in self.E}
                for (e1, e2), y in self.Y.items():
                    inv[e1].append(y); inv[e2].append(y)
                for e in self.E:
                    u = self.new()
                    self.cl.append([u] + inv[e])
                    excess += [u] * (self.t - 1)
            lits = excess + [H[v] for v in range(1, n + 1)] * self.t
            self.at_most(lits, self.t * (3 * n - 3) - len(self.E))
            self.euler_mode = 1  # also add the per-layer bounds below
        if self.euler_mode == 2:
            # U_e: e is crossed by no edge (forced true when all its Y are false). Such an edge can
            # be added to every layer, so layer_k + {uncrossed e} is still plane:
            #   #(C(e,k) or U_e) + h <= 3n-3.
            inv = {e: [] for e in self.E}
            for (e1, e2), y in self.Y.items():
                inv[e1].append(y); inv[e2].append(y)
            U = {}
            for e in self.E:
                U[e] = self.new()
                self.cl.append([U[e]] + inv[e])
            self.U = U
        for k in range(self.t):
            if self.euler_mode == 2:
                lits = []
                for e in self.E:
                    w = self.new()
                    self.cl += [[-self.C[(e, k)], w], [-self.U[e], w]]
                    lits.append(w)
            else:
                lits = [self.C[(e, k)] for e in self.E]
            self.at_most(lits + [H[v] for v in range(1, n + 1)], bound)

    def sweep(self, B, C):
        """Linear-functional sweep normal to hull edge BC (redundant, valid for point sets).
        L(x,y): d(x) > d(y), d = distance from line BC (tilted so d(B) < d(C) < all others).
        - total order (antisymmetric by construction, transitive);
        - convexity: y inside triangle abc -> d(y) < max(d(a),d(b),d(c));
        - segments: pq crosses xy -> max(d(p),d(q)) > min(d(x),d(y)) and vice versa."""
        n = self.n
        Lv = {}
        for x, y in itertools.combinations(range(1, n + 1), 2):
            Lv[(x, y)] = self.new()
        L = lambda x, y: Lv[(x, y)] if x < y else -Lv[(y, x)]
        for x, y, z in itertools.permutations(range(1, n + 1), 3):
            self.cl.append([-L(x, y), -L(y, z), L(x, z)])
        for x in range(1, n + 1):
            if x not in (B, C):
                self.cl += [[L(x, B)], [L(x, C)]]
        self.cl.append([L(C, B)])
        for (v, (a, b, c)), z in self.Zin.items():
            self.cl.append([-z, L(a, v), L(b, v), L(c, v)])
        for ((p, q), (x, y)), yv in self.Y.items():
            self.cl.append([-yv, L(p, x), L(p, y), L(q, x), L(q, y)])
            self.cl.append([-yv, L(x, p), L(x, q), L(y, p), L(y, q)])
        return L

    def at_most(self, lits, K):
        """Sequential counter (Sinz), at most K of lits true."""
        m = len(lits)
        if K >= m:
            return
        S = [[self.new() for _ in range(K)] for _ in range(m)]
        for i, x in enumerate(lits):
            self.cl.append([-x, S[i][0]])
            if i:
                for j in range(K):
                    self.cl.append([-S[i - 1][j], S[i][j]])
                for j in range(1, K):
                    self.cl.append([-x, -S[i - 1][j - 1], S[i][j]])
                self.cl.append([-x, -S[i - 1][K - 1]])
            else:
                for j in range(1, K):
                    self.cl.append([-S[0][j]])


def write_dimacs(enc, path):
    with open(path, "w") as f:
        f.write(f"p cnf {enc.nv} {len(enc.cl)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in enc.cl))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("t", type=int)
    ap.add_argument("out")
    ap.add_argument("--nop1", action="store_true")
    ap.add_argument("--nolex", action="store_true")
    ap.add_argument("--nocol", action="store_true")
    ap.add_argument("--euler", type=int, default=0, help="0 off, 1 layer+hull, 2 also uncrossed edges")
    ap.add_argument("--ydef", default="full")
    ap.add_argument("--tri", action="store_true", help="layers = covering triangulations (maximal)")
    a = ap.parse_args()
    e = Enc(a.n, a.t, p1=not a.nop1, lex=not a.nolex, col=not a.nocol, euler=a.euler, ydef=a.ydef, tri=a.tri)
    write_dimacs(e, a.out)
    print(f"n={a.n} t={a.t} vars={e.nv} clauses={len(e.cl)}")

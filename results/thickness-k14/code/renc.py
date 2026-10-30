"""Referee's independent encoder: "K_n has a straight-line drawing with t crossing-free layers".

Written from the definitions only (the authors' gen.py was NOT read before writing this).

Model.  Points 1..n in general position, labelled by increasing x (rotate so x-coordinates
are distinct).  X(a,b,c), a<b<c, is true iff (a,b,c) is counterclockwise.  For an arbitrary
triple, O(p,q,r) = X(sorted) xor parity(permutation).

Constraints (each is satisfied by the assignment induced by ANY real drawing):
  S  signotope axioms: for a<b<c<d the sequence X(abc),X(abd),X(acd),X(bcd) has at most one
     sign change (Felsner-Weil; necessary for x-sorted point sets, since moving a point d
     along ... -- see REFEREE.md for the 4-point argument).
  K  crossing: for each 4-set and each of its 3 pairings {pq, rs}: Y <-> [O(pqr)!=O(pqs)] and
     [O(rsp)!=O(rsq)] (exact segment-crossing test in general position);  Y -> the two edges
     do not share a layer.
  L  every edge in >= 1 layer (and <= 1, harmless).
  N  (optional, --nf) normal form X(1,a,b)=true for all 1<a<b (projective normalisation
     around a hull vertex; argument in REFEREE.md).
  P  (optional, --prec) layer value-precedence over edges in lex order.
  E  (optional, --euler) per-layer Euler bound: for each layer k,  [full=True: H, U, In are
     full equivalences, also valid for real drawings; plus the summed bound and units H1,H2,Hn]
        #{e : e in layer k or e uncrossed} + #{hull vertices} <= 3n-3,
     using one-sided indicator variables (H_v forced true when v lies in no triangle,
     U_e forced true when e crosses no edge, W_ek >= C_ek or U_e).
"""
import itertools, sys, json, argparse


class Enc:
    def __init__(self, n, t, nf=True, prec=True, euler=True, card="sinz", full=True):
        self.n, self.t, self.full = n, t, full
        self.nv = 0
        self.cl = []
        self.X = {}
        for tr in itertools.combinations(range(1, n + 1), 3):
            self.X[tr] = self.new()
        self.edges = list(itertools.combinations(range(1, n + 1), 2))
        self.C = {(e, k): self.new() for e in self.edges for k in range(t)}
        self.build(nf, prec, euler, card)

    def new(self):
        self.nv += 1
        return self.nv

    def O(self, p, q, r):
        """literal for 'p,q,r counterclockwise'"""
        s = sorted((p, q, r))
        perm = [s.index(p), s.index(q), s.index(r)]
        inv = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        v = self.X[tuple(s)]
        return v if inv % 2 == 0 else -v

    def xor(self, a, b):
        z = self.new()
        self.cl += [[-z, a, b], [-z, -a, -b], [z, -a, b], [z, a, -b]]
        return z

    def build(self, nf, prec, euler, card):
        n, t, X, C, cl = self.n, self.t, self.X, self.C, self.cl
        # S: signotope axioms
        for a, b, c, d in itertools.combinations(range(1, n + 1), 4):
            seq = [X[(a, b, c)], X[(a, b, d)], X[(a, c, d)], X[(b, c, d)]]
            for i, j, k in itertools.combinations(range(4), 3):
                cl.append([-seq[i], seq[j], -seq[k]])
                cl.append([seq[i], -seq[j], seq[k]])
        # L
        for e in self.edges:
            cl.append([C[(e, k)] for k in range(t)])
            for k1, k2 in itertools.combinations(range(t), 2):
                cl.append([-C[(e, k1)], -C[(e, k2)]])
        # K
        self.Y = {}
        cross_of = {e: [] for e in self.edges}
        for quad in itertools.combinations(range(1, n + 1), 4):
            a, b, c, d = quad
            for (p, q), (r, s) in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
                s1 = self.xor(self.O(p, q, r), self.O(p, q, s))
                s2 = self.xor(self.O(r, s, p), self.O(r, s, q))
                y = self.new()
                cl += [[-y, s1], [-y, s2], [y, -s1, -s2]]
                self.Y[((p, q), (r, s))] = y
                cross_of[(p, q)].append(y)
                cross_of[(r, s)].append(y)
                for k in range(t):
                    cl.append([-y, -C[((p, q), k)], -C[((r, s), k)]])
        # N
        if nf:
            for a, b in itertools.combinations(range(2, n + 1), 2):
                cl.append([X[(1, a, b)]])
        # P: value precedence of layers along the lex edge order
        if prec:
            E = self.edges
            cl.append([C[(E[0], 0)]])
            for i, e in enumerate(E):
                for k in range(1, t):
                    # layer k used at e  =>  layer k-1 used at some earlier edge
                    cl.append([-C[(e, k)]] + [C[(f, k - 1)] for f in E[:i]])
        # E: per-layer Euler bound with uncrossed edges
        if euler:
            H = {}
            for v in range(1, n + 1):
                H[v] = self.new()
                ins = []
                others = [u for u in range(1, n + 1) if u != v]
                for a, b, c in itertools.combinations(others, 3):
                    z = self.new()  # z -> v strictly inside triangle abc
                    o = self.O(a, b, c)
                    for (p, q) in ((a, b), (b, c), (c, a)):
                        w = self.O(p, q, v)
                        cl += [[-z, -o, w], [-z, o, -w]]
                    if self.full:  # z <- pattern, and H_v -> not z  (H_v <-> v in no triangle)
                        ws = [self.O(p, q, v) for (p, q) in ((a, b), (b, c), (c, a))]
                        cl += [[z, -o] + [-w for w in ws], [z, o] + ws, [-H[v], -z]]
                    ins.append(z)
                cl.append([H[v]] + ins)
            U = {}
            for e in self.edges:
                U[e] = self.new()
                cl.append([U[e]] + cross_of[e])
                if self.full:  # U_e <-> e crosses no edge
                    cl += [[-U[e], -y] for y in cross_of[e]]
            self.H, self.U = H, U
            for k in range(t):
                lits = []
                for e in self.edges:
                    w = self.new()
                    cl += [[w, -C[(e, k)]], [w, -U[e]]]
                    lits.append(w)
                lits += [H[v] for v in range(1, n + 1)]
                self.atmost(lits, 3 * n - 3)
            # summed form of the per-layer bounds (each edge in exactly one layer):
            #   sum_k |L_k u U| = |E| + (t-1)|U|  =>  (t-1)|U| + t*h <= t(3n-3) - |E|
            lits = [U[e] for e in self.edges] * (t - 1) + [H[v] for v in range(1, n + 1)] * t
            self.atmost(lits, t * (3 * n - 3) - len(self.edges))
            if nf:
                # 1, 2, n are hull vertices: 1 by N; for 2 and n, "v inside abc" contradicts N or
                # the signotope axiom on {v,a,b,c} (see REFEREE.md); true for real drawings anyway
                for v in (1, 2, n):
                    cl.append([H[v]])

    def atmost(self, lits, K):
        """Sinz sequential counter: s[i][j] >= [at least j+1 of lits[0..i] true]; forbid K+1."""
        m, cl = len(lits), self.cl
        if K >= m:
            return
        s = [[self.new() for _ in range(K)] for _ in range(m)]
        for i, x in enumerate(lits):
            cl.append([-x, s[i][0]])
            if i > 0:
                for j in range(K):
                    cl.append([-s[i - 1][j], s[i][j]])
                for j in range(1, K):
                    cl.append([-x, -s[i - 1][j - 1], s[i][j]])
                cl.append([-x, -s[i - 1][K - 1]])

    def dimacs(self, extra=()):
        out = [f"p cnf {self.nv} {len(self.cl) + len(extra)}\n"]
        out += [" ".join(map(str, c)) + " 0\n" for c in self.cl]
        out += [" ".join(map(str, c)) + " 0\n" for c in extra]
        return "".join(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int); ap.add_argument("t", type=int); ap.add_argument("out")
    ap.add_argument("--no-nf", action="store_true"); ap.add_argument("--no-prec", action="store_true")
    ap.add_argument("--no-euler", action="store_true"); ap.add_argument("--card", default="sinz")
    a = ap.parse_args()
    e = Enc(a.n, a.t, nf=not a.no_nf, prec=not a.no_prec, euler=not a.no_euler, card=a.card)
    with open(a.out, "w", newline="\n") as f:
        f.write(e.dimacs())
    print(a.out, e.nv, len(e.cl))

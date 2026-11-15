"""Projective spaces PG(n-1, q) built from scratch (stdlib only).

points : normalised vectors (first nonzero coordinate = 1), lexicographic order
lines  : all 2-dimensional subspaces, each a sorted tuple of q+1 point indices
planes : all hyperplanes (n = 4: planes), each a sorted tuple of point indices

Lines are found by *pair covering*: take any pair of points not yet on a
found line, generate the q+1 points of the span, mark all pairs covered.
The construction is verified (see `verify`): every pair of distinct points
is on exactly one line, line/point counts match the Gaussian binomials, etc.
"""
import itertools
from ff import GF


class PG:
    def __init__(self, n, q, F=None):
        self.n, self.q = n, q
        self.F = F or GF(q)
        F = self.F
        pts = []
        for v in itertools.product(range(q), repeat=n):
            nz = next((c for c in v if c), None)
            if nz == 1:
                pts.append(v)
        self.points = pts
        self.pidx = {v: i for i, v in enumerate(pts)}
        self.N = len(pts)
        assert self.N == (q ** n - 1) // (q - 1)
        self._build_lines()
        if n == 4:
            self._build_planes()
        elif n == 3:
            self.planes = self.lines

    # ---------------------------------------------------------------- basics
    def normalize(self, v):
        F = self.F
        nz = next((c for c in v if c), 0)
        if nz == 0:
            raise ValueError("zero vector")
        if nz == 1:
            return tuple(v)
        iv = F.inv[nz]
        return tuple(F.mul[iv][c] for c in v)

    def lincomb(self, a, u, b, v):
        F = self.F
        return tuple(F.add[F.mul[a][x]][F.mul[b][y]] for x, y in zip(u, v))

    def span_line(self, i, j):
        """indices of the q+1 points of the line through points i and j (i != j)."""
        u, v = self.points[i], self.points[j]
        res = [self.pidx[self.normalize(self.lincomb(1, u, b, v))] for b in range(self.q)]
        res.append(j)          # the point v itself ( = 0*u + 1*v )
        return tuple(sorted(res))

    def _build_lines(self):
        N, q = self.N, self.q
        covered = bytearray(N * N)
        lines = []
        for i in range(N):
            for j in range(i + 1, N):
                if covered[i * N + j]:
                    continue
                L = self.span_line(i, j)
                assert len(set(L)) == q + 1
                for a in L:
                    for b in L:
                        if a != b:
                            covered[a * N + b] = 1
                lines.append(L)
        self.lines = lines
        self.lidx = {L: k for k, L in enumerate(lines)}
        thr = [[] for _ in range(N)]
        for k, L in enumerate(lines):
            for p in L:
                thr[p].append(k)
        self.lines_through = thr

    def _build_planes(self):
        q, F = self.q, self.F
        planes = []
        for a in self.points:   # dual coordinates, also normalised vectors
            pl = tuple(i for i, v in enumerate(self.points)
                       if self._dot(a, v) == 0)
            planes.append(pl)
        self.planes = planes
        self.plane_normals = list(self.points)

    def _dot(self, a, v):
        F = self.F
        s = 0
        for x, y in zip(a, v):
            s = F.add[s][F.mul[x][y]]
        return s

    # ---------------------------------------------------------- verification
    def verify(self):
        n, q, N = self.n, self.q, self.N
        gb = lambda m, k: self._gauss(m, k, q)
        nl = gb(n, 2)
        assert len(self.lines) == nl, (len(self.lines), nl)
        for L in self.lines:
            assert len(L) == q + 1 and len(set(L)) == q + 1
        # every pair of distinct points on exactly one line
        cnt = bytearray(N * N)
        for L in self.lines:
            for a in L:
                for b in L:
                    if a < b:
                        assert cnt[a * N + b] == 0
                        cnt[a * N + b] = 1
        assert sum(cnt) == N * (N - 1) // 2
        r = (q ** (n - 1) - 1) // (q - 1)
        for p in range(N):
            assert len(self.lines_through[p]) == r
        # every line is closed under taking spans of two of its points
        for L in self.lines[:: max(1, len(self.lines) // 200)]:
            assert self.span_line(L[0], L[1]) == L
            assert self.span_line(L[1], L[-1]) == L
        if n == 4:
            assert len(self.planes) == N
            for pl in self.planes:
                assert len(pl) == (q ** 3 - 1) // (q - 1)
            # every line is in 1 plane or q+1 planes-through: count planes containing a line
            pl_sets = [set(p) for p in self.planes]
            for L in self.lines[:: max(1, len(self.lines) // 100)]:
                c = sum(1 for s in pl_sets if all(x in s for x in L))
                assert c == q + 1
        return True

    @staticmethod
    def _gauss(n, k, q):
        num, den = 1, 1
        for i in range(k):
            num *= q ** (n - i) - 1
            den *= q ** (i + 1) - 1
        return num // den

    # ------------------------------------------------------------ group action
    def perm_of_matrix(self, M, frob=0):
        """Permutation of point indices induced by P -> normalize(M * frob^f(P))
        (M is a list of rows; column vectors).  frob = number of times the
        Frobenius x -> x^p is applied to the coordinates first."""
        F, n = self.F, self.n
        perm = []
        for v in self.points:
            w = v
            for _ in range(frob):
                w = tuple(F.frob[c] for c in w)
            out = []
            for row in M:
                s = 0
                for a, x in zip(row, w):
                    s = F.add[s][F.mul[a][x]]
                out.append(s)
            perm.append(self.pidx[self.normalize(out)])
        assert sorted(perm) == list(range(self.N))
        return perm


if __name__ == "__main__":
    import sys, time
    for (n, q) in [(3, 2), (3, 3), (3, 9), (4, 2), (4, 3), (4, 4), (4, 5), (4, 7), (4, 8), (4, 9)]:
        t = time.time()
        P = PG(n, q)
        P.verify()
        print("PG(%d,%d): %d points, %d lines, %d planes   verified (%.1fs)" %
              (n - 1, q, P.N, len(P.lines), len(P.planes), time.time() - t))

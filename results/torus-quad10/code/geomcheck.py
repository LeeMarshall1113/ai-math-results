"""Independent exact checker: is (vertices, quadrilateral faces) an embedded polyhedral torus in R^3
with planar strictly convex faces and no coplanar edge-adjacent faces?  Standard library only."""
from fractions import Fraction as Q
import itertools

def sub(u, v): return tuple(a - b for a, b in zip(u, v))
def dot(u, v): return sum(a * b for a, b in zip(u, v))
def cross(u, v): return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])

def check_torus(verts, faces):
    X = [tuple(Q(c) for c in p) for p in verts]
    F = [tuple(f) for f in faces]
    rep = {"valid": False}
    V = len(X)
    if len(set(X)) != V: rep["err"] = "repeated vertex coordinates"; return rep
    if any(not (0 <= u < V) for f in F for u in f): rep["err"] = "face label out of range"; return rep
    # --- combinatorics
    edge_faces = {}
    for i, f in enumerate(F):
        if len(f) != 4 or len(set(f)) != 4: rep["err"] = f"face {i} not a 4-cycle"; return rep
        for k in range(4):
            edge_faces.setdefault(frozenset((f[k], f[(k+1) % 4])), []).append(i)
    if any(len(v) != 2 for v in edge_faces.values()): rep["err"] = "edge not in exactly 2 faces"; return rep
    for v in range(V):
        adj = {}
        for f in F:
            if v in f:
                k = f.index(v); a, c = f[k-1], f[(k+1) % 4]
                adj.setdefault(a, []).append(c); adj.setdefault(c, []).append(a)
        if not adj or any(len(t) != 2 for t in adj.values()): rep["err"] = f"bad link at {v}"; return rep
        s = next(iter(adj)); prev, cur, n = None, s, 0
        while True:
            nx = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]; prev, cur = cur, nx; n += 1
            if cur == s: break
        if n != len(adj): rep["err"] = f"link at {v} not one cycle"; return rep
    chi_e = V - len(edge_faces) + len(F)
    # orientability
    orient = {0: F[0]}; stack = [0]
    while stack:
        i = stack.pop(); f = orient[i]
        D = {(f[k], f[(k+1) % 4]) for k in range(4)}
        for j, g in enumerate(F):
            if j in orient: continue
            Dg = {(g[k], g[(k+1) % 4]) for k in range(4)}
            if any((b, a) in Dg for a, b in D): orient[j] = g; stack.append(j)
            elif any((a, b) in Dg for a, b in D): orient[j] = g[::-1]; stack.append(j)
    if len(orient) != len(F): rep["err"] = "disconnected"; return rep
    darts = [(o[k], o[(k+1) % 4]) for o in orient.values() for k in range(4)]
    if len(set(darts)) != len(darts): rep["err"] = "non-orientable"; return rep
    if chi_e != 0: rep["err"] = f"Euler characteristic {chi_e}"; return rep
    # --- geometry of faces
    planes = []
    for i, f in enumerate(F):
        P = [X[u] for u in f]
        n = cross(sub(P[1], P[0]), sub(P[2], P[0]))
        if n == (0, 0, 0) or dot(n, sub(P[3], P[0])) != 0: rep["err"] = f"face {i} not planar"; return rep
        turns = [dot(n, cross(sub(P[(k+1) % 4], P[k]), sub(P[(k+2) % 4], P[(k+1) % 4]))) for k in range(4)]
        if not (all(t > 0 for t in turns) or all(t < 0 for t in turns)):
            rep["err"] = f"face {i} not strictly convex"; return rep
        planes.append((n, dot(n, P[0])))
    for e, (i, j) in edge_faces.items():
        if cross(planes[i][0], planes[j][0]) == (0, 0, 0): rep["err"] = f"faces {i},{j} coplanar"; return rep
    # --- pairwise intersections
    counts = {"edge": 0, "vertex": 0, "empty": 0}
    for i, j in itertools.combinations(range(len(F)), 2):
        common = set(F[i]) & set(F[j])
        (n1, c1), (n2, c2) = planes[i], planes[j]
        d = cross(n1, n2)
        if d == (0, 0, 0):
            if c1 * sum(n2) != c2 * sum(n1) and not _parallel_distinct(n1, c1, n2, c2):
                pass
            if not _parallel_distinct(n1, c1, n2, c2):
                rep["err"] = f"faces {i},{j} coplanar (not handled)"; return rep
            if common: rep["err"] = f"faces {i},{j} parallel but share vertices"; return rep
            counts["empty"] += 1; continue
        def interval(f, n, c):
            ts = []
            P = [X[u] for u in f]
            g = [dot(n, p) - c for p in P]
            for k in range(4):
                a, b = P[k], P[(k+1) % 4]; ga, gb = g[k], g[(k+1) % 4]
                if ga == 0: ts.append(dot(d, a))
                if ga * gb < 0:
                    lam = ga / (ga - gb); x = tuple(a[m] + lam * (b[m] - a[m]) for m in range(3))
                    ts.append(dot(d, x))
            return (min(ts), max(ts)) if ts else None
        I1, I2 = interval(F[i], n2, c2), interval(F[j], n1, c1)
        ov = None
        if I1 and I2:
            lo, hi = max(I1[0], I2[0]), min(I1[1], I2[1])
            if lo <= hi: ov = (lo, hi)
        if len(common) == 2:
            a, b = common
            for f in (F[i], F[j]):   # referee fix: the shared pair must be an edge of BOTH faces
                k = f.index(a)
                if b not in (f[(k + 1) % 4], f[(k - 1) % 4]):
                    rep["err"] = f"faces {i},{j} share a diagonal pair {a},{b}"; return rep
            ta, tb = dot(d, X[a]), dot(d, X[b])
            exp = (min(ta, tb), max(ta, tb)); kind = "edge"
        elif len(common) == 1:
            v = next(iter(common)); t = dot(d, X[v]); exp = (t, t); kind = "vertex"
        elif not common:
            exp = None; kind = "empty"
        else:
            rep["err"] = f"faces {i},{j} share {len(common)} vertices"; return rep
        if ov != exp: rep["err"] = f"faces {i},{j}: intersection {ov} != expected {kind} {exp}"; return rep
        counts[kind] += 1
    rep.update(valid=True, V=V, E=len(edge_faces), F=len(F), pairs=counts)
    return rep

def _parallel_distinct(n1, c1, n2, c2):
    # n2 = k n1; planes distinct iff c2 != k c1
    for m in range(3):
        if n1[m] != 0:
            k = n2[m] / n1[m]; return c2 != k * c1
    return True

if __name__ == "__main__":
    import json, sys
    W = json.load(open(sys.argv[1]))
    print(check_torus(W["vertices"], W["faces"]))

"""Second, independently written encoding of the core conditions (GP + PL + AD) for {4,4}_(3,1).
Differences from encode.py: variables z_T ("chi(T) != 0") and s_T ("chi(T) > 0 if nonzero");
GP encoded by brute-force enumeration of all 3^6 sign assignments of the six brackets of each
relation (no auxiliary product variables); AD uses a different spanning triple of F; faces are
generated from lattice coordinates rather than the closed formula.
usage: encode2.py out.cnf [--drop=0,3]
--drop=i,j,...: omit the PL and AD constraints of faces Q_i (and every AD pair involving them);
GP is kept on all 10 points.  --drop=0,3 gives the 8-face core (Q0, Q3 share the edge {0,1})."""
import itertools, sys
DROP = ()
if len(sys.argv) > 2 and sys.argv[2].startswith("--drop="):
    DROP = tuple(int(i) for i in sys.argv.pop(2)[len("--drop="):].split(","))

# faces from Z^2 / <(10,0),(3,1)> : (x,y) -> (x - 3y) mod 10
lab = lambda x, y: (x - 3 * y) % 10
F = [(lab(x, 0), lab(x + 1, 0), lab(x + 1, 1), lab(x, 1)) for x in range(10)]
DROPPED = {F[i] for i in DROP}
NPT = 10
W = None
if len(sys.argv) > 2:   # soundness test mode: faces + coordinates of a known torus from JSON
    import json
    W = json.load(open(sys.argv[2])); F = [tuple(f) for f in W["faces"]]; NPT = len(W["vertices"])
quads = sorted(itertools.combinations(range(NPT), 4))
Z = {T: 1 + k for k, T in enumerate(quads)}            # nonzero
S = {T: 1 + len(quads) + k for k, T in enumerate(quads)}  # sign (meaningful only if nonzero)
NV = 2 * len(quads)
clauses = []

def parity(T):
    return sum(T[i] > T[j] for i in range(4) for j in range(i + 1, 4)) & 1

def cond(T, v):
    """list of literals whose conjunction says bracket(T) == v, v in {-1,0,1}"""
    key = tuple(sorted(T))
    if v == 0: return [-Z[key]]
    sgn = v if parity(T) == 0 else -v
    return [Z[key], S[key] if sgn > 0 else -S[key]]

def forbid(conj):
    clauses.append(sorted({-l for l in conj}))

# GP: for distinct a,b and 4-set {y1<y2<y3<y4} disjoint:
#   [ab y1 y2][ab y3 y4] - [ab y1 y3][ab y2 y4] + [ab y1 y4][ab y2 y3] = 0
for a, b in itertools.combinations(range(NPT), 2):
    rest = [u for u in range(NPT) if u not in (a, b)]
    for y1, y2, y3, y4 in itertools.combinations(rest, 4):
        Ts = [(a, b, y1, y2), (a, b, y3, y4), (a, b, y1, y3), (a, b, y2, y4), (a, b, y1, y4), (a, b, y2, y3)]
        for vals in itertools.product((-1, 0, 1), repeat=6):
            t = [vals[0] * vals[1], -vals[2] * vals[3], vals[4] * vals[5]]
            if (max(t) > 0 and min(t) >= 0) or (min(t) < 0 and max(t) <= 0):
                conj = []
                for T, v in zip(Ts, vals): conj += cond(T, v)
                forbid(conj)
# PL: every face is coplanar
for f in F:
    if f in DROPPED: continue
    forbid(cond(f, 1)); forbid(cond(f, -1))
# AD: for faces f, g sharing an edge {p,q}: the two vertices of g off the edge lie strictly on the
# same side of plane(f).  Spanning triple of f: q, p and the vertex of f adjacent to q (not p).
for f in F:
    for g in F:
        if f == g or f in DROPPED or g in DROPPED: continue
        sh = set(f) & set(g)
        if len(sh) != 2: continue
        k = [i for i in range(4) if f[i] in sh and f[(i + 1) % 4] in sh][0]
        p, q = f[k], f[(k + 1) % 4]
        r = f[(k + 2) % 4]
        x, y = [u for u in g if u not in sh]
        for vx in (-1, 0, 1):
            for vy in (-1, 0, 1):
                if vx * vy <= 0:
                    forbid(cond((q, p, r, x), vx) + cond((q, p, r, y), vy))

with open(sys.argv[1], "w", newline="\n") as fh:
    fh.write(f"p cnf {NV} {len(clauses)}\n")
    for c in clauses: fh.write(" ".join(map(str, c)) + " 0\n")
if W is not None:
    from fractions import Fraction as Fr
    X = [[Fr(c) for c in p] for p in W["vertices"]]
    def br(T):
        p = [X[t] for t in T]; u, v, w = [[p[k][m] - p[0][m] for m in range(3)] for k in (1, 2, 3)]
        d = u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0])
        return (d > 0) - (d < 0)
    val = {}
    for T in quads:
        b = br(T); val[Z[T]] = b != 0; val[S[T]] = b > 0
    bad = sum(1 for c in clauses if not any(val[abs(l)] == (l > 0) for l in c))
    print("soundness: clauses violated by the true chirotope:", bad)
print(sys.argv[1], "vars", NV, "clauses", len(clauses), "faces", [f for f in F if f not in DROPPED])

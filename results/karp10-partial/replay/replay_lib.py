"""replay_lib.py -- self-contained (standard library only) checker for the Karp n=10 exclusion certificates.

Nothing here uses nauty, GAP, z3 or the tables of the search pipeline: groups are closed by brute force,
orbits are union-find, graph (non-)embeddings are plain backtracking, isomorphism = embedding with equal edge
count.  Everything a third party needs to replay a certificate is in this file + replay_cert.py.

Conventions: vertices 0..9; edge (i,j), i<j, index = position in lexicographic order; a graph is a 45-bit
integer mask; a permutation g is a tuple with g[i] = image of i.
"""
import itertools, math, sys

# The checks below are written as `assert` statements.  Under `python -O` (or PYTHONOPTIMIZE) they would be
# silently skipped and every certificate would "pass", so refuse to run in that mode.
if not __debug__:
    raise SystemExit("replay_lib: assertions are disabled (python -O / PYTHONOPTIMIZE); the replay would be vacuous")


def is_prime(p):
    return isinstance(p, int) and p >= 2 and all(p % d for d in range(2, int(p ** 0.5) + 1))


N = 10
EDGES = [(i, j) for i in range(N) for j in range(i + 1, N)]
M = len(EDGES)
EIDX = {e: k for k, e in enumerate(EDGES)}
FULL = (1 << M) - 1


def eid(i, j):
    return EIDX[(i, j) if i < j else (j, i)]


def mask_of(edges):
    m = 0
    for (i, j) in edges:
        m |= 1 << eid(i, j)
    return m


def edges_of(mask):
    return [EDGES[k] for k in range(M) if (mask >> k) & 1]


def popcount(x):
    return bin(x).count("1")


def complement(mask):
    return FULL ^ mask


def rows_of(mask):
    rows = [0] * N
    for (i, j) in edges_of(mask):
        rows[i] |= 1 << j
        rows[j] |= 1 << i
    return rows


# ---------------------------------------------------------------- permutation groups (brute force)
def compose(a, b):
    """(a*b)(i) = a(b(i))"""
    return tuple(a[b[i]] for i in range(N))


def inverse(a):
    r = [0] * N
    for i, x in enumerate(a):
        r[x] = i
    return tuple(r)


IDENT = tuple(range(N))


def closure(gens):
    gens = [tuple(g) for g in gens]
    for g in gens:
        assert sorted(g) == list(range(N)), "not a permutation of 0..9"
    seen = {IDENT}
    frontier = [IDENT]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = compose(g, x)
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
        if len(seen) > 400000:
            raise RuntimeError("group too large for brute-force closure")
    return seen


def edge_orbits(gens):
    parent = list(range(M))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in gens:
        for k, (i, j) in enumerate(EDGES):
            a, b = find(k), find(eid(g[i], g[j]))
            if a != b:
                parent[a] = b
    orbs = {}
    for k in range(M):
        orbs[find(k)] = orbs.get(find(k), 0) | (1 << k)
    return sorted(orbs.values(), key=lambda m: m & -m)


def check_G1(gens, order, p, ngens, g):
    """Oliver class G^1 witness: Gamma = <gens> has the stated order; N = <ngens> is a normal p-subgroup;
    Gamma = <g> N (so Gamma/N is cyclic)."""
    assert is_prime(p), f"p = {p} is not a prime"
    G = closure(gens)
    assert len(G) == order, f"order mismatch {len(G)} != {order}"
    Nn = closure(ngens)
    assert Nn <= G, "N not inside Gamma"
    n = len(Nn)
    # |N| must be a power of the prime p
    t = n
    while t % p == 0:
        t //= p
    assert t == 1, f"|N|={n} is not a power of {p}"
    for s in gens:
        si = inverse(s)
        for x in Nn:
            assert compose(compose(s, x), si) in Nn, "N not normal"
    H = closure(list(ngens) + [g])
    assert H == G, "<g>N != Gamma"
    return len(G), n


# ---------------------------------------------------------------- embeddings
def find_embedding(gmask, hmask):
    """permutation phi of 0..9 with phi(E(G)) subset of E(H), or None."""
    if popcount(gmask) > popcount(hmask):
        return None
    grow, hrow = rows_of(gmask), rows_of(hmask)
    gdeg = [popcount(r) for r in grow]
    hdeg = [popcount(r) for r in hrow]
    if any(a > b for a, b in zip(sorted(gdeg, reverse=True), sorted(hdeg, reverse=True))):
        return None
    verts = [v for v in range(N) if gdeg[v] > 0]
    order, chosen = [], 0
    while len(order) < len(verts):
        best = max((v for v in verts if not (chosen >> v) & 1), key=lambda v: (popcount(grow[v] & chosen), gdeg[v]))
        order.append(best)
        chosen |= 1 << best
    phi = {}

    def rec(idx, used):
        if idx == len(order):
            return True
        u = order[idx]
        cand = (1 << N) - 1
        for w, img in phi.items():
            if (grow[u] >> w) & 1:
                cand &= hrow[img]
        cand &= ~used
        for v in range(N):
            if (cand >> v) & 1 and hdeg[v] >= gdeg[u]:
                phi[u] = v
                if rec(idx + 1, used | (1 << v)):
                    return True
                del phi[u]
        return False

    if not rec(0, 0):
        return None
    full = [None] * N
    for u, v in phi.items():
        full[u] = v
    rest = iter([v for v in range(N) if v not in phi.values()])
    for u in range(N):
        if full[u] is None:
            full[u] = next(rest)
    return tuple(full)


def check_embedding(phi, gmask, hmask):
    if sorted(phi) != list(range(N)):
        return False
    return all((hmask >> eid(phi[i], phi[j])) & 1 for (i, j) in edges_of(gmask))


def embeds(g, h):
    return find_embedding(g, h) is not None


def isomorphic(a, b):
    return popcount(a) == popcount(b) and embeds(a, b)


# ---------------------------------------------------------------- the 12 special graphs (Angel-Borja)
def cycle_edges(vs):
    return [(vs[t], vs[(t + 1) % len(vs)]) for t in range(len(vs))]


def special_graphs():
    U, P = list(range(5)), list(range(5, 10))
    G = {}
    G["G1"] = mask_of(cycle_edges(U))
    G["G2"] = mask_of(list(itertools.combinations(U, 2)))
    G["G3"] = mask_of(cycle_edges(U) + cycle_edges(P))
    G["G4"] = mask_of(cycle_edges(P) + list(itertools.combinations(U, 2)))
    G["G5"] = mask_of(list(itertools.combinations(U, 2)) + list(itertools.combinations(P, 2)))
    for k in ["G1", "G2", "G3", "G4", "G5"]:
        G[k + "bar"] = complement(G[k])
    return G


SPEC10 = ["G1", "G2", "G3", "G4", "G5", "G5bar", "G4bar", "G3bar", "G2bar", "G1bar"]

# Angel-Borja's nine order ideals as printed in the paper (Section 5)
AB_IDEALS = {
    1: {"G5bar"},
    2: {"G1", "G4bar", "G5bar"},
    3: {"G1", "G3", "G5bar"},
    4: {"G1", "G2", "G2bar", "G4bar", "G5bar"},
    5: {"G1", "G2", "G3", "G4", "G5"},
    6: {"G1", "G3", "G3bar", "G4bar", "G5bar"},
    7: {"G1", "G2", "G3", "G4", "G2bar", "G4bar", "G5bar"},
    8: {"G1", "G2", "G3", "G2bar", "G3bar", "G4bar", "G5bar"},
    9: {"G1", "G2", "G3", "G4", "G1bar", "G2bar", "G3bar", "G4bar", "G5bar"},
}


def derive_types():
    """Recompute the types: Sylow-5 subgroup S=<(01234),(56789)> of S10, chi(Delta^S)=1."""
    G = special_graphs()
    S = [(1, 2, 3, 4, 0, 5, 6, 7, 8, 9), (0, 1, 2, 3, 4, 6, 7, 8, 9, 5)]
    orbs = edge_orbits(S)
    assert len(orbs) == 5 and sorted(popcount(o) for o in orbs) == [5, 5, 5, 5, 25]
    weight = {n: 0 for n in SPEC10}
    allowed = dict(G)
    allowed["empty"] = 0
    allowed["K10"] = FULL
    unions = []
    for r in range(1, 6):
        for T in itertools.combinations(range(5), r):
            m = 0
            for t in T:
                m |= orbs[t]
            unions.append((r, m))
    for r, m in unions:
        hit = [n for n in allowed if isomorphic(allowed[n], m)]
        assert len(hit) == 1, f"union is not one of the 12 special classes: {hit}"
        if hit[0] in weight:
            weight[hit[0]] += (-1) ** (r - 1)
    below = {a: {b for b in SPEC10 if b != a and embeds(G[b], G[a])} for a in SPEC10}
    types = []
    for sub in range(1 << 10):
        I = {SPEC10[t] for t in range(10) if (sub >> t) & 1}
        if all(below[a] <= I for a in I) and sum(weight[a] for a in I) == 1:
            types.append(frozenset(I))
    return weight, below, types


# ---------------------------------------------------------------- statuses of a graph for a type
def status(mask, members, G):
    """('IN', S, phi) if mask embeds in a member S; ('OUT', S, phi) if a non-member S (or K10) embeds in mask; else ('FREE',)."""
    for S in members:
        phi = find_embedding(mask, G[S])
        if phi is not None:
            assert check_embedding(phi, mask, G[S])
            return ("IN", S, phi)
    for S in [s for s in SPEC10 if s not in members] + ["K10"]:
        sm = FULL if S == "K10" else G[S]
        phi = find_embedding(sm, mask)
        if phi is not None:
            assert check_embedding(phi, sm, mask)
            return ("OUT", S, phi)
    return ("FREE",)

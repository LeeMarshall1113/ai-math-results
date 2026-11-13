"""indep_graph.py -- a SECOND, independent implementation of graph isomorphism and embedding for the Karp n=10 replay.

Nothing here imports replay_lib (the cross-test script imports both and compares).  Graphs on 10 vertices are given as adjacency lists of bit rows
(row[v] = bitmask of neighbours), built from an edge list, so even the edge-mask encoding of replay_lib is not used.

Different algorithms from replay_lib:
  * isomorphism: CANONICAL FORM (colour refinement + individualisation search, minimum adjacency string over the leaves, with pruning by automorphisms
    that were found and verified); two graphs are isomorphic iff their canonical forms are equal.   replay_lib: backtracking embedding + equal edge count.
  * embedding G -> H (an injective vertex map sending edges to edges, i.e. G is a subgraph of H, spanning 10 vertices): PACKING formulation, G placed
    against the complement of H in natural vertex order with forward checking on the remaining candidate sets.   replay_lib: degree-ordered backtracking.
  * a third check with z3 is in crosscheck_independent.py (sampled).
"""
import itertools

N = 10


def rows_from_edges(edges):
    rows = [0] * N
    for a, b in edges:
        assert a != b
        rows[a] |= 1 << b
        rows[b] |= 1 << a
    return rows


def edges_from_rows(rows):
    return [(a, b) for a in range(N) for b in range(a + 1, N) if (rows[a] >> b) & 1]


def complement_rows(rows):
    full = (1 << N) - 1
    return [(~rows[v]) & full & ~(1 << v) for v in range(N)]


def n_edges(rows):
    return sum(bin(r).count("1") for r in rows) // 2


def is_connected(rows):
    seen = 1
    stack = [0]
    while stack:
        v = stack.pop()
        nb = rows[v] & ~seen
        seen |= nb
        for u in range(N):
            if (nb >> u) & 1:
                stack.append(u)
    return seen == (1 << N) - 1


# ------------------------------------------------------------------ canonical form
def _refine(rows, col):
    """colour refinement to an equitable partition; colours are ranks 0..k-1 assigned by sorting invariant signatures"""
    while True:
        sig = [(col[v], tuple(sorted(col[u] for u in range(N) if (rows[v] >> u) & 1))) for v in range(N)]
        order = sorted(set(sig))
        new = [order.index(s) for s in sig]
        if len(set(new)) == len(set(col)):
            return new
        col = new


def canonical(rows):
    """canonical form: the minimum, over the leaves of the individualisation tree, of the adjacency matrix read in colour order, as a tuple of 10 ints"""
    best = [None, None]            # [string, permutation (vertex -> label)]
    autos = []                     # verified automorphisms (as lists)

    def string_of(lab):
        inv = [0] * N
        for v in range(N):
            inv[lab[v]] = v
        return tuple(sum(((rows[inv[i]] >> inv[j]) & 1) << j for j in range(N)) for i in range(N))

    def is_auto(p):
        return all(((rows[a] >> b) & 1) == ((rows[p[a]] >> p[b]) & 1) for a in range(N) for b in range(N))

    def search(col, path):
        col = _refine(rows, col)
        k = len(set(col))
        if k == N:
            s = string_of(col)
            if best[0] is None or s < best[0]:
                best[0], best[1] = s, col[:]
            elif s == best[0]:
                # col and best[1] label the same graph identically: sigma = inv(best) o col is an automorphism
                inv = [0] * N
                for v in range(N):
                    inv[best[1][v]] = v
                sigma = [inv[col[v]] for v in range(N)]
                assert is_auto(sigma)
                autos.append(sigma)
            return
        # first smallest non-singleton cell (smallest size, then smallest colour): invariant choice
        sizes = {}
        for c in col:
            sizes[c] = sizes.get(c, 0) + 1
        cell_col = min((c for c in sizes if sizes[c] > 1), key=lambda c: (sizes[c], c))
        cell = [v for v in range(N) if col[v] == cell_col]
        tried = []
        for w in cell:
            # prune w if an automorphism fixing the whole path pointwise maps an already tried vertex to w
            if tried:
                fix = [a for a in autos if all(a[x] == x for x in path)]
                orb = {t: t for t in tried}
                changed = True
                reach = set(tried)
                while changed:
                    changed = False
                    for a in fix:
                        for x in list(reach):
                            if a[x] not in reach:
                                reach.add(a[x])
                                changed = True
                if w in reach:
                    continue
            tried.append(w)
            ncol = [2 * col[v] + (0 if v == w else 1) for v in range(N)]
            search(ncol, path + [w])

    search([0] * N, [])
    return best[0]


def isomorphic(r1, r2):
    return n_edges(r1) == n_edges(r2) and canonical(r1) == canonical(r2)


# ------------------------------------------------------------------ embedding (packing formulation)
def _order(rg):
    """static placement order of the source vertices: breadth-first from a vertex of maximum degree (components in turn), so that each vertex
    is placed next to already placed neighbours whenever possible"""
    deg = [bin(r).count("1") for r in rg]
    order, seen = [], 0
    while len(order) < N:
        start = max((v for v in range(N) if not (seen >> v) & 1), key=lambda v: (deg[v], -v))
        queue = [start]
        seen |= 1 << start
        while queue:
            v = queue.pop(0)
            order.append(v)
            for u in sorted((u for u in range(N) if (rg[v] >> u) & 1 and not (seen >> u) & 1), key=lambda u: (-deg[u], u)):
                seen |= 1 << u
                queue.append(u)
    return order


def embedding(rg, rh):
    """packing search: returns a witness list phi (phi[v] = image of v) with phi(edges of G) inside edges of H, or None"""
    if n_edges(rg) > n_edges(rh):
        return None
    order = _order(rg)
    pos = {v: i for i, v in enumerate(order)}
    dg = [bin(r).count("1") for r in rg]
    dh = [bin(r).count("1") for r in rh]
    phi = [None] * N

    def rec(i, used):
        if i == N:
            return True
        v = order[i]
        cand = ((1 << N) - 1) & ~used
        for u in range(N):
            if phi[u] is not None and (rg[v] >> u) & 1:
                cand &= rh[phi[u]]
        for t in range(N):
            if (cand >> t) & 1 and dh[t] >= dg[v]:
                phi[v] = t
                # forward check: every unplaced vertex keeps at least one candidate
                nu = used | (1 << t)
                ok = True
                for j in range(i + 1, N):
                    w = order[j]
                    c = ((1 << N) - 1) & ~nu
                    for u in range(N):
                        if phi[u] is not None and (rg[w] >> u) & 1:
                            c &= rh[phi[u]]
                    if c == 0:
                        ok = False
                        break
                if ok and rec(i + 1, nu):
                    return True
                phi[v] = None
        return False

    return list(phi) if rec(0, 0) else None


def embeds(rg, rh):
    return embedding(rg, rh) is not None


def check_witness(phi, rg, rh):
    return sorted(phi) == list(range(N)) and all(((rh[phi[a]] >> phi[b]) & 1) for a in range(N) for b in range(N) if (rg[a] >> b) & 1)

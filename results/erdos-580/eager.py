"""Eager avoidance-clause generation (all placements of a rooted core), used to
(a) diff against Zeraoulia's published CNFs and (b) build independent CNFs for his 4 cores."""
import sys
from lks import edge_vars, ev, bfs_order


def placements(cv, ce, marked, m):
    """Yield the edge-set clause of every injective phi: V(F)->V(H) with marked -> L and
    every core edge landing on an L-L or L-S pair."""
    var = edge_vars(m)
    order = bfs_order(cv, ce)
    cadj = {u: set() for u in cv}
    for a, b in ce:
        cadj[a].add(b); cadj[b].add(a)
    mk = set(marked)
    phi = {}
    used = set()
    out = set()

    def rec(i):
        if i == len(order):
            out.add(frozenset(-ev(var, phi[a], phi[b]) for a, b in ce))
            return
        u = order[i]
        for h in range(2 * m):
            if h in used or (u in mk and h >= m):
                continue
            if any(w in phi and ev(var, phi[w], h) is None for w in cadj[u]):
                continue  # would need an S-S edge
            phi[u] = h; used.add(h)
            rec(i + 1)
            del phi[u]; used.discard(h)
    rec(0)
    return out


# Zeraoulia's four rooted cores, as (core vertices, core edges, marked), on m = 9
P5 = ([0, 1, 2, 3, 4], [(0, 1), (1, 2), (2, 3), (3, 4)])
ZCORES = {
    "path_support_014": (P5[0], P5[1], [0, 1, 4]),
    "path_support_0124": (P5[0], P5[1], [0, 1, 2, 4]),
}


def spider_core(arms):
    """Leaf-deleted spider with given arm lengths: centre 0, arm i of length a_i loses its
    end leaf, so it contributes a_i - 1 core vertices; parents = ends of shortened arms."""
    cv, ce, mk = [0], [], []
    nxt = 1
    for a in arms:
        prev = 0
        for _ in range(a - 1):
            cv.append(nxt); ce.append((prev, nxt)); prev = nxt; nxt += 1
        mk.append(prev)
    return cv, ce, sorted(set(mk))


ZCORES["spider_224_core"] = spider_core((2, 2, 4))
ZCORES["spider_233_core"] = spider_core((2, 3, 3))


def automorphisms(cv, ce):
    """All permutations of cv preserving the edge set ce (brute force; cores have <= 8 vertices)."""
    from itertools import permutations
    E = {frozenset(e) for e in ce}
    out = []
    for p in permutations(cv):
        s = dict(zip(cv, p))
        if all(frozenset((s[a], s[b])) in E for a, b in ce):
            out.append(s)
    return out


def placements_stream(cv, ce, mk, m, inS=()):
    """Like placements() but streams each distinct clause exactly once without a global set:
    a placement phi is emitted only if it is lexicographically least (on cv order) among the
    valid placements phi o sigma, sigma in Aut(F), which are exactly those with the same
    edge set."""
    var = edge_vars(m)
    order = bfs_order(cv, ce)
    cadj = {u: set() for u in cv}
    for a, b in ce:
        cadj[a].add(b); cadj[b].add(a)
    mks = set(mk)
    ins = set(inS)
    auts = [s for s in automorphisms(cv, ce) if any(s[u] != u for u in cv)]
    key = list(cv)
    phi = {}
    used = set()

    def minimal():
        base = [phi[u] for u in key]
        for s in auts:
            alt = [phi[s[u]] for u in key]
            if any(alt[j] >= m for j, u in enumerate(key) if u in mks):
                continue                      # phi o sigma not valid (marked vertex into S)
            if any(alt[j] < m for j, u in enumerate(key) if u in ins):
                continue                      # phi o sigma not in the generated family
            if alt < base:
                return False
        return True

    def rec(i):
        if i == len(order):
            if minimal():
                yield sorted(-ev(var, phi[a], phi[b]) for a, b in ce)
            return
        u = order[i]
        for h in range(2 * m):
            if h in used or (u in mks and h >= m) or (u in ins and h < m):
                continue
            if any(w in phi and ev(var, phi[w], h) is None for w in cadj[u]):
                continue
            phi[u] = h; used.add(h)
            yield from rec(i + 1)
            del phi[u]; used.discard(h)
    yield from rec(0)


def placements_soft(cv, ce, mk, m, d):
    """Soft-mark placements: marked vertices may go anywhere; a marked vertex placed at s in S
    contributes the literal -d_s.  Only placements with at least one marked vertex in S are
    emitted (the rest are exactly the rooted-core clauses).  Dedupe by mark-preserving
    automorphisms (identical clauses); occasional duplicates are harmless."""
    var = edge_vars(m)
    order = bfs_order(cv, ce)
    cadj = {u: set() for u in cv}
    for a, b in ce:
        cadj[a].add(b); cadj[b].add(a)
    mks = set(mk)
    auts = [s for s in automorphisms(cv, ce)
            if any(s[u] != u for u in cv) and {s[u] for u in mks} == mks]
    key = list(cv)
    phi, used = {}, set()

    def rec(i):
        if i == len(order):
            if not any(phi[p] >= m for p in mks):
                return
            base = [phi[u] for u in key]
            if any([phi[s[u]] for u in key] < base for s in auts):
                return
            yield sorted([-ev(var, phi[a], phi[b]) for a, b in ce] +
                         [-d[phi[p]] for p in mks if phi[p] >= m])
            return
        u = order[i]
        for h in range(2 * m):
            if h in used:
                continue
            if any(w in phi and ev(var, phi[w], h) is None for w in cadj[u]):
                continue
            phi[u] = h; used.add(h)
            yield from rec(i + 1)
            del phi[u]; used.discard(h)
    yield from rec(0)

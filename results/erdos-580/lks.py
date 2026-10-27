"""Trusted core for the Erdos #580 (vertex form) finite check.  Pure Python, no deps.

Reduced host H_m (n = 2m), from the edge-minimal-counterexample lemma (REPORT.md sec 2):
  L = {0..m-1}, S = {m..2m-1}; S independent; every x in L has deg >= m
  (we only impose >= m, a relaxation of 'exactly m', so UNSAT is a stronger statement);
  every s in S has deg <= m-1, i.e. is not adjacent to all of L.

Variables: one per possible edge {u,v} with u<v, u in L (L-L and L-S pairs), numbered
1.. in lexicographic order of (u,v).  The base CNF is auxiliary-variable-free:
  deg(x) >= m  <=>  every set of (2m-1)-m+1 = m potential neighbours of x contains an edge.

Avoidance clauses: for a tree T on m vertices let F = T - leaves(T) (the core) and
P = parents of the leaves.  A clause  OR_{ab in E(F)} -x_{phi(a)phi(b)}  for an injective
map phi: V(F) -> V(H) with phi(P) subset of L and phi(E(F)) inside L-L or L-S pairs.
(Any such clause is implied by "H contains no copy of F with P mapped into L", which by the
leaf-extension lemma is implied by "H contains no copy of T".)
"""
from itertools import combinations


# ---------------------------------------------------------------- trees
def canon(adj):
    """Canonical string of a free tree given as adjacency dict (AHU on the centre)."""
    n = len(adj)
    if n == 1:
        return "()"
    deg = {v: len(adj[v]) for v in adj}
    layer = [v for v in adj if deg[v] <= 1]
    left = n
    removed = set()
    while left > 2:
        left -= len(layer)
        nxt = []
        for v in layer:
            removed.add(v)
            for w in adj[v]:
                if w not in removed:
                    deg[w] -= 1
                    if deg[w] == 1:
                        nxt.append(w)
        layer = nxt
    centres = [v for v in adj if v not in removed]

    def enc(v, p):
        return "(" + "".join(sorted(enc(w, v) for w in adj[v] if w != p)) + ")"
    if len(centres) == 1:
        return enc(centres[0], None)
    a, b = centres
    return "|".join(sorted([enc(a, b), enc(b, a)]))


def trees(m):
    """All non-isomorphic trees on m vertices, as sorted edge lists on 0..m-1."""
    cur = {"()": [[]]}
    cur = {canon({0: []}): []}
    level = {canon({0: []}): []}
    for k in range(2, m + 1):
        nxt = {}
        for edges in level.values():
            for v in range(k - 1):
                e2 = edges + [(v, k - 1)]
                adj = {i: [] for i in range(k)}
                for a, b in e2:
                    adj[a].append(b); adj[b].append(a)
                c = canon(adj)
                if c not in nxt:
                    nxt[c] = e2
        level = nxt
    return [level[c] for c in sorted(level)]


def core(edges, m):
    """(core vertices, core edges, marked=parents of leaves) for an m-vertex tree."""
    adj = {i: set() for i in range(m)}
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    leaves = {v for v in adj if len(adj[v]) == 1}
    cv = sorted(set(adj) - leaves)
    ce = sorted((a, b) for a, b in edges if a in cv and b in cv)
    marked = sorted({w for v in leaves for w in adj[v]} - leaves)
    return cv, ce, marked


def g6(edges, m):
    """graph6 string (for cross-reference with NetworkX / Zeraoulia's table)."""
    bits = []
    E = {tuple(sorted(e)) for e in edges}
    for j in range(1, m):
        for i in range(j):
            bits.append(1 if (i, j) in E else 0)
    while len(bits) % 6:
        bits.append(0)
    s = chr(m + 63)
    for i in range(0, len(bits), 6):
        s += chr(63 + int("".join(map(str, bits[i:i + 6])), 2))
    return s


# ---------------------------------------------------------------- host encoding
def edge_vars(m):
    var = {}
    k = 0
    for u in range(m):
        for v in range(u + 1, 2 * m):
            k += 1
            var[(u, v)] = k
    return var


def ev(var, a, b):
    return var.get((a, b) if a < b else (b, a))


def base_clauses(m, var):
    out = []
    for x in range(m):
        nb = [v for v in range(2 * m) if v != x]          # 2m-1 potential neighbours
        lits = [ev(var, x, v) for v in nb]
        for sub in combinations(lits, len(lits) - m + 1):  # deg >= m
            out.append(list(sub))
    for s in range(m, 2 * m):                              # deg(s) <= m-1
        out.append([-ev(var, l, s) for l in range(m)])
    return out


def is_avoidance_clause(clause, cv, ce, marked, m, var):
    """Trusted check: does `clause` equal {-x_e : e in phi(E(F))} for some valid phi?
    Brute force over injective phi consistent with the clause's edge set."""
    inv = {v: e for e, v in var.items()}
    if any(l >= 0 or -l not in inv for l in clause):
        return False
    E = {inv[-l] for l in clause}
    if len(E) != len(ce) or len(set(clause)) != len(clause):
        return False
    verts = sorted({u for e in E for u in e})
    if len(verts) != len(cv):
        return False
    adjE = {u: set() for u in verts}
    for a, b in E:
        adjE[a].add(b); adjE[b].add(a)
    # find a bijection cv -> verts mapping ce onto E with marked -> L
    order = bfs_order(cv, ce)
    cadj = {u: set() for u in cv}
    for a, b in ce:
        cadj[a].add(b); cadj[b].add(a)
    mk = set(marked)
    phi = {}
    used = set()

    def rec(i):
        if i == len(order):
            return True
        u = order[i]
        for h in verts:
            if h in used or (u in mk and h >= m):
                continue
            ok = all((w not in phi) or (phi[w] in adjE[h]) for w in cadj[u])
            if not ok:
                continue
            # degree must match exactly since |E| = |ce| and bijective
            if len(adjE[h]) != len(cadj[u]):
                continue
            phi[u] = h; used.add(h)
            if rec(i + 1):
                return True
            del phi[u]; used.discard(h)
        return False
    return rec(0)


def bfs_order(cv, ce):
    if not cv:
        return []
    adj = {u: [] for u in cv}
    for a, b in ce:
        adj[a].append(b); adj[b].append(a)
    root = max(cv, key=lambda u: len(adj[u]))
    seen = [root]
    for u in seen:
        for w in sorted(adj[u], key=lambda w: -len(adj[w])):
            if w not in seen:
                seen.append(w)
    return seen


def write_dimacs(path, nvars, clauses, comments=()):
    with open(path, "w", newline="\n") as f:
        for c in comments:
            f.write("c " + c + "\n")
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for cl in clauses:
            f.write(" ".join(map(str, cl)) + " 0\n")


def read_dimacs(path):
    nv = nc = None
    cls = []
    cur = []
    with open(path) as f:
        for line in f:
            if line.startswith("c"):
                continue
            if line.startswith("p"):
                _, _, nv, nc = line.split()
                nv, nc = int(nv), int(nc)
                continue
            for t in line.split():
                x = int(t)
                if x == 0:
                    cls.append(cur); cur = []
                else:
                    cur.append(x)
    return nv, nc, cls


# ---------------------------------------------------------------- counter base (v2)
def _counter(lits, k, nxt):
    """Unary counter with full equivalences: r[i][j] <=> 'at least j of lits[:i] true',
    for j = 1..k.  Returns (clauses, r_last, next_free_var); r_last[j] is a literal or
    True/False.  Every aux variable is functionally determined by lits."""
    T, F = True, False
    cls = []
    prev = [T] + [F] * k
    for i, x in enumerate(lits, 1):
        cur = [T]
        for j in range(1, k + 1):
            a, b = prev[j], prev[j - 1]          # r = a or (x and b)
            if a is T:
                cur.append(T); continue
            if a is F and b is F:
                cur.append(F); continue
            r = nxt; nxt += 1
            # r -> a or x ; r -> a or b ; a -> r ; (x and b) -> r
            def lit(v):
                return v
            ca = [] if a is F else [a]
            if b is T:
                cls.append([-r] + ca + [x])
                if a is not F: cls.append([-a, r])
                cls.append([-x, r])
            else:
                cls.append([-r] + ca + [x])
                cls.append([-r] + ca + [b])
                if a is not F: cls.append([-a, r])
                cls.append([-x, -b, r])
            cur.append(r)
        prev = cur
    return cls, prev, nxt


def base_clauses_counter(m, var, exact=True):
    """deg(x) = m (or >= m) for x in L via _counter; S-degree <= m-1.  Returns (clauses, nvars)."""
    nxt = len(var) + 1
    out = []
    for x in range(m):
        lits = [ev(var, x, v) for v in range(2 * m) if v != x]
        cls, r, nxt = _counter(lits, m + 1, nxt)
        out += cls
        out.append([r[m]])                        # at least m
        if exact:
            out.append([-r[m + 1]])               # not at least m+1
    for s in range(m, 2 * m):
        out.append([-ev(var, l, s) for l in range(m)])
    return out, nxt - 1


def core_keep(edges, m, keep):
    """Generalised core: delete the leaves of T except those whose parent is in `keep`.
    Marked = parents of the deleted leaves.  The leaf-extension lemma holds verbatim for any
    set of deleted leaves (r = number deleted, F' has m-r vertices; a marked image in L has
    >= m - (m-r-1) = r+1 unused neighbours; Hall).  keep = () gives core(); keeping every
    parent gives the whole tree with nothing marked."""
    adj = {i: set() for i in range(m)}
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    leaves = {v for v in adj if len(adj[v]) == 1}
    par = {v: next(iter(adj[v])) for v in leaves}
    dele = {v for v in leaves if par[v] not in set(keep)}
    cv = sorted(set(adj) - dele)
    ce = sorted(tuple(sorted((a, b))) for a, b in edges if a in cv and b in cv)
    marked = sorted({par[v] for v in dele})
    return cv, ce, marked


def variant_core(edges, m, spec):
    """spec: '' (rooted core), 'full', or 'keep=p1.p2...'"""
    if not spec:
        return core(edges, m)
    if spec == "full":
        return list(range(m)), sorted(tuple(sorted(e)) for e in edges), []
    assert spec.startswith("keep=")
    return core_keep(edges, m, [int(x) for x in spec[5:].split(".") if x])


def components(edges, m, spec):
    """spec = '+'-joined components; each is '' or 'core' (rooted core), 'full', 'keep=p.q'
    or 'keepS=p.q' (same F' as keep=, but only placements with the kept parents in S are
    generated -- those with a kept parent in L are subsumed by rooted-core clauses).
    Returns a list of (cv, ce, marked, inS)."""
    out = []
    for c in (spec or "core").split("+"):
        if c in ("", "core"):
            out.append(core(edges, m) + ([],))
        elif c == "variants":   # every leaf-keep variant (validity only; used by audit_cnf)
            from itertools import combinations as _cb
            P = core(edges, m)[2]
            for r in range(1, len(P) + 1):
                for keep in _cb(P, r):
                    out.append(core_keep(edges, m, keep) + ("novars",))
        elif c == "soft":
            out.append(core(edges, m) + ("soft",))
        elif c.startswith("keepS="):
            keep = [int(x) for x in c[6:].split(".") if x]
            out.append(core_keep(edges, m, keep) + (keep,))
        else:
            out.append(variant_core(edges, m, c) + ([],))
    return out


# ---------------------------------------------------------------- soft marks (v3)
# Soft leaf-extension lemma: if phi embeds the core F (m-r vertices) and every marked vertex p
# has deg(phi(p)) >= m-1, phi extends to T: phi(p) has >= (m-1)-(m-r-1) = r unused neighbours,
# so any k <= r leaves have >= r >= k candidates (Hall).  L-vertices have degree m, so the
# condition only bites for marked vertices placed in S.  d_s <=> deg(s) >= m-1 (s in S).
def base_clauses_soft(m, var):
    base, nv = base_clauses_counter(m, var, exact=True)
    d = {s: nv + 1 + (s - m) for s in range(m, 2 * m)}
    for s in range(m, 2 * m):
        xs = [ev(var, l, s) for l in range(m)]
        for a, b in combinations(xs, 2):          # d_s -> at most one non-neighbour in L
            base.append([-d[s], a, b])
        for j in range(m):                        # at most one non-neighbour -> d_s
            base.append([d[s]] + [-xs[k] for k in range(m) if k != j])
    return base, nv + m, d


def is_soft_clause(clause, cv, ce, marked, m, var, d):
    """clause == {-x_e : e in phi(E(F))} + {-d_phi(p) : p marked, phi(p) in S} for some phi."""
    inv = {v: e for e, v in var.items()}
    dinv = {v: s for s, v in d.items()}
    if any(l >= 0 for l in clause) or len(set(clause)) != len(clause):
        return False
    el = [-l for l in clause if -l in inv]
    dl = {dinv[-l] for l in clause if -l in dinv}
    if len(el) + len(dl) != len(clause) or len(el) != len(ce):
        return False
    E = {inv[k] for k in el}
    verts = sorted({u for e in E for u in e})
    if len(verts) != len(cv):
        return False
    adjE = {u: set() for u in verts}
    for a, b in E:
        adjE[a].add(b); adjE[b].add(a)
    order = bfs_order(cv, ce)
    cadj = {u: set() for u in cv}
    for a, b in ce:
        cadj[a].add(b); cadj[b].add(a)
    mk = set(marked)
    phi, used = {}, set()

    def rec(i):
        if i == len(order):
            return {phi[p] for p in mk if phi[p] >= m} == dl
        u = order[i]
        for h in verts:
            if h in used or len(adjE[h]) != len(cadj[u]):
                continue
            if any(w in phi and phi[w] not in adjE[h] for w in cadj[u]):
                continue
            phi[u] = h; used.add(h)
            if rec(i + 1):
                return True
            del phi[u]; used.discard(h)
        return False
    return rec(0)

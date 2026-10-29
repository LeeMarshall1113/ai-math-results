"""Method B (independent of the lattice/universal-cover argument): exhaustive backtracking over
labelled face sets.  By Lemma 1 every vertex has degree 4 and its 9 star vertices are distinct,
so WLOG vertex 0 has neighbours 1,2,3,4 in cyclic order, faces (0,1,5,2),(0,2,6,3),(0,3,7,4),
(0,4,8,1), and 9 is the remaining vertex.  The other 6 faces are 4-cycles on {1..9}.
We require: every edge in exactly two faces, every vertex link one cycle, conditions C1-C3,
Euler characteristic 0; then split orientable (torus) / non-orientable (Klein bottle) and
classify up to isomorphism with canon.py.
Also reports, as a sanity check, the counts when C3 is dropped (only C1,C2) -- still with the
same root labelling, so that count is only a cross-check, not a classification."""
import itertools, sys
from canon import canon_form

V = 10
base = [(0, 1, 5, 2), (0, 2, 6, 3), (0, 3, 7, 4), (0, 4, 8, 1)]
def E(f): return [frozenset((f[i], f[(i + 1) % 4])) for i in range(4)]
# all quads on {1..9} up to rotation/reversal
quads = []
for S in itertools.combinations(range(1, 10), 4):
    a = S[0]
    for p in itertools.permutations(S[1:]):
        if p[0] < p[2]:
            quads.append((a,) + p)
by_edge = {}
for q in quads:
    for e in E(q): by_edge.setdefault(e, []).append(q)

def c2_ok(f, faces):
    ef = set(E(f))
    for g in faces:
        s = set(f) & set(g)
        if len(s) > 2: return False
        if len(s) == 2:
            e = frozenset(s)
            if e not in ef or e not in set(E(g)): return False
    return True

def link_ok(faces, v):
    # link of v: edges (a,c) for each face (.., a, v, c ..); must be a single cycle
    adj = {}
    for f in faces:
        if v in f:
            i = f.index(v); a, c = f[(i - 1) % 4], f[(i + 1) % 4]
            adj.setdefault(a, []).append(c); adj.setdefault(c, []).append(a)
    if any(len(x) != 2 for x in adj.values()): return False
    start = next(iter(adj)); prev, cur, n = None, start, 0
    while True:
        nx = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        prev, cur = cur, nx; n += 1
        if cur == start: break
    return n == len(adj)

def c3_ok(faces):
    for v in range(V):
        star = set(); d = 0
        for f in faces:
            if v in f: star |= set(f); d += 1
        if len(star) != 2 * d + 1: return False
    return True

sols = []
def rec(faces, ecount, deg):
    open_e = [e for e, c in ecount.items() if c == 1]
    if not open_e:
        if len(faces) == 10: sols.append(list(faces))
        return
    if len(faces) == 10: return
    e = min(open_e, key=lambda e: (len(by_edge.get(e, [])), sorted(e)))
    for q in by_edge.get(e, []):
        if q in faces: continue
        if any(ecount.get(x, 0) >= 2 for x in E(q)): continue
        new_edges = [x for x in E(q) if x not in ecount]
        d2 = dict(deg)
        ok = True
        for x in new_edges:
            for u in x:
                d2[u] = d2.get(u, 0) + 1
                if d2[u] > 4: ok = False
        if not ok or not c2_ok(q, faces): continue
        ec2 = dict(ecount)
        for x in E(q): ec2[x] = ec2.get(x, 0) + 1
        rec(faces + [q], ec2, d2)

ec, dg = {}, {}
for f in base:
    for x in E(f):
        if x not in ec:
            for u in x: dg[u] = dg.get(u, 0) + 1
        ec[x] = ec.get(x, 0) + 1
rec(list(base), ec, dg)
print("labelled completions (C1,C2, degree<=4, closed):", len(sols))
good = [s for s in sols if all(link_ok(s, v) for v in range(V))]
print("  with every link a single cycle:", len(good))
good3 = [s for s in good if c3_ok(s)]
print("  and C3:", len(good3))
classes, nonor = {}, 0
for s in good3:
    nv = len({u for f in s for u in f}); ne = len({x for f in s for x in E(f)})
    assert nv - ne + 10 == 0
    try:
        k = canon_form(s)
    except AssertionError:
        nonor += 1; continue
    classes.setdefault(k, s)
print("  orientable:", len(good3) - nonor, " non-orientable:", nonor)
print("orientable isomorphism classes:", len(classes))
from enum_lattice import quotient_map
ref = canon_form(quotient_map(((10, 0), (3, 1))))
for k, s in classes.items():
    print("  ", s, " == {4,4}_(3,1):", k == ref)

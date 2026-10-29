"""Canonical form of an orientable quadrangulation given as a list of cyclic 4-tuples
(simple graph assumed).  Canonical = lexicographic minimum, over all darts and both
orientations, of the BFS rotation-system code.  Two maps get the same form iff they are
isomorphic (allowing reflection)."""
from collections import deque

def orient(faces):
    faces = [tuple(f) for f in faces]
    darts = lambda f: [(f[i], f[(i + 1) % 4]) for i in range(4)]
    out = [None] * len(faces); out[0] = faces[0]
    q = deque([0])
    while q:
        i = q.popleft()
        D = set(darts(out[i]))
        for j, g in enumerate(faces):
            if out[j] is not None: continue
            dg = set(darts(g))
            rg = set(darts(g[::-1]))
            if any((b, a) in dg for (a, b) in D):
                out[j] = g; q.append(j)
            elif any((b, a) in rg for (a, b) in D):
                out[j] = g[::-1]; q.append(j)
    assert all(o is not None for o in out), "disconnected"
    # check coherence
    S = set()
    for f in out:
        for d in darts(f):
            assert d not in S, "non-orientable or not a surface"
            S.add(d)
    return out

def rotations(of):
    succ = {}
    for f in of:
        for i in range(4):
            u, a, c = f[i], f[(i + 1) % 4], f[(i - 1) % 4]
            succ[(u, c)] = a  # around u: c -> a
    return succ

def code_from(succ, v, w):
    lab = {v: 0, w: 1}; order = [v, w]; start = {v: w, w: v}
    code = []
    k = 0
    while k < len(order):
        u = order[k]; k += 1
        s = start[u]
        n = s; row = []
        while True:
            if n not in lab:
                lab[n] = len(order); order.append(n); start[n] = u
            row.append(lab[n])
            n = succ[(u, n)]
            if n == s: break
        code.append(tuple(row))
    return tuple(code)

def canon_form(faces):
    of = orient(faces)
    best = None
    for o in (of, [f[::-1] for f in of]):
        succ = rotations(o)
        for (u, n) in succ:
            c = code_from(succ, u, n)
            if best is None or c < best: best = c
    return best

"""Symmetry breaking by the normaliser N = N_G(H) of the prescribed group H (GAP output
gap/normalizers.txt, permutations of OUR 820 point indices).

N permutes the H-orbits of points.  If C_1..C_m are the N-orbits on H-orbits, listed in a chosen
order, then for every H-invariant solution K there is a first index j such that K contains an
H-orbit of C_j; applying an element of N we may assume K contains the representative R_j of C_j
and no orbit of C_1..C_{j-1}.  The m cases  (z_{R_j} = 1, z_O = 0 for O in C_1 u ... u C_{j-1})
therefore cover all H-invariant solutions (soundness does not depend on N being the full
normaliser: any subgroup of N_G(H) would do).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load_normalizers(path=None):
    path = path or os.path.join(HERE, "..", "gap", "normalizers.txt")
    res = {}
    cur = None
    for line in open(path):
        line = line.strip()
        if line.startswith("CLASS "):
            parts = line.split()
            cur = parts[1]
            res[cur] = dict(order=int(parts[3]), gens=[])
        elif line:
            res[cur]["gens"].append([int(x) - 1 for x in line.split(",")])
    return res


def classes_on_orbits(orbs, oid, gens):
    """N-orbits on the H-orbits (orbit indices), N given by permutations of the points."""
    n = len(orbs)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in gens:
        for o, pts in enumerate(orbs):
            a, b = find(o), find(oid[g[pts[0]]])
            if a != b:
                parent[a] = b
    cl = {}
    for o in range(n):
        cl.setdefault(find(o), []).append(o)
    return sorted(cl.values(), key=lambda c: c[0])


def cases(orbs, oid, gens, order=None):
    """list of (class_members, representative_orbit, forbidden_orbits)"""
    cls = classes_on_orbits(orbs, oid, gens)
    if order is not None:
        cls = [cls[i] for i in order]
    res = []
    forb = []
    for c in cls:
        res.append((c, c[0], list(forb)))
        forb = forb + c
    return res

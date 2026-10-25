"""Shared helpers: 3-graphs on <= 7 vertices as bitmasks over the 35 triples of range(7)."""
import itertools, pickle, sys
from k4enum import TRI, TIDX

TT = {}
for t in TRI:
    for p in itertools.permutations(t):
        TT[p] = 1 << TIDX[t]

def relabel(m, verts):
    out = 0
    for t in itertools.combinations(range(len(verts)), 3):
        if m & TT[(verts[t[0]], verts[t[1]], verts[t[2]])]:
            out |= 1 << TIDX[t]
    return out

def flag_key(m, typ, ext):
    best = None
    for p in itertools.permutations(ext):
        k = relabel(m, list(typ) + list(p))
        if best is None or k < best:
            best = k
    return best

cl = pickle.load(open("classes.pkl", "rb"))

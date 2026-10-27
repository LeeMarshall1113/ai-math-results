"""Semantic check of Zeraoulia's base (totalizer) clauses.

Claim to verify: projected onto edge variables, his base clauses are EXACTLY
  deg(x) = 9 for x in L  and  each s in S not adjacent to all of L.
Method: split the base into components by auxiliary variables; every component must touch
the 17 edge variables of a single L-vertex x (or be an S-degree clause).  For each component,
enumerate all 2^17 assignments of those edges and decide satisfiability of the component
with a SAT solver under assumptions; SAT must hold iff exactly 9 are true.
Run on the CNF given as argv[1] with its varmap (from diff_zenodo.py) in argv[2]."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from itertools import product
from pysat.solvers import Solver
from lks import read_dimacs, edge_vars

m = 9
cnf, vm = sys.argv[1], sys.argv[2]
nv, nc, cls = read_dimacs(cnf)
evs = list(map(int, open(vm).read().split()))
emap = {v: i + 1 for i, v in enumerate(evs)}          # his var -> our edge var
inv = {k: e for e, k in edge_vars(m).items()}         # our edge var -> (u,v)
base = cls[:3735]
# all clauses after the base must be pure avoidance clauses (checked by diff_zenodo.py)
assert all(all(l < 0 and -l in emap for l in c) for c in cls[3735:])
# union-find over aux vars
par = {}
def f(x):
    while par.setdefault(x, x) != x:
        par[x] = par[par[x]]; x = par[x]
    return x
for c in base:
    aux = [abs(l) for l in c if abs(l) not in emap]
    for a in aux[1:]:
        par[f(a)] = f(aux[0])
comps = {}
pure = []
for c in base:
    aux = [abs(l) for l in c if abs(l) not in emap]
    if not aux:
        pure.append(c)
    else:
        comps.setdefault(f(aux[0]), []).append(c)
print("components", len(comps), "edge-only base clauses", len(pure))
# edge-only base clauses must be the 9 S-degree clauses
sdeg = {frozenset(-edge_vars(m)[(l, s)] for l in range(m)) for s in range(m, 2 * m)}
assert {frozenset(-emap[-l] for l in c) for c in pure} == sdeg and len(pure) == 9
byx = {}
for k, cc in comps.items():
    E = sorted({abs(l) for c in cc for l in c if abs(l) in emap})
    verts = [inv[emap[e]] for e in E]
    common = set(verts[0])
    for uv in verts:
        common &= set(uv)
    assert len(common) == 1 and len(E) == 2 * m - 1, (len(E), common)
    byx.setdefault(common.pop(), []).append((E, cc))
seen_x = set(byx)
for x, parts in sorted(byx.items()):
    E = parts[0][0]
    assert all(p[0] == E for p in parts)
    cc = [c for p in parts for c in p[1]]
    bad = 0
    with Solver(name="minisat22", bootstrap_with=cc) as s:
        for bits in product((0, 1), repeat=len(E)):
            ass = [e if b else -e for e, b in zip(E, bits)]
            sat = s.solve(assumptions=ass)
            if sat != (sum(bits) == m):
                bad += 1
    print(f"vertex {x}: {len(parts)} components, {len(cc)} clauses, 2^{len(E)} assignments, mismatches {bad}")
    assert bad == 0
assert seen_x == set(range(m))
print("OK: base == (deg = 9 on L) & (S-degree <= 8), exactly, on edge variables")

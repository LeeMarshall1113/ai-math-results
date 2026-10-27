"""UNTRUSTED producer for trees whose rooted-core CNF is SAT.
Start from base + all rooted-core placement clauses; while SAT, take the model host H, find
the smallest leaf-keep variant F' of T (lks.core_keep; 'full' = T itself) that embeds in H with
its marked vertices in L, and add clauses for up to K of its placements plus R random images
under permutations of L and of S.  Every added clause is a placement clause of some variant,
which audit_cnf.py re-checks; the final UNSAT is certified by LRAT + cake_lpr as usual.
If some H contains no copy of T at all, it is printed as a COUNTEREXAMPLE.

usage: python hybrid.py m index extra_out [K R]
"""
import sys, os, time, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pysat.solvers import Solver
from lks import trees, core, core_keep, edge_vars, base_clauses_counter
from eager import placements_stream
from cegar import find_embeddings

m, ti, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
K = int(sys.argv[4]) if len(sys.argv) > 4 else 4
R = int(sys.argv[5]) if len(sys.argv) > 5 else 50
T = trees(m)[ti]
cv, ce, mk = core(T, m)
var = edge_vars(m)
P = mk
variants = []
for r in range(1, len(P) + 1):
    for keep in itertools.combinations(P, r):
        variants.append(core_keep(T, m, keep))
variants.sort(key=lambda v: (len(v[0]), -len(v[2])))
s = Solver(name="cadical195")
for c in base_clauses_counter(m, var, exact=True)[0]:
    s.add_clause(c)
n0 = 0
for c in placements_stream(cv, ce, mk, m):
    s.add_clause(c); n0 += 1
rng = random.Random(1)
extra, seen = [], set()
t0 = time.time()
it = 0
while s.solve():
    it += 1
    mdl = s.get_model()
    H = {v: set() for v in range(2 * m)}
    for (a, b), k in var.items():
        if mdl[k - 1] > 0:
            H[a].add(b); H[b].add(a)
    embs = None
    for j, (vcv, vce, vmk) in enumerate(variants):
        last = j == len(variants) - 1          # the full tree: search exhaustively
        embs = find_embeddings(H, m, vcv, vce, vmk, K, rng, budget=None if last else 20000)
        if embs:
            break
    if not embs:
        print("COUNTEREXAMPLE host", {v: sorted(H[v]) for v in H}, flush=True)
        sys.exit(3)
    imgs = list(embs)
    for _ in range(R):
        pl = list(range(m)); rng.shuffle(pl)
        ps = list(range(m, 2 * m)); rng.shuffle(ps)
        sig = pl + ps
        imgs += [frozenset(tuple(sorted((sig[a], sig[b]))) for a, b in E) for E in embs]
    for E in imgs:
        cl = tuple(sorted(-var[e] for e in E))
        if cl not in seen:
            seen.add(cl); extra.append(cl); s.add_clause(list(cl))
    if it % 20 == 0 or it < 5:
        print(f"  it {it} extra {len(extra)} last variant size {len(vcv)} t {time.time()-t0:.0f}s", flush=True)
with open(outp, "w", newline="\n") as f:
    for c in extra:
        f.write(" ".join(map(str, c)) + " 0\n")
print(f"hybrid m={m} tree {ti}: UNSAT after {it} iterations, {n0} core + {len(extra)} extra clauses,"
      f" {time.time()-t0:.0f}s")

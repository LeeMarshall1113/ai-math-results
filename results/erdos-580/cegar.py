"""UNTRUSTED producer: lazily finds a small set of avoidance clauses whose conjunction with
the base is UNSAT, for one tree.  Its output is re-checked by audit.py (each clause must be
a genuine rooted-core placement) and the UNSAT claim by CaDiCaL LRAT + cake_lpr.

usage: python cegar.py m tree_index outdir [K]
"""
import sys, os, json, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from lks import trees, core, edge_vars, ev, bfs_order, g6, canon, base_clauses_counter


def find_embeddings(H, m, cv, ce, marked, K, rng, full=False, budget=None):
    """Up to K placements (as frozensets of host edges) of the core in host H
    (adjacency sets), marked vertices into L.  If full, cv/ce is the whole tree, no marks."""
    order = bfs_order(cv, ce)
    cadj = {u: set() for u in cv}
    for a, b in ce:
        cadj[a].add(b); cadj[b].add(a)
    parent = {}
    for i, u in enumerate(order):
        for w in order[:i]:
            if w in cadj[u]:
                parent[u] = w; break
    mk = set(marked)
    found = set()
    phi = {}
    used = set()
    verts = list(range(2 * m))

    nodes = [0]

    def rec(i):
        nodes[0] += 1
        if len(found) >= K or (budget and nodes[0] > budget):
            return
        if i == len(order):
            found.add(frozenset(tuple(sorted((phi[a], phi[b]))) for a, b in ce))
            return
        u = order[i]
        cand = list(H[phi[parent[u]]]) if u in parent else verts[:]
        rng.shuffle(cand)
        for h in cand:
            if h in used or (u in mk and h >= m):
                continue
            if any(w in phi and phi[w] not in H[h] for w in cadj[u]):
                continue
            phi[u] = h; used.add(h)
            rec(i + 1)
            del phi[u]; used.discard(h)
            if len(found) >= K or (budget and nodes[0] > budget):
                return
    rec(0)
    return found


def run(m, ti, outdir, K=8, exact=False, seed=1, timeout=None, R=0, verbose=False):
    T = trees(m)[ti]
    cv, ce, marked = core(T, m)
    var = edge_vars(m)
    nE = len(var)
    name = f"m{m}_t{ti:03d}"
    info = dict(m=m, tree=ti, edges=T, g6=g6(T, m), core_vertices=cv, core_edges=ce,
                marked=marked)
    if not ce:
        info["status"] = "trivial (core has no edge: star; any L-vertex has degree m >= m-1)"
        return info, []
    s = Solver(name="cadical195")
    for c in base_clauses_counter(m, var, exact=exact)[0]:
        s.add_clause(c)
    rng = random.Random(seed)
    clauses = []
    seen = set()
    it = 0
    t0 = time.time()
    while s.solve():
        it += 1
        mdl = s.get_model()
        H = {v: set() for v in range(2 * m)}
        for (a, b), k in var.items():
            if mdl[k - 1] > 0:
                H[a].add(b); H[b].add(a)
        embs = find_embeddings(H, m, cv, ce, marked, K, rng)
        if not embs:
            full = find_embeddings(H, m, list(range(m)), sorted(T), [], 1, rng)
            if not full:
                info["status"] = "COUNTEREXAMPLE"
                info["host"] = {v: sorted(H[v]) for v in H}
                return info, clauses
            info["status"] = "core lemma insufficient (host has T but no rooted core)"
            info["host"] = {v: sorted(H[v]) for v in H}
            return info, clauses
        imgs = list(embs)
        for _ in range(R):
            pl = list(range(m)); rng.shuffle(pl)
            ps = list(range(m, 2 * m)); rng.shuffle(ps)
            sig = pl + ps
            imgs += [frozenset(tuple(sorted((sig[a], sig[b]))) for a, b in E) for E in embs]
        if verbose and (it < 20 or it % 50 == 0):
            print(f"  it {it} clauses {len(clauses)} t {time.time()-t0:.0f}s", flush=True)
        for E in imgs:
            cl = sorted(-var[e] for e in E)
            key = tuple(cl)
            if key in seen:
                continue
            seen.add(key)
            clauses.append(cl)
            s.add_clause(cl)
        if timeout and time.time() - t0 > timeout:
            info["status"] = "timeout"
            return info, clauses
    info["status"] = "UNSAT"
    info["iterations"] = it
    info["clauses"] = len(clauses)
    info["seconds"] = round(time.time() - t0, 1)
    return info, clauses


if __name__ == "__main__":
    m, ti, outdir = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    K = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    R = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    exact = not (len(sys.argv) > 6 and sys.argv[6] == "atleast")
    os.makedirs(outdir, exist_ok=True)
    info, cls = run(m, ti, outdir, K, exact, R=R, verbose=True)
    name = f"m{m}_t{ti:03d}"
    with open(os.path.join(outdir, name + ".av"), "w") as f:
        for c in cls:
            f.write(" ".join(map(str, c)) + "\n")
    with open(os.path.join(outdir, name + ".json"), "w") as f:
        json.dump(info, f)
    print(name, info["status"], info.get("iterations"), info.get("clauses"), info.get("seconds"))

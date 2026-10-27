"""Diff Zeraoulia's 4 published CNFs against our independent regeneration.
Avoidance clauses: compared as sets of literal-sets under the lexicographic edge-variable map.
Base clauses: his totalizer base is checked semantically (see REPORT.md) -- here we only
report which clauses use auxiliary variables."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lks import edge_vars, read_dimacs, core
from eager import placements, ZCORES
Z = "zenodo/"
m = 9
var = edge_vars(m)
nE = len(var)
for name, (cv, ce, mk) in ZCORES.items():
    nv, nc, cls = read_dimacs(Z + name + ".cnf")
    mine = placements(cv, ce, mk, m)
    edge_only = [c for c in cls if all(abs(l) <= nE for l in c)]
    aux = [c for c in cls if any(abs(l) > nE for l in c)]
    theirs_av = {frozenset(c) for c in edge_only if all(l < 0 for l in c) and len(c) == len(ce)}
    other = [c for c in edge_only if not (all(l < 0 for l in c) and len(c) == len(ce))]
    print(f"{name}: vars {nv} clauses {nc}; aux-using {len(aux)}; edge-only {len(edge_only)}"
          f" of which avoidance-shaped {len(theirs_av)}, other {len(other)}")
    print(f"   ours {len(mine)}  equal={mine == theirs_av}  ours-theirs={len(mine-theirs_av)}"
          f"  theirs-ours={len(theirs_av-mine)}")
    print("   other edge-only clauses:", other[:12], "...")

# ---- variable map: his edge variables are interleaved with totalizer auxiliaries.
# Hypothesis: the k-th smallest variable occurring in his avoidance + S-degree clauses is our
# edge variable k (both lexicographic in (u,v)).  Test it on all four instances.
print("--- under sorted-order variable map ---")
for name, (cv, ce, mk) in ZCORES.items():
    nv, nc, cls = read_dimacs(Z + name + ".cnf")
    neg = [c for c in cls if all(l < 0 for l in c) and (len(c) == len(ce) or len(c) == m)]
    evs = sorted({-l for c in neg for l in c})
    assert len(evs) == nE, len(evs)
    mp = {v: i + 1 for i, v in enumerate(evs)}
    av = [c for c in cls if all(l < 0 and -l in mp for l in c)]
    th = {frozenset(-mp[-l] for l in c) for c in av if len(c) == len(ce)}
    sdeg = {frozenset(-mp[-l] for l in c) for c in av if len(c) == m}
    mine = placements(cv, ce, mk, m)
    mysdeg = {frozenset(-var[(l, s)] for l in range(m)) for s in range(m, 2 * m)}
    print(f"{name}: avoidance equal={th == mine} ({len(th)} vs {len(mine)}),"
          f" S-degree clauses equal={sdeg == mysdeg}, his edge vars = {evs[:4]}..{evs[-2:]}")
    # also save the map for the base-semantics check
    with open(f"work/{name}.varmap", "w") as f:
        f.write(" ".join(str(v) for v in evs) + "\n")

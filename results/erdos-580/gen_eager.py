"""Write the eager CNF for tree `index` on m vertices:
counter base (deg = m on L, S-degree <= m-1) + one clause per placement of the rooted core
(lks.core: leaves deleted, parents -> L), each distinct clause once (streamed, constant memory).
usage: python gen_eager.py m index out.cnf [full | keep=p1.p2]
(full: clauses for placements of the whole tree T, no leaf lemma, no marks)
"""
import sys, os, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lks import trees, core, edge_vars, base_clauses_counter, base_clauses_soft, g6, components
from eager import placements_stream, placements_soft

m, ti, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
T = trees(m)[ti]
MAXC = int(os.environ.get("MAX_CLAUSES", "40000000"))
SPEC = sys.argv[4] if len(sys.argv) > 4 else ""
FULL = SPEC == "full"
comps = components(T, m, SPEC)
assert all(c[1] for c in comps), "star: trivial"
var = edge_vars(m)
SOFT = any(c[3] == "soft" for c in comps)
if SOFT:
    base, nv, D = base_clauses_soft(m, var)
else:
    base, nv = base_clauses_counter(m, var, exact=True)
body = out + ".body"
n = 0
with open(body, "w", newline="\n") as f:
    for cv, ce, mk, ins in comps:
      if ins == "novars":
        continue                 # 'variants' only declares what the extra clauses may be
      gen = placements_soft(cv, ce, mk, m, D) if ins == "soft" else placements_stream(cv, ce, mk, m, ins)
      for c in gen:
        f.write(" ".join(map(str, c)) + " 0\n")
        n += 1
        if n > MAXC:
            break
      if n > MAXC:
        break
    if os.environ.get("EXTRA_CLAUSES"):       # e.g. hybrid.py output, re-checked by audit_cnf
        for line in open(os.environ["EXTRA_CLAUSES"]):
            if line.strip():
                f.write(line if line.endswith("\n") else line + "\n"); n += 1
if n > MAXC:
    os.remove(body)
    sys.exit(f"more than {MAXC} placement clauses: skipped")
with open(out, "w", newline="\n") as f:
    f.write(f"c Erdos 580 vertex form, n={2*m}, tree #{ti} g6={g6(T, m)} (eager, {SPEC or "rooted core"})\n")
    for cv, ce, mk, ins in comps:
        f.write(f"c component: vertices {cv} edges {ce} marked(->L) {mk} forced(->S) {ins}\n")
    f.write(f"c {len(base)} base clauses, {n} placement clauses\n")
    f.write(f"p cnf {nv} {len(base) + n}\n")
    for c in base:
        f.write(" ".join(map(str, c)) + " 0\n")
    with open(body) as b:
        shutil.copyfileobj(b, f, 1 << 22)
os.remove(body)
print(f"{out}: {len(base)} base + {n} placement clauses, {nv} vars")

"""Independent re-check of a generated certificate CNF (does not use eager.py):
the header and the first clauses must equal lks.base_clauses_counter exactly, and every
remaining clause must be a genuine rooted-core placement clause (lks.is_avoidance_clause).
usage: python audit_cnf.py m index file.cnf"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lks import trees, core, components, edge_vars, base_clauses_counter, base_clauses_soft, is_soft_clause, is_avoidance_clause, read_dimacs

m, ti, path = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
T = trees(m)[ti]
comps = components(T, m, sys.argv[4] if len(sys.argv) > 4 else "")
# a clause is accepted if it is a placement clause of some component's F' (where the forced-
# into-S restriction is irrelevant for validity: any placement of F' with marks in L is valid)
bylen = {}
for cv, ce, mk, ins in comps:
    if ins != "soft":
        bylen.setdefault(len(ce), []).append((cv, ce, mk))
soft = [(cv, ce, mk) for cv, ce, mk, ins in comps if ins == "soft"]
var = edge_vars(m)
SOFT = any(c[3] == "soft" for c in comps)
if SOFT:
    base, nv, D = base_clauses_soft(m, var)
else:
    base, nv = base_clauses_counter(m, var, exact=True)
hv, hc, cls = read_dimacs(path)
assert hv == nv and hc == len(cls), "header"
assert cls[:len(base)] == base, "base mismatch"
cache = {}
for c in cls[len(base):]:
    if not (any(is_avoidance_clause(c, cv, ce, mk, m, var) for cv, ce, mk in bylen.get(len(c), []))
            or any(is_soft_clause(c, cv, ce, mk, m, var, D) for cv, ce, mk in soft)):
        print("REJECT", c); sys.exit(1)
print(f"AUDIT OK {os.path.basename(path)}: base {len(base)} exact, {len(cls)-len(base)} placement clauses valid")

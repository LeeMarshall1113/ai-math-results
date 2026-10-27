"""Untrusted helper: for a tree whose rooted-core CNF is SAT, try keeping the leaves of
subsets of parents (smallest resulting F' first) until the eager CNF is UNSAT within a
short CaDiCaL limit.  Prints the spec to pass to run_certs.py.
usage: python find_variant.py m index [max_clauses]"""
import sys, os, subprocess, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lks import trees, core, core_keep
from eager import placements_stream
CAD = "cadical"
m, i = int(sys.argv[1]), int(sys.argv[2])
cap = int(sys.argv[3]) if len(sys.argv) > 3 else 3_000_000
T = trees(m)[i]
_, _, P = core(T, m)
cands = []
for r in range(1, len(P) + 1):
    for keep in itertools.combinations(P, r):
        cv, ce, mk = core_keep(T, m, keep)
        cands.append((len(cv), -len(mk), keep))
cands.sort()
os.makedirs(os.path.join(HERE, "work", "fv"), exist_ok=True)
for size, _, keep in cands:
    spec = "keep=" + ".".join(map(str, keep))
    cnf = os.path.join(HERE, "work", "fv", f"m{m}_t{i}.cnf")
    cv, ce, mk = core_keep(T, m, keep)
    n = 0
    for _ in placements_stream(cv, ce, mk, m):
        n += 1
        if n > cap:
            break
    if n > cap:
        print(spec, "too many clauses", flush=True); continue
    subprocess.run([sys.executable, os.path.join(HERE, "gen_eager.py"), str(m), str(i), cnf, spec],
                   capture_output=True)
    try:
        rc = subprocess.run([CAD, "-q", cnf], capture_output=True, text=True, timeout=60).returncode
    except subprocess.TimeoutExpired:
        rc = "timeout"
    print(spec, "size", size, "marked", mk, "clauses", n, "rc", rc, flush=True)
    if rc == 20:
        print("FOUND", f"{m}:{i}@{spec}"); break

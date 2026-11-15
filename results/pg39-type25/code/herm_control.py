"""Positive control: the orbit-SAT pipeline (allowed counts {1,4,10}, |K| = 280, invariant under an element of order 7)
finds a set with the fingerprint of the Hermitian surface H(3,9); compared with the explicit surface."""
import json, os, sys, time
from collections import Counter
import groups
from orbsat import OrbitInstance, run_solver, check_K
from controls import hermitian_points, line_stats, plane_stats
from pg import PG

F, P, gens, reps = groups.build_all()
M, f = reps["7_singer3"]
inst = OrbitInstance(P, [P.perm_of_matrix(M, f)])
cnf = inst.build_cnf(280, (1, 4, 10), amo=True)
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cnf", "herm_control_7.cnf")
cnf.to_dimacs(path)
t = time.time()
st, model, secs, info = run_solver(path, "cadical", 300)
tr = set(m for m in model if m > 0)
K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
ok, bad = check_K(P, K, (1, 4, 10))
H = hermitian_points(P, 3)
res = dict(status=st, solver_seconds=round(secs, 2), conflicts=info.get("conflicts"), size=len(K), exact_line_check=ok,
           lines_found=dict(line_stats(P, K)), planes_found=dict(plane_stats(P, K)),
           lines_explicit_H=dict(line_stats(P, H)), planes_explicit_H=dict(plane_stats(P, H)))
print(json.dumps(res))
assert dict(line_stats(P, K)) == dict(line_stats(P, H)) and dict(plane_stats(P, K)) == dict(plane_stats(P, H))
json.dump(dict(res, K=K), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "hermitian_control.json"), "w"))
print("fingerprint of the SAT-found set equals that of the explicit Hermitian surface")

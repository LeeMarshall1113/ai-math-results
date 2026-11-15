"""Class 2a (homology g = diag(1,1,1,-1), axis plane pi0 = {x4 = 0}, centre C = (0,0,0,1)).
K n pi0 = S0 is a planar (2,5)-set of size m in {26,35}; with c = [C in K] the sizes force
 (c=0,m=26,k=260) (c=0,m=35,k=287) (c=1,m=26,k=287) (c=1,m=35,k=260).
The centraliser induces PGammaL(3,9) on pi0, so WLOG S0 is one of the class representatives of
results/planar_classes.json.  One SAT instance per (class, c)."""
import json, os, sys, time
import groups
from orbsat import OrbitInstance, run_solver, check_K
from pg import PG
HERE = os.path.dirname(os.path.abspath(__file__))
F, P, gens, reps = groups.build_all()
M, f = reps["2a_homology_3+1"]
inst = OrbitInstance(P, [P.perm_of_matrix(M, f)])
P2 = PG(3, 9, F)
classes = json.load(open(os.path.join(HERE, "..", "results", sys.argv[2] if len(sys.argv) > 2 else "planar_classes.json")))
timeout = int(sys.argv[1]) if len(sys.argv) > 1 else 600
pi0 = {}   # P2 index -> P index of (x1,x2,x3,0)
for i, v in enumerate(P2.points):
    pi0[i] = P.pidx[P.normalize((v[0], v[1], v[2], 0))]
C = P.pidx[(0, 0, 0, 1)]
oid = inst.oid
cases = []
for ci, cl in enumerate(classes):
    m = cl["k"]
    for c in (0, 1):
        k = {(0, 26): 260, (0, 35): 287, (1, 26): 287, (1, 35): 260}[(c, m)]
        cases.append((ci, m, c, k, cl))
only = None
if len(sys.argv) > 3:
    only = set(int(x) for x in sys.argv[3].split(","))
for n, (ci, m, c, k, cl) in enumerate(cases):
    if only is not None and n not in only:
        continue
    S0 = set(pi0[p] for p in cl["rep"])
    units = []
    for i in range(P2.N):
        units.append((oid[pi0[i]], 1 if pi0[i] in S0 else 0))
    units.append((oid[C], c))
    cnf = inst.build_cnf(k, (2, 5), amo=True, units=units)
    path = os.path.join(HERE, "..", "cnf", "2a_c%d_m%d_cls%d_k%d.cnf" % (c, m, ci, k))
    cnf.to_dimacs(path, comments=["class 2a k %d c %d m %d class %d aut %d" % (k, c, m, ci, cl["aut_order"])])
    st, model, t, info = run_solver(path, "cadical", timeout)
    rec = dict(case=n, cls_index=ci, m=m, c=c, k=k, aut=cl["aut_order"], status=st, seconds=round(t, 1), **info)
    print(json.dumps(rec), flush=True)
    if st == "SATISFIABLE":
        tr = set(x for x in model if x > 0)
        K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
        ok, bad = check_K(P, K, (2, 5))
        print("  MODEL |K| = %d exact line check %s" % (len(K), ok), flush=True)
        json.dump(K, open(os.path.join(HERE, "..", "results", "MODEL_2a_case%d.json" % n), "w"))

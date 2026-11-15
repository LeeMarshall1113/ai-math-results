"""Class 2b (g = diag(1,1,-1,-1), fixed lines l1 = P(<e1,e2>), l2 = P(<e3,e4>)).
F1 = K n l1, F2 = K n l2, sizes (|F1|,|F2|,k) in {(5,2,287),(2,5,287)~swap,(5,5,260),(2,2,260)} (pencil counts, see NOTES).
GL(2,9) x GL(2,9) acts independently on l1, l2 (it centralises g), so (F1,F2) may be replaced by representatives of
PGL(2,9)-orbits on subsets of PG(1,9) separately.  One SAT instance per (F1-class, F2-class, k)."""
import itertools, json, os, sys, time
import groups
from orbsat import OrbitInstance, run_solver, check_K
from pg import PG
from ff import GF
HERE = os.path.dirname(os.path.abspath(__file__))
F, P, gens, reps = groups.build_all()
M, f = reps["2b_homology_2+2"]
inst = OrbitInstance(P, [P.perm_of_matrix(M, f)])
P1 = PG(2, 9, F)       # PG(1,9)
omega = F.prim
pgl_gens = [P1.perm_of_matrix([[1, 1], [0, 1]]), P1.perm_of_matrix([[omega, 0], [0, 1]]), P1.perm_of_matrix([[0, 1], [1, 0]])]
def orbits_pgl(sz):
    subsets = [tuple(c) for c in itertools.combinations(range(10), sz)]
    seen = {}
    reps_ = []
    for s in subsets:
        if s in seen: continue
        orb = {s}; fr = [s]
        while fr:
            nx = []
            for t in fr:
                for g in pgl_gens:
                    u = tuple(sorted(g[i] for i in t))
                    if u not in orb: orb.add(u); nx.append(u)
            fr = nx
        for t in orb: seen[t] = len(reps_)
        reps_.append((s, len(orb)))
    return reps_
o2, o5 = orbits_pgl(2), orbits_pgl(5)
print("PGL(2,9)-orbits on 2-subsets: %s ; on 5-subsets: %s" % ([r[1] for r in o2], [r[1] for r in o5]), flush=True)
# point index on l1 / l2 of PG(3,9) from PG(1,9) index
def on_l1(i):
    v = P1.points[i]; return P.pidx[P.normalize((v[0], v[1], 0, 0))]
def on_l2(i):
    v = P1.points[i]; return P.pidx[P.normalize((0, 0, v[0], v[1]))]
oid = inst.oid
cases = []
for (s1, n1) in o5:
    for (s2, n2) in o2:
        cases.append((5, s1, 2, s2, 287))
for (s1, n1) in o5:
    for (s2, n2) in o5:
        cases.append((5, s1, 5, s2, 260))
for (s1, n1) in o2:
    for (s2, n2) in o2:
        cases.append((2, s1, 2, s2, 260))
timeout = int(sys.argv[1]) if len(sys.argv) > 1 else 300
for n, (a1, s1, a2, s2, k) in enumerate(cases):
    S1 = set(on_l1(i) for i in s1); S2 = set(on_l2(i) for i in s2)
    units = [(oid[on_l1(i)], 1 if on_l1(i) in S1 else 0) for i in range(10)] + [(oid[on_l2(i)], 1 if on_l2(i) in S2 else 0) for i in range(10)]
    cnf = inst.build_cnf(k, (2, 5), amo=True, units=units)
    path = os.path.join(HERE, "..", "cnf", "2b_case%d_k%d.cnf" % (n, k))
    cnf.to_dimacs(path, comments=["class 2b case %d |F1|=%d |F2|=%d k %d F1 %s F2 %s" % (n, a1, a2, k, s1, s2)])
    st, model, t, info = run_solver(path, "cadical", timeout)
    print(json.dumps(dict(case=n, F1=a1, F2=a2, k=k, status=st, seconds=round(t, 1), **info)), flush=True)
    if st == "SATISFIABLE":
        tr = set(x for x in model if x > 0)
        K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
        ok, bad = check_K(P, K, (2, 5))
        print("  MODEL |K| = %d exact line check %s" % (len(K), ok), flush=True)
        json.dump(K, open(os.path.join(HERE, "..", "results", "MODEL_2b_case%d.json" % n), "w"))

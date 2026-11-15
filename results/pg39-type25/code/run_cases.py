"""Run a class with the normaliser case split:  python run_cases.py CLASS K [--timeout S] [--lrat] [--only j]
Cases are listed in the order given by symbreak.cases() with the N-classes sorted so that classes of
'fixed' objects come first (index order as printed)."""
import argparse, json, os, sys, time
import groups, symbreak
from orbsat import OrbitInstance, run_solver, check_K
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cls"); ap.add_argument("k", type=int)
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--order", default="")     # comma list of class indices
    ap.add_argument("--only", type=int, default=-1)
    ap.add_argument("--lrat", action="store_true")
    a = ap.parse_args()
    F, P, gens, reps = groups.build_all()
    N = symbreak.load_normalizers()
    M, f = reps[a.cls]
    inst = OrbitInstance(P, [P.perm_of_matrix(M, f)])
    order = [int(x) for x in a.order.split(",")] if a.order else None
    cs = symbreak.cases(inst.orbs, inst.oid, N[a.cls]["gens"], order)
    print(a.cls, "k=%d" % a.k, "cases:", [(len(c[0]), c[1], len(c[2])) for c in cs], flush=True)
    for j, (members, rep, forb) in enumerate(cs):
        if a.only >= 0 and j != a.only: continue
        units = [(rep, 1)] + [(o, 0) for o in forb]
        cnf = inst.build_cnf(a.k, (2, 5), amo=True, units=units)
        tag = "%s_k%d_case%d" % (a.cls.replace("+", "p"), a.k, j)
        cnfp = os.path.join(HERE, "..", "cnf", tag + ".cnf")
        cnf.to_dimacs(cnfp, comments=["class %s k %d case %d rep %d forbidden %d" % (a.cls, a.k, j, rep, len(forb))])
        lr = os.path.join(HERE, "..", "certs", tag + ".lrat") if a.lrat else None
        st, model, t, info = run_solver(cnfp, "cadical", a.timeout, lr)
        rec = dict(cls=a.cls, k=a.k, case=j, rep=rep, nforb=len(forb), status=st, seconds=round(t, 1), vars=cnf.nvars, clauses=len(cnf.clauses), **info)
        if st == "SATISFIABLE":
            tr = set(m for m in model if m > 0)
            K = sorted(p for i, v in enumerate(inst.z) if v in tr for p in inst.orbs[i])
            ok, bad = check_K(P, K, (2, 5))
            rec.update(model_size=len(K), line_check=ok, K=K)
        print(json.dumps({x: y for x, y in rec.items() if x != "K"}), flush=True)
        open(os.path.join(HERE, "..", "results", "case_runs.jsonl"), "a").write(json.dumps(rec) + "\n")

if __name__ == "__main__":
    main()

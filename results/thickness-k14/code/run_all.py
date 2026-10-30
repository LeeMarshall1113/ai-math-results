"""Production cube-and-conquer for "K_14 has geometric thickness 3?" with per-leaf LRAT + cake_lpr.

Cases (p1 normal form: vertices 1, 2, n are always hull vertices, and hull order = label
order); S = the other hull vertices among 3..n-1:
  hS  : |S| <= 2, exact hull units, plus sweep lemmas for every hull edge (gen.Enc.sweep);
  ge6 : a 3-subset T of 3..n-1 forced to be hull (covers every |S| >= 3), no sweeps.
Every node (case, cube over orientation variables) is solved by CaDiCaL with an LRAT proof and a
wall-clock limit; UNSAT proofs are checked by $LRAT_CHECK (cake_lpr) and then deleted (SHA-256 kept). On a
timeout the node is split on the next k free orientation variables into all 2^k sign patterns.
A SAT node writes its model. Records go to $TK_RES/tree.jsonl (append-only, resumable).
cover_check.py verifies completeness; replay.py re-certifies the leaves of a finished tree.
Solver and checker come from the environment, see tools.py.

usage: python run_all.py n t [--workers 2] [--timeout 300] [--k 4] [--only CASE]
"""
import argparse, hashlib, itertools, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor
from gen import Enc
import tools

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.environ.get("TK_WORK", os.path.join(HERE, "work")))
RES = os.path.abspath(os.environ.get("TK_RES", os.path.join(HERE, "results")))


def case_list(n):
    mid = list(range(3, n))
    cs = [("h3", ())]
    cs += [("h4", (v,)) for v in mid]
    cs += [("h5", T) for T in itertools.combinations(mid, 2)]
    cs += [("ge6", T) for T in itertools.combinations(mid, 3)]
    return cs


def case_clauses(n, t, kind, S):
    e = Enc(n, t, euler=3)
    for v in (1, 2, n):
        e.cl.append([e.H[v]])
    if kind == "ge6":
        for v in S:
            e.cl.append([e.H[v]])
    else:
        for v in range(3, n):
            e.cl.append([e.H[v] if v in S else -e.H[v]])
        hv = sorted(set(S) | {1, 2, n})
        for i in range(len(hv) if not os.environ.get("TK_NOSWEEP") else 0):
            e.sweep(hv[i], hv[(i + 1) % len(hv)])
    return e.nv, e.cl, e


sha = tools.sha_file


def cname(case, cube):
    kind, S = case
    c = "".join(("p" if l > 0 else "m") + str(abs(l)) for l in cube) or "root"
    return f"{kind}_{'-'.join(map(str, S)) or 'none'}_{c}"


class Runner:
    def __init__(self, n, t, timeout, k):
        self.n, self.t, self.timeout, self.k = n, t, timeout, k
        X = {tr: None for tr in itertools.combinations(range(1, n + 1), 3)}
        self.Xid = {tr: i + 1 for i, tr in enumerate(X)}  # gen.Enc numbering (X first)
        self.cand = [self.Xid[tr] for tr in sorted((tr for tr in X if tr[0] != 1),
                                                   key=lambda tr: (tr[2] - tr[0], tr))]
        self.cache = {}
        os.makedirs(WORK, exist_ok=True); os.makedirs(RES, exist_ok=True)
        open(os.path.join(RES, "tree.jsonl"), "a").close()
        self.log = open(os.path.join(RES, "tree.jsonl"), "a")
        self.done = set()
        p = os.path.join(RES, "tree.jsonl")
        for line in open(p):
            r = json.loads(line)
            self.done.add((r["kind"], tuple(r["S"]), tuple(r["cube"])))

    def record(self, **r):
        self.log.write(json.dumps(r) + "\n"); self.log.flush()

    def solve(self, case, cube, nv, cl):
        kind, S = case
        name = hashlib.sha1(cname(case, cube).encode()).hexdigest()[:16]
        cnf = os.path.join(WORK, name + ".cnf"); lrat = os.path.join(WORK, name + ".lrat")
        with open(cnf, "w", newline="\n") as f:
            f.write(f"p cnf {nv} {len(cl) + len(cube)}\n")
            f.write("".join(" ".join(map(str, c)) + " 0\n" for c in cl))
            f.write("".join(f"{l} 0\n" for l in cube))
        t0 = time.time()
        # CaDiCaL 3.0.1 aborts (rc 3) inside congruence closure on some leaves; CADICAL_OPTS defaults
        # to --congruence=false. Any other rc is retried without inprocessing, then fatal.
        for extra in ([], ["--inprocessing=false"]):
            rc, stdout = tools.solve(cnf, lrat, timeout=self.timeout, extra=extra)
            if rc in (10, 20, None):
                break
            print("solver rc", rc, "on", cname(case, cube), "options", extra, flush=True)
        else:
            raise SystemExit(f"solver failed rc={rc} on {cname(case, cube)}")
        dt = time.time() - t0
        base = dict(kind=kind, S=list(S), cube=cube)
        if rc == 20:
            h1 = sha(cnf)
            ok, out = tools.check_lrat(cnf, lrat)
            h2, size = sha(lrat), os.path.getsize(lrat)  # hash of the file the checker read
            self.record(**base, status="UNSAT", cake=ok, solve_s=round(dt, 2), cnf_sha=h1, lrat_sha=h2, lrat_bytes=size)
            os.remove(lrat); os.remove(cnf)
            if not ok:
                print("CAKE FAILED", cname(case, cube), out[-300:], flush=True)
            return "UNSAT"
        if rc == 10:
            model = os.path.join(RES, cname(case, cube)[:60] + "_" + name + ".model")
            open(model, "w").write(stdout)
            self.record(**base, status="SAT", solve_s=round(dt, 2), model=os.path.basename(model))
            os.remove(cnf)
            if os.path.exists(lrat): os.remove(lrat)
            print("SAT FOUND", cname(case, cube), flush=True)
            return "SAT"
        for p in (cnf, lrat):
            if os.path.exists(p): os.remove(p)
        free = [v for v in self.cand if v not in set(abs(l) for l in cube)][: self.k]
        self.record(**base, status="SPLIT", vars=free, solve_s=round(dt, 2))
        return ("SPLIT", free)

    def node(self, case, cube, pool):
        key = (case[0], tuple(case[1]), tuple(cube))
        nv, cl, _ = self.cache[case]
        if key in self.done:  # resume: re-derive split vars from the log is unnecessary; skip solved
            return []
        res = self.solve(case, cube, nv, cl)
        if isinstance(res, tuple):
            free = res[1]
            kids = []
            for signs in itertools.product((1, -1), repeat=len(free)):
                kids.append(cube + [s * v for s, v in zip(signs, free)])
            return kids
        return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int); ap.add_argument("t", type=int)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    R = Runner(a.n, a.t, a.timeout, a.k)
    cases = case_list(a.n)
    if a.only:
        cases = [c for c in cases if c[0] == a.only]
    # cheap cases first (ge6, h5, h4), h3 last
    order = {"ge6": 0, "h5": 1, "h4": 2, "h3": 3}
    cases.sort(key=lambda c: order[c[0]])
    # resume: find pending nodes = SPLIT children without a record
    split = {}
    for line in open(os.path.join(RES, "tree.jsonl")):
        r = json.loads(line)
        if r["status"] == "SPLIT":
            split[(r["kind"], tuple(r["S"]), tuple(r["cube"]))] = r["vars"]
    with ThreadPoolExecutor(a.workers) as pool:
        for case in cases:
            R.cache[case] = case_clauses(a.n, a.t, *case)
            frontier = [[]]
            while frontier:
                todo = []
                for cube in frontier:
                    key = (case[0], tuple(case[1]), tuple(cube))
                    if key in split:  # resumed internal node
                        free = split[key]
                        todo_k = [cube + [s * v for s, v in zip(signs, free)]
                                  for signs in itertools.product((1, -1), repeat=len(free))]
                        frontier_extra = todo_k
                        todo.append(("expand", frontier_extra))
                    else:
                        todo.append(("solve", cube))
                nxt = []
                futs = [pool.submit(R.node, case, c, pool) for kind, c in todo if kind == "solve"]
                for kind, c in todo:
                    if kind == "expand":
                        nxt += c
                for f in futs:
                    nxt += f.result()
                frontier = nxt
            print(time.strftime("%H:%M:%S"), "case done", case, flush=True)
            del R.cache[case]
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()

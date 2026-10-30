"""Referee cube-and-conquer over renc.py's CNF.  Adaptive binary tree: each node = base CNF +
unit cube; solve with CaDiCaL (LRAT) for <= T s; on UNSAT check the LRAT with $LRAT_CHECK (cake_lpr)
and delete it, keeping SHA-256; on timeout split on the next variable of a fixed list (both
signs).  Records go to <out>/tree.jsonl.  SAT anywhere -> stop and save the model.
Solver and checker come from the environment (tools.py); scratch files go to $TK_WORK.

usage: cnc.py n t outdir [T] [workers]
"""
import os, sys, json, time, hashlib, threading, itertools
from concurrent.futures import ThreadPoolExecutor
from renc import Enc
import tools

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.path.abspath(os.environ.get("TK_WORK", os.path.join(HERE, "..", "work", "referee")))


sha_file = tools.sha_file


def split_vars(enc):
    """orientation variables of triples not containing point 1 (fixed by the normal form),
    short spans first"""
    return [enc.X[tr] for tr in sorted((tr for tr in enc.X if tr[0] != 1), key=lambda tr: (tr[2] - tr[0], tr))]


def hull_roots(enc):
    """Partition of all assignments of H_3..H_{n-1} (a binary tree over H_3, H_4, ... in order,
    where a branch stops as soon as it has 3 positive literals): leaves = complete patterns with
    <= 2 positives, and prefixes ending in their 3rd positive."""
    hv = [enc.H[v] for v in range(3, enc.n)]
    out = []
    def rec(i, cube, pos):
        if pos == 3 or i == len(hv):
            out.append(cube); return
        rec(i + 1, cube + [hv[i]], pos + 1)
        rec(i + 1, cube + [-hv[i]], pos)
    rec(0, [], 0)
    return out


class Run:
    def __init__(self, n, t, out, T, workers):
        self.enc = Enc(n, t)
        self.base = self.enc.dimacs()
        self.base_body = self.base.split("\n", 1)[1]
        self.ncl = len(self.enc.cl)
        self.base_sha = hashlib.sha256(self.base.encode()).hexdigest()
        self.vars = split_vars(self.enc)
        self.out, self.T, self.K = out, T, 4
        os.makedirs(out, exist_ok=True); os.makedirs(SCR, exist_ok=True)
        self.log = open(os.path.join(out, "tree.jsonl"), "a")
        self.lock = threading.Lock()
        self.pool = ThreadPoolExecutor(workers)
        self.pending = 0
        self.done = threading.Event()
        self.sat = None
        self.record(kind="meta", n=n, t=t, nv=self.enc.nv, ncl=self.ncl, base_sha=self.base_sha,
                    split_vars=self.vars[:60], T=T, K=4)

    def record(self, **r):
        with self.lock:
            self.log.write(json.dumps(r) + "\n"); self.log.flush()

    def submit(self, cube):
        with self.lock:
            self.pending += 1
        self.pool.submit(self.node, cube)

    def node(self, cube):
        try:
            self._node(cube)
        except Exception as ex:
            self.record(kind="error", cube=cube, err=repr(ex))
            print("ERROR", cube, ex, flush=True)
        finally:
            with self.lock:
                self.pending -= 1
                if self.pending == 0:
                    self.done.set()

    def _node(self, cube):
        if self.sat is not None:
            return
        name = hashlib.sha1(json.dumps(cube).encode()).hexdigest()[:16]
        cnf = os.path.join(SCR, name + ".cnf"); lrat = os.path.join(SCR, name + ".lrat")
        with open(cnf, "w", newline="\n") as f:
            f.write(f"p cnf {self.enc.nv} {self.ncl + len(cube)}\n")
            f.write(self.base_body)
            f.write("".join(f"{l} 0\n" for l in cube))
        t0 = time.time()
        rc, stdout = tools.solve(cnf, lrat, timeout=self.T)  # CADICAL_OPTS default: --congruence=false
        dt = round(time.time() - t0, 1)
        if rc == 20:
            h1 = sha_file(cnf)
            ok, o = tools.check_lrat(cnf, lrat)
            h2, sz = sha_file(lrat), os.path.getsize(lrat)  # hash of the file the checker read
            self.record(kind="leaf", cube=cube, status="UNSAT", cake=ok, solve_s=dt, cnf_sha=h1,
                        lrat_sha=h2, lrat_bytes=sz)
            print("UNSAT", len(cube), dt, "cake", ok, sz >> 20, "MB", flush=True)
            os.remove(lrat); os.remove(cnf)
            return
        if os.path.exists(lrat):
            os.remove(lrat)
        if rc == 10:
            model = [int(x) for ln in stdout.splitlines() if ln.startswith("v") for x in ln[1:].split()]
            self.sat = model
            self.record(kind="leaf", cube=cube, status="SAT", solve_s=dt, model=model)
            print("SAT!", cube, flush=True)
            os.remove(cnf)
            return
        os.remove(cnf)
        if rc is not None:
            self.record(kind="error", cube=cube, rc=rc)
            print("solver rc", rc, cube, flush=True)
            return
        used = {abs(l) for l in cube}
        vs = [x for x in self.vars if x not in used][:self.K]
        self.record(kind="split", cube=cube, vars=vs, solve_s=dt)
        print("SPLIT", len(cube), "on", vs, flush=True)
        for signs in itertools.product((1, -1), repeat=len(vs)):
            self.submit(cube + [sg * v for sg, v in zip(signs, vs)])

    def go(self, roots=([],)):
        for c in roots:
            self.submit(list(c))
        self.done.wait()
        self.pool.shutdown()


if __name__ == "__main__":
    n, t, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    T = int(sys.argv[4]) if len(sys.argv) > 4 else 300
    w = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    R = Run(n, t, out, T, w)
    roots = hull_roots(R.enc)
    R.record(kind="roots", roots=roots)
    print(len(roots), "hull-pattern roots", flush=True)
    R.go(roots)
    print("FINISHED", "SAT" if R.sat else "all leaves done", flush=True)

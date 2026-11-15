"""Orbit-collapsed SAT instances for 'K is invariant under the group H and meets every
line of PG(3,9) in a number of points from the allowed set'.

  python orbsat.py CLASS K [--allowed 2,5] [--noamo] [--lrat PATH] [--timeout S] [--solver cadical|kissat]
  CLASS is one of the representative names of groups.build_all() or 'trivial'.

Variables z_O (one per point orbit O of H).  For every line orbit the multiset of orbit
multiplicities {O: |l n O|} is collected (identical multisets are encoded once):
    sum_O |l n O| z_O  in  allowed.
|K| is constrained to equal K (sum of orbit sizes).
"""
import argparse
import os
import re
import subprocess
import sys
import time
from collections import Counter

from cnfenc import CNF, sum_in_set, sum_equals
from pg import PG
import groups

CADICAL = r"C:\Users\user\src\build-cadical\cadical.exe"
KISSAT = r"C:\Users\user\src\build-kissat\rel\kissat.exe"


def orbits_of(n, perms):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in perms:
        for i in range(n):
            a, b = find(i), find(g[i])
            if a != b:
                parent[a] = b
    groups_ = {}
    for i in range(n):
        groups_.setdefault(find(i), []).append(i)
    orbs = sorted(groups_.values(), key=lambda o: o[0])
    oid = [0] * n
    for k, o in enumerate(orbs):
        for p in o:
            oid[p] = k
    return orbs, oid


def line_perm(P, g):
    return [P.lidx[tuple(sorted(g[p] for p in L))] for L in P.lines]


class OrbitInstance:
    def __init__(self, P, perms):
        self.P = P
        self.perms = perms
        self.orbs, self.oid = orbits_of(P.N, perms)
        lperms = [line_perm(P, g) for g in perms]
        self.lorbs, self.loid = orbits_of(len(P.lines), lperms)
        # signatures
        sig = {}
        for lo in self.lorbs:
            L = P.lines[lo[0]]
            s = tuple(sorted(Counter(self.oid[p] for p in L).items()))
            sig.setdefault(s, []).append(lo[0])
        self.signatures = sig       # signature -> list of representative lines

    def summary(self):
        return dict(points=self.P.N, point_orbits=len(self.orbs),
                    fixed_points=sum(1 for o in self.orbs if len(o) == 1),
                    lines=len(self.P.lines), line_orbits=len(self.lorbs),
                    distinct_constraints=len(self.signatures),
                    orbit_sizes=dict(Counter(len(o) for o in self.orbs)))

    def plane_data(self):
        """plane orbits: signature (point orbit multiplicities) for each plane orbit, and for each
        point orbit the multiset of line signatures through its representative"""
        if hasattr(self, "_pd"):
            return self._pd
        P = self.P
        plidx = {tuple(sorted(pl)): i for i, pl in enumerate(P.planes)}
        pperms = [[plidx[tuple(sorted(g[p] for p in pl))] for pl in P.planes] for g in self.perms]
        porbs, _ = orbits_of(len(P.planes), pperms)
        psig = {}
        for po in porbs:
            pl = P.planes[po[0]]
            s = tuple(sorted(Counter(self.oid[p] for p in pl).items()))
            psig.setdefault(s, []).append(po[0])
        # line signature of every line
        line_sig = {}
        for s, reps in self.signatures.items():
            pass
        sig_of_line = {}
        for lo in self.lorbs:
            L = P.lines[lo[0]]
            s = tuple(sorted(Counter(self.oid[p] for p in L).items()))
            for l in lo:
                sig_of_line[l] = s
        through = []
        for o in self.orbs:
            p = o[0]
            through.append(Counter(sig_of_line[l] for l in P.lines_through[p]))
        self._pd = (psig, through)
        return self._pd

    def build_cnf(self, k, allowed=(2, 5), amo=True, k_set=None, units=(), local=False, planes=False,
                  plane_sizes=(26, 35)):
        """units: iterable of (orbit_index, 0/1) forced values.
        local : add, for every point orbit, 'number of 5-secants through its representative = base + 30 z_O'
        planes: add, for every plane orbit, 'plane section size in plane_sizes'"""
        cnf = CNF()
        z = [cnf.new_var() for _ in self.orbs]
        self.z = z
        yvar = {}
        for s in sorted(self.signatures):
            lits = [z[o] for o, m in s]
            ws = [m for o, m in s]
            fin = sum_in_set(cnf, lits, ws, set(allowed), amo=amo)
            yvar[s] = fin.get(5)
        self.yvar = yvar
        if k is not None:
            sum_equals(cnf, z, [len(o) for o in self.orbs], k, amo=False)
        if local:
            psig, through = self.plane_data()
            base = (k - 182) // 3
            assert (k - 182) % 3 == 0
            for oi, cnt in enumerate(through):
                sigs = sorted(cnt)
                lits = [yvar[sg] for sg in sigs]
                ws = [cnt[sg] for sg in sigs]
                fin = sum_in_set(cnf, lits, ws, {base, base + 30}, amo=False)
                # z_O <-> (sum == base + 30)
                v_in, v_out = fin.get(base + 30), fin.get(base)
                if v_in is not None:
                    cnf.add(-z[oi], v_in)
                    cnf.add(z[oi], -v_in)
                if v_out is not None:
                    cnf.add(z[oi], v_out)
        if planes:
            psig, through = self.plane_data()
            for sg in sorted(psig):
                lits = [z[o] for o, m in sg]
                ws = [m for o, m in sg]
                sum_in_set(cnf, lits, ws, set(plane_sizes), amo=False)
        for o, v in units:
            cnf.add(z[o] if v else -z[o])
        return cnf


def check_K(P, K, allowed=(2, 5)):
    S = set(K)
    bad = [i for i, L in enumerate(P.lines) if sum(1 for p in L if p in S) not in allowed]
    return not bad, bad


def parse_cadical(out):
    status = None
    model = []
    for line in out.splitlines():
        if line.startswith("s "):
            status = line[2:].strip()
        elif line.startswith("v "):
            model.extend(int(x) for x in line[2:].split())
    return status, [m for m in model if m != 0]


WSL_CADICAL = "/home/user/ext/cadical/build/cadical"
WSL_CAKE = "/home/user/lrchk/bin/cake_lpr"


def wsl_path(p):
    """Windows path -> WSL path"""
    p = os.path.abspath(p).replace("\\", "/")
    if len(p) > 1 and p[1] == ":":
        p = "/mnt/" + p[0].lower() + p[2:]
    return p


_WSL_BAD = {"flag": False}


def run_solver(cnf_path, solver="cadical", timeout=600, lrat=None, extra=(), retries=1):
    """CaDiCaL 3.0.1.  Preferred: the WSL build.  Fallback when WSL is unavailable (E_UNEXPECTED):
    the Windows build with --congruence=false (its default configuration aborts on some CNFs).
    returns (status, model, seconds, info); status in SATISFIABLE / UNSATISFIABLE / TIMEOUT / ERROR"""
    assert solver == "cadical"
    t = time.time()
    use_win = _WSL_BAD["flag"]
    for attempt in range(retries + 2):
        if not use_win:
            args = ["wsl", "-d", "Ubuntu-24.04", "-e", WSL_CADICAL, "-t", str(int(timeout))]
            if lrat:
                args += ["--lrat=true", "--binary=false", "--congruence=false"]
            args += list(extra) + [wsl_path(cnf_path)]
            if lrat:
                args += [wsl_path(lrat)]
        else:
            args = [CADICAL, "--congruence=false"]
            if lrat:
                args += ["--lrat=true", "--binary=false"]
            args += list(extra) + [cnf_path]
            if lrat:
                args += [lrat]
        try:
            r = subprocess.run(args, capture_output=True, text=True, timeout=timeout + 60)
            out = r.stdout.replace(chr(0), "")
            err = r.stderr.replace(chr(0), "")
        except subprocess.TimeoutExpired:
            return "TIMEOUT", [], time.time() - t, {"backend": "win" if use_win else "wsl"}
        status, model = parse_cadical(out)
        info = {"backend": "win" if use_win else "wsl"}
        m = re.search(r"^c conflicts:\s+(\d+)", out, re.M)
        if m:
            info["conflicts"] = int(m.group(1))
        elapsed = time.time() - t
        if status is None:
            if r.returncode == 10:
                status = "SATISFIABLE"
            elif r.returncode == 20:
                status = "UNSATISFIABLE"
            elif r.returncode == 0 and elapsed >= 0.9 * timeout:
                status = "TIMEOUT"
            elif "E_UNEXPECTED" in err or "E_UNEXPECTED" in out:
                _WSL_BAD["flag"] = True
                use_win = True
                continue
            else:
                info["rc"] = r.returncode
                info["err"] = (err or out)[-200:].strip()
                if attempt < retries + 1:
                    time.sleep(3)
                    continue
                status = "ERROR"
        return status, model, time.time() - t, info
    return "ERROR", [], time.time() - t, info


def get_perms(P, class_name, F=None):
    F_, P_, gens, reps = groups.build_all()
    if class_name == "trivial":
        return []
    assert class_name in reps, (class_name, list(reps))
    return [P.perm_of_matrix(*reps[class_name][:1], frob=reps[class_name][1])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cls")
    ap.add_argument("k", type=int)
    ap.add_argument("--allowed", default="2,5")
    ap.add_argument("--noamo", action="store_true")
    ap.add_argument("--lrat")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--solver", default="cadical")
    ap.add_argument("--cnf")
    args = ap.parse_args()
    allowed = tuple(int(x) for x in args.allowed.split(","))
    P = PG(4, 9)
    perms = get_perms(P, args.cls)
    inst = OrbitInstance(P, perms)
    print(args.cls, inst.summary())
    cnf = inst.build_cnf(args.k, allowed, amo=not args.noamo)
    path = args.cnf or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cnf", "%s_k%d.cnf" % (args.cls.replace("+", "p").replace("(", "").replace(")", "").replace(",", "_"), args.k))
    cnf.to_dimacs(path, comments=["class %s k %d allowed %s" % (args.cls, args.k, allowed)])
    print("cnf: %d vars %d clauses -> %s" % (cnf.nvars, len(cnf.clauses), path))
    status, model, t, info = run_solver(path, args.solver, args.timeout, args.lrat)
    print("status", status, "time %.1fs" % t, info)
    if status == "SATISFIABLE":
        true_orbits = [i for i, v in enumerate(inst.z) if v in set(m for m in model if m > 0)]
        K = sorted(p for i in true_orbits for p in inst.orbs[i])
        ok, bad = check_K(P, K, allowed)
        print("model: |K| = %d, exact line check %s" % (len(K), ok))
        print("K =", K)


if __name__ == "__main__":
    main()

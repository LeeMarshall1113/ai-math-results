#!/usr/bin/env python3
"""Independent replay of the group-class exclusions (stdlib only; imports nothing from ../code).

  python replay.py CLASS K [--cake /path/to/cake_lpr] [--units "orbit:val,orbit:val"]

Rebuilds from scratch
  * GF(9) as pairs a+b*i (i^2 = -1, element code a+3b),
  * PG(3,9): points = normalised vectors (lexicographic order of codes), lines = row-reduced
    2x4 matrices (a DIFFERENT enumeration from the main code, which uses pair covering),
  * the permutation of the class representative given in ../certs/classes.json,
  * point orbits (ordered by smallest member), line orbits, and for every line orbit the multiset
    {point orbit : multiplicity},
  * the CNF  'sum_O mult_O z_O in {2,5} for every distinct multiset, and sum_O |O| z_O = K'
    with the same exact encoder (copied here; its semantics are tested exhaustively below),
and checks that the DIMACS file produced has the same SHA-256 as ../certs/<tag>.cnf recorded in
../certs/manifest.jsonl.  If --cake is given, it also runs cake_lpr on (cnf, lrat).
"""
import hashlib
import itertools
import json
import os
import subprocess
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")

# ------------------------------------------------------------------ GF(9) = GF(3)[i]/(i^2+1)
def e_add(x, y):
    return ((x % 3 + y % 3) % 3) + 3 * (((x // 3) + (y // 3)) % 3)


def e_mul(x, y):
    a, b, c, d = x % 3, x // 3, y % 3, y // 3
    return ((a * c - b * d) % 3) + 3 * ((a * d + b * c) % 3)


def e_frob(x):            # x -> x^3 : a+bi -> a-bi
    return (x % 3) + 3 * ((-(x // 3)) % 3)


E_INV = {x: next(y for y in range(1, 9) if e_mul(x, y) == 1) for x in range(1, 9)}


def decode(pair):
    return pair[0] + 3 * pair[1]


def normalize(v):
    for c in v:
        if c:
            iv = E_INV[c]
            return tuple(e_mul(iv, x) for x in v)
    raise ValueError


def build_points():
    pts = []
    for v in itertools.product(range(9), repeat=4):
        nz = next((c for c in v if c), None)
        if nz == 1:
            pts.append(v)
    return pts


def build_lines_rref(points):
    """all 2-dim subspaces via reduced row echelon forms; each line = sorted tuple of point indices"""
    pidx = {p: i for i, p in enumerate(points)}
    lines = []
    for c1, c2 in itertools.combinations(range(4), 2):
        # row1: 1 at c1, 0 at columns < c1 and at c2, free at columns > c1 except c2
        free1 = [j for j in range(c1 + 1, 4) if j != c2]
        free2 = [j for j in range(c2 + 1, 4)]
        for f1 in itertools.product(range(9), repeat=len(free1)):
            for f2 in itertools.product(range(9), repeat=len(free2)):
                r1 = [0] * 4
                r2 = [0] * 4
                r1[c1] = 1
                r2[c2] = 1
                for j, x in zip(free1, f1):
                    r1[j] = x
                for j, x in zip(free2, f2):
                    r2[j] = x
                members = set()
                for a, b in [(1, t) for t in range(9)] + [(0, 1)]:
                    v = tuple(e_add(e_mul(a, x), e_mul(b, y)) for x, y in zip(r1, r2))
                    members.add(pidx[normalize(v)])
                assert len(members) == 10
                lines.append(tuple(sorted(members)))
    assert len(lines) == len(set(lines)) == 7462
    return lines


def point_perm(points, pidx, M, frob):
    perm = []
    for v in points:
        w = tuple(e_frob(c) for c in v) if frob else v
        out = []
        for row in M:
            s = 0
            for a, x in zip(row, w):
                s = e_add(s, e_mul(a, x))
            out.append(s)
        perm.append(pidx[normalize(tuple(out))])
    assert sorted(perm) == list(range(len(points)))
    return perm


def orbits_of(n, perms):
    seen = [-1] * n
    orbs = []
    for s in range(n):
        if seen[s] >= 0:
            continue
        stack, orb = [s], []
        seen[s] = len(orbs)
        while stack:
            x = stack.pop()
            orb.append(x)
            for g in perms:
                y = g[x]
                if seen[y] < 0:
                    seen[y] = len(orbs)
                    stack.append(y)
        orbs.append(sorted(orb))
    return orbs, seen


# ------------------------------------------------------------------ CNF encoder (same spec as code/cnfenc.py)
class CNF:
    def __init__(self):
        self.nvars = 0
        self.clauses = []

    def new_var(self):
        self.nvars += 1
        return self.nvars

    def add(self, *lits):
        self.clauses.append(tuple(lits))

    def dimacs(self, comments=()):
        out = []
        for c in comments:
            out.append("c %s\n" % c)
        out.append("p cnf %d %d\n" % (self.nvars, len(self.clauses)))
        for cl in self.clauses:
            out.append(" ".join(map(str, cl)) + " 0\n")
        return "".join(out)


def sum_in_set(cnf, lits, weights, allowed, amo=True):
    n = len(lits)
    allowed = set(allowed)
    top = max(allowed)
    suffix = [None] * (n + 1)
    suffix[n] = {0}
    for i in range(n - 1, -1, -1):
        w = weights[i]
        suffix[i] = suffix[i + 1] | {s + w for s in suffix[i + 1] if s + w <= top}
    layer = [None] * (n + 1)
    layer[0] = {0} if any((a - 0) in suffix[0] for a in allowed) else set()
    for i in range(n):
        w = weights[i]
        nxt = set()
        for s in layer[i]:
            for t in (s, s + w):
                if t <= top and any((a - t) in suffix[i + 1] for a in allowed):
                    nxt.add(t)
        layer[i + 1] = nxt
    layer[n] = {s for s in layer[n] if s in allowed}
    var = [dict() for _ in range(n + 1)]
    for i in range(n + 1):
        for s in sorted(layer[i]):
            var[i][s] = cnf.new_var()
    if not layer[0]:
        cnf.add()
        return {}
    cnf.add(var[0][0])
    for i in range(n):
        z, w = lits[i], weights[i]
        for s, v in var[i].items():
            t0 = var[i + 1].get(s)
            t1 = var[i + 1].get(s + w)
            cnf.add(*([-v, z] + ([t0] if t0 else [])))
            cnf.add(*([-v, -z] + ([t1] if t1 else [])))
        for t, v in var[i + 1].items():
            p0 = var[i].get(t)
            p1 = var[i].get(t - w)
            cnf.add(*([-v, z] + ([p0] if p0 else [])))
            cnf.add(*([-v, -z] + ([p1] if p1 else [])))
        if amo and i + 1 < n + 1:
            vs = list(var[i + 1].values())
            for a in range(len(vs)):
                for b in range(a + 1, len(vs)):
                    cnf.add(-vs[a], -vs[b])
    return dict(var[n])


def eval_ext(cnf, ins, assignment):
    """is there an extension of the input assignment satisfying all clauses?  (tiny DPLL)"""
    clauses = cnf.clauses

    def solve(val):
        val = dict(val)
        changed = True
        while changed:
            changed = False
            for cl in clauses:
                un, sat = [], False
                for l in cl:
                    v = abs(l)
                    if v in val:
                        if (val[v] and l > 0) or ((not val[v]) and l < 0):
                            sat = True
                            break
                    else:
                        un.append(l)
                if sat:
                    continue
                if not un:
                    return False
                if len(un) == 1:
                    val[abs(un[0])] = un[0] > 0
                    changed = True
        for cl in clauses:
            for l in cl:
                if abs(l) not in val:
                    for b in (True, False):
                        v2 = dict(val)
                        v2[abs(l)] = b
                        if solve(v2):
                            return True
                    return False
        return True

    return solve(dict(zip(ins, assignment)))


def test_shape(weights, allowed=(2, 5)):
    """exhaustive semantic test of the encoder on one weight vector"""
    cnf = CNF()
    ins = [cnf.new_var() for _ in weights]
    sum_in_set(cnf, ins, weights, set(allowed), amo=True)
    for a in itertools.product((False, True), repeat=len(weights)):
        want = sum(w for w, b in zip(weights, a) if b) in allowed
        if eval_ext(cnf, ins, a) != want:
            return False
    return True


# ------------------------------------------------------------------ main
def build_cnf_for_class(name, k, points, lines, pidx, cls, units=(), allowed=(2, 5), header=None):
    M = [[decode(e) for e in row] for row in cls["matrix"]]
    perm = point_perm(points, pidx, M, cls["frob"])
    orbs, oid = orbits_of(len(points), [perm])
    lidx = {L: i for i, L in enumerate(lines)}
    lperm = [lidx[tuple(sorted(perm[p] for p in L))] for L in lines]
    lorbs, _ = orbits_of(len(lines), [lperm])
    sigs = {}
    for lo in lorbs:
        L = lines[lo[0]]
        s = tuple(sorted(Counter(oid[p] for p in L).items()))
        sigs.setdefault(s, []).append(lo[0])
    cnf = CNF()
    z = [cnf.new_var() for _ in orbs]
    shapes = set()
    for s in sorted(sigs):
        lits = [z[o] for o, m in s]
        ws = [m for o, m in s]
        shapes.add(tuple(sorted(ws)))
        sum_in_set(cnf, lits, ws, set(allowed), amo=True)
    cnf.clauses  # total size
    sum_in_set(cnf, z, [len(o) for o in orbs], {k}, amo=False)
    for o, v in units:
        cnf.add(z[o] if v else -z[o])
    info = dict(point_orbits=len(orbs), line_orbits=len(lorbs), distinct=len(sigs),
                orbit_sizes=dict(Counter(len(o) for o in orbs)), shapes=sorted(shapes))
    comment = header or ("class %s k %d allowed 2,5 amo=True" % (name, k))
    return cnf, comment, info


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("cls")
    ap.add_argument("k", type=int)
    ap.add_argument("--cake")
    ap.add_argument("--no-shape-test", action="store_true")
    a = ap.parse_args()
    classes = json.load(open(os.path.join(CERTS, "classes.json")))
    points = build_points()
    pidx = {p: i for i, p in enumerate(points)}
    lines = build_lines_rref(points)
    print("PG(3,9): %d points, %d lines (row-reduced enumeration)" % (len(points), len(lines)))
    cnf, comment, info = build_cnf_for_class(a.cls, a.k, points, lines, pidx, classes[a.cls])
    print("orbits:", info["point_orbits"], "line orbits:", info["line_orbits"], "distinct constraints:", info["distinct"],
          "orbit sizes:", info["orbit_sizes"])
    if not a.no_shape_test:
        bad = [s for s in info["shapes"] if not test_shape(list(s))]
        print("encoder exhaustive test on %d constraint shapes: %s" % (len(info["shapes"]), "OK" if not bad else "FAIL %s" % bad))
        assert not bad
    text = cnf.dimacs([comment])
    h = hashlib.sha256(text.encode()).hexdigest()
    tag = "%s_k%d" % (a.cls.replace("(", "").replace(")", "").replace(",", "_").replace("+", "p"), a.k)
    rec = None
    for line in open(os.path.join(CERTS, "manifest.jsonl")):
        r = json.loads(line)
        if r["tag"] == tag:
            rec = r
    assert rec is not None, "no manifest record for " + tag
    print("regenerated CNF: %d vars %d clauses sha256 %s" % (cnf.nvars, len(cnf.clauses), h))
    print("manifest sha256 :", rec["cnf_sha256"], "->", "MATCH" if h == rec["cnf_sha256"] else "MISMATCH")
    assert h == rec["cnf_sha256"]
    lr = os.path.join(CERTS, tag + ".lrat")
    if os.path.exists(lr):
        lh = hashlib.sha256(open(lr, "rb").read()).hexdigest()
        print("LRAT sha256 matches manifest:", lh == rec.get("lrat_sha256"))
    if a.cake:
        cn = os.path.join(CERTS, tag + ".cnf")
        r = subprocess.run([a.cake, "--CML_HEAP_SIZE=1500", "--CML_STACK_SIZE=500", cn, lr], capture_output=True, text=True)
        print("cake_lpr:", (r.stdout + r.stderr).strip().splitlines()[-1:])


if __name__ == "__main__":
    main()

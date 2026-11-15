#!/usr/bin/env python3
"""Stdlib replay of the star cover: rebuild EVERY CNF of the proof chain from scratch and compare the SHA-256 with certs/manifest.jsonl.

  python rebuild_all.py                    everything (about 30-60 minutes, 4 GB disk for temporary bodies, < 1 GB RAM)
  python rebuild_all.py --only planar,base,cubes,children,cover1,cover2     a subset (comma list)
  python rebuild_all.py --cases A0,C3      restrict the star cases (base, cubes, children, cover1)
  python rebuild_all.py --selftest         first run the exhaustive semantic test of the encoder (small shapes)

Imports nothing from ../code (only replay.py and star_setup.py of this folder, stdlib).  Independent re-implementation of: GF(9) (a+bi, i^2 = -1), PG(3,9) with its
lines in pair-covering order (= sorted by the two smallest points), the star set-up of NOTES 5.1 (P0, Sigma representative, first non-collinear triple, units), the exact
DP/BDD encoder (replay.sum_in_set, exhaustively tested by --selftest), and for every CNF the exact byte stream, so that a hash equality means the same file.
Inputs that are DATA and not rebuilt here (checked elsewhere): results/planar_classes.json (class representatives of the planar (2,5)-sets, certified by the
planar-completeness CNFs below + nauty), certs/cubes_*.json (the homothety-orbit representatives of the level-1 sections), certs/split_*.json (level-2 children),
results/planar_k*_R*.json (the enumerated frame solutions, certified complete by the planar-completeness CNFs).

Items and manifest tags:
  planar : planar_complete_k{26,35}_R{0,1}   (frame enumeration of planar (2,5)-sets + blocking clauses of all solutions)
  base   : full_{A,C}_c<i>.cnf   (10 base CNFs; hashes in results/modal_base_sha256.json)
  cubes  : <case><i>_c<jj>      (174 cube CNFs = base + 92 unit clauses)       tags of certs/manifest.jsonl
  children: <cube>_p<plane>_k<n> (level-2 children = cube CNF + 80 unit clauses)
  cover1 : cover1_<case><i>      (level-1 completeness: plane CNF + blocking clauses of the 8 homothety images of the representatives)
  cover2 : cover2_<cube>         (level-2 completeness: plane CNF with the cube section + blocking clauses of the children)
"""
import hashlib
import itertools
import json
import os
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import replay as R            # noqa: E402  (GF(9), points, lines, exact DP encoder)
import star_setup as SS       # noqa: E402  (points of the quotient plane, noncollinear_triple)

CERTS = os.path.join(HERE, "..", "certs")
RES = os.path.join(HERE, "..", "results")


# ------------------------------------------------------------------ small helpers
class StreamCNF:
    """clauses are streamed to a body file (no clause list in memory); same interface as replay.CNF for the encoder"""

    def __init__(self, path):
        self.fh = open(path, "w", newline="\n")
        self.nvars = 0
        self.nclauses = 0

    def new_var(self):
        self.nvars += 1
        return self.nvars

    def add(self, *lits):
        self.fh.write((" ".join(map(str, lits)) + " 0\n") if lits else "0\n")
        self.nclauses += 1

    def close(self):
        self.fh.close()


def sum_equals(cnf, lits, weights, value, amo=False):
    return R.sum_in_set(cnf, lits, weights, {value}, amo=amo)


def sha_text_body(header, body_path, tail=""):
    h = hashlib.sha256(header.encode())
    with open(body_path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 24), b""):
            h.update(ch)
    h.update(tail.encode())
    return h.hexdigest()


def header_text(comments, nv, nc):
    return "".join("c %s\n" % c for c in comments) + "p cnf %d %d\n" % (nv, nc)


def lit(p, v):
    return (p + 1) if v else -(p + 1)


class Geometry:
    def __init__(self):
        t = time.time()
        self.points4 = R.build_points()
        self.pidx4 = {p: i for i, p in enumerate(self.points4)}
        self.lines4 = sorted(R.build_lines_rref(self.points4))              # pair-covering order = sorted by the two smallest points
        for a, b in zip(self.lines4, self.lines4[1:]):
            assert (a[0], a[1]) < (b[0], b[1])
        self.points3 = SS.build_points3()
        self.pidx3 = {p: i for i, p in enumerate(self.points3)}
        self.lines3 = sorted(SS.lines3(self.points3))                       # quotient lines in the order of the main code's P2.lines
        assert len(self.lines3) == 91
        self.P0 = self.pidx4[(0, 0, 0, 1)]
        self.proj = {q: self.pidx3[R.normalize(v[:3])] for q, v in enumerate(self.points4) if q != self.P0}
        self.zero = {d: self.pidx4[R.normalize(self.points3[d] + (0,))] for d in range(91)}
        self.lines_thru = [[] for _ in range(len(self.points4))]
        for j, L in enumerate(self.lines4):
            for p in L:
                self.lines_thru[p].append(j)
        # lines through P0 in ascending index, with their direction (= direction of the smallest non-P0 point)
        self.line_of = {}
        for j in self.lines_thru[self.P0]:
            q = next(p for p in self.lines4[j] if p != self.P0)
            self.line_of[self.proj[q]] = j
        assert len(self.line_of) == 91
        self.line_pts = {}
        for L2, Lq in enumerate(self.lines3):
            s = set(Lq)
            self.line_pts[L2] = sorted([self.P0] + [q for q, d in self.proj.items() if d in s])
        self.homs = [[self.pidx4[R.normalize((v[0], v[1], v[2], R.e_mul(lam, v[3])))] for v in self.points4] for lam in range(1, 9)]
        print("geometry built (%.1fs): %d points, %d lines" % (time.time() - t, len(self.points4), len(self.lines4)), flush=True)


class Case:
    """star set-up of one top case (NOTES 5.1), re-derived"""

    def __init__(self, G, case, ci, classes):
        self.G, self.case, self.ci = G, case, ci
        cl = [c for c in classes if c["k"] == (35 if case == "C" else 26)][ci]
        self.aut = cl["aut_order"]
        self.Sigma = set(cl["rep"])
        self.k = 287 if case == "A" else 260
        self.tri = SS.noncollinear_triple(self.Sigma, G.lines3)
        self.pt_units = {G.P0: 1}
        for v in self.tri:
            self.pt_units[G.zero[v]] = 1
        self.line_units = [(j, 0 if d in self.Sigma else 1) for d, j in G.line_of.items()]
        self.line_y = dict(self.line_units)
        self.comment = "trivial group case %s class %d (aut %d) k %d" % (case, ci, self.aut, self.k)

    def plane_size(self, L2):
        return 26 if len(self.Sigma & set(self.G.lines3[L2])) == 5 else 35

    # ---------------------------------------------------------------- base CNF (body to a temp file)
    def build_base(self, body_path):
        G = self.G
        cnf = StreamCNF(body_path)
        x = [cnf.new_var() for _ in range(len(G.points4))]
        y = []
        for L in G.lines4:
            fin = R.sum_in_set(cnf, [x[p] for p in L], [1] * len(L), {2, 5}, amo=True)
            y.append(fin[5])
        sum_equals(cnf, x, [1] * len(G.points4), self.k, amo=False)
        for p, v in self.pt_units.items():
            cnf.add(x[p] if v else -x[p])
        for j, v in self.line_units:
            cnf.add(y[j] if v else -y[j])
        cnf.close()
        self.nv, self.nc = cnf.nvars, cnf.nclauses

    # ---------------------------------------------------------------- planar CNF of a plane through P0 (Setup.plane_cnf)
    def plane_cnf(self, L2, fixed):
        G = self.G
        pts = G.line_pts[L2]
        ps = set(pts)
        cnf = R.CNF()
        x = {p: cnf.new_var() for p in pts}
        lines = [j for j, L in enumerate(G.lines4) if all(p in ps for p in L)]
        assert len(lines) == 91
        y = {}
        for j in lines:
            fin = R.sum_in_set(cnf, [x[p] for p in G.lines4[j]], [1] * 10, {2, 5}, amo=True)
            y[j] = fin[5]
        m = self.plane_size(L2)
        sum_equals(cnf, [x[p] for p in pts], [1] * 91, m, amo=False)
        d = (m - 20) // 3
        for p in pts:
            thr = [y[j] for j in lines if p in G.lines4[j]]
            fin = R.sum_in_set(cnf, thr, [1] * len(thr), {d, d + 3}, amo=False)
            cnf.add(-x[p], fin[d + 3])
            cnf.add(x[p], -fin[d + 3])
        for j in lines:
            if j in self.line_y:
                cnf.add(y[j] if self.line_y[j] else -y[j])
        for p, v in {**self.pt_units, **fixed}.items():
            if p in x:
                cnf.add(x[p] if v else -x[p])
        return cnf, x, pts


def load_manifest():
    man = {}
    for l in open(os.path.join(CERTS, "manifest.jsonl")):
        try:
            r = json.loads(l)
        except Exception:
            continue
        if "cnf_sha256" in r and "tag" in r:
            man.setdefault(r["tag"], set()).add(r["cnf_sha256"])
    return man


class Report:
    def __init__(self, man):
        self.man = man
        self.ok = 0
        self.bad = []
        self.missing = []
        self.checked = set()

    def check(self, tag, h, what=""):
        self.checked.add(tag)
        if tag not in self.man:
            self.missing.append(tag)
            return
        if self.man[tag] == {h}:
            self.ok += 1
        else:
            self.bad.append((tag, h, sorted(self.man[tag])))
            print("MISMATCH %s %s: rebuilt %s, manifest %s" % (what, tag, h, sorted(self.man[tag])), flush=True)


# ------------------------------------------------------------------ planar completeness CNFs (frame enumeration of planar (2,5)-sets)
def planar_items():
    t_start = time.time()
    # PG(2,9) as 3-vectors, lines pair-covering order, l0 = {z = 0}
    pts = SS.build_points3()
    pidx = {p: i for i, p in enumerate(pts)}
    lines = sorted(SS.lines3(pts))
    thru = [[] for _ in pts]
    for j, L in enumerate(lines):
        for p in L:
            thru[p].append(j)
    l0 = [i for i, v in enumerate(pts) if v[2] == 0]
    assert len(l0) == 10
    # PGammaL(2,9) on the 10 points of PG(1,9): closure of generators, orbits on 5-subsets (reps = first member in lexicographic order)
    p1 = [v for v in itertools.product(range(9), repeat=2) if next((c for c in v if c), None) == 1]
    p1idx = {p: i for i, p in enumerate(p1)}

    def perm2(M, frob):
        out = []
        for v in p1:
            w = tuple(R.e_frob(c) for c in v) if frob else v
            u0 = R.e_add(R.e_mul(M[0][0], w[0]), R.e_mul(M[0][1], w[1]))
            u1 = R.e_add(R.e_mul(M[1][0], w[0]), R.e_mul(M[1][1], w[1]))
            out.append(p1idx[R.normalize((u0, u1))])
        return tuple(out)
    prim = 4                      # code of 1+i, a primitive element of GF(9)
    pw, e = [], 1
    for _ in range(8):
        e = R.e_mul(e, prim)
        pw.append(e)
    assert len(set(pw)) == 8 and pw[-1] == 1, pw
    gens = [perm2([[1, 1], [0, 1]], 0), perm2([[prim, 0], [0, 1]], 0), perm2([[0, 1], [1, 0]], 0), perm2([[1, 0], [0, 1]], 1)]
    ident = tuple(range(10))
    seen, frontier = {ident}, [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                c = tuple(h[g[i]] for i in range(10))
                if c not in seen:
                    seen.add(c)
                    nxt.append(c)
        frontier = nxt
    assert len(seen) == 1440, len(seen)
    orbit_of, reps = {}, []
    for s in itertools.combinations(range(10), 5):
        if s in orbit_of:
            continue
        orb = set(tuple(sorted(g[i] for i in s)) for g in seen)
        for t in orb:
            orbit_of[t] = len(reps)
        reps.append((s, len(orb)))
    assert [r[1] for r in reps] == [180, 72], reps
    origin = pidx[(0, 0, 1)]
    results = []
    for ri, (R1, osz) in enumerate(reps):
        Rpts = set(pidx[R.normalize((p1[i][0], p1[i][1], 0))] for i in R1)
        # stabiliser H_R: diag(A,1) frob^e preserving R  (A in GL(2,9)), in the enumeration order of the main code
        H = []
        for a, b, c, d in itertools.product(range(9), repeat=4):
            det = R.e_add(R.e_mul(a, d), R.e_mul(R.e_mul(b, c), NEG1))
            if det == 0:
                continue
            for e in (0, 1):
                perm = []
                for v in pts:
                    w = tuple(R.e_frob(cc) for cc in v) if e else v
                    o0 = R.e_add(R.e_mul(a, w[0]), R.e_mul(b, w[1]))
                    o1 = R.e_add(R.e_mul(c, w[0]), R.e_mul(d, w[1]))
                    perm.append(pidx[R.normalize((o0, o1, w[2]))])
                if set(perm[p] for p in Rpts) == Rpts:
                    H.append(perm)
        results.append((ri, osz, Rpts, H))
    print("planar frame data built (%.1fs): |H_R| = %s" % (time.time() - t_start, [len(h[3]) for h in results]), flush=True)
    return pts, pidx, lines, thru, l0, origin, results


NEG1 = 2          # -1 in GF(3) is the element 2 (code a+3b with a = 2, b = 0)


def planar_cnf(k, pts, pidx, lines, thru, l0, origin, Rpts, H, sols):
    cnf = R.CNF()
    N = len(pts)
    x = [cnf.new_var() for _ in range(N)]
    y = []
    for L in lines:
        fin = R.sum_in_set(cnf, [x[p] for p in L], [1] * len(L), {2, 5}, amo=True)
        y.append(fin[5])
    sum_equals(cnf, x, [1] * N, k, amo=False)
    d = (k - 20) // 3
    for p in range(N):
        t = [y[j] for j in thru[p]]
        fin = R.sum_in_set(cnf, t, [1] * len(t), {d, d + 3}, amo=False)
        cnf.add(-x[p], fin[d + 3])
        cnf.add(x[p], -fin[d + 3])
    l0s = set(l0)
    for p in l0:
        cnf.add(x[p] if p in Rpts else -x[p])
    cnf.add(x[origin])
    order = [p for p in range(N) if p not in l0s]
    ident = list(range(N))
    for h in H:
        if h != ident:
            prev = None
            for p in order:
                a, b = x[p], x[h[p]]
                if a == b:
                    continue
                if prev is None:
                    cnf.add(-a, b)
                else:
                    cnf.add(-prev, -a, b)
                e = cnf.new_var()
                if prev is None:
                    cnf.add(a, b, e)
                    cnf.add(-a, -b, e)
                else:
                    cnf.add(-prev, a, b, e)
                    cnf.add(-prev, -a, -b, e)
                prev = e
    for K in sols:
        Ks = set(K)
        cnf.add(*[(-x[p] if p in Ks else x[p]) for p in order])
    return cnf


def run_planar(rep):
    pts, pidx, lines, thru, l0, origin, results = planar_items()
    for ri, osz, Rpts, H in results:
        for k in (26, 35):
            d = json.load(open(os.path.join(RES, "planar_k%d_R%d.json" % (k, ri))))
            assert sorted(Rpts) == d["R"], "frame R differs from the stored one"
            cnf = planar_cnf(k, pts, pidx, lines, thru, l0, origin, Rpts, H, d["solutions"])
            txt = cnf.dimacs(["planar completeness k %d R#%d blocked %d solutions" % (k, ri, len(d["solutions"]))])
            rep.check("planar_complete_k%d_R%d" % (k, ri), hashlib.sha256(txt.encode()).hexdigest(), "planar")
            print("planar_complete_k%d_R%d rebuilt (%d vars, %d clauses)" % (k, ri, cnf.nvars, len(cnf.clauses)), flush=True)


# ------------------------------------------------------------------ star cover
def cube_units(C, cube_json, i):
    """unit literals of cube i: pt_units first (their order), then the section points ascending that are not in pt_units"""
    bu = cube_json["base_units"]
    assert [int(p) for p in bu] == list(C.pt_units)
    cu = cube_json["cubes"][i]
    assert cu["index"] == i
    units = [lit(int(p), v) for p, v in bu.items()]
    sec = cu["section"]
    for p, v in sec:
        if str(p) not in bu:
            units.append(lit(p, v))
    return units, {p: v for p, v in sec}


def run_star(rep, G, classes, cases, only, tmp):
    split_files = {}
    for fn in os.listdir(CERTS):
        if fn.startswith("split_") and fn.endswith(".json"):
            split_files[fn[6:-5]] = os.path.join(CERTS, fn)
    base_hashes = json.load(open(os.path.join(RES, "modal_base_sha256.json")))
    for case, ci in cases:
        C = Case(G, case, ci, classes)
        name = "full_%s_c%d" % (case, ci)
        t = time.time()
        body = os.path.join(tmp, name + ".body")
        C.build_base(body)
        h_base = sha_text_body(header_text([C.comment], C.nv, C.nc), body)
        if "base" in only:
            if base_hashes.get(name) == h_base:
                rep.ok += 1
            else:
                rep.bad.append((name, h_base, base_hashes.get(name)))
                print("MISMATCH base %s: rebuilt %s, recorded %s" % (name, h_base, base_hashes.get(name)), flush=True)
        print("%s base CNF rebuilt: %d vars, %d clauses, sha256 %s (%.0fs)" % (name, C.nv, C.nc, h_base[:16], time.time() - t), flush=True)
        cj = json.load(open(os.path.join(CERTS, "cubes_%s%d.json" % (case, ci))))
        # header of cubes: 92 units; children: 92 + 80
        by_header = {}

        def cnf_hash(units):
            n = len(units)
            if n not in by_header:
                h = hashlib.sha256(header_text([C.comment], C.nv, C.nc + n).encode())
                with open(body, "rb") as fh:
                    for ch in iter(lambda: fh.read(1 << 24), b""):
                        h.update(ch)
                by_header[n] = h
            h = by_header[n].copy()
            h.update("".join("%d 0\n" % u for u in units).encode())
            return h.hexdigest()
        # cubes
        n_cubes = 0
        for i, cu in enumerate(cj["cubes"]):
            units, secd = cube_units(C, cj, i)
            tag = "%s%d_c%02d" % (case, ci, i)
            n_cubes += 1
            if "cubes" in only:
                rep.check(tag, cnf_hash(units), "cube")
            # children of a split cube
            for sid, sp_path in split_files.items():
                if sid == tag and "children" in only:
                    sp = json.load(open(sp_path))
                    for kk, ch in enumerate(sp["children"]):
                        rep.check("%s_p%d_k%03d" % (tag, sp["plane"], kk), cnf_hash(units + ch), "child")
        # level 1 completeness
        if "cover1" in only:
            L = cj["plane"]
            blocked = set()
            for cu in cj["cubes"]:
                key = sorted(p for p, v in cu["section"] if v)
                for hperm in G.homs:
                    blocked.add(tuple(sorted(hperm[p] for p in key)))
            assert len(blocked) == cj["sections"], (len(blocked), cj["sections"])
            cnf, x, pts = C.plane_cnf(L, {})
            free = [p for p in pts if p not in C.pt_units]
            for b in sorted(blocked):
                bs = set(b)
                cnf.add(*[(-x[p] if p in bs else x[p]) for p in free])
            txt = cnf.dimacs(["level-1 completeness case %s class %d plane %d: %d sections = homothety orbits of %d cube representatives" % (case, ci, L, len(blocked), len(cj["cubes"]))])
            rep.check("cover1_%s%d" % (case, ci), hashlib.sha256(txt.encode()).hexdigest(), "cover1")
            print("cover1_%s%d rebuilt: %d sections from %d representatives" % (case, ci, len(blocked), len(cj["cubes"])), flush=True)
        # level 2 completeness
        if "cover2" in only:
            for sid, sp_path in sorted(split_files.items()):
                if not sid.startswith("%s%d_" % (case, ci)):
                    continue
                i = int(sid.split("_c")[1])
                units, secd = cube_units(C, cj, i)
                sp = json.load(open(sp_path))
                fixed = {abs(u) - 1: int(u > 0) for u in units}
                cnf, x, pts = C.plane_cnf(sp["plane"], fixed)
                free = [p for p in pts if p not in fixed and p not in C.pt_units]
                blocked = set()
                for ch in sp["children"]:
                    assert sorted(abs(u) - 1 for u in ch) == free
                    blocked.add(tuple(sorted(abs(u) - 1 for u in ch if u > 0)))
                assert len(blocked) == sp["n"]
                for b in sorted(blocked):
                    bs = set(b)
                    cnf.add(*[(-x[p] if p in bs else x[p]) for p in free])
                txt = cnf.dimacs(["level-2 completeness cube %s plane %d: %d children" % (sid, sp["plane"], len(blocked))])
                rep.check("cover2_%s" % sid, hashlib.sha256(txt.encode()).hexdigest(), "cover2")
                print("cover2_%s rebuilt: %d children" % (sid, len(blocked)), flush=True)
        os.remove(body)


def main():
    args = sys.argv[1:]
    only = set((args[args.index("--only") + 1] if "--only" in args else "planar,base,cubes,children,cover1,cover2").split(","))
    sel = set(args[args.index("--cases") + 1].split(",")) if "--cases" in args else None
    if "--selftest" in args:
        for w in ([1] * 4, [1] * 7, [1, 1, 1, 1, 1, 1, 1, 1]):
            assert R.test_shape(w, (2, 5)), w
        assert R.test_shape([1] * 6, (1, 4))
        print("encoder selftest OK")
    man = load_manifest()
    rep = Report(man)
    classes = json.load(open(os.path.join(RES, "planar_classes.json")))
    t0 = time.time()
    if "planar" in only:
        run_planar(rep)
    star_items = {"base", "cubes", "children", "cover1", "cover2"} & only
    if star_items:
        G = Geometry()
        cases = [("A", i) for i in range(3)] + [("C", i) for i in range(7)]
        if sel:
            cases = [c for c in cases if "%s%d" % c in sel]
        tmp = tempfile.mkdtemp(prefix="pg39_replay_")
        run_star(rep, G, classes, cases, only, tmp)
        try:
            os.rmdir(tmp)
        except OSError:
            pass
    import re
    pat = re.compile(r"^([AC]\d_c\d\d(_p\d+_k\d+)?|cover[12]_.*|planar_complete_.*)$")
    unchecked = sorted(t for t in man if pat.match(t) and t not in rep.checked) if (only >= {"planar", "base", "cubes", "children", "cover1", "cover2"} and not sel) else []
    if unchecked:
        print("manifest tags of the cover chain that were NOT rebuilt:", unchecked[:20])
    print("SUMMARY: %d CNFs rebuilt and equal to the recorded hash, %d MISMATCHES, %d tags not found in the manifest (%.0fs)" % (rep.ok, len(rep.bad), len(rep.missing), time.time() - t0))
    if rep.missing:
        print("not in manifest:", rep.missing[:12], "..." if len(rep.missing) > 12 else "")
    sys.exit(1 if (rep.bad or unchecked) else 0)


if __name__ == "__main__":
    main()

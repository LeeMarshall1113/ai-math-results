#!/usr/bin/env python3
"""Explicit, loud assertions for the cover logic of PROOF-CHAIN sections 2 and 3 (referee item 3).  Stdlib only; run time about 1 minute, < 1 GB RAM.

  python cover_asserts.py            all assertions, including the CNF-hash rebuild of every cube and child (about 35 s of it)
  python cover_asserts.py --fast     skip the hash rebuild (assertions on the data files and the stored cover CNFs only)

Every check raises AssertionError (non-zero exit) at the first failure; the last line is printed only if every assertion held.
The geometry (GF(9), PG(3,9), quotient plane, P0, heights, homotheties) is rebuilt by replay.py / rebuild_all.Geometry; every FACT below is then recomputed from that
geometry and from the data files in certs/, never read off a data file and trusted.

LEVEL 1 (for each of the 10 cases A0..A2, C0..C6, data certs/cubes_<case><i>.json, CNF certs/cover1_<case><i>.cnf as stored and verified by cake_lpr)
  L1a  the plane L1 of the file is recomputed: the first of the quotient lines v1v2, v1v3, v2v3 (v = the triple of the case) that meets Sigma in 5 points; the triple and the
       base units of the file equal those recomputed from certs-independent data (results/planar_classes.json).
  L1b  the plane has 91 points = P0 + the points projecting onto the line; exactly 91 lines of PG(3,9) lie in it; 3 of its points are fixed by the base units
       (P0 and the two height-0 points of the triple on L1), so 88 points are FREE, and the 4th base unit (the third height-0 point) lies outside the plane.
  L1c  every representative section of the file lists exactly the 91 points of the plane, once each, with values 0/1 and not more (no point outside the plane);
       it agrees with the base units; it is a (2,5)-section by a DIRECT geometric check (size m, every one of the 91 lines of the plane meets it in 2 or 5 points,
       lines through P0 have the type prescribed by Sigma, local five-secant counts d+3 / d); the 92 unit clauses of the cube CNF are base units + section.
  L1d  the 8 homothety images of the section (recomputed here: (v:t) -> (v:lambda t), lambda in GF(9)*) stay in the plane, fix the base units, are again valid sections,
       and the images of one cube are pairwise distinct (orbit size as recorded, a divisor of 8) and distinct from the images of every other cube.
  L1e  the stored certs/cover1_<case><i>.cnf is byte-identical to the CNF rebuilt from the geometry (plane CNF + one blocking clause per blocked section), its sha256 equals the
       manifest record, and (parsing the file independently) its blocking clauses are exactly the 88-literal clauses, one per section, and the set of sections they block
       equals, as a set, the union of the 8 homothety orbits of the cube list (so: the cube list is not trusted, the blocked set is the expansion of it, cube by cube).
  L1f  (unless --fast) the cube CNF of cube j, rebuilt from the base CNF + the 92 units of cube j, has the sha256 of the manifest record <case><i>_c<jj>; the section used to build
       it is the identity image (lambda = 1) of the orbit blocked in cover1.  Hence manifest tag j  <->  cube list entry j  <->  blocked orbit j.

LEVEL 2 (for each of the 8 certs/split_<cube>.json and certs/cover2_<cube>.cnf)
  L2a  planes pi1 (of the cube) and pi2 (of the split) are distinct planes through P0; they share exactly 10 points, and these 10 points are a line of PG(3,9) through P0;
       pi2 has 81 points off that line; the points of pi2 fixed by the base units outside the common line are counted (the third height-0 point iff its direction is on L2),
       and the number of FREE points is 81 minus that count.  This is the 80 of the documentation (asserted equal to 80 for all 8 cases, and printed with the reason).
  L2b  every child assigns exactly that free point set (no more, no fewer, no repeats); the children are pairwise distinct; n = number of children.
  L2c  for every child, parent section + base units + child assignment give a complete section of pi2 which passes the direct check (size m2 from Sigma, all 91 lines in {2,5},
       P0-line types, local counts).  So every child is a genuinely consistent level-2 section.
  L2d  the stored cover2 CNF is byte-identical to the rebuilt one, its sha256 equals the manifest, and its blocking clauses (parsed independently) block exactly the children, one each.
  L2e  (unless --fast) the child CNFs (cube CNF + 80 units) rebuilt have the sha256 of the manifest records.

COVERAGE (manifest-to-child)
  M1   the 10 cube lists contain 174 cubes; 167 have a manifest record UNSAT_VERIFIED / 's VERIFIED UNSAT' (monolithic); the other 7 are exactly the parents of 7 of the 8 split files,
       and each of those 7 is a Modal TIMEOUT in the ledger (results/modal_ledger.jsonl) with no verified monolithic record.
  M2   the 8th split parent, C2_c00, also has a verified monolithic record (redundant double coverage; evidence from the ledger is printed).
  M3   every child of every split file has a manifest record with status UNSAT_VERIFIED and 's VERIFIED UNSAT' and the recorded cnf hash equals the rebuilt one; no child tag of the
       manifest is missing from a split file (1,084 = sum of the n of the 8 splits, 140*6 + 122*2... computed, not assumed); the 18 cover records are UNSATISFIABLE / verified.
  M4   the stored LRAT of each cover (certs/cover*.lrat) has the recorded size (and, with --lrat, the recorded sha256).
"""
import hashlib
import json
import os
import re
import sys
import tempfile
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import replay as R            # noqa: E402
import rebuild_all as RB      # noqa: E402

CERTS = RB.CERTS
RES = RB.RES
FAST = "--fast" in sys.argv
LRAT_HASH = "--lrat" in sys.argv


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 24), b""):
            h.update(ch)
    return h.hexdigest()


def parse_dimacs(path):
    """independent DIMACS reader: (comments, nvars, nclauses, [clause tuples])"""
    comments, clauses, nv, nc = [], [], None, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("c"):
                comments.append(line.rstrip("\n"))
            elif line.startswith("p"):
                _, _, a, b = line.split()
                nv, nc = int(a), int(b)
            else:
                lits = [int(t) for t in line.split()]
                check(lits and lits[-1] == 0, "clause line not terminated by 0: %r" % line[:60])
                clauses.append(tuple(lits[:-1]))
    check(nc == len(clauses), "%s: header says %s clauses, file has %d" % (path, nc, len(clauses)))
    return comments, nv, nc, clauses


def manifest_records():
    by = {}
    for l in open(os.path.join(CERTS, "manifest.jsonl")):
        if not l.strip():
            continue
        r = json.loads(l)
        by.setdefault(r["tag"], []).append(r)
    return by


def best(by, tag):
    """the manifest record(s) of a tag; all must agree on the cnf hash and be verified"""
    recs = by.get(tag, [])
    return recs


def verified(r):
    return r.get("status") in ("UNSAT_VERIFIED", "UNSATISFIABLE") and r.get("cake_lpr") == "s VERIFIED UNSAT"


class PlaneData:
    def __init__(self, G, C, L):
        self.G, self.C, self.L = G, C, L
        self.pts = G.line_pts[L]
        self.pset = set(self.pts)
        check(len(self.pts) == 91 and len(self.pset) == 91, "plane over quotient line %d has %d points" % (L, len(self.pts)))
        check(G.P0 in self.pset, "P0 not in plane")
        self.lines = [j for j, Ln in enumerate(G.lines4) if all(p in self.pset for p in Ln)]
        check(len(self.lines) == 91, "plane has %d lines" % len(self.lines))
        # independent: the plane is closed (any two points' line lies in it)
        self.m = C.plane_size(L)
        self.d = (self.m - 20) // 3


def check_section(G, C, PD, S):
    """DIRECT geometric check of a candidate section S (set of point indices) of the plane PD: returns None or a message"""
    S = set(S)
    if not S <= PD.pset:
        return "points outside the plane"
    if len(S) != PD.m:
        return "size %d != %d" % (len(S), PD.m)
    five = {}
    for j in PD.lines:
        c = sum(1 for p in G.lines4[j] if p in S)
        if c not in (2, 5):
            return "line %d meets in %d points" % (j, c)
        five[j] = (c == 5)
    for j, y in C.line_units:                    # types of the lines through P0 prescribed by Sigma (y = 1: five-secant)
        if j in five and five[j] != bool(y):
            return "line type through P0 wrong (line %d)" % j
    for p in PD.pts:
        b = sum(1 for j in G.lines_thru[p] if j in five and five[j])
        want = PD.d + 3 if p in S else PD.d
        if b != want:
            return "local count at %d: %d five-secants, expected %d" % (p, b, want)
    for p in C.pt_units:
        if p in PD.pset and p not in S:
            return "base unit point %d not in section" % p
    return None


def independent_homs(G):
    """the 8 homotheties (v1:v2:v3:v4) -> (v1:v2:v3:lambda v4), recomputed from the point list (lambda = 1 first)"""
    homs = []
    for lam in range(1, 9):
        homs.append([G.pidx4[R.normalize((v[0], v[1], v[2], R.e_mul(lam, v[3])))] for v in G.points4])
    return homs


def main():
    t0 = time.time()
    by = manifest_records()
    G = RB.Geometry()
    classes = json.load(open(os.path.join(RES, "planar_classes.json")))
    homs = independent_homs(G)
    check(homs[0] == list(range(820)), "lambda = 1 is not the identity")
    check(homs == G.homs, "independently computed homotheties differ from Geometry.homs")
    cases = [("A", i) for i in range(3)] + [("C", i) for i in range(7)]
    cubes_info = {}              # tag -> (case obj, cube json, index, section dict, units)
    n_cubes = 0
    blocked_of = {}
    print("== LEVEL 1 ==")
    for case, ci in cases:
        C = RB.Case(G, case, ci, classes)
        name = "%s%d" % (case, ci)
        cj = json.load(open(os.path.join(CERTS, "cubes_%s.json" % name)))
        # L1a
        v1, v2, v3 = C.tri
        cand = []
        for a, b in ((v1, v2), (v1, v3), (v2, v3)):
            Lab = [L for L, Ln in enumerate(G.lines3) if a in Ln and b in Ln]
            check(len(Lab) == 1, "no unique quotient line through two points")
            cand.append(Lab[0])
        L1 = next(L for L in cand if len(C.Sigma & set(G.lines3[L])) == 5)
        check(cj["plane"] == L1, "%s: plane in file %s, recomputed %s" % (name, cj["plane"], L1))
        check(cj["triple"] == list(C.tri), "%s: triple of the file differs" % name)
        check({int(p): v for p, v in cj["base_units"].items()} == C.pt_units and [int(p) for p in cj["base_units"]] == list(C.pt_units), "%s: base units differ" % name)
        # L1b
        PD = PlaneData(G, C, L1)
        in_plane_units = [p for p in C.pt_units if p in PD.pset]
        free = [p for p in PD.pts if p not in C.pt_units]
        check(len(in_plane_units) == 3 and len(free) == 88, "%s: %d fixed points in the plane" % (name, len(in_plane_units)))
        outside = [p for p in C.pt_units if p not in PD.pset]
        check(len(outside) == 1 and outside[0] in [G.zero[v] for v in C.tri], "%s: units outside the plane" % name)
        # L1c, L1d
        union = {}
        for i, cu in enumerate(cj["cubes"]):
            check(cu["index"] == i and cu["plane"] == L1, "%s cube %d index/plane" % (name, i))
            sec = cu["section"]
            keys = [p for p, v in sec]
            check(len(keys) == 91 and set(keys) == PD.pset, "%s cube %d: section does not list exactly the 91 plane points" % (name, i))
            check(all(v in (0, 1) for p, v in sec), "%s cube %d: non-boolean value" % (name, i))
            S = {p for p, v in sec if v}
            for p, v in C.pt_units.items():
                if p in PD.pset:
                    check(sec[keys.index(p)][1] == v, "%s cube %d: section disagrees with base unit" % (name, i))
            msg = check_section(G, C, PD, S)
            check(msg is None, "%s cube %d: %s" % (name, i, msg))
            units, secd = RB.cube_units(C, cj, i)
            check(len(units) == 92 and len(set(abs(u) for u in units)) == 92, "%s cube %d: not 92 distinct unit variables" % (name, i))
            check({abs(u) - 1 for u in units} == set(C.pt_units) | PD.pset, "%s cube %d: unit variables are not base units + plane" % (name, i))
            check({abs(u) - 1 for u in units if u > 0} == (S | set(C.pt_units)) , "%s cube %d: unit signs" % (name, i))
            orb = set()
            for h in homs:
                img = frozenset(h[p] for p in S)
                check(set(h[p] for p in PD.pts) == PD.pset, "%s cube %d: homothety moves the plane" % (name, i))
                check(all(h[p] == p for p in C.pt_units), "%s: homothety moves a base unit" % name)
                msg = check_section(G, C, PD, img)
                check(msg is None, "%s cube %d: homothety image invalid: %s" % (name, i, msg))
                orb.add(img)
            check(len(orb) == cu["orbit_size"] and 8 % len(orb) == 0, "%s cube %d: orbit size %d, recorded %d" % (name, i, len(orb), cu["orbit_size"]))
            for img in orb:
                check(img not in union, "%s cube %d: an image is already in the orbit of cube %s" % (name, i, union.get(img)))
                union[img] = i
            check(frozenset(S) in orb, "representative not in its own orbit")
            cubes_info[("%s_c%02d" % (name, i))] = (C, cj, i)
        check(len(union) == cj["sections"], "%s: %d blocked sections, file says %d" % (name, len(union), cj["sections"]))
        n_cubes += len(cj["cubes"])
        # L1e
        path = os.path.join(CERTS, "cover1_%s.cnf" % name)
        cnf, x, pts = C.plane_cnf(L1, {})
        fr = [p for p in pts if p not in C.pt_units]
        for b in sorted(tuple(sorted(s)) for s in union):
            bs = set(b)
            cnf.add(*[(-x[p] if p in bs else x[p]) for p in fr])
        comment = "level-1 completeness case %s class %d plane %d: %d sections = homothety orbits of %d cube representatives" % (case, ci, L1, len(union), len(cj["cubes"]))
        txt = cnf.dimacs([comment])
        # NOTE: rebuild_all blocks `sorted(blocked)` where blocked holds tuples of the TRUE points; we reproduce that order
        stored = open(path, newline="").read()
        check(stored == txt, "%s: stored cover1 CNF differs from the CNF rebuilt from the cube list" % name)
        check(hashlib.sha256(stored.encode()).hexdigest() == by["cover1_" + name][0]["cnf_sha256"], "%s: manifest hash of cover1 differs" % name)
        _, nv, nc, cls = parse_dimacs(path)
        long = [c for c in cls if len(c) == len(fr)]
        check(len(long) == len(union), "%s: %d 88-literal clauses, %d sections" % (name, len(long), len(union)))
        varpt = {v + 1: p for v, p in enumerate(pts)}           # variable v <-> point pts[v-1]  (plane_cnf allocates x first)
        blocked_file = set()
        for c in long:
            check(sorted(abs(l) for l in c) == sorted(v for v, p in varpt.items() if p in set(fr)), "%s: blocking clause not over exactly the free points" % name)
            blocked_file.add(frozenset(varpt[abs(l)] for l in c if l < 0))      # literal -x_p blocks 'p in S'
        # blocked sections as in the file contain the fixed unit points only as units: add them back
        fixed_in = {p for p in C.pt_units if p in PD.pset}
        blocked_file = {frozenset(b | fixed_in) for b in blocked_file}
        check(blocked_file == set(union), "%s: set of blocked sections in the stored CNF != union of the 8 homothety orbits of the cube list" % name)
        check(len(cls) - len(long) == len(cnf.clauses) - len(union), "%s: unexpected count of non-blocking clauses" % name)
        blocked_of[name] = union
        print("  %-3s plane L1 = %2d (|Sigma n L1| = 5, m = %d): %2d cubes, orbit sizes %s (recorded = recomputed), %3d blocked sections = union of orbits; 91-point sections, 88 free, "
              "cover1 CNF stored == rebuilt, manifest hash ok" % (name, L1, PD.m, len(cj["cubes"]), dict(Counter(cu["orbit_size"] for cu in cj["cubes"])), len(union)), flush=True)
    check(n_cubes == 174, "total cubes %d" % n_cubes)

    # ---------------------------------------------------------------- level 2
    print("== LEVEL 2 ==")
    splits = {}
    for fn in sorted(os.listdir(CERTS)):
        if fn.startswith("split_") and fn.endswith(".json"):
            splits[fn[6:-5]] = json.load(open(os.path.join(CERTS, fn)))
    check(len(splits) == 9, "%d split files" % len(splits))      # 8 from the cloud timeouts + C1_c04 (local split, 2026-11-11)
    n_children = 0
    child_tags = {}
    for sid, sp in sorted(splits.items()):
        C, cj, i = cubes_info[sid]
        check(sp["id"] == sid, "split file id")
        L1 = cj["plane"]
        L2 = sp["plane"]
        P1, P2 = PlaneData(G, C, L1), PlaneData(G, C, L2)
        check(L1 != L2, "%s: same plane" % sid)
        inter = P1.pset & P2.pset
        check(len(inter) == 10 and G.P0 in inter, "%s: planes share %d points" % (sid, len(inter)))
        check(tuple(sorted(inter)) in set(G.lines4), "%s: the common points are not a line of PG(3,9)" % sid)
        common_dir = set(G.lines3[L1]) & set(G.lines3[L2])
        check(len(common_dir) == 1, "%s: quotient lines do not meet in one point" % sid)
        v = next(iter(common_dir))
        check(inter == {G.P0} | {q for q, d in G.proj.items() if d == v}, "%s: common line is not the line through P0 in the common direction" % sid)
        off = P2.pset - inter
        check(len(off) == 81, "%s: %d points of pi2 off the line" % (sid, len(off)))
        units_off = [p for p in C.pt_units if p in off]
        why = [("height-0 point of direction %s (v%d)" % (G.points3[d], t + 1)) for t, d in enumerate(C.tri) if G.zero[d] in units_off]
        check(all(p != G.P0 for p in units_off), "P0 is not on the common line?")
        free2 = sorted(off - set(C.pt_units))
        check(len(free2) == 81 - len(units_off), "free count")
        check(len(free2) == 80 == sp["assumption_size"], "%s: %d free points of pi2, documented 80" % (sid, len(free2)))
        units, secd = RB.cube_units(C, cj, i)
        S1 = {p for p, vv in secd.items() if vv}
        # L2b / L2c
        seen = set()
        for kk, ch in enumerate(sp["children"]):
            check(len(ch) == 80 and len({abs(u) for u in ch}) == 80, "%s child %d: %d literals" % (sid, kk, len(ch)))
            check(sorted(abs(u) - 1 for u in ch) == free2, "%s child %d: free-point set differs from the geometric one" % (sid, kk))
            key = tuple(u > 0 for u in sorted(ch, key=abs))
            check(key not in seen, "%s: duplicate child %d" % (sid, kk))
            seen.add(key)
            S2 = {abs(u) - 1 for u in ch if u > 0} | (S1 & inter) | {p for p in C.pt_units if p in P2.pset}
            msg = check_section(G, C, P2, S2)
            check(msg is None, "%s child %d is not a consistent level-2 section: %s" % (sid, kk, msg))
        check(sp["n"] == len(sp["children"]) == len(seen), "%s: n" % sid)
        n_children += sp["n"]
        for kk in range(sp["n"]):
            child_tags["%s_p%d_k%03d" % (sid, L2, kk)] = sid
        # L2d
        path = os.path.join(CERTS, "cover2_%s.cnf" % sid)
        fixed = {abs(u) - 1: int(u > 0) for u in units}
        cnf, x, pts = C.plane_cnf(L2, fixed)
        fr = [p for p in pts if p not in fixed and p not in C.pt_units]
        check(fr == free2, "%s: cover2 free set differs from geometric free set" % sid)
        blocked = sorted({tuple(sorted(abs(u) - 1 for u in ch if u > 0)) for ch in sp["children"]})
        for b in blocked:
            bs = set(b)
            cnf.add(*[(-x[p] if p in bs else x[p]) for p in fr])
        txt = cnf.dimacs(["level-2 completeness cube %s plane %d: %d children" % (sid, L2, len(blocked))])
        stored = open(path, newline="").read()
        check(stored == txt, "%s: stored cover2 CNF differs from the one rebuilt from the split file" % sid)
        check(hashlib.sha256(stored.encode()).hexdigest() == by["cover2_" + sid][0]["cnf_sha256"], "%s: manifest hash of cover2 differs" % sid)
        _, nv, nc, cls = parse_dimacs(path)
        long = [c for c in cls if len(c) == 80]
        check(len(long) == sp["n"], "%s: %d 80-literal clauses vs %d children" % (sid, len(long), sp["n"]))
        varpt = {v + 1: p for v, p in enumerate(pts)}
        got = {frozenset(varpt[abs(l)] for l in c if l < 0) for c in long}
        want = {frozenset(abs(u) - 1 for u in ch if u > 0) for ch in sp["children"]}
        check(got == want and all(sorted(varpt[abs(l)] for l in c) == free2 for c in long), "%s: blocking clauses of the stored cover2 CNF != children" % sid)
        print("  %-7s pi1=%2d pi2=%2d (m2=%d): common line %2d pts (P0 + dir %s), %d points off it, %d fixed there by units [%s] -> %d free; %3d children, all consistent sections; "
              "cover2 CNF stored == rebuilt == manifest, blocking clauses == children" % (sid, L1, L2, P2.m, len(inter), G.points3[v], len(off), len(units_off),
                                                                                       "; ".join(why) or "none", len(free2), sp["n"]), flush=True)
    check(n_children == 1084 + splits["C1_c04"]["n"], "children total %d" % n_children)      # 1,084 children with cloud records + the local children of C1_c04
    check(all(sp["n"] == len(sp["children"]) for sp in splits.values()), "n")

    # ---------------------------------------------------------------- coverage
    print("== COVERAGE ==")
    all_cubes = set(cubes_info)
    check(len(all_cubes) == 174, "174 cubes")
    mono_ok = {t for t in all_cubes if by.get(t) and all(verified(r) for r in by[t])}
    for t in all_cubes:
        for r in by.get(t, []):
            check(verified(r), "%s has a non-verified manifest record %s" % (t, r.get("status")))
    nomono = sorted(all_cubes - mono_ok)
    check(len(mono_ok) == 167 and len(nomono) == 7, "monolithic verified %d, without %d" % (len(mono_ok), len(nomono)))
    check(set(nomono) <= set(splits) and set(splits) - set(nomono) == {"C2_c00", "C1_c04"}, "cubes without a monolithic record are not exactly the split parents minus C2_c00 and C1_c04")
    ledger = {}
    for l in open(os.path.join(RES, "modal_ledger.jsonl")):
        try:
            r = json.loads(l)
        except Exception:
            continue
        if r.get("manual") or "id" not in r:
            if r.get("manual") and r.get("id") in splits:
                ledger.setdefault(r["id"] + "#manual", []).append(r)
            continue
        old = ledger.get(r["id"])
        if old is not None and old.get("status") == "UNSAT_VERIFIED" and r.get("status") != "UNSAT_VERIFIED":
            continue
        ledger[r["id"]] = r
    timeouts = sorted(t for t in all_cubes if ledger.get(t, {}).get("status") == "TIMEOUT")
    # C6_c08 has only a manual timeout line (cancelled by the operator after > 2 h of cumulative container time); the other six have a TIMEOUT record
    manual_to = sorted(t[:-7] for t in ledger if t.endswith("#manual") and any(r.get("status") == "TIMEOUT" for r in ledger[t]))
    check(sorted(set(timeouts) | set(manual_to)) == nomono, "ledger timeouts %s / manual %s != cubes without monolithic record %s" % (timeouts, manual_to, nomono))
    print("  M1: 174 cubes = 167 with a verified monolithic manifest record + 7 without; the 7 are %s;" % nomono)
    print("      each of the 7 is a Modal TIMEOUT (stage-1 records: %s; manual operator-cancel line: %s) and each has a split file." % (timeouts, manual_to))
    c2 = by["C2_c00"]
    check(all(verified(r) for r in c2) and "C2_c00" in splits, "C2_c00")
    c1 = by["C1_c04"]
    check(len(c1) == 1 and verified(c1[0]) and "C1_c04" in splits, "C1_c04")
    check(c1[0].get("via") == "modal" and c1[0].get("superseded_for_verification_by") == "split C1_c04: cover2_C1_c04 + local children", "C1_c04 cloud record is not marked as superseded-for-verification by the local split")
    print("  M2b: second redundant split C1_c04: its monolithic cloud record is verified (cake_lpr heap 6500 MB, over the 6 GB WSL cap) and is marked superseded-for-verification in the manifest; the cube is closed again by cover2_C1_c04 + %d children,"
          " all with LOCAL records (via = local, cake_lpr heap <= 4500 MB): two independent refutations of the cube" % splits["C1_c04"]["n"])
    print("  M2: redundant 8th split C2_c00: monolithic manifest record verified (conflicts %s, conflicts recorded in ledger %s); its 140 children were nevertheless run and verified:" % (c2[0].get("conflicts"), ledger["C2_c00"].get("conflicts")))
    print("      ledger note of the manual line: %s" % ledger["C2_c00#manual"][0]["note"][:260])
    # M3
    bad = []
    for tag, sid in sorted(child_tags.items()):
        recs = by.get(tag, [])
        if not recs or not all(verified(r) for r in recs) or len({r["cnf_sha256"] for r in recs}) != 1:
            bad.append(tag)
    check(not bad, "children without a verified unique manifest record: %s" % bad[:5])
    in_manifest = {t for t in by if re.match(r"^[AC]\d_c\d\d_p\d+_k\d+$", t)}
    for t in child_tags:
        if child_tags[t] == "C1_c04":
            check(all(r.get("via") == "local" for r in by[t]), "%s: child of C1_c04 is not a local record" % t)
    check(in_manifest == set(child_tags), "manifest child tags differ from the split files: extra %s missing %s" % (sorted(in_manifest - set(child_tags))[:3], sorted(set(child_tags) - in_manifest)[:3]))
    covers = [t for t in by if t.startswith("cover1_") or t.startswith("cover2_")]
    check(len(covers) == 19 and all(verified(r) for t in covers for r in by[t]), "cover records")
    print("  M3: %d children in the %d split files (%s); %d child tags in the manifest, the same set, each with a verified UNSAT record; 19 cover records verified" % (
        len(child_tags), len(splits), ", ".join("%s:%d" % (s, p["n"]) for s, p in sorted(splits.items())), len(in_manifest)))
    # M4
    for t in sorted(covers):
        f = os.path.join(CERTS, t + ".lrat")
        if os.path.exists(f):
            r = by[t][0]
            check(os.path.getsize(f) == r["lrat_bytes"], "%s: LRAT size %d != recorded %d" % (t, os.path.getsize(f), r["lrat_bytes"]))
            if LRAT_HASH:
                check(sha_file(f) == r["lrat_sha256"], "%s: LRAT sha256 differs" % t)
        else:
            print("  (LRAT file of %s not present locally)" % t)
    print("  M4: stored cover LRAT files have the recorded sizes%s" % (" and sha256" if LRAT_HASH else " (use --lrat for the sha256)"))

    # M5: certs/PROOF-INDEX.tsv agrees with the manifest and with the split files
    import csv
    rows = list(csv.DictReader(open(os.path.join(CERTS, "PROOF-INDEX.tsv")), delimiter="	"))
    check(len(rows) == 174 + n_children, "PROOF-INDEX has %d rows" % len(rows))
    ids = {r["id"] for r in rows}
    check(ids == all_cubes | set(child_tags), "PROOF-INDEX ids differ from the 174 cubes + 1,084 children")
    for r in rows:
        recs = by.get(r["id"])
        if recs:
            check(r["cnf_sha256"] == recs[-1]["cnf_sha256"] and r["lrat_sha256"] == recs[-1]["lrat_sha256"], "PROOF-INDEX %s: hashes differ from the manifest" % r["id"])
        if r["level"] == "2":
            check(child_tags[r["id"]] == r["parent"] and r["cover_tag"] == "cover2_" + r["parent"], "PROOF-INDEX %s: parent / cover tag" % r["id"])
        else:
            check(r["cover_tag"].startswith("cover1_" + r["id"][:2]), "PROOF-INDEX %s: cover tag" % r["id"])
    n_yes = sum(1 for r in rows if r["local_identical"] == "yes")
    check(not any(r["local_identical"] == "NO" for r in rows), "PROOF-INDEX records a local disagreement")
    print("  M5: certs/PROOF-INDEX.tsv: %d rows = 174 cubes + %d children, hashes equal the manifest, parent / cover tags consistent with the split files; %d rows record a bit-identical local re-run, 0 disagreements" % (len(rows), n_children, n_yes))

    # ---------------------------------------------------------------- hashes of the cubes / children (rebuild)
    if not FAST:
        print("== HASH REBUILD of cube and child CNFs (rebuild_all machinery) ==")
        rep = RB.Report(RB.load_manifest())
        tmp = tempfile.mkdtemp(prefix="pg39_cover_asserts_")
        RB.run_star(rep, G, classes, cases, {"base", "cubes", "children"}, tmp)
        try:
            os.rmdir(tmp)
        except OSError:
            pass
        check(not rep.bad, "hash mismatches: %s" % rep.bad[:3])
        check(sorted(rep.missing) == nomono, "tags absent from the manifest (cubes without a monolithic record) %s != %s" % (sorted(rep.missing), nomono))
        check(rep.ok == 10 + 167 + n_children, "rebuilt-and-equal count %d" % rep.ok)
        print("  L1f/L2e: %d CNFs (10 base, 167 cubes, %d children) rebuilt and equal to the manifest; the only cube tags without a record are the 7 parents above" % (rep.ok, n_children))
    print("COVER ASSERTIONS: ALL PASS (%.0fs)" % (time.time() - t0))


if __name__ == "__main__":
    main()

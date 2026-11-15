"""Benchmark CNFs for the speed check (runs inside a lease).

Cubes (all already timed with CaDiCaL --congruence=false, no proof; the sampler is deterministic so they are reproduced exactly):
  b_A0s0 : case A class 0, seed 1000, 1st random section        365 s, 2.31 M conflicts
  b_A0s1 : case A class 0, seed 1100, 1st                        704 s, 3.46 M
  b_A1s1 : case A class 1, seed 1001, 2nd                        352 s, 2.43 M
  b_C0s3 : case C class 0, seed 1300, 1st                        601 s, 5.45 M
  b_A1s0 : case A class 1, seed 1001, 1st                       1957 s, 19.7 M   (heavy)
  b_C2s0 : case C class 2, seed 1002, 1st                       2272 s, 13.8 M   (heavy)
Variants of the base CNF (for b_A1s1 and b_A0s0): v_noamo (no at-most-one per DP layer), v_noksum (no |K| = k sum), v_planes0 (exact plane sums
for the 91 planes through P0: 26 resp. 35 points, implied by the parameter table), v_loc91 (local counts 'five-secants through p = 35 + 30 x_p' at the
91 points of another plane through P0; implied by the line counts and |K| = k).  Cover CNF cov4 : 4 cubes of class A0 with selectors.
"""
import json
import os
import random
import sys

import full
from cnfenc import sum_equals, sum_in_set
from trivial_tree import Setup, write_cube_cnf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
CNF = os.path.join(ROOT, "cnf")


def cube_lits(S, ci, seed_off, nth, cap=20000):
    v1, v2, v3 = S.tri
    Ls = [S.quotient_line_through(v1, v2), S.quotient_line_through(v1, v3), S.quotient_line_through(v2, v3)]
    Ls = [L for L in Ls if S.plane_size(L) == 26] + [L for L in Ls if S.plane_size(L) != 26]
    secs = S.sections(Ls[0], {}, cap=cap)
    rng = random.Random(1000 + ci + 100 * seed_off)
    sec = None
    for _ in range(nth + 1):
        sec = rng.choice(secs)
    return [((p + 1) if v else -(p + 1)) for p, v in {**S.pt_units, **sec}.items()], Ls, sec


def variant(path, S, amo=True, ksum=True, plane_lines=(), local_pts=()):
    P, k = S.P, S.k
    cnf = full.CNFFile(path)
    x = [cnf.new_var() for _ in range(P.N)]
    y = []
    for L in P.lines:
        fin = sum_in_set(cnf, [x[p] for p in L], [1] * 10, {2, 5}, amo=amo)
        y.append(fin[5])
    if ksum:
        sum_equals(cnf, x, [1] * P.N, k, amo=False)
    for p, v in S.pt_units.items():
        cnf.add(x[p] if v else -x[p])
    for j, v in S.line_y.items():
        cnf.add(y[j] if v else -y[j])
    for L2 in plane_lines:
        pts = S.line_pts[L2]
        sum_equals(cnf, [x[p] for p in pts], [1] * len(pts), S.plane_size(L2), amo=False)
    base = (k - 182) // 3
    for p in local_pts:
        lits = [y[j] for j in P.lines_through[p]]
        fin = sum_in_set(cnf, lits, [1] * len(lits), {base, base + 30}, amo=False)
        cnf.add(-x[p], fin[base + 30])
        cnf.add(x[p], -fin[base + 30])
        cnf.add(x[p], fin[base])
    nv, nc = cnf.nvars, cnf.nclauses
    cnf.close(["variant amo=%s ksum=%s planes=%d local=%d" % (amo, ksum, len(plane_lines), len(local_pts))])
    return nv, nc


def cover_cnf(base, cubes_lits, out):
    """base CNF + selector variable per cube (-sel v units) + (OR sel)"""
    with open(base) as fh:
        line = fh.readline()
        comments = []
        while line.startswith("c"):
            comments.append(line)
            line = fh.readline()
        parts = line.split()
        nv, nc = int(parts[2]), int(parts[3])
    extra = []
    sels = []
    for j, lits in enumerate(cubes_lits):
        s = nv + 1 + j
        sels.append(s)
        for u in lits:
            extra.append("%d %d 0\n" % (-s, u))
    extra.append(" ".join(map(str, sels)) + " 0\n")
    with open(base) as fh, open(out, "w", newline="\n") as o:
        line = fh.readline()
        while line.startswith("c"):
            o.write(line)
            line = fh.readline()
        o.write("p cnf %d %d\n" % (nv + len(sels), nc + len(extra)))
        while True:
            ch = fh.read(1 << 22)
            if not ch:
                break
            o.write(ch)
        o.writelines(extra)


def main():
    out = {}
    specs = {"b_A0s0": ("A", 0, 0, 0), "b_A0s1": ("A", 0, 1, 0), "b_A1s1": ("A", 1, 0, 1), "b_C0s3": ("C", 0, 3, 0),
             "b_A1s0": ("A", 1, 0, 0), "b_C2s0": ("C", 2, 0, 0)}
    setups = {}
    for name, (case, ci, off, nth) in specs.items():
        if (case, ci) not in setups:
            S = Setup(case, ci)
            base = os.path.join(CNF, "full_%s_c%d.cnf" % (case, ci))
            if not os.path.exists(base):
                full.make_case(case, ci, base)
            setups[(case, ci)] = (S, base)
        S, base = setups[(case, ci)]
        lits, Ls, sec = cube_lits(S, ci, off, nth)
        write_cube_cnf(base, lits, os.path.join(CNF, name + ".cnf"))
        out[name] = dict(case=case, cls=ci, units=len(lits))
        print("made", name, flush=True)
    # variants on class A1 (cube b_A1s1) and A0 (b_A0s0)
    for tag, (case, ci, off, nth, name) in {"A1": ("A", 1, 0, 1, "b_A1s1"), "A0": ("A", 0, 0, 0, "b_A0s0")}.items():
        S, base = setups[(case, ci)]
        lits, Ls, sec = cube_lits(S, ci, off, nth)
        planes0 = [L for L in range(len(S.P2.lines))]          # all 91 planes through P0
        other = [L for L in range(len(S.P2.lines)) if L != Ls[0]][0]
        loc91 = [p for p in S.line_pts[other] if p != S.P0]
        for vname, kw in {"noamo": dict(amo=False), "noksum": dict(ksum=False), "planes0": dict(plane_lines=planes0), "loc91": dict(local_pts=loc91)}.items():
            vb = os.path.join(CNF, "vbase_%s_%s.cnf" % (tag, vname))
            nv, nc = variant(vb, S, **kw)
            write_cube_cnf(vb, lits, os.path.join(CNF, "v_%s_%s.cnf" % (name[2:], vname)))
            os.remove(vb)
            print("variant", tag, vname, nv, nc, flush=True)
    # cover CNF: 4 cubes of class A0 (the two timed ones + the first two orbit representatives)
    S, base = setups[("A", 0)]
    import cube_run
    L, reps, nsec = cube_run.cubes_for(S)
    lit_sets = [cube_lits(S, 0, 0, 0)[0], cube_lits(S, 0, 1, 0)[0]]
    for sec, orb in reps[:2]:
        lit_sets.append([((p + 1) if v else -(p + 1)) for p, v in {**S.pt_units, **sec}.items()])
        write_cube_cnf(base, lit_sets[-1], os.path.join(CNF, "b_A0r%d.cnf" % (len(lit_sets) - 3)))
    cover_cnf(base, lit_sets, os.path.join(CNF, "cov4_A0.cnf"))
    print("cover cnf made", flush=True)
    json.dump(out, open(os.path.join(ROOT, "results", "bench_cubes.json"), "w"))


if __name__ == "__main__":
    main()

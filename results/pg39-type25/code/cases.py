"""Case descriptors for the classes that need a case split (2a, 2b, 2d).

A case is a dict  {name, cls, k, ins:[points forced into K], outs:[points forced out of K]}  (point indices of our PG(3,9)).
The CNF of a case is the orbit CNF of the class, |K| = k, plus one unit clause per orbit touched by ins/outs
(units sorted by orbit index; see orbsat.OrbitInstance.build_cnf and replay/replay_cases.py).

Justification of each case list (the 'WLOG' arguments) is in NOTES.md, section 4.
"""
import itertools
import json
import os

import groups
from ff import GF
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))


def load_planar_classes():
    return json.load(open(os.path.join(HERE, "..", "results", "planar_classes.json")))


def cases_2a(P=None):
    """homology diag(1,1,1,-1): axis plane pi0 = {x4 = 0}, centre C = (0,0,0,1).  S0 = K n pi0 is a planar (2,5)-set of size m;
    (c = [C in K], m) determine k.  One case per (planar class, c)."""
    F = GF(9)
    P = P or PG(4, 9, F)
    P2 = PG(3, 9, F)
    pi0 = [P.pidx[P.normalize((v[0], v[1], v[2], 0))] for v in P2.points]
    C = P.pidx[(0, 0, 0, 1)]
    out = []
    for ci, cl in enumerate(load_planar_classes()):
        m = cl["k"]
        S0 = set(pi0[p] for p in cl["rep"])
        for c in (0, 1):
            k = {(0, 26): 260, (0, 35): 287, (1, 26): 287, (1, 35): 260}[(c, m)]
            ins = sorted(S0 | ({C} if c else set()))
            outs = sorted((set(pi0) - S0) | (set() if c else {C}))
            out.append(dict(name="2a_m%d_cls%d_c%d" % (m, ci, c), cls="2a_homology_3+1", k=k, ins=ins, outs=outs,
                            note="planar class %d (|Aut|=%d)" % (ci, cl["aut_order"])))
    return out


def pgl2_orbits_5():
    F = GF(9)
    P1 = PG(2, 9, F)
    omega = F.prim
    gens = [P1.perm_of_matrix([[1, 1], [0, 1]]), P1.perm_of_matrix([[omega, 0], [0, 1]]), P1.perm_of_matrix([[0, 1], [1, 0]])]

    def orbit_reps(sz):
        seen, reps = set(), []
        for s in itertools.combinations(range(10), sz):
            if s in seen:
                continue
            orb, fr = {s}, [s]
            while fr:
                nx = []
                for t in fr:
                    for g in gens:
                        u = tuple(sorted(g[i] for i in t))
                        if u not in orb:
                            orb.add(u)
                            nx.append(u)
                fr = nx
            seen |= orb
            reps.append(s)
        return reps
    return P1, orbit_reps(2), orbit_reps(5)


def cases_2b(P=None):
    """g = diag(1,1,-1,-1): fixed lines l1 = P(<e1,e2>), l2 = P(<e3,e4>); (|F1|,|F2|,k) in {(5,2,287),(5,5,260),(2,2,260)}
    ((2,5,287) is the image of (5,2,287) under the swap l1 <-> l2).  F_i up to PGL(2,9) acting on each line separately."""
    F = GF(9)
    P = P or PG(4, 9, F)
    P1, r2, r5 = pgl2_orbits_5()

    def on_l1(i):
        v = P1.points[i]
        return P.pidx[P.normalize((v[0], v[1], 0, 0))]

    def on_l2(i):
        v = P1.points[i]
        return P.pidx[P.normalize((0, 0, v[0], v[1]))]
    out = []
    combos = [(5, 2, 287)] + [(5, 5, 260)] + [(2, 2, 260)]
    for (a1, a2, k) in combos:
        for s1 in (r5 if a1 == 5 else r2):
            for s2 in (r5 if a2 == 5 else r2):
                ins = sorted([on_l1(i) for i in s1] + [on_l2(i) for i in s2])
                allf = [on_l1(i) for i in range(10)] + [on_l2(i) for i in range(10)]
                outs = sorted(set(allf) - set(ins))
                out.append(dict(name="2b_F%d_%d_k%d_%s_%s" % (a1, a2, k, "".join(map(str, s1)), "".join(map(str, s2))),
                                cls="2b_homology_2+2", k=k, ins=ins, outs=outs, note="F1=%s F2=%s" % (s1, s2)))
    return out


def cases_2d(P=None):
    F = GF(9)
    P = P or PG(4, 9, F)
    sets = json.load(open(os.path.join(HERE, "..", "results", "baer_fixed_sets.json")))
    by = {}
    for S in sets:
        by.setdefault(len(S), []).append(S)
    allB = sorted(set(p for S in sets for p in S) | set())
    B = [i for i, v in enumerate(P.points) if all(c in (0, 1, 2) for c in v)]
    out = []
    for name, S, k in (("S_empty", by[0][0], 260), ("S_point", by[1][0], 287), ("S_27", by[27][0], 287)):
        ins = sorted(S)
        outs = sorted(set(B) - set(S))
        out.append(dict(name="2d_" + name, cls="2d_baer", k=k, ins=ins, outs=outs, note="fixed set size %d" % len(S)))
    return out


FAMILIES = {"2a": cases_2a, "2b": cases_2b, "2d": cases_2d}


def build_case_cnf(inst, case, orbit_of_point):
    units = {}
    for p in case["ins"]:
        units[orbit_of_point[p]] = 1
    for p in case["outs"]:
        o = orbit_of_point[p]
        assert units.get(o, 0) != 1, "point forced in and out through one orbit"
        units[o] = 0
    return inst.build_cnf(case["k"], (2, 5), amo=True, units=sorted(units.items()))

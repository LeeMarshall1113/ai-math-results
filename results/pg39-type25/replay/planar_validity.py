#!/usr/bin/env python3
"""Owner audit of the certified planar solution lists (referee item 4).  Stdlib only, self-contained (own GF(9), own PG(2,9), own group construction; the only
import of the proof chain is OPTIONAL and at the very end, to compare the group of rebuild_all with the independent one).  Run time about 1 minute, < 300 MB.

  python planar_validity.py

For each of the four lists results/planar_k{26,35}_R{0,1}.json (the data whose completeness is certified by certs/planar_complete_*.cnf/.lrat) it checks EXHAUSTIVELY:

 (a) every listed solution is a set of exactly k points of PG(2,9) meeting EVERY one of the 91 lines in 2 or 5 points (direct geometric check, no CNF), and the
     local five-secant counts d+3 (point in the set) / d (point outside), d = (k-20)/3, hold (a consequence the CNF uses as extra clauses);
 (b) the normalisation of the frame: K meets the line l0 = {z = 0} exactly in the 5-set R of the file, and the origin (0,0,1) is in K;
 (c) H_R is constructed HERE from its definition  H_R = { v -> diag(A,1) frob^e (v) : A in GL(2,9), e in {0,1}, the induced map on l0 maps R onto R }
     (11,520 candidate maps, as permutations of the 91 points; distinct permutations only), its order is reported and compared with (i) the order stored in the file
     and (ii) the order forced by orbit-stabiliser, 8 * |PGammaL(2,9)| / |orbit of R| (the 8 = the kernel of the action on l0, the homologies diag(c,c,1));
     it is checked to be a group (closed under composition, contains the identity) fixing l0, the origin and R.  Then each solution x is checked to satisfy
     x <=_lex x o h for every h in H_R (the lex-leader clause set of rebuild_all.planar_cnf: positions = the 81 points off l0 in increasing index, 0 < 1),
     i.e. to be the lexicographically least member of its H_R orbit, and (as a consequence, also checked explicitly through canonical forms) no two listed
     solutions lie in one H_R orbit;
 (d) the R choices cover all orbits: PGammaL(2,9) is generated here (1,440 distinct permutations of the 10 points of l0, = 2 * 720), it has exactly two orbits on the 252 5-subsets,
     of sizes 180 and 72, and the R's of the files are one member of each.  (Together with: PGammaL(3,9) is transitive on lines, hand proof in PROOF-CHAIN 1.3 step 1,
     and the stabiliser of l0 and of the origin in PGammaL(3,9) is exactly {diag(A,c) frob^e} ~ H = {diag(A,1) frob^e}: a matrix fixing the point (0,0,1) has third
     column (0,0,c), a matrix preserving z = 0 has last row (0,0,c), hence block diagonal.)
 (e) the counting identity of planar_count.py is recomputed with THIS group: |Omega0|, N(R), N_total(m) = 91 (180 N(R0) + 72 N(R1)) / x5(m) and compared with the class list
     (sum |PGammaL(3,9)| / |Aut_i| from results/planar_classes.json).

Nothing is sampled: every solution of every list is checked against every element of H_R.
"""
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ---------------------------------------------------------------- GF(9) = GF(3)[i], i^2 = -1; element a + b i has code a + 3 b
ADD = [[0] * 9 for _ in range(9)]
MUL = [[0] * 9 for _ in range(9)]
for x in range(9):
    for y in range(9):
        a, b, c, d = x % 3, x // 3, y % 3, y // 3
        ADD[x][y] = (a + c) % 3 + 3 * ((b + d) % 3)
        MUL[x][y] = (a * c - b * d) % 3 + 3 * ((a * d + b * c) % 3)
FROB = [(x % 3) + 3 * ((-(x // 3)) % 3) for x in range(9)]            # x -> x^3 = conjugation
for x in range(9):
    cube = MUL[MUL[x][x]][x]                                          # x^3 computed directly as a cross-check of the table
    check(FROB[x] == cube, "frobenius table")
INV = [0] + [next(y for y in range(1, 9) if MUL[x][y] == 1) for x in range(1, 9)]
check(all(MUL[x][INV[x]] == 1 for x in range(1, 9)), "inverse")
NEG = [next(y for y in range(9) if ADD[x][y] == 0) for x in range(9)]


def normalize(v):
    for c in v:
        if c:
            iv = INV[c]
            return tuple(MUL[iv][t] for t in v)
    raise ValueError("zero vector")


PTS = [v for v in itertools.product(range(9), repeat=3) if next((c for c in v if c), None) == 1]
PIDX = {v: i for i, v in enumerate(PTS)}
check(len(PTS) == 91, "PG(2,9) has 91 points")


def dot(a, v):
    s = 0
    for x, y in zip(a, v):
        s = ADD[s][MUL[x][y]]
    return s


LINES = [tuple(i for i, v in enumerate(PTS) if dot(a, v) == 0) for a in PTS]      # line a: {x : a.x = 0}
check(all(len(L) == 10 for L in LINES) and len(set(LINES)) == 91, "lines")
THRU = [[j for j, L in enumerate(LINES) if p in L] for p in range(91)]
check(all(len(t) == 10 for t in THRU), "10 lines through a point")
L0 = [i for i, v in enumerate(PTS) if v[2] == 0]
ORIGIN = PIDX[(0, 0, 1)]
check(len(L0) == 10, "l0")
L0SET = set(L0)
ORDER = [p for p in range(91) if p not in L0SET]             # lex positions (points off l0, increasing index), as in rebuild_all.planar_cnf


def act(A, e, v):
    """the collineation  v -> diag(A,1) . frob^e(v)  on a vector"""
    w = tuple(FROB[c] for c in v) if e else v
    o0 = ADD[MUL[A[0]][w[0]]][MUL[A[1]][w[1]]]
    o1 = ADD[MUL[A[2]][w[0]]][MUL[A[3]][w[1]]]
    return normalize((o0, o1, w[2]))


def perm_of(A, e):
    return tuple(PIDX[act(A, e, v)] for v in PTS)


def det(A):
    return ADD[MUL[A[0]][A[3]]][NEG[MUL[A[1]][A[2]]]]


# ---------------------------------------------------------------- (c)/(d): the groups
def build_groups():
    t = time.time()
    perms = {}                                  # distinct permutations of the 91 points -> some (A,e)
    on_l0 = {}                                  # restriction to l0, as a tuple over the sorted l0 indices
    n_el = 0
    for A in itertools.product(range(9), repeat=4):
        if det(A) == 0:
            continue
        for e in (0, 1):
            n_el += 1
            p = perm_of(A, e)
            perms[p] = (A, e)
    check(n_el == 5760 * 2, "|GL(2,9)| = %d, expected 5760" % (n_el // 2))
    print("group of maps diag(A,1) frob^e: %d pairs (A,e) -> %d distinct permutations of the 91 points (%.1fs)" % (n_el, len(perms), time.time() - t))
    check(len(perms) == n_el, "the maps diag(A,1) frob^e act faithfully on PG(2,9)")
    pos = {p: i for i, p in enumerate(L0)}
    res = {}
    for p in perms:
        res[p] = tuple(pos[p[q]] for q in L0)                                        # permutation of the 10 points of l0
    pgaml2 = set(res.values())
    check(len(pgaml2) == 2 * 9 * (81 - 1) == 1440, "induced group on l0 has order %d, |PGammaL(2,9)| = 1440" % len(pgaml2))
    # orbits on the 252 five-subsets of l0
    orb_of, reps = {}, []
    for s in itertools.combinations(range(10), 5):
        if s in orb_of:
            continue
        orb = {tuple(sorted(g[i] for i in s)) for g in pgaml2}
        for u in orb:
            orb_of[u] = len(reps)
        reps.append((s, len(orb)))
    check(sorted(r[1] for r in reps) == [72, 180] and sum(r[1] for r in reps) == 252, "orbits on 5-subsets: %s" % reps)
    print("PGammaL(2,9) induced on l0: %d permutations; orbits on the 252 five-subsets of l0 have sizes %s" % (len(pgaml2), [r[1] for r in reps]))
    return perms, res, orb_of, reps


def stabiliser(perms, res, Rl0):
    """H_R as a list of permutations of the 91 points: elements whose action on l0 maps R (positions in l0) onto R"""
    H = [p for p in perms if {res[p][i] for i in Rl0} == set(Rl0)]
    return H


def compose(g, h):                               # (g o h)(p) = g[h[p]]
    return tuple(g[h[p]] for p in range(91))


def lexcmp_ok(Sset, h):
    """x <=_lex x o h  on the positions ORDER, with x the indicator of S, (x o h)[p] = x[h[p]]; 0 < 1"""
    for p in ORDER:
        a = 1 if p in Sset else 0
        b = 1 if h[p] in Sset else 0
        if a != b:
            return a < b
    return True


def main():
    t0 = time.time()
    perms, res, orb_of, reps = build_groups()
    identity = tuple(range(91))
    G_total = 2 * ((9 ** 3 - 1) * (9 ** 3 - 9) * (9 ** 3 - 81) // 8)                 # |PGammaL(3,9)| = 2 |PGL(3,9)|
    classes = json.load(open(os.path.join(RES, "planar_classes.json")))
    summary = {}
    Hs = {}
    for ri in (0, 1):
        d0 = json.load(open(os.path.join(RES, "planar_k26_R%d.json" % ri)))
        Rpts = d0["R"]
        pos = {p: i for i, p in enumerate(L0)}
        check(all(p in L0SET for p in Rpts) and len(Rpts) == 5, "R of the file is a 5-subset of l0")
        Rl0 = tuple(sorted(pos[p] for p in Rpts))
        osz = len({tuple(sorted(g[i] for i in Rl0)) for g in set(res.values())})
        check(osz == d0["orbit_size"], "orbit size of R%d: computed %d, file %d" % (ri, osz, d0["orbit_size"]))
        H = stabiliser(perms, res, Rl0)
        theory = 8 * 1440 // osz
        check(len(H) == theory == d0["Hsize"], "H_R%d order: constructed %d, orbit-stabiliser %d, file %d" % (ri, len(H), theory, d0["Hsize"]))
        Hset = set(H)
        check(identity in Hset, "identity in H")
        check(all(h[ORIGIN] == ORIGIN and {h[p] for p in L0} == L0SET and {h[p] for p in Rpts} == set(Rpts) for h in H), "H_R fixes origin, l0, R")
        check(all(compose(g, h) in Hset for g in H for h in H), "H_R%d is not closed under composition" % ri)
        print("H_R%d: orbit of R under PGammaL(2,9) has %d members; |H_R| = %d (constructed) = 8*1440/%d (orbit-stabiliser) = %d (file); group closure checked on all %d products" % (
            ri, osz, len(H), osz, d0["Hsize"], len(H) ** 2))
        Hs[ri] = H
        # cross-check against rebuild_all's own construction (optional import, same chain)
    try:
        sys.path.insert(0, HERE)
        import rebuild_all as RB
        pts_, pidx_, lines_, thru_, l0_, origin_, results_ = RB.planar_items()
        check(pts_ == PTS and l0_ == L0 and origin_ == ORIGIN, "point order differs from rebuild_all")
        for ri, osz, Rpts, H in results_:
            same = {tuple(h) for h in H} == set(Hs[ri])
            check(same, "H_R%d of rebuild_all differs from the independently constructed group" % ri)
            check(len(H) == len(Hs[ri]), "size")
        print("cross-check: the groups H_R0, H_R1 used by rebuild_all.planar_cnf are IDENTICAL (as sets of permutations) to the ones constructed here")
    except ImportError:
        print("(rebuild_all not importable; cross-check skipped)")
    # lists
    counts = {}
    for m in (26, 35):
        d_ = (m - 20) // 3
        omega = {}
        for ri in (0, 1):
            f = os.path.join(RES, "planar_k%d_R%d.json" % (m, ri))
            dd = json.load(open(f))
            H = Hs[ri]
            Rset = set(dd["R"])
            sols = dd["solutions"]
            canon = set()
            om0 = 0
            n_nonid_stab = 0
            for si, sol in enumerate(sols):
                S = set(sol)
                check(len(sol) == len(S) == m and list(sol) == sorted(sol) and all(0 <= p < 91 for p in sol), "k%d R%d sol %d: not a sorted %d-set of points" % (m, ri, si, m))
                # (a)
                for j, L in enumerate(LINES):
                    c = sum(1 for p in L if p in S)
                    check(c in (2, 5), "k%d R%d sol %d: line %d meets in %d points" % (m, ri, si, j, c))
                five = [sum(1 for p in L if p in S) == 5 for L in LINES]
                for p in range(91):
                    b = sum(1 for j in THRU[p] if five[j])
                    check(b == (d_ + 3 if p in S else d_), "k%d R%d sol %d: local count at %d" % (m, ri, si, p))
                # (b)
                check(S & L0SET == Rset, "k%d R%d sol %d: K n l0 != R" % (m, ri, si))
                check(ORIGIN in S, "k%d R%d sol %d: origin not in K" % (m, ri, si))
                # (c)
                for h in H:
                    check(lexcmp_ok(S, h), "k%d R%d sol %d is not lex-least in its H_R orbit" % (m, ri, si))
                imgs = {frozenset(p for p in range(91) if h[p] in S) for h in H}            # the orbit {x o h} as sets
                stab = sum(1 for h in H if all((h[p] in S) == (p in S) for p in range(91)))
                check(len(imgs) * stab == len(H), "orbit-stabiliser inconsistency")
                key = min(tuple(1 if h[p] in S else 0 for p in ORDER) for h in H)
                check(key == tuple(1 if p in S else 0 for p in ORDER), "canonical form (min over H_R) is not the listed solution itself")
                check(key not in canon, "k%d R%d: two listed solutions in the same H_R orbit" % (m, ri))
                canon.add(key)
                om0 += len(H) // stab
                n_nonid_stab += (stab > 1)
            check(len(canon) == len(sols) == len({tuple(s) for s in sols}), "duplicates in the list")
            omega[ri] = om0
            print("k=%d R#%d: %d solutions: all (2,5)-sets (91 lines each), K n l0 = R, origin in K, lex-least in the H_R orbit (|H_R| = %d, checked against every element), pairwise distinct orbits; "
                  "%d have a non-trivial stabiliser; |Omega0| = %d" % (m, ri, len(sols), len(H), n_nonid_stab, om0))
        # (e)
        x5 = (10 * m - 182) // 3
        N = {ri: 81 * omega[ri] // (m - 5) for ri in (0, 1)}
        check(all((81 * omega[ri]) % (m - 5) == 0 for ri in (0, 1)), "divisibility")
        sizes = {ri: json.load(open(os.path.join(RES, "planar_k%d_R%d.json" % (m, ri))))["orbit_size"] for ri in (0, 1)}
        tot = sizes[0] * N[0] + sizes[1] * N[1]
        check((91 * tot) % x5 == 0, "divisibility 2")
        n_total = 91 * tot // x5
        cl = [c for c in classes if c["k"] == m]
        n_cls = sum(G_total // c["aut_order"] for c in cl)
        print("k=%d: N(R0), N(R1) = %d, %d;  N_total = %d  (frame count, own group);  class list sum |PGammaL(3,9)|/|Aut_i| = %d  (|G| = %d, |Aut| = %s)  ->  %s" % (
            m, N[0], N[1], n_total, n_cls, G_total, [c["aut_order"] for c in cl], "EQUAL" if n_total == n_cls else "DIFFERENT"))
        check(n_total == n_cls, "counting identity fails for m = %d" % m)
        counts[m] = n_total
    print("PLANAR VALIDITY: ALL CHECKS PASS  (N_total(26) = %d, N_total(35) = %d; %.0fs)" % (counts[26], counts[35], time.time() - t0))


if __name__ == "__main__":
    main()

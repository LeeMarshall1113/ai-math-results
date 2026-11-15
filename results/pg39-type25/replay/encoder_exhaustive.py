#!/usr/bin/env python3
"""Exhaustive / differential semantic tests of the encoders of the proof chain (referee item 5; companion of ../ENCODER-SEMANTICS.md).  Stdlib only.

  python encoder_exhaustive.py [--quick] [--only E1,E2,...]

The encoder under test is the REAL one: replay.sum_in_set (the encoder whose output bytes are hashed in the manifest), the real rebuild_all.planar_cnf for the real planar
shape, and the real lex-leader loop (text-compared with rebuild_all.py).  Method (see ENCODER-SEMANTICS.md, section 3): for an input vector z the CNF has an extension iff
unit propagation from the inputs (plus the unit clauses) ends WITHOUT conflict in a total assignment that satisfies every clause; a propagation conflict refutes every extension
(unit propagation is sound), so "valid input => extension" and "invalid input => none" are decided per input, with nothing assumed.  A third kind of result,
'no conflict but not total', would be reported as a failure of the test (it never happens).

  E1  production line shape: 10 inputs of weight 1, allowed {2,5}, amo=True and amo=False.  ALL 1,024 inputs; additionally the CNF is model-COUNTED by an independent DPLL counter
      (no input fixed): the number of models equals 45 + 252 = 297, so each valid input has exactly one extension and there are no spurious models.
  E2  production local-count shapes: 10 inputs, allowed {2,5} (m=26) and {5,8} (m=35), amo=False; and sum shapes.  ALL 1,024 inputs each; model count = number of valid inputs.
  E3  random shapes (weights 1..3, random allowed sets, amo on/off), n <= 12: ALL 2^n inputs each (300 shapes); and the unit-weight shapes n = 14 (targets {5}, {7}, {2,5}, {5,8}) and
      n = 16 ({2,5}): ALL 2^n inputs.
  E4  real planar shape, q = 9: rebuild_all.planar_cnf (91 point variables, 91 lines {2,5}, sum = m, local counts, origin unit; no lex, no blocking) for m = 26 and 35, (+) the
      1,345 real solutions of results/planar_k*_R*.json, all of them: must extend; (-) for each of them a random swap (one point out, one in) and a one-point deletion: UP
      must conflict EXACTLY when the DIRECT geometric check says the set is not a valid (2,5)-set containing the origin (a swap may produce another genuine set; then UP must accept).
  E5  generic plane pipeline on other planes (replica of planar_cnf, checked clause-for-clause identical to planar_cnf on PG(2,9)): PG(2,3) type (1,4), m = 4: ALL 8,192 subsets;
      PG(2,4) type (0,2), m = 6 (hyperovals) and type (3,5), m = 15: ALL C(21,6) and C(21,15) subsets, and all subsets of size 5, 7, 14, 16.  CNF-accepted sets == brute-force sets.
  E6  lex-leader gadget (the loop of rebuild_all.planar_cnf, text-compared): random permutations / groups on 10 and 12 points, ALL 2^n inputs: accepted inputs ==
      { x : x <=_lex x o h for every h }.
  E7  real sizes of the plain sum constraints: [1]*91 = 26, 35 and [1]*820 = 260, 287: random valid inputs (must extend, total assignment satisfies all clauses) and inputs of
      weight target +-1, +-2 and 0, all (must conflict).  SAMPLED, not exhaustive (the sets are far too large: C(91,26) ~ 1e23); the argument is the proof of section 2.
"""
import itertools
import json
import os
import random
import sys
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import replay as R            # noqa: E402
import rebuild_all as RB      # noqa: E402

RES = os.path.join(HERE, "..", "results")
QUICK = "--quick" in sys.argv
ONLY = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else None


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ---------------------------------------------------------------- unit-propagation engine
class UP:
    def __init__(self, nvars, clauses):
        self.nv = nvars
        self.cl = [tuple(c) for c in clauses]
        self.occ = defaultdict(list)
        for i, c in enumerate(self.cl):
            for l in c:
                self.occ[l].append(i)
        self.base = None
        self.base_conflict = self._init()

    def _init(self):
        self.val = [0] * (self.nv + 1)                 # 0 unassigned, 1 true, -1 false
        self.nfalse = [0] * len(self.cl)
        self.sat = [False] * len(self.cl)
        conflict = False
        queue = []
        for i, c in enumerate(self.cl):
            if len(c) == 0:
                conflict = True
            elif len(c) == 1:
                queue.append(c[0])
        if not conflict:
            conflict = not self._propagate(queue)
        self.base = (self.val[:], self.nfalse[:], self.sat[:])
        return conflict

    def _propagate(self, queue):
        val, nfalse, sat, cl, occ = self.val, self.nfalse, self.sat, self.cl, self.occ
        while queue:
            l = queue.pop()
            v = abs(l)
            want = 1 if l > 0 else -1
            if val[v] == want:
                continue
            if val[v] == -want:
                return False
            val[v] = want
            for ci in occ[l]:
                sat[ci] = True
            for ci in occ[-l]:
                nfalse[ci] += 1
                if sat[ci]:
                    continue
                c = cl[ci]
                k = len(c) - nfalse[ci]
                if k == 0:
                    return False
                if k == 1:
                    for u in c:
                        if val[abs(u)] == 0:
                            queue.append(u)
                            break
        return True

    def run(self, assumptions):
        """returns (conflict?, total assignment as list or None).  assumptions: list of literals"""
        if self.base_conflict:
            return True, None
        self.val, self.nfalse, self.sat = self.base[0][:], self.base[1][:], self.base[2][:]
        if not self._propagate(list(assumptions)):
            return True, None
        return False, self.val

    def total_and_model(self, val):
        if any(val[v] == 0 for v in range(1, self.nv + 1)):
            return False
        for c in self.cl:
            if not any((val[abs(l)] == 1) == (l > 0) for l in c):
                return False
        return True


def decide(up, inputs, bits, search=False):
    """True iff the CNF (with its units) has an extension of the input assignment `bits`.
    Unit propagation alone decides it for the DP encoder (conflict, or a total model).  With search=True (used for the lex gadget, whose auxiliary flags are NOT functions of x)
    a complete DPLL search over the remaining variables is run when propagation leaves variables open."""
    assum = [(v if b else -v) for v, b in zip(inputs, bits)]
    conflict, val = up.run(assum)
    if conflict:
        return False
    if up.total_and_model(val):
        return True
    if not search:
        raise AssertionError("unit propagation neither conflicted nor produced a total model (inputs %s)" % (bits,))
    state = (up.val[:], up.nfalse[:], up.sat[:])

    def rec(state):
        v0, nf, st = state
        free = next((v for v in range(1, up.nv + 1) if v0[v] == 0), None)
        if free is None:
            return up.total_and_model(v0)
        for lit in (free, -free):
            up.val, up.nfalse, up.sat = v0[:], nf[:], st[:]
            if up._propagate([lit]) and rec((up.val, up.nfalse, up.sat)):
                return True
        return False
    return rec(state)


def count_models(nv, clauses, limit=10 ** 6):
    """independent DPLL model counter over ALL variables (no input fixed)"""
    up = UP(nv, clauses)
    if up.base_conflict:
        return 0
    cnt = 0
    base = up.base
    # recursive DPLL using the same propagation engine, restoring states by copying (small instances only)
    def rec(state):
        nonlocal cnt
        val, nfalse, sat = state
        free = next((v for v in range(1, nv + 1) if val[v] == 0), None)
        if free is None:
            cnt += 1
            check(cnt <= limit, "model count limit")
            return
        for lit in (free, -free):
            up.val, up.nfalse, up.sat = val[:], nfalse[:], sat[:]
            if up._propagate([lit]):
                rec((up.val, up.nfalse, up.sat))
    rec(base)
    return cnt


def build(shape_weights, allowed, amo):
    cnf = R.CNF()
    ins = [cnf.new_var() for _ in shape_weights]
    R.sum_in_set(cnf, ins, shape_weights, set(allowed), amo=amo)
    return cnf, ins


def exhaustive_shape(weights, allowed, amo, count=False):
    cnf, ins = build(weights, allowed, amo)
    up = UP(cnf.nvars, cnf.clauses)
    n = len(weights)
    nvalid = 0
    for bits in itertools.product((False, True), repeat=n):
        want = sum(w for w, b in zip(weights, bits) if b) in allowed
        got = decide(up, ins, bits)
        check(got == want, "shape %s allowed %s amo %s input %s: expected %s, CNF says %s" % (weights, sorted(allowed), amo, bits, want, got))
        nvalid += want
    if count:
        c = count_models(cnf.nvars, cnf.clauses)
        check(c == nvalid, "model count %d != valid inputs %d (spurious or non-unique auxiliary assignments) for %s %s" % (c, nvalid, weights, sorted(allowed)))
    return nvalid, cnf.nvars, len(cnf.clauses)


def want_test(name):
    return ONLY is None or name in ONLY


# ---------------------------------------------------------------- E1 / E2
def e1():
    print("E1/E2: production shapes, all 2^10 inputs, plus independent model counting")
    for allowed, amo, label in (({2, 5}, True, "line shape {2,5} amo=True"), ({2, 5}, False, "line shape {2,5} amo=False"),
                                ({2, 5}, False, "local count m=26 {d, d+3} = {2,5}"), ({5, 8}, False, "local count m=35 {d, d+3} = {5,8}"),
                                ({26}, False, "x"), ):
        if label == "x":
            continue
        nv, nvars, ncl = exhaustive_shape([1] * 10, allowed, amo, count=True)
        print("  %-38s: %4d valid inputs of 1024 accepted, all others rejected; models of the CNF counted = %d (= valid inputs); %d variables, %d clauses" % (label, nv, nv, nvars, ncl), flush=True)
        if allowed == {2, 5}:
            check(nv == 45 + 252, "297")


def e3():
    rnd = random.Random(20261111)
    print("E3: random shapes and unit-weight shapes up to n = 16, ALL inputs")
    t = time.time()
    n_shapes = 0
    tot = 0
    for _ in range(100 if QUICK else 300):
        n = rnd.randint(1, 12 if not QUICK else 9)
        weights = [rnd.randint(1, 3) for _ in range(n)]
        allowed = set(rnd.sample(range(0, sum(weights) + 2), rnd.randint(1, 3)))
        amo = rnd.random() < 0.5
        nv, _, _ = exhaustive_shape(weights, allowed, amo, count=(n <= 7))
        n_shapes += 1
        tot += 2 ** n
    print("  %d random weighted shapes (n <= %d, weights 1..3, 1-3 allowed sums, amo on/off): all %d inputs correct; model count = valid inputs for the n <= 7 ones (%.0fs)" % (n_shapes, 12 if not QUICK else 9, tot, time.time() - t), flush=True)
    for n, allowed in ((14, {5}), (14, {7}), (14, {2, 5}), (14, {5, 8}), (16, {2, 5})):
        if QUICK and n > 12:
            continue
        t = time.time()
        nv, nvars, ncl = exhaustive_shape([1] * n, allowed, False)
        print("  unit weights n=%d allowed %s, amo=False: all %d inputs correct (%d valid; %d vars, %d clauses; %.0fs)" % (n, sorted(allowed), 2 ** n, nv, nvars, ncl, time.time() - t), flush=True)


# ---------------------------------------------------------------- E4: real planar shape
def planes_data():
    pts, pidx, lines, thru, l0, origin, results = RB.planar_items()
    return pts, pidx, lines, thru, l0, origin


def direct_valid(S, lines, thru, origin, m):
    if len(S) != m or origin not in S:
        return False
    for L in lines:
        if sum(1 for p in L if p in S) not in (2, 5):
            return False
    return True


def e4():
    print("E4: real planar CNF (rebuild_all.planar_cnf), q = 9")
    pts, pidx, lines, thru, l0, origin = planes_data()
    rnd = random.Random(7)
    for m in (26, 35):
        cnf = RB.planar_cnf(m, pts, pidx, lines, thru, [], origin, set(), [], [])
        up = UP(cnf.nvars, cnf.clauses)
        check(not up.base_conflict, "base conflict")
        ins = list(range(1, 92))                      # x[p] = variable p+1 (planar_cnf allocates the point variables first)
        sols = []
        for ri in (0, 1):
            sols += json.load(open(os.path.join(RES, "planar_k%d_R%d.json" % (m, ri))))["solutions"]
        if QUICK:
            sols = sols[:60]
        t = time.time()
        pos_ok = neg_checked = neg_genuine = 0
        for sol in sols:
            S = set(sol)
            check(decide(up, ins, [p in S for p in range(91)]), "m=%d: a real solution is rejected by the CNF" % m)
            pos_ok += 1
            for kind in ("swap", "swap", "delete", "add"):
                T = set(S)
                if kind == "swap":
                    T.remove(rnd.choice(sorted(T)))
                    T.add(rnd.choice([p for p in range(91) if p not in S]))
                elif kind == "delete":
                    T.remove(rnd.choice(sorted(T)))
                else:
                    T.add(rnd.choice([p for p in range(91) if p not in S]))
                want = direct_valid(T, lines, thru, origin, m)
                got = decide(up, ins, [p in T for p in range(91)])
                check(got == want, "m=%d %s: CNF says %s, direct check says %s, set %s" % (m, kind, got, want, sorted(T)))
                neg_checked += 1
                neg_genuine += want
        print("  m=%d (%d vars, %d clauses): %d real solutions all accepted; %d perturbed sets (swap/delete/add) all decided as the direct geometric check decides (%d of them were genuine (2,5)-sets and accepted) (%.0fs)" % (
            m, cnf.nvars, len(cnf.clauses), pos_ok, neg_checked, neg_genuine, time.time() - t), flush=True)


# ---------------------------------------------------------------- E5: replica pipeline on other planes
def singer_plane(q):
    n = q * q + q + 1
    for D in itertools.combinations(range(n), q + 1):
        if D[0] != 0:
            break
        diffs = [(a - b) % n for a in D for b in D if a != b]
        if len(set(diffs)) == len(diffs):
            lines = [tuple(sorted((d + s) % n for d in D)) for s in range(n)]
            return n, lines
    raise RuntimeError


def generic_plane_cnf(n, lines, a, b, m, q, origin=None):
    """replica of rebuild_all.planar_cnf without frame/lex/blocking, for line numbers {a,b}: variables x[p] first, y_l = [line meets in b]; local counts"""
    thru = [[j for j, L in enumerate(lines) if p in L] for p in range(n)]
    cnf = R.CNF()
    x = [cnf.new_var() for _ in range(n)]
    y = []
    for L in lines:
        fin = R.sum_in_set(cnf, [x[p] for p in L], [1] * len(L), {a, b}, amo=True)
        y.append(fin[b])
    R.sum_in_set(cnf, x, [1] * n, {m}, amo=False)
    # local counts: through p in S: v_in = (m-1-(q+1)(a-1))/(b-a) lines of type b ; p not in S: v_out = (m-(q+1)a)/(b-a); v_in = v_out + q/(b-a)
    num = m - (q + 1) * a
    check(num % (b - a) == 0 and q % (b - a) == 0, "integrality of the local counts")
    d = num // (b - a)
    dl = d + q // (b - a)
    for p in range(n):
        t = [y[j] for j in thru[p]]
        fin = R.sum_in_set(cnf, t, [1] * len(t), {d, dl}, amo=False)
        cnf.add(-x[p], fin[dl])
        cnf.add(x[p], -fin[dl])
    return cnf, x, thru


def e5():
    rnd = random.Random(5)
    print("E5: the plane pipeline (line shapes, size, local counts) on PG(2,3) and PG(2,4): brute force vs CNF")
    # replica identical to planar_cnf on PG(2,9)?
    pts, pidx, lines9, thru9, l0, origin = planes_data()
    ref = RB.planar_cnf(26, pts, pidx, lines9, thru9, [], origin, set(), [], [])
    rep, _, _ = generic_plane_cnf(91, lines9, 2, 5, 26, 9)
    rep.add(1 + origin)                               # planar_cnf ends with the origin unit
    check(rep.clauses == ref.clauses and rep.nvars == ref.nvars, "replica differs from rebuild_all.planar_cnf")
    print("  replica clause stream == rebuild_all.planar_cnf (PG(2,9), m=26, %d clauses): identical" % len(ref.clauses), flush=True)
    for q, a, b, m, label in ((3, 1, 4, 4, "PG(2,3) type (1,4), m=4 (the lines)"), (4, 0, 2, 6, "PG(2,4) type (0,2), m=6 (hyperovals)"), (4, 3, 5, 15, "PG(2,4) type (3,5), m=15 (complements of hyperovals)")):
        n, lines = singer_plane(q)
        cnf, x, thru = generic_plane_cnf(n, lines, a, b, m, q)
        up = UP(cnf.nvars, cnf.clauses)
        sizes = list(range(n + 1)) if q == 3 else [m]                       # PG(2,3): every subset; PG(2,4): every subset of size m ...
        t = time.time()
        brute = 0
        tested = 0
        subsets = [S for sz in sizes for S in itertools.combinations(range(n), sz)]
        if q == 4:                                                          # ... and 3000 random subsets of each size m-1, m+1 (the neighbouring sizes)
            for sz in (m - 1, m + 1):
                subsets += [tuple(sorted(rnd.sample(range(n), sz))) for _ in range(300 if QUICK else 3000)]
        for S in subsets:
            Ss = set(S)
            want = len(S) == m and all(sum(1 for p in L if p in Ss) in (a, b) for L in lines)
            got = decide(up, x, [p in Ss for p in range(n)])
            check(got == want, "%s: set %s brute force %s CNF %s" % (label, S, want, got))
            brute += want
            tested += 1
        print("  %-52s: %d subsets tested (all subsets of size %s, plus random neighbours for q = 4): CNF == brute force for every one; %d genuine sets (%.0fs)" % (label, tested, sizes if q == 4 else "0..13", brute, time.time() - t), flush=True)
        if q == 3:
            check(tested == 2 ** 13, "all subsets")
        if q == 4:
            check(brute == 168, "number of hyperovals / complements in PG(2,4) must be 168, got %d" % brute)


# ---------------------------------------------------------------- E6: lex gadget
LEX_SRC = """    for h in H:
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
"""


def lex_gadget(cnf, x, H, order, ident):
    """VERBATIM copy of the lex-leader loop of rebuild_all.planar_cnf (the text is compared with the source below)"""
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


def e6():
    print("E6: lex-leader clause gadget")
    src = open(os.path.join(HERE, "rebuild_all.py")).read()
    check(LEX_SRC in src, "the lex loop of rebuild_all.py is not the text tested here")
    rnd = random.Random(3)
    n_inst = 0
    for trial in range(30 if QUICK else 120):
        n = rnd.choice([8, 10, 12])
        order = sorted(rnd.sample(range(n), rnd.randint(n - 3, n)))             # some points not compared (like the points of l0)
        H = []
        for _ in range(rnd.randint(1, 4)):
            p = list(range(n))
            rnd.shuffle(p)
            if rnd.random() < 0.5:                                              # make it sparse (many fixed points), as most elements of H_R are not
                p = list(range(n))
                for _ in range(rnd.randint(1, 3)):
                    i, j = rnd.sample(range(n), 2)
                    p[i], p[j] = p[j], p[i]
            H.append(p)
        cnf = R.CNF()
        x = [cnf.new_var() for _ in range(n)]
        lex_gadget(cnf, x, H, order, list(range(n)))
        up = UP(cnf.nvars, cnf.clauses)
        for bits in itertools.product((False, True), repeat=n):
            want = all(next((bits[p] < bits[h[p]] for p in order if bits[p] != bits[h[p]]), True) for h in H)
            got = decide(up, x, bits, search=True)
            check(got == want, "lex gadget: input %s perms %s order %s: expected %s got %s" % (bits, H, order, want, got))
        n_inst += 1
    print("  %d random instances (n in {8,10,12}, 1-4 permutations, partial compare order), ALL 2^n inputs: accepted inputs == {x : x <=_lex x o h for all h}" % n_inst, flush=True)


# ---------------------------------------------------------------- E7: real sizes, sampled
def e7():
    print("E7: plain sum constraints at the real sizes (SAMPLED)")
    rnd = random.Random(11)
    shapes = [(91, 26), (91, 35)] + ([] if QUICK else [(820, 260), (820, 287)])
    for n, k in shapes:
        cnf = R.CNF()
        ins = [cnf.new_var() for _ in range(n)]
        R.sum_in_set(cnf, ins, [1] * n, {k}, amo=False)
        t = time.time()
        up = UP(cnf.nvars, cnf.clauses)
        trials = 100 if n == 91 else 3
        pos = neg = 0
        for _ in range(trials):
            S = set(rnd.sample(range(n), k))
            check(decide(up, ins, [p in S for p in range(n)]), "valid input rejected")
            pos += 1
            for w in (k - 2, k - 1, k + 1, k + 2):
                T = set(rnd.sample(range(n), w))
                check(not decide(up, ins, [p in T for p in range(n)]), "invalid input accepted (weight %d)" % w)
                neg += 1
        for w in (0, n):
            T = set(range(w))
            check(not decide(up, ins, [p in T for p in range(n)]), "extreme input accepted")
            neg += 1
        print("  [1]*%d = %d (%d vars, %d clauses): %d random inputs of weight %d extend (unique total model satisfying every clause), %d inputs of weight %d+-1,+-2, 0, %d rejected (%.0fs)" % (
            n, k, cnf.nvars, len(cnf.clauses), pos, k, neg, k, n, time.time() - t), flush=True)
        del up


def main():
    t0 = time.time()
    for name, fn in (("E1", e1), ("E3", e3), ("E4", e4), ("E5", e5), ("E6", e6), ("E7", e7)):
        if want_test(name) or (name == "E1" and want_test("E2")):
            fn()
    print("ENCODER TESTS: ALL PASS (%.0fs)%s" % (time.time() - t0, "  [--quick]" if QUICK else ""))


if __name__ == "__main__":
    main()

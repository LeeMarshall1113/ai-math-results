"""CNF builder and exact encoders for weighted-sum constraints (stdlib only).

CNF            : clause list with fresh-variable allocation and DIMACS output
sum_in_set     : sum_i w_i * [lit_i] in A   (DP/BDD encoding, one variable per
                 reachable-and-co-reachable (layer, partial sum) node, forward and
                 backward channeling clauses, pairwise AMO per layer optional).
                 Returns the dict {final_sum: state variable} so callers can use
                 "the sum equals s" literals (e.g. y_l = [sum = 5]).
sum_equals     : special case A = {k}.

The encoders are validated exhaustively in `selftest()` : for random weight
vectors the projection of the CNF models on the input literals is compared with
the intended predicate (all 2^n input assignments, one SAT-free unit-propagation
completeness check through brute-force over the auxiliary variables is avoided
by extending each input assignment through forward propagation and checking
that no clause is violated and that the instance is UNSAT when the predicate is
false -- see the test code).
"""
import itertools
import random
import sys


class CNF:
    def __init__(self):
        self.nvars = 0
        self.clauses = []

    def new_var(self):
        self.nvars += 1
        return self.nvars

    def add(self, *lits):
        self.clauses.append(tuple(lits))

    def add_clause(self, lits):
        self.clauses.append(tuple(lits))

    def to_dimacs(self, path, comments=()):
        with open(path, "w", newline="\n") as fh:
            for c in comments:
                fh.write("c %s\n" % c)
            fh.write("p cnf %d %d\n" % (self.nvars, len(self.clauses)))
            for cl in self.clauses:
                fh.write(" ".join(map(str, cl)) + " 0\n")


def sum_in_set(cnf, lits, weights, allowed, amo=True):
    """Encode  sum_i weights[i]*lits[i]  in `allowed`  (exactly, in both directions).
    Returns {s: var} for the final-layer states (the variable is true iff the sum is s)."""
    n = len(lits)
    allowed = set(allowed)
    top = max(allowed)
    # co-reachability: suffix[i] = sums achievable by items i..n-1 (capped at top)
    suffix = [None] * (n + 1)
    suffix[n] = {0}
    for i in range(n - 1, -1, -1):
        w = weights[i]
        suffix[i] = suffix[i + 1] | {s + w for s in suffix[i + 1] if s + w <= top}
    # forward reachability with co-reachability pruning
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
    if not layer[0]:       # infeasible constraint
        cnf.add()
        return {}
    cnf.add(var[0][0])
    for i in range(n):
        z, w = lits[i], weights[i]
        for s, v in var[i].items():
            # forward
            t0 = var[i + 1].get(s)
            t1 = var[i + 1].get(s + w)
            cnf.add(*([-v, z] + ([t0] if t0 else [])))
            cnf.add(*([-v, -z] + ([t1] if t1 else [])))
        for t, v in var[i + 1].items():
            # backward
            p0 = var[i].get(t)
            p1 = var[i].get(t - w)
            cnf.add(*([-v, z] + ([p0] if p0 else [])))
            cnf.add(*([-v, -z] + ([p1] if p1 else [])))
        if amo and i + 1 < n + 1:
            vs = list(var[i + 1].values())
            for a in range(len(vs)):
                for b in range(a + 1, len(vs)):
                    cnf.add(-vs[a], -vs[b])
    # last layer: at least one allowed state must hold (follows from forward clauses)
    return dict(var[n])


def sum_equals(cnf, lits, weights, value, amo=False):
    return sum_in_set(cnf, lits, weights, {value}, amo=amo)


# -------------------------------------------------------------------------------
def _eval_cnf_with_inputs(cnf, input_vars, assignment):
    """Complete the auxiliary variables by unit propagation + exhaustive search of what is
    left (tiny); return True iff some extension satisfies all clauses."""
    val = {}
    for v, b in zip(input_vars, assignment):
        val[v] = b
    clauses = cnf.clauses
    # simple DPLL on the remaining variables
    def solve(val):
        # unit propagate
        changed = True
        val = dict(val)
        while changed:
            changed = False
            for cl in clauses:
                unassigned = []
                sat = False
                for l in cl:
                    v = abs(l)
                    if v in val:
                        if (val[v] and l > 0) or ((not val[v]) and l < 0):
                            sat = True
                            break
                    else:
                        unassigned.append(l)
                if sat:
                    continue
                if not unassigned:
                    return False
                if len(unassigned) == 1:
                    l = unassigned[0]
                    val[abs(l)] = l > 0
                    changed = True
        # pick a free variable
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
    return solve(val)


def selftest(trials=60, seed=1):
    rnd = random.Random(seed)
    for t in range(trials):
        n = rnd.randint(1, 7)
        weights = [rnd.randint(1, 4) for _ in range(n)]
        allowed = set(rnd.sample(range(0, sum(weights) + 2), rnd.randint(1, 3)))
        amo = rnd.random() < 0.5
        cnf = CNF()
        ins = [cnf.new_var() for _ in range(n)]
        sum_in_set(cnf, ins, weights, allowed, amo=amo)
        for assignment in itertools.product((False, True), repeat=n):
            s = sum(w for w, b in zip(weights, assignment) if b)
            want = s in allowed
            got = _eval_cnf_with_inputs(cnf, ins, assignment)
            assert want == got, (weights, allowed, assignment, s, want, got, amo)
    return True


if __name__ == "__main__":
    print("selftest", selftest())

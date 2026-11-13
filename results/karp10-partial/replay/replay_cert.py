"""replay_cert.py -- verify a Karp n=10 exclusion certificate from scratch (standard library only).

usage:  python replay_cert.py certs/type4.json [-v]

A certificate states: for the given Angel-Borja type T, the listed fixed-point equations (for subgroups Gamma
of S10 in Oliver's class G^1, on Delta_P, or on the link / deletion of a Gamma-fixed edge) cannot all hold
for any monotone P of type T.  The script re-derives every ingredient:
  1. the nine types (Sylow-5 equation) -- replay_lib.derive_types();
  2. Gamma in G^1:   Gamma = <gens>, N = <ngens> normal p-subgroup, Gamma = <g> N;
  3. the edge orbits of Gamma, the unions of orbits, their isomorphism classes;
  4. each union's status for type T: IN (embeds in a member of the type: witness permutation checked),
     OUT (contains a non-member special graph or K10: witness checked) or FREE;
  5. the Euler sums  sum_{T != 0} (-1)^{|T|-1} x[union]  (fix), or over unions containing/avoiding the fixed
     edge e (lk / del);  exact:  chi = 1;
  6. exhaustive search over all 0/1 values of the FREE classes, respecting monotonicity
     (x[a] <= x[b] whenever b embeds in a), for an assignment satisfying every listed equation.
It prints NO SOLUTION (certificate valid) or lists a solution (certificate invalid).
"""
import json, sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from replay_lib import *


def parity_sign(r):
    return 1 if r % 2 == 1 else -1



def solve_dfs(free, imps, comp, forced, limit=3):
    """DFS with fixpoint propagation.  Variables = free classes (0/1); implications x[a] <= x[b];
    equations const + sum(w x) = 1.  Returns up to `limit` solutions (list of dicts)."""
    idx = {c: i for i, c in enumerate(free)}
    nv = len(free)
    assert len(idx) == nv, "duplicate class in the free list"
    assert all(c in idx for c in forced), "a forced assumption is not a free class (invariant I1 violated)"
    assert all(a in idx and b in idx for a, b in imps), "an implication mentions a non-free class"
    assert all(c in idx for _, _, coef in comp for c in coef), "an equation mentions a non-free class"
    ims = [(idx[a], idx[b]) for a, b in imps]
    eqs = [(const, [(idx[c], w) for c, w in coef.items()]) for _, const, coef in comp]
    sols = []

    def propagate(val):
        changed = True
        while changed:
            changed = False
            for a, b in ims:                      # x[a] <= x[b]
                if val[a] == 1:
                    if val[b] == 0:
                        return False
                    if val[b] == -1:
                        val[b] = 1
                        changed = True
                if val[b] == 0:
                    if val[a] == 1:
                        return False
                    if val[a] == -1:
                        val[a] = 0
                        changed = True
            for const, terms in eqs:
                lo = hi = const
                for u, w in terms:
                    if val[u] == -1:
                        if w > 0:
                            hi += w
                        else:
                            lo += w
                    else:
                        lo += w * val[u]
                        hi += w * val[u]
                if lo > 1 or hi < 1:
                    return False
                if lo == 1 or hi == 1:
                    for u, w in terms:
                        if val[u] == -1:
                            val[u] = (0 if w > 0 else 1) if lo == 1 else (1 if w > 0 else 0)
                            changed = True
        return True

    def rec(val):
        if len(sols) >= limit:
            return
        if not propagate(val):
            return
        try:
            v = val.index(-1)
        except ValueError:
            if all(const + sum(w * val[u] for u, w in terms) == 1 for const, terms in eqs):
                sols.append({free[i]: val[i] for i in range(nv)})
            return
        for choice in (0, 1):
            v2 = val[:]
            v2[v] = choice
            rec(v2)

    val0 = [-1] * nv
    for c, w in forced.items():
        val0[idx[c]] = w
    rec(val0)
    return sols


def clique_union(a, b):
    """the graph K_a u K_b on the ten vertices 0..9 (a + b = 10), built from its definition"""
    assert a >= 1 and b >= 1 and a + b == N
    return mask_of([(i, j) for i in range(a) for j in range(i + 1, a)] + [(i, j) for i in range(a, N) for j in range(i + 1, N)])


def main(path, verbose=False, declared=None):
    """declared = (a, b) from the command line flag --assume-out-cliques=a,b: the identity of the single assume_out graph (K_a u K_b), which the
    replay then CHECKS by an isomorphism test.  A certificate with assumptions and no declaration (or a declaration without assumptions) is refused."""
    cert = json.load(open(path))
    T = cert["type"]
    print(f"certificate {path}: type {T}; claim: {cert['claim']}")
    weight, below, types = derive_types()
    ab = {k: frozenset(v) for k, v in AB_IDEALS.items()}
    assert set(types) == set(ab.values()), "derived types differ from Angel-Borja"
    members = sorted(ab[T], key=SPEC10.index)
    print(f"  [1] nine types re-derived (Sylow-5, exact weights {weight}); type {T} = {members}")
    G = special_graphs()

    groups = {}
    for name, gd in cert["groups"].items():
        gens = [tuple(x) for x in gd["gens"]]
        order, n = check_G1(gens, gd["order"], gd["p"], [tuple(x) for x in gd["N_gens"]], tuple(gd["g"]))
        orbs = edge_orbits(gens)
        groups[name] = dict(gens=gens, orbs=orbs)
        print(f"  [2] {name}: |Gamma|={order}, normal {gd['p']}-subgroup N of order {n}, Gamma=<g>N  => Gamma in G^1;"
              f" edge orbits {[popcount(o) for o in orbs]}")

    classes = []     # list of dict(mask, status, witness)
    frozen = {"on": False}
    def classify(mask):
        for i, c in enumerate(classes):
            if isomorphic(c["mask"], mask):
                return i
        assert not frozen["on"], "invariant I1 violated: a new class would be created after the free list was built"
        st = status(mask, members, G)
        classes.append(dict(mask=mask, status=st))
        return len(classes) - 1

    eqs = []
    for eq in cert["equations"]:
        gname, kind = eq["group"], eq["kind"]
        gr = groups[gname]
        orbs = gr["orbs"]
        k = len(orbs)
        terms = {}      # class idx -> weight
        if kind == "fix":
            idxs = list(range(k))
            e_mask = 0
        else:
            ea, eb = eq["fixed_edge"]
            e_mask = 1 << eid(ea, eb)
            assert e_mask in orbs, f"{(ea, eb)} is not a Gamma-fixed edge of {gname}"
            idxs = [i for i in range(k) if orbs[i] != e_mask]
        for r in range(1, len(idxs) + 1):
            for sub in itertools.combinations(idxs, r):
                m = 0
                for i in sub:
                    m |= orbs[i]
                if kind == "lk":
                    m |= e_mask
                c = classify(m)
                terms[c] = terms.get(c, 0) + parity_sign(r)
        eqs.append((f"{gname}:{kind}" + (f"@{tuple(eq['fixed_edge'])}" if kind != "fix" else ""), terms))

    # ---- INVARIANT I1 (see DEPENDENCIES.md / the referee reply): every class that any equation or any assumption mentions is classified
    # BEFORE the list of free classes is built, and `classify` is frozen afterwards; so `free` is exactly the set of FREE classes of the final
    # class list, and every forced assumption is an element of it (or has been recorded as a contradiction).
    n_after_equations = len(classes)
    n_assump = len(cert.get("assume_out", [])) + len(cert.get("assume_in", []))
    if n_assump:
        assert declared is not None, "the certificate has assumptions but the identity of the assumed graph was not declared (--assume-out-cliques=a,b)"
        assert len(cert.get("assume_out", [])) == 1 and not cert.get("assume_in"), "expected exactly one assume_out graph and no assume_in"
        gm0 = mask_of([tuple(x) for x in cert["assume_out"][0]])
        assert isomorphic(gm0, clique_union(*declared)), f"the assumed graph is not K_{declared[0]} u K_{declared[1]}"
        print(f"  [assume-identity] the assumed graph (edges {popcount(gm0)}) is isomorphic to K_{declared[0]} u K_{declared[1]} (built from its definition)")
    else:
        assert declared is None, "an assumption identity was declared but the certificate has no assumptions"
    assumed = []     # (kind, want, class index), classified now, resolved below
    for kind, edges in [("out", e) for e in cert.get("assume_out", [])] + [("in", e) for e in cert.get("assume_in", [])]:
        gm = mask_of([tuple(x) for x in edges])
        assumed.append((kind, 0 if kind == "out" else 1, gm, classify(gm)))
    frozen["on"] = True
    n_classes_at_freeze = len(classes)

    print(f"  [3,4] {len(classes)} distinct union classes (equations and assumptions); statuses:")
    for i, c in enumerate(classes):
        st = c["status"]
        extra = f" via {st[1]}" if st[0] != "FREE" else ""
        if verbose or st[0] == "FREE":
            print(f"       class {i}: {popcount(c['mask'])} edges, {st[0]}{extra}")
    free = [i for i, c in enumerate(classes) if c["status"][0] == "FREE"]
    print(f"      free classes: {free}")

    # monotone implications among free classes: x[a] <= x[b] when b embeds in a, a != b
    imps = [(a, b) for a in free for b in free if a != b and popcount(classes[b]['mask']) < popcount(classes[a]['mask'])
            and embeds(classes[b]["mask"], classes[a]["mask"])]
    print(f"      implications x[a] <= x[b] among free classes: {imps}")

    # equation constants / coefficients
    comp = []
    for name, terms in eqs:
        const = sum(w for c, w in terms.items() if classes[c]["status"][0] == "IN")
        coef = {c: w for c, w in terms.items() if classes[c]["status"][0] == "FREE" and w != 0}
        comp.append((name, const, coef))
        print(f"  [5] {name}: forced-IN part = {const}, free terms = {coef}; equation: {const} + sum(free) = 1")

    # constraints: assumptions on specific graphs (given explicitly in the certificate)
    free_set = set(free)
    forced = {}
    contradiction = False
    for kind, want, gm, c in assumed:
        st_c = classes[c]["status"][0]
        if st_c == "FREE":
            assert c in free_set, "invariant I1 violated: a forced assumption is not in the free list"
            if forced.get(c, want) != want:
                print("  two assumptions put the same class both in and out of P: inconsistent")
                contradiction = True
            forced[c] = want
        else:
            val = 1 if st_c == "IN" else 0
            if val != want:
                print(f"  assumption on a {st_c} graph contradicts the type itself: trivially inconsistent")
                contradiction = True
        print(f"  [assume] graph with {popcount(gm)} edges assumed {'NOT in P' if want == 0 else 'in P'} (status for the type: {st_c}; class {c}, "
              f"{'already produced by an equation' if c < n_after_equations else 'new, created by the assumption'})")
    assert len(classes) == n_classes_at_freeze, "invariant I1 violated: a class was created after the free list"
    if contradiction:
        sols = []
        free = []
    else:
        sols = None
    if sols is None:
        sols = solve_dfs(free, imps, comp, forced)
        assert len(classes) == n_classes_at_freeze
    if sols:
        print("  RESULT: solutions exist -> certificate does NOT exclude the type:", sols[:3])
        return 1
    if forced:
        print(f"  [6] complete propagating DFS over {len(free)} free classes: NO SOLUTION  ==> type {T} together with the stated assumption(s) is impossible: "
              f"the assumed membership statement is FALSE for every non-evasive P of type {T} (certificate valid)")
    else:
        print(f"  [6] complete propagating DFS over {len(free)} free classes: NO SOLUTION  ==> type {T} is excluded (certificate valid)")
    return 0


if __name__ == "__main__":
    decl = None
    for a in sys.argv[2:]:
        if a.startswith("--assume-out-cliques="):
            decl = tuple(int(x) for x in a.split("=", 1)[1].split(","))
        else:
            assert a == "-v", f"unknown argument {a}"
    sys.exit(main(sys.argv[1], "-v" in sys.argv, decl))

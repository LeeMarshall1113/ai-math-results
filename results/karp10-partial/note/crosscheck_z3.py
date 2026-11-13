"""crosscheck_z3.py -- second opinion on the exhaustive search in replay_cert.py (needs z3; NOT part of the standard-library replay).

For every certificate in ../certs that lists groups and equations (types 1-4, 6-9 and the four Lemma A certificates) this script
  1. rebuilds the 0/1 system again, with its own code, from the group generators (the primitives closure / edge_orbits / embeds /
     isomorphic of replay_lib are reused; the table of classes, statuses, equations and the search are written separately),
  2. imposes ALL monotonicity constraints x[a] <= x[b] (b a subgraph of a) between every pair of classes, not only between free
     classes as replay_cert.py does,
  3. asks z3 whether the system has a 0/1 solution:  certificate valid  <=>  UNSAT,
  4. runs a negative control: the same groups and equations against type 5 (and, for the Lemma A certificates, the same system
     without the assumption on K_a u K_{10-a}); these must be satisfiable, otherwise the solver step proves nothing.

usage:  python note/crosscheck_z3.py            (from the karp10 directory; needs z3-solver; under a minute)
"""
import glob
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "replay"))
from replay_lib import *  # noqa: E402,F401,F403
import z3  # noqa: E402

G = special_graphs()


def build(cert, type_override=None):
    T = type_override or cert["type"]
    members = AB_IDEALS[T]
    groups = {}
    for name, gd in cert["groups"].items():
        gens = [tuple(x) for x in gd["gens"]]
        check_G1(gens, gd["order"], gd["p"], [tuple(x) for x in gd["N_gens"]], tuple(gd["g"]))
        groups[name] = (gens, edge_orbits(gens))
    classes = []          # representative masks
    status = []           # 1 = in, 0 = out, None = free

    def cls(mask):
        for i, r in enumerate(classes):
            if popcount(r) == popcount(mask) and embeds(r, mask):
                return i
        st = None
        if any(embeds(mask, G[s]) for s in members):
            st = 1
        elif any(embeds(G[s], mask) for s in SPEC10 if s not in members) or popcount(mask) == M:
            st = 0
        classes.append(mask)
        status.append(st)
        return len(classes) - 1

    eqs = []
    for eq in cert["equations"]:
        gens, orbs = groups[eq["group"]]
        k = len(orbs)
        if eq["kind"] == "fix":
            idx, e_mask = list(range(k)), 0
        else:
            ea, eb = eq["fixed_edge"]
            e_mask = 1 << eid(ea, eb)
            assert e_mask in orbs
            idx = [i for i in range(k) if orbs[i] != e_mask]
        coef = {}
        for r in range(1, len(idx) + 1):
            for sub in itertools.combinations(idx, r):
                m = 0
                for i in sub:
                    m |= orbs[i]
                if eq["kind"] == "lk":
                    m |= e_mask
                c = cls(m)
                coef[c] = coef.get(c, 0) + (1 if r % 2 == 1 else -1)
        eqs.append(coef)
    assumed = []
    for key, val in (("assume_out", 0), ("assume_in", 1)):
        for edges in cert.get(key, []):
            assumed.append((cls(mask_of([tuple(x) for x in edges])), val))
    return classes, status, eqs, assumed


def solve(classes, status, eqs, assumed, use_assumptions=True):
    n = len(classes)
    x = [z3.Int(f"x{i}") for i in range(n)]
    s = z3.Solver()
    for i in range(n):
        s.add(x[i] >= 0, x[i] <= 1)
        if status[i] is not None:
            s.add(x[i] == status[i])
    for a in range(n):
        for b in range(n):
            if a != b and popcount(classes[b]) < popcount(classes[a]) and embeds(classes[b], classes[a]):
                s.add(x[a] <= x[b])                       # b is a subgraph of a
    for coef in eqs:
        s.add(z3.Sum([c * x[i] for i, c in coef.items()]) == 1)
    if use_assumptions:
        for i, v in assumed:
            s.add(x[i] == v)
    return s.check()


def main():
    t0 = time.time()
    bad = 0
    for f in sorted(glob.glob(os.path.join(ROOT, "certs", "*.json"))):
        cert = json.load(open(f))
        if "groups" not in cert:
            continue
        name = os.path.basename(f)
        classes, status, eqs, assumed = build(cert)
        r = solve(classes, status, eqs, assumed)
        nfree = sum(1 for v in status if v is None)
        # negative control
        if assumed:
            ctrl = solve(classes, status, eqs, assumed, use_assumptions=False)
            ctrl_txt = f"control (no assumption): {ctrl}"
            ok_ctrl = (ctrl == z3.sat)
        else:
            c2, s2, e2, a2 = build(cert, type_override=5)
            ctrl = solve(c2, s2, e2, a2)
            ctrl_txt = f"control (same equations against type 5): {ctrl}"
            ok_ctrl = (ctrl == z3.sat)
        ok = (r == z3.unsat) and ok_ctrl
        bad += 0 if ok else 1
        print(f"{name:26s} classes={len(classes):4d} free={nfree:3d}  z3: {r}   {ctrl_txt}   {'OK' if ok else 'PROBLEM'}", flush=True)
    print(f"z3 cross-check of all Euler-type certificates: {'ALL CONFIRMED' if bad == 0 else 'PROBLEMS: %d' % bad}  ({time.time() - t0:.0f} s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

"""check_z5_absolute.py [certs/lemma_z5_euler.json] -- the EXACT acyclicity test behind Lemma Z5 (standard library only; run by replay_all.py).

replay_smith.py uses the RELATIVE-HOMOLOGY FILTER  H_s(K, K_D) = H~_{s-1}(K_D)  to discard down-sets X of the CC-class poset.  That filter is only
NECESSARY for F_p-acyclicity of K = K_D u {faces whose union class is in X} (equal dimensions in the long exact sequence of the pair do not force the
connecting maps to be isomorphisms; example: K_D = four points, K = a triangle boundary plus the fourth point), so its survivors are a SUPERSET of the
admissible patterns, which is all that the soundness of the certificate needs.

This script removes the doubt in the other direction.  For EVERY down-set X (2,539 of them) it computes the reduced F_p-homology of the ABSOLUTE
complex K (augmented, with the empty face) and keeps X iff it vanishes.  It then
  1. asserts that the set of absolutely acyclic X equals the set of survivors of the relative filter (so for this complex the filter is exact);
  2. applies the Euler equations (class G^1, verified) of the extra groups of the certificate to the absolutely acyclic X;
  3. checks every claimed-in / claimed-out graph of the certificate against the resulting patterns;
  4. asserts the counts pinned in the certificate key "expect".
So "exactly 16 patterns, 8 after two Euler equations, nine graphs in and nine out" is replayed from absolute acyclicity, not from the filter.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from replay_lib import *  # noqa: F401,F403
from replay_smith import connected, homology_by_size  # noqa: E402


def main(path):
    cert = json.load(open(path))
    p = cert["prime"]
    assert is_prime(p), f"prime = {p} is not a prime"
    gens = [tuple(g) for g in cert["gens"]]
    n_grp = len(closure(gens))
    t = n_grp
    while t % p == 0:
        t //= p
    assert t == 1 and n_grp > 1, f"group order {n_grp} is not a power of {p}"
    orbs = edge_orbits(gens)
    k = len(orbs)
    unions = {}
    for T in range(1 << k):
        m = 0
        for i in range(k):
            if (T >> i) & 1:
                m |= orbs[i]
        unions[T] = m
    status = {T: ("D" if not connected(m) else ("N" if not connected(complement(m)) else "CC")) for T, m in unions.items() if T}
    classes, cls_of = [], {}
    for T, s in status.items():
        if s != "CC":
            continue
        for ci, rep in enumerate(classes):
            if isomorphic(rep, unions[T]):
                cls_of[T] = ci
                break
        else:
            classes.append(unions[T])
            cls_of[T] = len(classes) - 1
    nC = len(classes)
    below = [[j for j in range(nC) if j != i and popcount(classes[j]) < popcount(classes[i]) and embeds(classes[j], classes[i])] for i in range(nC)]
    dfaces = [T for T, s in status.items() if s == "D"]
    fbc = {}
    for T, c in cls_of.items():
        fbc.setdefault(c, []).append(T)
    order = sorted(range(nC), key=lambda i: popcount(classes[i]))
    good_abs, good_rel, total = [], [], [0]
    dd = {0: [0]}
    for T in dfaces:
        dd.setdefault(popcount(T), []).append(T)
    expected = {s + 1: v for s, v in homology_by_size(dd, p).items() if v}      # the filter: H_s(relative) must equal reduced H_{s-1}(K_D)

    def rec(pos, cur):
        if pos == nC:
            total[0] += 1
            faces = {0: [0]}
            for T in dfaces:
                faces.setdefault(popcount(T), []).append(T)
            for c in cur:
                for T in fbc[c]:
                    faces.setdefault(popcount(T), []).append(T)
            H = {s: v for s, v in homology_by_size(faces, p).items() if v}
            if not H:
                good_abs.append(frozenset(cur))
            rel = {}
            for c in cur:
                for T in fbc[c]:
                    rel.setdefault(popcount(T), []).append(T)
            Hr = {s: v for s, v in homology_by_size(rel, p).items() if v} if rel else {}
            if Hr == expected:
                good_rel.append(frozenset(cur))
            return
        i = order[pos]
        rec(pos + 1, cur)
        if all(j in cur for j in below[i]):
            rec(pos + 1, cur | {i})

    rec(0, frozenset())
    print(f"  [1] {nC} CC classes; {total[0]} down-sets; {len(good_abs)} have vanishing ABSOLUTE reduced F_{p}-homology; {len(good_rel)} pass the relative filter")
    assert set(good_abs) == set(good_rel), "absolute acyclicity and the relative filter select different patterns"
    print("      the two criteria select the same patterns (for this complex the necessary filter is exact)")
    exp = cert["expect"]       # REQUIRED: the pinned counts are part of the certificate; deleting them is an error, not a relaxation
    assert ("after_euler" in exp) == bool(cert.get("euler")), "the key expect.after_euler must be present exactly when the certificate has Euler groups"
    assert "downsets" in exp and "relative_survivors" in exp
    if "downsets" in exp:
        assert total[0] == exp["downsets"], f"expected {exp['downsets']} down-sets, found {total[0]}"
    if "relative_survivors" in exp:
        assert len(good_abs) == exp["relative_survivors"], f"expected {exp['relative_survivors']} acyclic patterns, found {len(good_abs)}"
    good = list(good_abs)
    euler_counts = []
    for ent in cert.get("euler", []):
        gam = [tuple(x) for x in ent["gens"]]
        check_G1(gam, ent["order"], ent["p"], [tuple(x) for x in ent["N_gens"]], tuple(ent["g"]))
        gorbs = edge_orbits(gam)
        gk = len(gorbs)
        gtab = {}
        for T in range(1, 1 << gk):
            m = 0
            for i in range(gk):
                if (T >> i) & 1:
                    m |= gorbs[i]
            if not connected(m):
                gtab[T] = ("D", None)
            elif not connected(complement(m)):
                gtab[T] = ("N", None)
            else:
                ci = next((c for c, rep in enumerate(classes) if isomorphic(rep, m)), None)
                assert ci is not None, "a CC union of the Euler group is not a union class of Q"
                gtab[T] = ("CC", ci)
        keep = []
        for pat in good:
            val = sum((-1) ** (popcount(T) - 1) * (1 if st == "D" else (0 if st == "N" else (1 if ci in pat else 0))) for T, (st, ci) in gtab.items())
            if val == 1:
                keep.append(pat)
        print(f"  [2] Euler equation of a group of order {ent['order']} ({gk} orbits, G^1 verified): {len(keep)} of {len(good)} absolutely acyclic patterns survive")
        euler_counts.append(len(keep))
        good = keep
    if "after_euler" in exp:
        assert euler_counts == exp["after_euler"], f"expected {exp['after_euler']}, found {euler_counts}"
    assert good
    ints = set.intersection(*[set(s) for s in good])
    print(f"  [3] classes in all {len(good)} patterns: {len(ints)} (edge counts {sorted(popcount(classes[i]) for i in ints)})")
    ok = True
    for name, key, want in (("IN", "claims_in", True), ("OUT", "claims_out", False)):
        for edges in cert.get(key, []):
            gm = mask_of([tuple(e) for e in edges])
            ci = next((c for c, rep in enumerate(classes) if isomorphic(rep, gm)), None)
            holds = ci is not None and all((ci in pat) == want for pat in good)
            ok = ok and holds
            print(f"  [4] {name}: graph with {popcount(gm)} edges (class {ci}) {'in all' if want else 'in none of'} the {len(good)} patterns? {holds}")
    print("  RESULT:", "ABSOLUTE ACYCLICITY CONFIRMS THE CERTIFICATE" if ok else "CERTIFICATE NOT CONFIRMED")
    return 0 if ok else 1


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "certs", "lemma_z5_euler.json")))

"""replay_nonselfdual.py <certs/lemma_z5_euler.json> -- standard-library-only check that a non-evasive monotone P of type 5 is NEVER self-dual (P != P*,
with P* = {G : complement(G) not in P}).

What is re-derived here (same code path as replay_smith.py, no data from the search pipeline):
  1. Q = <gens> is a p-group (its order is a power of the prime p: checked), with its edge orbits;
  2. the 32 isomorphism classes of CC unions of Q-orbits (connected graphs with connected complement);
  3. the complement map on the 32 classes: a fixed-point-free involution (no graph on 10 vertices is self-complementary, since 45 is odd);
  4. every down-set X of the CC-class poset that passes the RELATIVE-HOMOLOGY FILTER of replay_smith.py (H_s(K, K_D) = H~_{s-1}(K_D), a NECESSARY
     condition for F_p-acyclicity of the fixed complex K, not an equivalent one: replay/check_z5_absolute.py tests the whole complex and selects
     the same patterns) and then the Euler equations of the extra Oliver groups of the certificate (each verified to be in G^1);
  5. the conclusion.  If P = P*, then for every class c:  c in X  iff  complement(c) not in X, so X contains exactly ONE class of each
     complementary pair and |X| = 32/2 = 16.  The trace of every admissible P is among the survivors of step 4 (a SUPERSET of the truly
     admissible patterns), and none of the survivors has 16 classes.  Hence P != P*.
     Consistency checks also done: no survivor is fixed by the duality map X -> {c : complement(c) not in X}, and the survivors are closed under it
     (P* is again non-evasive of type 5).
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from replay_lib import *            # noqa
from replay_smith import connected, homology_by_size


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
    status_ = {}
    for T, m in unions.items():
        if T:
            status_[T] = "D" if not connected(m) else ("N" if not connected(complement(m)) else "CC")
    classes, cls_of = [], {}
    for T, s in status_.items():
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
    dfaces = {0: [0]}
    for T, s in status_.items():
        if s == "D":
            dfaces.setdefault(popcount(T), []).append(T)
    hD = {s: v for s, v in homology_by_size(dfaces, p).items() if v}
    expected = {s + 1: v for s, v in hD.items()}
    faces_by_class = {}
    for T, c in cls_of.items():
        faces_by_class.setdefault(c, []).append(T)
    order = sorted(range(nC), key=lambda i: popcount(classes[i]))
    good = []

    def rec(pos, cur):
        if pos == nC:
            faces = {}
            for c in cur:
                for T in faces_by_class[c]:
                    faces.setdefault(popcount(T), []).append(T)
            H = {s: v for s, v in homology_by_size(faces, p).items() if v} if faces else {}
            if H == expected:
                good.append(frozenset(cur))
            return
        i = order[pos]
        rec(pos + 1, cur)
        if all(j in cur for j in below[i]):
            rec(pos + 1, cur | {i})

    rec(0, frozenset())
    print(f"  [1] Q has order {n_grp} (a {p}-group); {nC} CC classes; {len(good)} down-sets pass the relative-homology filter (necessary for F_{p}-acyclicity)")
    exp = cert["expect"]       # REQUIRED: the pinned counts are part of the certificate; deleting them is an error, not a relaxation
    assert ("after_euler" in exp) == bool(cert.get("euler")), "the key expect.after_euler must be present exactly when the certificate has Euler groups"
    assert "relative_survivors" in exp
    if "relative_survivors" in exp:
        assert len(good) == exp["relative_survivors"], f"expected {exp['relative_survivors']} survivors, found {len(good)}"
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
                ci = next(c for c, rep in enumerate(classes) if isomorphic(rep, m))
                gtab[T] = ("CC", ci)
        keep = []
        for pat in good:
            val = sum((-1) ** (popcount(T) - 1) * (1 if st == "D" else (0 if st == "N" else (1 if ci in pat else 0))) for T, (st, ci) in gtab.items())
            if val == 1:
                keep.append(pat)
        print(f"  [2] Euler equation of a group of order {ent['order']} ({gk} orbits, G^1 verified): {len(keep)} of {len(good)} patterns survive")
        euler_counts.append(len(keep))
        good = keep
    if "after_euler" in exp:
        assert euler_counts == exp["after_euler"], f"expected {exp['after_euler']}, found {euler_counts}"
    # complement map on the classes
    cmap = {}
    for c in range(nC):
        cm = complement(classes[c])
        cmap[c] = next(d for d in range(nC) if isomorphic(classes[d], cm))
    assert all(cmap[cmap[c]] == c for c in range(nC)) and all(cmap[c] != c for c in range(nC))
    assert nC % 2 == 0
    print(f"  [3] the complement map is a fixed-point-free involution on the {nC} classes (no self-complementary graph on 10 vertices: 45 is odd)")
    sizes = sorted(len(p_) for p_ in good)
    half = nC // 2
    print(f"  [4] a self-dual trace would contain exactly one class of each complementary pair, i.e. exactly {half} of the {nC} classes;")
    print(f"      the {len(good)} survivors have sizes {sizes}")
    size_ok = half not in sizes
    # consistency checks: duality map
    dual = lambda pat: frozenset(c for c in range(nC) if cmap[c] not in pat)
    good_set = set(good)
    selfdual = 0
    for pat in good:
        d = dual(pat)
        assert d in good_set, "the surviving patterns are not closed under duality"
        if d == pat:
            selfdual += 1
    print(f"  [5] consistency: the survivors are closed under duality (P* is again non-evasive of type 5); survivors fixed by duality: {selfdual}")
    ok = size_ok and selfdual == 0
    print("  RESULT:", "NO SURVIVING PATTERN IS SELF-DUAL: a non-evasive type-5 property is never self-dual" if ok else "a self-dual pattern exists")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

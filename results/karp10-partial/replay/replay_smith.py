"""replay_smith.py -- verify a 'Smith-local' certificate from scratch (standard library only).

Certificate: a p-subgroup Q of S10 (generators), the prime p, a list of graphs claimed to lie in P and a list claimed NOT to lie in P.

Premises (cited, not re-proved here):
  (S1) a non-evasive nontrivial monotone P on 10 vertices of Angel-Borja type 5 satisfies D <= P <= D*
       (certified by certs/type5_force_K*_K*.json), so P = D u X with X a down-set of connected and co-connected (CC) graphs;
  (S2) Delta_P is Z-acyclic (Kahn-Saks-Sturtevant), hence by Smith theory the fixed complex Delta_P^Q (vertices = Q-orbits of edges,
       faces = sets of orbits whose union graph is in P) is F_p-acyclic for the p-group Q.
Computation (all re-done here): orbits of Q; the union graph of every set of orbits; its status D / N / CC; isomorphism classes of the CC
unions (embedding + equal edge count); the embedding order among CC classes; the F_p-homology of the D-part K_D; then EXHAUSTIVE
enumeration of every down-set X of the CC classes (monotone in the embedding order) and the RELATIVE-HOMOLOGY FILTER: the homology of the relative
complex of the X-faces must equal the reduced homology of K_D shifted by one, H_s(K, K_D) = H~_{s-1}(K_D) for all s.
  This filter is NECESSARY for F_p-acyclicity of K = K_D u {X-faces} (long exact sequence of the pair), but NOT sufficient in general: equal
  dimensions do not force the connecting maps to be isomorphisms (example: K_D = four points, K = a triangle boundary plus the fourth point
  has H_1(K, K_D) = F^3 = H~_0(K_D) and is not acyclic).  The surviving X are therefore a SUPERSET of the patterns that really give an acyclic K.
The certificate is valid iff every surviving X contains every claimed-in class and no claimed-out class; this is sound for any superset of the
admissible patterns, so the necessary filter suffices.  The exact statement "these are exactly the patterns with F_p-acyclic K" is replayed
separately by replay/check_z5_absolute.py (homology of the whole complex), which replay_all.py runs and requires to select the same set.

Optional certificate key "expect" = {"downsets": N, "relative_survivors": M, "after_euler": [m1, m2, ...]} pins the counts printed below;
the replay asserts them, so the counts quoted in the notes are replayed, not just reported.

usage: python replay_smith.py certs/lemma_z5.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from replay_lib import *  # noqa: E402,F401,F403


def connected(mask):
    rows = rows_of(mask)
    seen, fr = 1, 1
    while fr:
        v = (fr & -fr).bit_length() - 1
        fr &= fr - 1
        nb = rows[v] & ~seen
        seen |= nb
        fr |= nb
    return seen == (1 << N) - 1


def rank_mod(rows, p):
    piv = {}
    rank = 0
    for r in rows:
        r = {c: v % p for c, v in r.items() if v % p}
        while r:
            c = min(r)
            if c in piv:
                pr = piv[c]
                f = r[c] * pow(pr[c], -1, p) % p
                for cc, vv in pr.items():
                    nv = (r.get(cc, 0) - f * vv) % p
                    if nv:
                        r[cc] = nv
                    else:
                        r.pop(cc, None)
            else:
                piv[c] = r
                rank += 1
                break
    return rank


def homology_by_size(faces_by_size, p):
    """dim H_s over F_p of the chain complex whose s-chains are spanned by the given faces (bitmasks of size s); faces whose boundary
    facet is missing from the family simply lose that term (relative boundary).  Face of size 0 (mask 0) = augmentation."""
    idx = {s: {t: i for i, t in enumerate(fs)} for s, fs in faces_by_size.items()}
    ranks = {}
    for s, fs in faces_by_size.items():
        if s == 0 or (s - 1) not in faces_by_size:
            ranks[s] = 0
            continue
        tgt = idx[s - 1]
        rows = []
        for t in fs:
            row = {}
            bits = [b for b in range(32) if (t >> b) & 1]
            for pos, b in enumerate(bits):
                j = tgt.get(t ^ (1 << b))
                if j is not None:
                    row[j] = (-1) ** pos
            rows.append(row)
        ranks[s] = rank_mod(rows, p)
    return {s: len(fs) - ranks[s] - ranks.get(s + 1, 0) for s, fs in faces_by_size.items()}


def main(path):
    cert = json.load(open(path))
    p = cert["prime"]
    assert is_prime(p), f"prime = {p} is not a prime"
    gens = [tuple(g) for g in cert["gens"]]
    G = closure(gens)
    n = len(G)
    t = n
    while t % p == 0:
        t //= p
    assert t == 1 and n > 1, f"group order {n} is not a power of {p}"
    orbs = edge_orbits(gens)
    k = len(orbs)
    print(f"certificate {path}: {cert['claim']}")
    print(f"  [1] Q = <gens> has order {n} (a {p}-group); {k} edge orbits of sizes {sorted(popcount(o) for o in orbs)}")
    unions = {}
    for T in range(1 << k):
        m = 0
        for i in range(k):
            if (T >> i) & 1:
                m |= orbs[i]
        unions[T] = m
    status = {}
    for T, m in unions.items():
        if T == 0:
            continue
        status[T] = "D" if not connected(m) else ("N" if not connected(complement(m)) else "CC")
    cnt = {s: sum(1 for v in status.values() if v == s) for s in ("D", "N", "CC")}
    print(f"  [2] {len(status)} nonempty unions: {cnt}")
    classes = []
    cls_of = {}
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
    print(f"  [3] {nC} isomorphism classes of CC unions")
    below = [[j for j in range(nC) if j != i and popcount(classes[j]) < popcount(classes[i]) and embeds(classes[j], classes[i])]
             for i in range(nC)]
    dfaces = {0: [0]}
    for T, s in status.items():
        if s == "D":
            dfaces.setdefault(popcount(T), []).append(T)
    hD = {s: v for s, v in homology_by_size(dfaces, p).items() if v}
    expected = {s + 1: v for s, v in hD.items()}
    print(f"  [4] reduced F_{p}-homology of the D-part K_D by face size: {hD}  =>  the relative X-complex must have homology {expected} (necessary condition for acyclicity of K)")
    faces_by_class = {}
    for T, c in cls_of.items():
        faces_by_class.setdefault(c, []).append(T)
    order = sorted(range(nC), key=lambda i: popcount(classes[i]))
    good = []
    total = [0]

    def rec(pos, cur):
        if pos == nC:
            total[0] += 1
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
    print(f"  [5] {total[0]} down-sets of the CC-class poset; {len(good)} pass the relative-homology filter (necessary for F_{p}-acyclicity of the fixed complex; "
          f"the exact test is replay/check_z5_absolute.py)")
    exp = cert["expect"]       # REQUIRED: the pinned counts are part of the certificate; deleting them is an error, not a relaxation
    assert ("after_euler" in exp) == bool(cert.get("euler")), "the key expect.after_euler must be present exactly when the certificate has Euler groups"
    assert "downsets" in exp and "relative_survivors" in exp
    if "downsets" in exp:
        assert total[0] == exp["downsets"], f"expected {exp['downsets']} down-sets, found {total[0]}"
    if "relative_survivors" in exp:
        assert len(good) == exp["relative_survivors"], f"expected {exp['relative_survivors']} survivors of the relative filter, found {len(good)}"
    assert good, "no pattern passes the necessary filter: type 5 would be excluded"
    euler_counts = []
    # optional: Euler (Oliver, class G^1) equations chi(Delta_P^Gamma) = 1 of further groups Gamma whose invariant CC graphs all lie among the CC
    # union classes of Q; they are decided by the pattern.  Each group is verified to be in G^1 from the witness (normal p-subgroup N, Gamma = <g>N).
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
                ci = None
                for c, rep in enumerate(classes):
                    if isomorphic(rep, m):
                        ci = c
                        break
                assert ci is not None, "a CC union of the Euler group is not a union class of Q: equation not decided by the pattern"
                gtab[T] = ("CC", ci)
        survivors = []
        for pat in good:
            val = 0
            for T, (st, ci) in gtab.items():
                v = 1 if st == "D" else (0 if st == "N" else (1 if ci in pat else 0))
                val += (-1) ** (popcount(T) - 1) * v
            if val == 1:
                survivors.append(pat)
        print(f"  [5b] Euler equation (fix) of a group of order {ent['order']} with {gk} edge orbits (G^1 verified): {len(survivors)} of {len(good)} patterns survive")
        euler_counts.append(len(survivors))
        good = survivors
        assert good, "no pattern survives: type 5 would be excluded"
    if "after_euler" in exp:
        assert euler_counts == exp["after_euler"], f"expected survivor counts {exp['after_euler']} after the Euler equations, found {euler_counts}"
    ok = True
    for name, key, want in (("IN", "claims_in", True), ("OUT", "claims_out", False)):
        for edges in cert[key]:
            gm = mask_of([tuple(e) for e in edges])
            ci = None
            for c, rep in enumerate(classes):
                if isomorphic(rep, gm):
                    ci = c
                    break
            if ci is None:
                print(f"  claimed graph with {len(edges)} edges is not a CC union class of Q")
                ok = False
                continue
            holds = all((ci in pat) == want for pat in good)
            print(f"  [6] {name}: graph with {popcount(gm)} edges (class {ci}): "
                  f"{'in all' if want else 'in none of'} the {len(good)} surviving patterns? {holds}")
            ok = ok and holds
    print("  RESULT:", "CERTIFICATE VALID" if ok else "CERTIFICATE INVALID")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

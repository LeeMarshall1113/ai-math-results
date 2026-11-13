"""z5_table.py -- the 32 classes and the 8 admissible patterns of Lemma Z5, for the appendix of the note.

Re-derives, with the standard-library replay code, exactly what replay_nonselfdual.py derives
(Q = <(01234)(56789)>, the 32 isomorphism classes of connected graphs with connected complement that are unions of
Q-orbits, the 16 patterns with F_5-acyclic fixed complex, the 8 that also satisfy the Euler equations of the two groups
in certs/lemma_z5_euler.json), then writes note/z5_patterns.json and note/z5_table.tex.

usage:  python note/z5_table.py          (from the karp10 directory; about 15 seconds)
"""
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "replay"))
from replay_lib import *  # noqa: E402,F401,F403
from replay_smith import connected, homology_by_size  # noqa: E402


def charpoly(rows):
    """characteristic polynomial coefficients of the adjacency matrix (Faddeev-LeVerrier, exact integers)"""
    n = N
    A = [[(rows[i] >> j) & 1 for j in range(n)] for i in range(n)]
    I = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    Mk = [[0] * n for _ in range(n)]
    c = [1]
    for k in range(1, n + 1):
        # M_k = A M_{k-1} + c_{n-k+1} I
        AM = [[sum(A[i][t] * Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        Mk = [[AM[i][j] + c[-1] * I[i][j] for j in range(n)] for i in range(n)]
        AMk = [[sum(A[i][t] * Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        tr = sum(AMk[i][i] for i in range(n))
        assert tr % k == 0
        c.append(-tr // k)
    return tuple(c)


def invariants(mask):
    rows = rows_of(mask)
    deg = [popcount(r) for r in rows]
    tri_v = [0] * N
    tri = 0
    for a, b, c in itertools.combinations(range(N), 3):
        if (rows[a] >> b) & 1 and (rows[a] >> c) & 1 and (rows[b] >> c) & 1:
            tri += 1
            tri_v[a] += 1
            tri_v[b] += 1
            tri_v[c] += 1
    quad = 0                                   # number of 4-cycles
    for a, b in itertools.combinations(range(N), 2):
        common = popcount(rows[a] & rows[b] & ~((1 << a) | (1 << b)))
        quad += common * (common - 1) // 2
    quad //= 2
    return dict(edges=popcount(mask), deg=sorted(deg, reverse=True), tri=tri, quad=quad,
                local=sorted(zip(deg, tri_v), reverse=True), cp=charpoly(rows))


def degstr(deg):
    """degree multiset as LaTeX math (exponent = multiplicity)"""
    from collections import Counter
    c = Counter(deg)
    return "$" + r"\,".join(f"{d}^{{{c[d]}}}" for d in sorted(c, reverse=True)) + "$"


def main():
    cert = json.load(open(os.path.join(ROOT, "certs", "lemma_z5_euler.json")))
    p = cert["prime"]
    assert is_prime(p)
    gens = [tuple(g) for g in cert["gens"]]
    orbs = edge_orbits(gens)
    k = len(orbs)
    unions = {}
    for T in range(1 << k):
        m = 0
        for i in range(k):
            if (T >> i) & 1:
                m |= orbs[i]
        unions[T] = m
    st = {T: ("D" if not connected(m) else ("N" if not connected(complement(m)) else "CC")) for T, m in unions.items() if T}
    classes, cls_of = [], {}
    for T, s in st.items():
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
    below = [[j for j in range(nC) if j != i and popcount(classes[j]) < popcount(classes[i]) and embeds(classes[j], classes[i])]
             for i in range(nC)]
    dfaces = {0: [0]}
    for T, s in st.items():
        if s == "D":
            dfaces.setdefault(popcount(T), []).append(T)
    hD = {s: v for s, v in homology_by_size(dfaces, p).items() if v}
    expected = {s + 1: v for s, v in hD.items()}
    fbc = {}
    for T, c in cls_of.items():
        fbc.setdefault(c, []).append(T)
    order = sorted(range(nC), key=lambda i: popcount(classes[i]))
    good = []

    def rec(pos, cur):
        if pos == nC:
            faces = {}
            for c in cur:
                for T in fbc[c]:
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
    n16 = len(good)
    for ent in cert["euler"]:
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
            val = sum((-1) ** (popcount(T) - 1) * (1 if s == "D" else (0 if s == "N" else (1 if ci in pat else 0)))
                      for T, (s, ci) in gtab.items())
            if val == 1:
                keep.append(pat)
        good = keep
    assert n16 == 16 and len(good) == 8, (n16, len(good))
    cmap = {c: next(d for d in range(nC) if isomorphic(classes[d], complement(classes[c]))) for c in range(nC)}
    assert all(cmap[cmap[c]] == c and cmap[c] != c for c in range(nC))
    dual = lambda pat: frozenset(c for c in range(nC) if cmap[c] not in pat)
    assert all(dual(pt) in set(good) and dual(pt) != pt for pt in good)

    # the nine graphs of the certificate, identified by isomorphism with explicit constructions
    def circ(S):
        return mask_of([(i, (i + s) % 10) for i in range(10) for s in S if i != (i + s) % 10])
    def cyc(vs):
        return [(vs[t], vs[(t + 1) % len(vs)]) for t in range(len(vs))]
    K5a = list(itertools.combinations(range(5), 2))
    K5b = list(itertools.combinations(range(5, 10), 2))
    known = {
        "C10": circ([1]),
        "C5 + pendant edges": mask_of(cyc([0, 1, 2, 3, 4]) + [(i, i + 5) for i in range(5)]),
        "Petersen": mask_of(cyc([0, 1, 2, 3, 4]) + [(i, i + 5) for i in range(5)] + [(5, 7), (7, 9), (9, 6), (6, 8), (8, 5)]),
        "prism C5 x K2": mask_of(cyc([0, 1, 2, 3, 4]) + cyc([5, 6, 7, 8, 9]) + [(i, i + 5) for i in range(5)]),
        "K5 + pendant edges": mask_of(K5a + [(i, i + 5) for i in range(5)]),
        "5-sun": mask_of(cyc([0, 1, 2, 3, 4]) + [(i, 5 + i) for i in range(5)] + [((i + 1) % 5, 5 + i) for i in range(5)]),
        "C10^2": circ([1, 2]),
        "(K5 u C5)+PM": mask_of(K5a + cyc([5, 6, 7, 8, 9]) + [(i, i + 5) for i in range(5)]),
        "K5 x K2": mask_of(K5a + K5b + [(i, i + 5) for i in range(5)]),
    }
    nine_ids = {}
    for nm, m in known.items():
        hit = [c for c in range(nC) if isomorphic(classes[c], m)]
        assert len(hit) == 1, (nm, hit)
        nine_ids[nm] = hit[0]
    # the nine graphs of the certificate are exactly these
    cert_in = [mask_of([tuple(e) for e in edges]) for edges in cert["claims_in"]]
    for m in cert_in:
        assert sum(1 for c in nine_ids.values() if isomorphic(classes[c], m)) == 1
    assert len(set(nine_ids.values())) == 9 and len(cert_in) == 9

    inv = [invariants(classes[c]) for c in range(nC)]
    key = lambda c: (inv[c]["edges"], [-d for d in inv[c]["deg"]], inv[c]["tri"], inv[c]["quad"], inv[c]["local"], inv[c]["cp"], classes[c])
    sort_ids = sorted(range(nC), key=key)
    label = {c: i + 1 for i, c in enumerate(sort_ids)}
    sig = {}
    for c in range(nC):
        s = (inv[c]["edges"], tuple(inv[c]["deg"]), inv[c]["tri"], inv[c]["quad"], tuple(inv[c]["local"]), inv[c]["cp"])
        assert s not in sig, f"invariants do not separate classes {sig[s]} and {c}: {s[:4]}"
        sig[s] = c
    nunions = {c: len(fbc[c]) for c in range(nC)}
    pats = sorted(good, key=lambda pt: (len(pt), sorted(label[c] for c in pt)))
    plab = {pt: i + 1 for i, pt in enumerate(pats)}
    dualpat = {plab[pt]: plab[dual(pt)] for pt in pats}
    out = dict(
        classes=[dict(label=label[c], edges=inv[c]["edges"], degrees=inv[c]["deg"], triangles=inv[c]["tri"], four_cycles=inv[c]["quad"],
                      complement=label[cmap[c]], orbit_unions=nunions[c], rep_edges=edges_of(classes[c]),
                      name=next((nm for nm, i in nine_ids.items() if i == c), None))
                 for c in sort_ids],
        patterns={plab[pt]: sorted(label[c] for c in pt) for pt in pats},
        dual_pairs=sorted({tuple(sorted((a, b))) for a, b in dualpat.items()}),
    )
    json.dump(out, open(os.path.join(HERE, "z5_patterns.json"), "w"), indent=1)
    # LaTeX table, two halves side by side (classes 1-16 and 17-32)
    def row(c):
        marks = " & ".join(r"$\bullet$" if c in pt else "" for pt in pats)
        nm = next((n for n, i in nine_ids.items() if i == c), None)
        tag = r"$^{*}$" if nm else ""
        return f"{label[c]}{tag} & {inv[c]['edges']} & {degstr(inv[c]['deg'])} & {inv[c]['tri']} & {inv[c]['quad']} & {label[cmap[c]]} & {marks}"
    half = nC // 2
    lines = []
    lines.append(r"\begin{tabular}{@{}rrlrrr@{\hspace{6pt}}cccccccc@{\hspace{14pt}}rrlrrr@{\hspace{6pt}}cccccccc@{}}")
    lines.append(r"\toprule")
    hdr = r"class & $e$ & degrees & $t$ & $q$ & cpl. & 1&2&3&4&5&6&7&8"
    lines.append(hdr + " & " + hdr + r" \\")
    lines.append(r"\midrule")
    for i in range(half):
        lines.append(row(sort_ids[i]) + " & " + row(sort_ids[i + half]) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    open(os.path.join(HERE, "z5_table.tex"), "w").write("\n".join(lines) + "\n")
    print("classes:", nC, "patterns:", [len(pt) for pt in pats], "dual pairs:", out["dual_pairs"])
    print("the nine graphs are classes:", {nm: label[i] for nm, i in nine_ids.items()})
    pat_nine = [all(i in pt for i in nine_ids.values()) for pt in pats]
    print("all nine in every pattern:", all(pat_nine))
    print("pattern memberships:", {plab[pt]: sorted(label[c] for c in pt) for pt in pats})
    print("wrote note/z5_patterns.json and note/z5_table.tex")


if __name__ == "__main__":
    main()

"""Standalone exact checker for a K4^- flag-algebra certificate (pure Python, integers/Fractions).

    python check_cert.py cert.json[.gz]

Conventions.  A 3-graph on vertex set {0..6} is an int whose bit i is the triple
TRI[i], TRI = list(itertools.combinations(range(7), 3)).  A sigma-flag on f vertices has
type vertices 0..s-1 (inducing the type exactly) and free vertices s..f-1.

What is checked.
 (1) Catalogue: every H in cert["graphs"] is K4^- -free (no 4 vertices spanning >= 3 triples),
     the H are pairwise non-isomorphic, and sum_H 7!/|Aut H| equals the number of labelled
     K4^- -free 3-graphs on 7 vertices, which is computed here independently.  So the list is
     a complete set of isomorphism-class representatives.
 (2) Blocks: each flag restricts to the block's type on 0..s-1, flags are distinct and in
     canonical form; Gram matrices are G = sum_k v_k v_k^T with integer v_k (PSD by construction).
 (3) For every H the exact value
        val(H) = e(H)/35 + sum_b c_b(H) / (DEN_b * M^2)
                 + sum_F mu_F * lin_F(H)/420
                 + sum_g [ p_g (cD_g(H) - t1 cG_g(H)) - q_g (cD_g(H) - t2 cG_g(H)) ] / 105
     is computed, where c_b(H) = sum over injective theta:[s]->V(H) inducing the type and over
     ordered splits (A,B) of the other 7-s vertices into two (7-s)/2-sets of
     G[flag(theta+A)][flag(theta+B)]  (0 if a flag is not listed), DEN_b = (7)_s * C(7-s,(7-s)/2);
     lin_F(H) = sum over (v, 3-set S) with (H[v+S],v) ~ F of (#{pairs P in R: vP in H} - 3[R in H]),
        R = V - v - S  (this is 420*([[d F]]_1 - e*[[F]]_1));
     cG_g(H) = #{(v, 4-set S): (H[v+S],v) ~ g},  cD_g(H) = same with the extra condition vR in H
        (R = the remaining pair)  (so cD/105 = [[d g]]_1 and cG/105 = [[g]]_1);
     p_g, q_g >= 0;  plus product rows  sum_g [pp_g (cE_g - t1 cP_g) - pq_g (cE_g - t2 cP_g)]/35 over
     4-vertex graphs g, cP_g(H) = #{4-sets U: H[U] ~ g}, cE_g = same with V-U an edge
     (cE/35 = [e * g], cP/35 = [g]; phi(e*g) = phi(e) phi(g) for every limit phi), pp, pq >= 0.
     The checker reports max_H val(H) and confirms it is <= cert["bound"].
Mathematical consequence (see REPORT.md): if phi is an extremal limit that is degree-regular
and phi(e) = pi lies in [t1, t2], then pi <= max_H val(H).  Without interval rows: pi <= max val.
"""
import gzip, itertools, json, sys, time
from fractions import Fraction
from math import comb, factorial

TRI = list(itertools.combinations(range(7), 3))
BIT = {}
for i, t in enumerate(TRI):
    for p in itertools.permutations(t):
        BIT[p] = 1 << i
QUAD = {}

def k4m_free(m, n):
    for q in itertools.combinations(range(n), 4):
        c = 0
        for t in itertools.combinations(q, 3):
            if m & BIT[t]:
                c += 1
        if c >= 3:
            return False
    return True

def relabel(m, verts):
    """induced subgraph on the ordered vertex list, relabelled 0..k-1"""
    out = 0
    k = len(verts)
    for i in range(k):
        for j in range(i + 1, k):
            for l in range(j + 1, k):
                if m & BIT[(verts[i], verts[j], verts[l])]:
                    out |= BIT[(i, j, l)]
    return out

def flag_key(m, typ, ext):
    return min(relabel(m, list(typ) + list(p)) for p in itertools.permutations(ext))

# ---------------- (1) catalogue ----------------
def canon_aut(m):
    n = 7
    deg = [0] * n; cod = [[0] * n for _ in range(n)]
    for t in TRI:
        if m & BIT[t]:
            for a in t:
                deg[a] += 1
            for a, b in itertools.combinations(t, 2):
                cod[a][b] += 1; cod[b][a] += 1
    inv = [(deg[v], tuple(sorted(cod[v][u] for u in range(n) if u != v))) for v in range(n)]
    keys = sorted(set(inv))
    cells = [[v for v in range(n) if inv[v] == k] for k in keys]
    best = None; cnt = 0
    for parts in itertools.product(*[list(itertools.permutations(c)) for c in cells]):
        order = [v for part in parts for v in part]
        r = relabel(m, order)
        if best is None or r < best:
            best, cnt = r, 1
        elif r == best:
            cnt += 1
    return best, cnt   # cnt = |Aut(m)| (orderings within the invariant cells giving the minimum)

def labelled_count7():
    """number of labelled K4^- -free 3-graphs on {0..6}, by extending labelled 6-vertex ones."""
    tri6 = [t for t in TRI if max(t) < 6]                       # 20 triples
    b6 = {t: 1 << i for i, t in enumerate(tri6)}
    quads6 = []
    for q in itertools.combinations(range(6), 4):
        quads6.append(sum(b6[t] for t in itertools.combinations(q, 3)))
    ind = bytearray(1 << 20)
    L6 = 0
    for m in range(1 << 20):
        ok = True
        for q in quads6:
            if (m & q).bit_count() >= 3:
                ok = False; break
        if ok:
            ind[m] = 1; L6 += 1
    g = list(ind)                                                # g[U] = #{G in L6 : G subset of U}
    for i in range(20):
        bit = 1 << i
        for U in range(1 << 20):
            if U & bit:
                g[U] += g[U ^ bit]
    pairs = list(itertools.combinations(range(6), 2))
    full = (1 << 20) - 1
    L7 = 0
    for link in range(1 << 15):
        ls = {pairs[i] for i in range(15) if link >> i & 1}
        bad = 0; ok = True
        for t in tri6:
            c = ((t[0], t[1]) in ls) + ((t[0], t[2]) in ls) + ((t[1], t[2]) in ls)
            if c == 3:
                ok = False; break
            if c == 2:
                bad |= b6[t]
        if ok:
            L7 += g[full ^ bad]
    return L6, L7

def check_catalogue(graphs):
    seen = set(); orbit = 0
    for m in graphs:
        assert k4m_free(m, 7), "catalogue graph contains K4^-"
        c, aut = canon_aut(m)
        assert c not in seen, "two catalogue graphs are isomorphic"
        seen.add(c); orbit += factorial(7) // aut
    L6, L7 = labelled_count7()
    print(f"  labelled K4^- -free: n=6 {L6}, n=7 {L7}; catalogue orbit sum {orbit}; classes {len(graphs)}")
    assert L6 == 33652
    assert orbit == L7, "catalogue incomplete"

# ---------------- (2)+(3) ----------------
def check(cert):
    t0 = time.time()
    graphs = cert["graphs"]
    print("(1) catalogue"); check_catalogue(graphs); print(f"  ok  {time.time()-t0:.0f}s")
    M = int(cert["scale"]); M2 = M * M
    blocks = []
    print("(2) blocks")
    for blk in cert["blocks"]:
        s = int(blk["s"]); typ = int(blk["type"]); f = (7 + s) // 2
        assert s in (1, 3, 5)
        assert k4m_free(typ, s) and relabel(typ, list(range(s))) == typ
        flags = [int(x) for x in blk["flags"]]
        assert len(set(flags)) == len(flags)
        for F in flags:
            assert relabel(F, list(range(s))) == typ, "flag does not extend the type"
            assert k4m_free(F, f)
            assert flag_key(F, list(range(s)), list(range(s, f))) == F, "flag not canonical"
        n = len(flags)
        G = [[0] * n for _ in range(n)]
        for v in blk["factors"]:
            v = [int(x) for x in v]; assert len(v) == n
            nz = [(i, x) for i, x in enumerate(v) if x]
            for i, a in nz:
                Gi = G[i]
                for j, b in nz:
                    Gi[j] += a * b
        den = factorial(7) // factorial(7 - s) * comb(7 - s, (7 - s) // 2)
        blocks.append((s, typ, {F: i for i, F in enumerate(flags)}, G, den))
    print(f"  {len(blocks)} blocks ok  {time.time()-t0:.0f}s")
    lin = cert.get("linear_stat") or {"flags": [], "mu": []}
    linF = {int(F): Fraction(mu) for F, mu in zip(lin["flags"], lin["mu"])}
    itv = cert.get("interval_stat")
    if itv:
        t1, t2 = Fraction(itv["t1"]), Fraction(itv["t2"])
        pg = {int(g): Fraction(p) for g, p in zip(itv.get("flags", []), itv.get("p", []))}
        qg = {int(g): Fraction(q) for g, q in zip(itv.get("flags", []), itv.get("q", []))}
        pp = {int(g): Fraction(p) for g, p in zip(itv.get("graphs4", []), itv.get("pp", []))}
        pq = {int(g): Fraction(q) for g, q in zip(itv.get("graphs4", []), itv.get("pq", []))}
        for d in (pg, qg, pp, pq):
            assert all(x >= 0 for x in d.values())
        for g in pg:
            assert relabel(g, [0]) == 0 and k4m_free(g, 5) and flag_key(g, [0], [1, 2, 3, 4]) == g
        for g in pp:
            assert k4m_free(g, 4) and g == min(relabel(g, list(p)) for p in itertools.permutations(range(4)))
    for F in linF:
        assert relabel(F, [0]) == 0 and k4m_free(F, 4)
    print("(3) coefficients")
    V = list(range(7))
    best = None; arg = None
    for h, H in enumerate(graphs):
        val = Fraction(bin(H).count("1"), 35)
        for s, typ, fidx, G, den in blocks:
            f = (7 + s) // 2; k = f - s
            acc = 0
            for theta in itertools.permutations(V, s):
                if relabel(H, list(theta)) != typ:
                    continue
                rest = [v for v in V if v not in theta]
                for A in itertools.combinations(rest, k):
                    B = [v for v in rest if v not in A]
                    ia = fidx.get(flag_key(H, theta, A)); ib = fidx.get(flag_key(H, theta, B))
                    if ia is not None and ib is not None:
                        acc += G[ia][ib]
            val += Fraction(acc, den * M2)
        if linF:
            for v in V:
                R0 = [x for x in V if x != v]
                for S in itertools.combinations(R0, 3):
                    F = flag_key(H, [v], S)
                    if F in linF:
                        R = [x for x in R0 if x not in S]
                        n1 = sum(1 for a, b in itertools.combinations(R, 2) if H & BIT[(v, a, b)])
                        n2 = 1 if H & BIT[tuple(R)] else 0
                        val += linF[F] * (n1 - 3 * n2) / 420
        if itv:
            for v in V:
                R0 = [x for x in V if x != v]
                for S in itertools.combinations(R0, 4):
                    g = flag_key(H, [v], S)
                    if g in pg:
                        a, b = [x for x in R0 if x not in S]
                        dd = 1 if H & BIT[(v, a, b)] else 0
                        val += (pg[g] * (dd - t1) - qg[g] * (dd - t2)) / 105
            for U in itertools.combinations(V, 4):
                r = relabel(H, list(U))
                g = min(relabel(r, list(p)) for p in itertools.permutations(range(4)))
                if g in pp:
                    T = tuple(x for x in V if x not in U)
                    ee = 1 if H & BIT[T] else 0
                    val += (pp[g] * (ee - t1) - pq[g] * (ee - t2)) / 35
        if best is None or val > best:
            best, arg = val, H
        if h % 1000 == 0:
            print(f"  {h}/{len(graphs)}  max so far {float(best):.10f}  {time.time()-t0:.0f}s", flush=True)
    claim = Fraction(cert["bound"])
    print(f"max_H val(H) = {best}  = {float(best):.12f}   (attained at H = {arg})")
    print(f"claimed bound {claim} = {float(claim):.12f}")
    assert best <= claim, "CLAIM FAILS"
    if itv:
        print(f"interval [{t1}, {t2}] = [{float(t1):.9f}, {float(t2):.9f}]:",
              "EXCLUDED (bound < t1)" if best < t1 else "not excluded")
    print(f"PASS  {time.time()-t0:.0f}s")
    return best

if __name__ == "__main__":
    p = sys.argv[1]
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        cert = json.load(fh)
    check(cert)

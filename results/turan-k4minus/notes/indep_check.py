"""Independent referee re-implementation of the K4^- certificate evaluation.

Different algorithm from check_cert.py: every statistic is computed by summing over all
7! orderings of V(H) (vectorised), dividing by the stabiliser of the split; flags are
identified via a lookup table of ALL labelled copies of each listed flag (no canonical form).

usage: python indep_check.py CERT [nH]   (nH = number of graphs to evaluate; default all)
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import gzip, json, itertools, random, time
from fractions import Fraction
from math import factorial, comb
import numpy as np

TRI = list(itertools.combinations(range(7), 3))
TI = {t: i for i, t in enumerate(TRI)}
PERMS = np.array(list(itertools.permutations(range(7))), dtype=np.int64)  # 5040x7


def tri_map(order):
    """order: array (P,7) of vertex orderings. Returns (P,35): for new triple t=(a,b,c) the
    index of the old triple {order[a],order[b],order[c]}."""
    P = order.shape[0]
    out = np.zeros((P, 35), dtype=np.int64)
    for t, (a, b, c) in enumerate(TRI):
        x = np.sort(order[:, [a, b, c]], axis=1)
        out[:, t] = [TI[tuple(r)] for r in x.tolist()]
    return out


def code(bits, f):
    """bits (P,35) bool of relabelled graph; code of induced graph on new vertices 0..f-1."""
    idx = [t for t, tr in enumerate(TRI) if tr[2] < f]
    w = np.array([1 << t for t in idx], dtype=np.int64)
    return bits[:, idx].astype(np.int64) @ w


def mask_relabel(m, verts):
    out = 0
    for t, (a, b, c) in enumerate(TRI):
        if c < len(verts):
            if m >> TI[tuple(sorted((verts[a], verts[b], verts[c])))] & 1:
                out |= 1 << t
    return out


def k4m_free(m, n):
    for q in itertools.combinations(range(n), 4):
        if sum(m >> TI[t] & 1 for t in itertools.combinations(q, 3)) >= 3:
            return False
    return True


def lookup_table(flags, s, f):
    """all labelled copies (permuting free vertices) -> flag index"""
    d = {}
    for i, F in enumerate(flags):
        for p in itertools.permutations(range(s, f)):
            c = mask_relabel(F, list(range(s)) + list(p))
            if c in d and d[c] != i:
                raise SystemExit("two listed flags are isomorphic")
            d[c] = i
    keys = np.array(sorted(d), dtype=np.int64)
    vals = np.array([d[k] for k in keys.tolist()], dtype=np.int64)
    return keys, vals


def find(keys, vals, c):
    pos = np.searchsorted(keys, c)
    pos = np.minimum(pos, len(keys) - 1)
    ok = keys[pos] == c
    return np.where(ok, vals[pos], -1)


def gram(factors, n):
    F = np.array([[int(x) for x in v] for v in factors], dtype=object).reshape(-1, n)
    Fi = np.array(F, dtype=np.int64) if len(F) else np.zeros((0, n), np.int64)
    assert np.all(np.abs(Fi) < 2 ** 30)
    hi = Fi >> 15; lo = Fi - (hi << 15)
    a = hi.T @ hi; b = hi.T @ lo; c = lo.T @ lo
    G = [[(int(a[i, j]) << 30) + ((int(b[i, j]) + int(b[j, i])) << 15) + int(c[i, j])
          for j in range(n)] for i in range(n)]
    return G


def main():
    cert = json.load(gzip.open(sys.argv[1], "rt"))
    graphs = [int(g) for g in cert["graphs"]]
    nH = int(sys.argv[2]) if len(sys.argv) > 2 else len(graphs)
    sel = list(range(len(graphs)))
    if nH < len(graphs):
        random.seed(12345); sel = sorted(random.sample(sel, nH))
    M2 = int(cert["scale"]) ** 2
    T0 = tri_map(PERMS)
    blocks = []
    for blk in cert["blocks"]:
        s = int(blk["s"]); typ = int(blk["type"]); f = (7 + s) // 2; k = f - s
        flags = [int(x) for x in blk["flags"]]
        for F in flags:
            assert mask_relabel(F, list(range(s))) == typ and k4m_free(F, f)
        keys, vals = lookup_table(flags, s, f)
        n = len(flags)
        G = gram(blk["factors"], n)
        full2 = np.concatenate([PERMS[:, :s], PERMS[:, f:], PERMS[:, s:f]], axis=1)
        blocks.append(dict(s=s, f=f, k=k, typc=typ, keys=keys, vals=vals, n=n, G=G,
                           T2=tri_map(full2)))
    itv = cert.get("interval_stat")
    if itv:
        t1, t2 = Fraction(itv["t1"]), Fraction(itv["t2"])
        gl = [int(g) for g in itv["flags"]]
        P = [Fraction(x) for x in itv["p"]]; Q = [Fraction(x) for x in itv["q"]]
        g4 = [int(g) for g in itv["graphs4"]]
        PP = [Fraction(x) for x in itv["pp"]]; PQ = [Fraction(x) for x in itv["pq"]]
        assert min(P + Q + PP + PQ) >= 0
        gkeys, gvals = lookup_table(gl, 1, 5)
        d4 = {}
        for i, g in enumerate(g4):
            for p in itertools.permutations(range(4)):
                d4[mask_relabel(g, list(p))] = i
        k4keys = np.array(sorted(d4), dtype=np.int64)
        k4vals = np.array([d4[x] for x in k4keys.tolist()], dtype=np.int64)
        pair_t = TI[(0, 5, 6)]      # root 0 with remaining pair (5,6)
        rest_t = TI[(4, 5, 6)]      # complement triple of U = {0,1,2,3}
    assert not cert.get("linear_stat"), "linear rows not implemented here"
    t0 = time.time(); best = None; vals_out = {}
    for cnt, h in enumerate(sel):
        H = graphs[h]
        hb = np.array([(H >> i) & 1 for i in range(35)], dtype=bool)
        bits = hb[T0]
        val = Fraction(bin(H).count("1"), 35)
        for b in blocks:
            s, f, k = b["s"], b["f"], b["k"]
            ct = code(bits, s)
            ok = ct == b["typc"]
            ia = find(b["keys"], b["vals"], code(bits, f))
            bits2 = hb[b["T2"]]
            ib = find(b["keys"], b["vals"], code(bits2, f))
            good = ok & (ia >= 0) & (ib >= 0)
            n = b["n"]
            cnts = np.bincount(ia[good] * n + ib[good], minlength=n * n)
            nz = np.nonzero(cnts)[0]
            acc = sum(int(cnts[z]) * b["G"][z // n][z % n] for z in nz.tolist())
            stab = factorial(k) ** 2
            den = factorial(7) // factorial(7 - s) * comb(7 - s, k)
            val += Fraction(acc, stab * den * M2)
        if itv:
            gi = find(gkeys, gvals, code(bits, 5))
            dd = bits[:, pair_t]
            for i in range(len(gl)):
                sel_i = gi == i
                cG = Fraction(int(sel_i.sum()), 48); cD = Fraction(int((sel_i & dd).sum()), 48)
                val += (P[i] * (cD - t1 * cG) - Q[i] * (cD - t2 * cG)) / 105
            ui = find(k4keys, k4vals, code(bits, 4))
            ee = bits[:, rest_t]
            for i in range(len(g4)):
                sel_i = ui == i
                cP = Fraction(int(sel_i.sum()), 144); cE = Fraction(int((sel_i & ee).sum()), 144)
                val += (PP[i] * (cE - t1 * cP) - PQ[i] * (cE - t2 * cP)) / 35
        vals_out[h] = val
        if best is None or val > best[0]:
            best = (val, H)
        if cnt % 500 == 0:
            print(f"{cnt}/{len(sel)} max {float(best[0]):.12f} {time.time()-t0:.0f}s", flush=True)
    print(f"evaluated {len(sel)} graphs; max val = {best[0]} = {float(best[0]):.13f} at H={best[1]}")
    print(f"claimed bound {cert['bound']} = {float(Fraction(cert['bound'])):.13f}")
    print("max <= claim:", best[0] <= Fraction(cert["bound"]))


if __name__ == "__main__":
    main()

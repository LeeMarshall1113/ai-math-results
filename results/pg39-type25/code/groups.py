"""PGammaL(4,9) acting on the 820 points of PG(3,9): generators and the
representatives of the (conjectured complete) list of classes of subgroups of
prime order.  Everything is built from our own GF(9) / PG(3,9) code.

Matrices act on column vectors, v -> M v, then normalise.  A semilinear
element is a pair (M, f) = P -> normalize(M * frob^f(P)).
"""
import itertools
from ff import GF
from pg import PG


def matmul(F, A, B):
    n = len(A)
    return [[_dot(F, [A[i][k] for k in range(n)], [B[k][j] for k in range(n)]) for j in range(n)]
            for i in range(n)]


def _dot(F, u, v):
    s = 0
    for x, y in zip(u, v):
        s = F.add[s][F.mul[x][y]]
    return s


def ident(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def matpow(F, A, e):
    R, B = ident(len(A)), A
    while e:
        if e & 1:
            R = matmul(F, R, B)
        B = matmul(F, B, B)
        e >>= 1
    return R


def blockdiag(*blocks):
    n = sum(len(b) for b in blocks)
    M = [[0] * n for _ in range(n)]
    o = 0
    for b in blocks:
        for i in range(len(b)):
            for j in range(len(b)):
                M[o + i][o + j] = b[i][j]
        o += len(b)
    return M


def companion(F, coeffs):
    """companion matrix of x^n + c_{n-1} x^{n-1} + ... + c_0, coeffs = [c_0..c_{n-1}]"""
    n = len(coeffs)
    M = [[0] * n for _ in range(n)]
    for i in range(1, n):
        M[i][i - 1] = 1
    for i in range(n):
        M[i][n - 1] = F.neg[coeffs[i]]
    return M


def mat_order(F, A, bound=10 ** 6):
    n = len(A)
    I = ident(n)
    R = A
    k = 1
    while R != I:
        R = matmul(F, R, A)
        k += 1
        if k > bound:
            return None
    return k


def is_scalar(A):
    n = len(A)
    return all(A[i][j] == 0 for i in range(n) for j in range(n) if i != j) and len({A[i][i] for i in range(n)}) == 1


def primitive_companion(F, n):
    """companion matrix of a monic degree-n polynomial over GF(q) whose root generates GF(q^n)*"""
    q = F.q
    target = q ** n - 1
    factors = [p for p in range(2, target + 1) if target % p == 0 and all(p % d for d in range(2, int(p ** .5) + 1))]
    for coeffs in itertools.product(range(q), repeat=n):
        if coeffs[0] == 0:
            continue
        C = companion(F, list(coeffs))
        if matpow(F, C, target) != ident(n):
            continue
        if all(matpow(F, C, target // p) != ident(n) for p in factors):
            return C
    raise RuntimeError("no primitive polynomial")


def order5_blocks(F):
    """the two classes of 2x2 matrices of order 5 over GF(9) (traces are the roots of t^2+t-1)"""
    q = F.q
    res = {}
    for a, b, c, d in itertools.product(range(q), repeat=4):
        M = [[a, b], [c, d]]
        if M == ident(2):
            continue
        if matpow(F, M, 5) != ident(2):
            continue
        tr = F.add[a][d]
        res.setdefault(tr, M)
    assert len(res) == 2, res
    return res


def build_all(q=9):
    F = GF(q)
    P = PG(4, q, F)
    omega = F.prim
    n = 4
    Id = ident(4)

    def E(i, j, a):
        M = ident(4)
        M[i][j] = a
        return M

    def diag(*d):
        M = [[0] * 4 for _ in range(4)]
        for i, x in enumerate(d):
            M[i][i] = x
        return M

    def permmat(images):    # e_i -> e_{images[i]}
        M = [[0] * 4 for _ in range(4)]
        for i, j in enumerate(images):
            M[j][i] = 1
        return M

    gens = {
        "cyc": (permmat([1, 2, 3, 0]), 0),
        "swap": (permmat([1, 0, 2, 3]), 0),
        "E12_1": (E(0, 1, 1), 0),
        "E12_w": (E(0, 1, omega), 0),
        "dw": (diag(omega, 1, 1, 1), 0),
        "frob": (Id, 1),
    }
    minus1 = F.neg[1]
    lam = omega      # a non-square of GF(9)*  (a generator is a non-square)
    B = [[0, lam], [1, 0]]     # B^2 = lam * I
    reps = {}
    # involutions
    reps["2a_homology_3+1"] = (diag(1, 1, 1, minus1), 0)
    reps["2b_homology_2+2"] = (diag(1, 1, minus1, minus1), 0)
    reps["2c_fixedpointfree"] = (blockdiag(B, B), 0)
    reps["2d_baer"] = (Id, 1)
    # order 3 (unipotent)
    reps["3a_J211"] = (E(0, 1, 1), 0)
    M22 = ident(4); M22[0][1] = 1; M22[2][3] = 1
    reps["3b_J22"] = (M22, 0)
    M31 = ident(4); M31[0][1] = 1; M31[1][2] = 1
    reps["3c_J31"] = (M31, 0)
    # order 5
    blocks = order5_blocks(F)
    tr = sorted(blocks)
    M1, M2 = blocks[tr[0]], blocks[tr[1]]
    reps["5a_scalar_GF81"] = (blockdiag(M1, M1), 0)
    reps["5b_fixline"] = (blockdiag(ident(2), M1), 0)
    reps["5c_mixed"] = (blockdiag(M1, M2), 0)
    # order 7, 13 : GL(3,9) Singer
    C3 = primitive_companion(F, 3)
    assert mat_order(F, C3) == 728
    reps["7_singer3"] = (blockdiag(matpow(F, C3, 104), [[1]]), 0)
    reps["13_singer3"] = (blockdiag(matpow(F, C3, 56), [[1]]), 0)
    # order 41 : GL(4,9) Singer
    C4 = primitive_companion(F, 4)
    assert mat_order(F, C4) == 6560
    reps["41_singer4"] = (matpow(F, C4, 160), 0)
    return F, P, gens, reps


def perm_of(P, elt):
    M, f = elt
    return P.perm_of_matrix(M, f)


def perm_order(perm):
    seen = [False] * len(perm)
    from math import gcd
    o = 1
    for i in range(len(perm)):
        if not seen[i]:
            l = 0
            j = i
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                l += 1
            o = o * l // gcd(o, l)
    return o


def cycle_type(perm):
    from collections import Counter
    seen = [False] * len(perm)
    cnt = Counter()
    for i in range(len(perm)):
        if not seen[i]:
            l = 0
            j = i
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                l += 1
            cnt[l] += 1
    return dict(sorted(cnt.items()))


def write_gap(path, P, gens, reps):
    with open(path, "w") as fh:
        fh.write("# generated by groups.py : permutations of the 820 points of PG(3,9) (1-based)\n")
        for name, elt in gens.items():
            perm = perm_of(P, elt)
            fh.write("gen_%s := PermList([%s]);;\n" % (name.replace("-", "_"), ",".join(str(x + 1) for x in perm)))
        fh.write("gens := [%s];;\n" % ",".join("gen_" + n for n in gens))
        fh.write("reps := [];;\nrepnames := [];;\n")
        for name, elt in reps.items():
            perm = perm_of(P, elt)
            fh.write("Add(repnames, \"%s\");; Add(reps, PermList([%s]));;\n" % (name, ",".join(str(x + 1) for x in perm)))


if __name__ == "__main__":
    F, P, gens, reps = build_all()
    for name, elt in reps.items():
        perm = perm_of(P, elt)
        ct = cycle_type(perm)
        print("%-20s order %-3d cycle type %s" % (name, perm_order(perm), ct))
    import sys, os
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gap", "pgammal49_data.g")
    write_gap(out, P, gens, reps)
    print("wrote", out)


def export_classes_json(path):
    """matrices of the 13 representatives in a field-independent format: entries are [a, b]
    meaning a + b*i in GF(9) = GF(3)[i]/(i^2+1) (our element code is a + 3b)."""
    F, P, gens, reps = build_all()
    assert F.q == 9 and F.modulus == [1, 0, 1]

    def enc(c):
        return [c % 3, c // 3]
    out = {}
    for name, (M, f) in reps.items():
        out[name] = dict(frob=f, matrix=[[enc(c) for c in row] for row in M])
    import json
    json.dump(out, open(path, "w"), indent=1)
    return out

#!/usr/bin/env python3
"""Independent computation (stdlib, no nauty) of |Aut| of the planar class representatives and of their pairwise inequivalence under PGammaL(3,9).

  python planar_aut.py [--procs 2]

PGammaL(3,9) = {M o frob^e} (M in PGL(3,9), e in {0,1}; the decomposition is unique).  For the FIXED standard frame F0 (which frob fixes pointwise) the map g = M o frob^e -> (g(F0), e) = (M(F0), e)
is a BIJECTION from PGammaL(3,9) onto the set of pairs (ordered projective frame F' = (b1,b2,b3,b4), no three collinear; e in {0,1}): for every such pair there is exactly one g = M o frob^e with
g(F0) = F' (M is determined by F' alone since PGL(3,9) is sharply transitive on ordered frames).  [Wording corrected 2026-11-11 after the referee's remark: this is a bijection of GROUP ELEMENTS with
pairs relative to F0, NOT an action of PGammaL(3,9) on frames alone; frob fixes F0, so the group does not act regularly on frames.  The counting below only uses the bijection.]  Hence for sets S1, S2 of PG(2,9):
    #{g : g(S1) = S2} = #{(F', e) : F' an ordered frame inside S2, g_{F',e}(S1) = S2}     (F0 chosen inside S1; any g with g(S1) = S2 maps F0 into S2).
|Aut(S)| = that number for S1 = S2 = S.  The script computes |Aut| of the 3 + 7 representatives of results/planar_classes.json (compare with nauty's aut_order), and, for the pairs of
representatives with equal |Aut|, searches for an isomorphism (none must exist).  Together with planar_count.py (the number of planar (2,5)-sets of size m counted without nauty equals
sum |PGammaL(3,9)| / |Aut_i|) this shows that the listed representatives are pairwise inequivalent and represent ALL classes.
"""
import json
import multiprocessing
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import replay as R   # noqa: E402
import star_setup as SS   # noqa: E402

RES = os.path.join(HERE, "..", "results")
ADD = [[R.e_add(a, b) for b in range(9)] for a in range(9)]
MUL = [[R.e_mul(a, b) for b in range(9)] for a in range(9)]
NEG = [MUL[a][2] for a in range(9)]
INV = [0] + [R.E_INV[a] for a in range(1, 9)]
FROB = [R.e_frob(a) for a in range(9)]
PTS = SS.build_points3()


def cof(x1, x2, x3, x4):
    return ADD[MUL[x1][x4]][NEG[MUL[x2][x3]]]


def inv3(m):
    (a, b, c), (d, e, f), (g, h, i) = m
    C11, C12, C13 = cof(e, f, h, i), NEG[cof(d, f, g, i)], cof(d, e, g, h)
    C21, C22, C23 = NEG[cof(b, c, h, i)], cof(a, c, g, i), NEG[cof(a, b, g, h)]
    C31, C32, C33 = cof(b, c, e, f), NEG[cof(a, c, d, f)], cof(a, b, d, e)
    det = ADD[ADD[MUL[a][C11]][MUL[b][C12]]][MUL[c][C13]]
    if det == 0:
        return None
    t = INV[det]
    return [[MUL[t][C11], MUL[t][C21], MUL[t][C31]], [MUL[t][C12], MUL[t][C22], MUL[t][C32]], [MUL[t][C13], MUL[t][C23], MUL[t][C33]]]


def matvec(m, v):
    return tuple(ADD[ADD[MUL[r[0]][v[0]]][MUL[r[1]][v[1]]]][MUL[r[2]][v[2]]] for r in m)


def matmul(A, B):
    cols = [(B[0][c], B[1][c], B[2][c]) for c in range(3)]
    return [[ADD[ADD[MUL[A[r][0]][col[0]]][MUL[A[r][1]][col[1]]]][MUL[A[r][2]][col[2]]] for col in cols] for r in range(3)]


def normalize(v):
    for c in v:
        if c:
            t = INV[c]
            return tuple(MUL[t][x] for x in v)
    return None


def columns(v1, v2, v3):
    return [[v1[0], v2[0], v3[0]], [v1[1], v2[1], v3[1]], [v1[2], v2[2], v3[2]]]


def frame_data(u1, u2, u3, u4):
    """Ainv with A = [c1 u1 | c2 u2 | c3 u3], u4 = sum c_i u_i ; None if (u1..u4) is not a frame"""
    Ui = inv3(columns(u1, u2, u3))
    if Ui is None:
        return None
    c = matvec(Ui, u4)
    if 0 in c:
        return None
    return inv3(columns(tuple(MUL[c[0]][x] for x in u1), tuple(MUL[c[1]][x] for x in u2), tuple(MUL[c[2]][x] for x in u3)))


def count_maps(S1, S2, first_only=False):
    """number of g in PGammaL(3,9) with g(S1) = S2 (S1, S2 lists of point indices)"""
    S1v = [PTS[p] for p in S1]
    S2v = [PTS[p] for p in S2]
    S2set = set(S2v)
    base = None
    n1 = len(S1v)
    for i in range(n1):
        for j in range(n1):
            if j == i:
                continue
            for k in range(n1):
                if k in (i, j):
                    continue
                for l in range(n1):
                    if l in (i, j, k):
                        continue
                    if frame_data(S1v[i], S1v[j], S1v[k], S1v[l]) is not None:
                        base = (S1v[i], S1v[j], S1v[k], S1v[l])
                        break
                if base:
                    break
            if base:
                break
        if base:
            break
    Ainv, S1e = [], []
    for e in (0, 1):
        fr = [tuple(FROB[x] for x in u) for u in base] if e else list(base)
        Ainv.append(frame_data(*fr))
        S1e.append([tuple(FROB[x] for x in s) for s in S1v] if e else S1v)
    count = 0
    n2 = len(S2v)
    for i in range(n2):
        vi = S2v[i]
        for j in range(n2):
            if j == i:
                continue
            vj = S2v[j]
            for k in range(n2):
                if k in (i, j):
                    continue
                vk = S2v[k]
                Vi = inv3(columns(vi, vj, vk))
                if Vi is None:
                    continue
                for l in range(n2):
                    if l in (i, j, k):
                        continue
                    d = matvec(Vi, S2v[l])
                    if 0 in d:
                        continue
                    B = columns(tuple(MUL[d[0]][x] for x in vi), tuple(MUL[d[1]][x] for x in vj), tuple(MUL[d[2]][x] for x in vk))
                    for e in (0, 1):
                        M = matmul(B, Ainv[e])
                        good = True
                        for w in S1e[e]:
                            if normalize(matvec(M, w)) not in S2set:
                                good = False
                                break
                        if good:
                            count += 1
                            if first_only:
                                return count
    return count


def aut_job(args):
    name, S = args
    t = time.time()
    n = count_maps(S, S)
    return name, n, time.time() - t


def main():
    nproc = int(sys.argv[sys.argv.index("--procs") + 1]) if "--procs" in sys.argv else 2
    classes = json.load(open(os.path.join(RES, "planar_classes.json")))
    reps = []
    for k in (26, 35):
        for i, c in enumerate([c for c in classes if c["k"] == k]):
            reps.append(("%s%d" % ("A" if k == 26 else "C", i), c["rep"], c["aut_order"]))
    with multiprocessing.Pool(nproc) as pool:
        res = pool.map(aut_job, [(n, S) for n, S, a in reps])
    ok = True
    for (name, S, nauty), (_, cnt, secs) in zip(reps, res):
        print("%s: |Aut| by the frame count = %d, nauty = %d (%.0fs) %s" % (name, cnt, nauty, secs, "EQUAL" if cnt == nauty else "DIFFERENT"), flush=True)
        ok &= (cnt == nauty)
    for pre in ("A", "C"):
        grp = [(n, S, a) for n, S, a in reps if n[0] == pre]
        for x in range(len(grp)):
            for y in range(x + 1, len(grp)):
                if grp[x][2] == grp[y][2]:
                    n = count_maps(grp[x][1], grp[y][1], first_only=True)
                    print("%s ~ %s (equal |Aut| = %d): isomorphism %s" % (grp[x][0], grp[y][0], grp[x][2], "FOUND (classes are equal!)" if n else "none exists (inequivalent)"), flush=True)
                    ok &= (n == 0)
    print("ALL |Aut| CONFIRMED AND ALL REPRESENTATIVES PAIRWISE INEQUIVALENT" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

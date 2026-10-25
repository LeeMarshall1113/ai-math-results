"""Level-1 case split (TRUSTED; used for k >= 8, where one whole-gate LRAT proof is too big).

    python level1.py K P            prints the number of nodes / cover leaves / dead leaves

Let T1 = {s in 1..H : (K+1)s >= P}, H = (P-1)/2, and for a class c in 2..H let
cov[c] = {s in T1 : (K+1) dist_P(s c) < P}.  Call Y, a subset of {2..H}, a COVER when
|Y| <= K-1 and the union of cov[c] over c in Y is all of T1.  The y-part of every model of the
gate CNF G(K, L, P) is a cover: clauses A1 say exactly 'cov of Y meets every s in T1' and
A2 says |Y| <= K-1 (gen_cnf.py, block A).

THEOREM.  leaves(K, P) returns pairs (I, E) of disjoint class sets such that every cover Y
satisfies  I <= Y  and  Y & E = {}  for at least one returned pair.  Hence, if the CNF
G /\\ {y_c : c in I} /\\ {-y_c : c in E} is unsatisfiable for every returned pair, G is
unsatisfiable.

PROOF (depth-first search, node = (I, E), invariant: every cover consistent with the node is
consistent with some returned leaf below it).  Let U = T1 minus the union of cov[I].
  * U empty: return (I, E) itself.
  * |I| = K-1: a consistent cover Y contains I and has |Y| <= K-1, so Y = I, but I does not
    cover U: no consistent cover.
  * Counting bound: a consistent cover adds at most r = K-1-|I| classes outside E u I, each
    meeting at most m = max_{c not in E u I} |cov[c] & U| elements of U; if r*m < |U| none
    covers U: no consistent cover.
  * Otherwise pick s in U (fewest candidates; any choice is valid).  A consistent cover
    contains some candidate c in cand = {c not in E : s in cov[c]} (classes of I miss s,
    since s is in U).  With cand = (c_1..c_n) in increasing order and j the least index with
    c_j in Y, Y is consistent with child j = (I + {c_j}, E + {c_1..c_{j-1}}).  If n = 0 there
    is no consistent cover.   QED
Deterministic: the leaf list (and its order) depends only on (K, P).
"""
import sys


def fold(r, n):
    r %= n
    return min(r, n - r)


def leaves(K, P, stats=None):
    H = (P - 1) // 2
    T1 = [s for s in range(1, H + 1) if (K + 1) * s >= P]
    bit = {s: 1 << i for i, s in enumerate(T1)}
    classes = list(range(2, H + 1))
    cov = {c: sum(bit[s] for s in T1 if (K + 1) * fold(s * c, P) < P) for c in classes}
    full = (1 << len(T1)) - 1
    out = []
    st = stats if stats is not None else {}
    st.update(nodes=0, dead=0)

    def dfs(I, E, covered):
        st['nodes'] += 1
        U = full & ~covered
        if U == 0:
            out.append((tuple(I), tuple(sorted(E))))
            return
        if len(I) == K - 1:
            st['dead'] += 1
            return
        free = [c for c in classes if c not in E and c not in I]
        r = K - 1 - len(I)
        if r * max((bin(cov[c] & U).count('1') for c in free), default=0) < bin(U).count('1'):
            st['dead'] += 1
            return
        best = None
        m = U
        while m:
            low = m & -m
            m ^= low
            cand = [c for c in free if cov[c] & low]
            if best is None or len(cand) < len(best):
                best = cand
                if not cand:
                    break
        if not best:
            st['dead'] += 1
            return
        Ej = set(E)
        for c in best:
            dfs(I + [c], Ej, covered | cov[c])
            Ej = Ej | {c}

    dfs([], set(), 0)
    st['leaves'] = len(out)
    return out


if __name__ == '__main__':
    K, P = int(sys.argv[1]), int(sys.argv[2])
    st = {}
    lv = leaves(K, P, st)
    print(K, P, st)

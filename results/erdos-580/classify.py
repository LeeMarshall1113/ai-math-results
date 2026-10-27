"""Classify the trees on m vertices (k = m-1 edges) by the published exact LKS results, for a
host of order n = 2m:
  D  diameter <= 5                                    Piguet-Stein 2008, Thm 2
  P  path, or path + star (centre on a path vertex)   Bazgan-Li-Wozniak 2000
  C  T(k, l, c): path of length k-l plus stars with a, b leaves (a+b = l) centred at two path
     vertices at distance c, with l >= c and (c even or l + c >= floor(n/2) + 1)
                                                      Piguet-Stein 2008, Thm 3 (as published)
  E  none of these -> needs the SAT certificate.
usage: python classify.py m
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lks import trees, canon, core, g6


def cn(edges, n):
    adj = {i: [] for i in range(n)}
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    return canon(adj)


def diameter(edges, n):
    adj = {i: [] for i in range(n)}
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    def far(s):
        d = {s: 0}; q = [s]
        for u in q:
            for w in adj[u]:
                if w not in d:
                    d[w] = d[u] + 1; q.append(w)
        u = max(d, key=d.get); return u, d[u]
    u, _ = far(0)
    return far(u)[1]


def classify(m, host_n=None):
    k = m - 1
    n = host_n or 2 * m
    P, C = {}, {}
    for l in range(0, k + 1):
        L = k - l                        # path length
        for v in range(L + 1):           # path + star with l leaves at v
            e = [(i, i + 1) for i in range(L)] + [(v, L + 1 + j) for j in range(l)]
            P.setdefault(cn(e, m), dict(path_length=L, centre=v, leaves=l))
        for c in range(1, L + 1):
            if not (l >= c and (c % 2 == 0 or l + c >= n // 2 + 1)):
                continue
            for p in range(0, L - c + 1):
                for a in range(0, l + 1):
                    e = [(i, i + 1) for i in range(L)]
                    nx = L + 1
                    for _ in range(a):
                        e.append((p, nx)); nx += 1
                    for _ in range(l - a):
                        e.append((p + c, nx)); nx += 1
                    C.setdefault(cn(e, m), dict(path_length=L, l=l, c=c, centres=[p, p + c],
                                                leaves=[a, l - a]))
    rows = []
    for i, e in enumerate(trees(m)):
        c = cn(e, m)
        d = diameter(e, m) if m > 1 else 0
        if d <= 5:
            r = ("D", {"diameter": d})
        elif c in P:
            r = ("P", P[c])
        elif c in C:
            r = ("C", C[c])
        else:
            r = ("E", {"diameter": d})
        cv, ce, mk = core(e, m)
        rows.append(dict(index=i, g6=g6(e, m), route=r[0], witness=r[1],
                         core_size=len(cv), marked=len(mk)))
    return rows


if __name__ == "__main__":
    m = int(sys.argv[1])
    rows = classify(m)
    from collections import Counter
    print("m", m, "trees", len(rows), Counter(r["route"] for r in rows))
    for r in rows:
        if r["route"] == "E":
            print(" E", r)
    with open(f"classification_m{m}.json", "w") as f:
        json.dump(rows, f, indent=1)

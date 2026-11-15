"""Stage 2 of the planar classification: canonical forms with nauty (pynauty), run in WSL:

   wsl -d Ubuntu-24.04 -e python3 planar_classify.py

Reads results/planar_k<k>_R<i>.json (all solutions of the stage-1 enumeration), builds the
coloured incidence graph of PG(2,9) (points / lines get different colours, so the automorphism
group of the graph is PGammaL(3,9)), colours K-points, non-K points, 5-secants, 2-secants
separately and groups the solutions by nauty's canonical certificate.
Writes results/planar_classes.json : one entry per isomorphism class with a representative,
the stabiliser order in PGammaL(3,9), and where it was found.
"""
import glob
import json
import os
import sys
from collections import Counter, defaultdict

import pynauty

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pg import PG

HERE = os.path.dirname(os.path.abspath(__file__))


def incidence_graph(P, K):
    S = set(K)
    n = P.N
    nl = len(P.lines)
    adj = {v: [] for v in range(n + nl)}
    cnt = []
    for j, L in enumerate(P.lines):
        c = sum(1 for p in L if p in S)
        cnt.append(c)
        for p in L:
            adj[p].append(n + j)
            adj[n + j].append(p)
    col = [set(), set(), set(), set()]
    for p in range(n):
        col[0 if p in S else 1].add(p)
    for j in range(nl):
        col[2 if cnt[j] == 5 else 3].add(n + j)
    col = [c for c in col if c]
    return pynauty.Graph(n + nl, directed=False, adjacency_dict=adj, vertex_coloring=col), cnt


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else "planar_k*_R*.json"
    outname = sys.argv[2] if len(sys.argv) > 2 else "planar_classes.json"
    P = PG(3, 9)
    classes = {}
    found_in = defaultdict(set)
    nsol = Counter()
    for path in sorted(glob.glob(os.path.join(HERE, "..", "results", pattern))):
        d = json.load(open(path))
        k = d["k"]
        tag = os.path.basename(path)
        for K in d["solutions"]:
            g, cnt = incidence_graph(P, K)
            cert = pynauty.certificate(g)
            key = (k, cert)
            if key not in classes:
                gens, g1, g2, orbs, norb = pynauty.autgrp(g)
                size = int(round(g1 * 10 ** g2))
                classes[key] = dict(k=k, rep=list(K), aut_order=size, lines=dict(Counter(cnt)))
            found_in[key].add(tag)
            nsol[key] += 1
    out = []
    for key, c in sorted(classes.items(), key=lambda kv: (kv[1]["k"], -kv[1]["aut_order"])):
        c = dict(c)
        c["found_in"] = sorted(found_in[key])
        c["solutions_seen"] = nsol[key]
        out.append(c)
    json.dump(out, open(os.path.join(HERE, "..", "results", outname), "w"))
    for k in (26, 35):
        cl = [c for c in out if c["k"] == k]
        print("k=%d : %d isomorphism classes under PGammaL(3,9)" % (k, len(cl)))
        for c in cl:
            print("   |Aut| = %-6d lines %s  found in %s (%d labelled solutions)" % (
                c["aut_order"], c["lines"], c["found_in"], c["solutions_seen"]))


if __name__ == "__main__":
    main()

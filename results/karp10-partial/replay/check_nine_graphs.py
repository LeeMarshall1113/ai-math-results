"""check_nine_graphs.py -- the nine graphs of Lemma Z5 are what the note says they are (standard library only; run by replay_all.py).

The certificates certs/lemma_z5.json and certs/lemma_z5_euler.json list the claimed-in and claimed-out graphs as raw edge lists.  The note names them
(C10, C5 with pendant edges, Petersen, prism, K5 with pendant edges, C5 with a triangle on each edge, the square of C10, (K5 u C5)+perfect matching,
K5 x K2).  Degree sequences do not determine a graph, so this script BUILDS each named graph from its definition and checks, by explicit isomorphism
tests (an embedding with equal edge count), that

  * the nine claimed-in graphs of lemma_z5_euler.json are pairwise non-isomorphic and match the nine constructions one to one;
  * the nine claimed-out graphs are exactly the complements of the nine constructions, one to one;
  * the seven claimed-in graphs of lemma_z5.json are seven of the nine constructions, and its claimed-out graphs their complements.
"""
import itertools
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from replay_lib import *  # noqa: F401,F403


def cyc(vs):
    return [(vs[t], vs[(t + 1) % len(vs)]) for t in range(len(vs))]


def circ(S):
    return mask_of([(i, (i + s) % 10) for i in range(10) for s in S])


K5a = list(itertools.combinations(range(5), 2))
K5b = list(itertools.combinations(range(5, 10), 2))
SPOKES = [(i, i + 5) for i in range(5)]
NINE = {
    "C10": circ([1]),
    "C5 with a pendant edge at each vertex": mask_of(cyc([0, 1, 2, 3, 4]) + SPOKES),
    "Petersen graph": mask_of(cyc([0, 1, 2, 3, 4]) + SPOKES + [(5, 7), (7, 9), (9, 6), (6, 8), (8, 5)]),
    "pentagonal prism C5 x K2": mask_of(cyc([0, 1, 2, 3, 4]) + cyc([5, 6, 7, 8, 9]) + SPOKES),
    "K5 with a pendant edge at each vertex": mask_of(K5a + SPOKES),
    "C5 with a triangle on each edge": mask_of(cyc([0, 1, 2, 3, 4]) + [(i, 5 + i) for i in range(5)] + [((i + 1) % 5, 5 + i) for i in range(5)]),
    "square of C10 (circulant +-1, +-2)": circ([1, 2]),
    "(K5 u C5) + perfect matching": mask_of(K5a + cyc([5, 6, 7, 8, 9]) + SPOKES),
    "K5 x K2": mask_of(K5a + K5b + SPOKES),
}


def degs(m):
    r = [popcount(x) for x in rows_of(m)]
    return " ".join(f"{d}^{r.count(d)}" for d in sorted(set(r), reverse=True))


def match(claims, targets):
    """returns the list of target indices matched by the claims; asserts each claim matches exactly one target"""
    hit = []
    for edges in claims:
        gm = mask_of([tuple(e) for e in edges])
        ids = [i for i, t in enumerate(targets) if isomorphic(t, gm)]
        assert len(ids) == 1, f"a claimed graph with {popcount(gm)} edges matches {len(ids)} of the constructions"
        hit.append(ids[0])
    return hit


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    names = list(NINE)
    ins = [NINE[n] for n in names]
    outs = [complement(m) for m in ins]
    for a, b in itertools.combinations(range(9), 2):
        assert not isomorphic(ins[a], ins[b]), f"{names[a]} and {names[b]} are isomorphic"
    print("  [1] the nine constructed graphs are pairwise non-isomorphic:")
    for n, m in zip(names, ins):
        print(f"      {n}: {popcount(m)} edges, degrees {degs(m)}")
    ok = True
    for fn, n_in in (("lemma_z5_euler.json", 9), ("lemma_z5.json", 7)):
        cert = json.load(open(os.path.join(here, "..", "certs", fn)))
        hin = match(cert["claims_in"], ins)
        hout = match(cert["claims_out"], outs)
        good = len(cert["claims_in"]) == n_in and len(set(hin)) == n_in and len(cert["claims_out"]) == n_in and len(set(hout)) == n_in and sorted(hin) == sorted(hout)
        if n_in == 9:
            good = good and sorted(hin) == list(range(9))
        print(f"  [2] {fn}: {len(cert['claims_in'])} claimed-in graphs match {len(set(hin))} distinct constructions, "
              f"{len(cert['claims_out'])} claimed-out graphs are the complements of {len(set(hout))} distinct constructions; same indices: {sorted(hin) == sorted(hout)}")
        ok = ok and good
    print("  RESULT:", "THE CLAIMED GRAPHS ARE EXACTLY THE NINE NAMED GRAPHS AND THEIR COMPLEMENTS" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

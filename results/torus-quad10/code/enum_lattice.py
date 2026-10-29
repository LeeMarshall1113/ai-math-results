"""Method A: every degree-4 quadrangulation of the torus is a quotient of the square
tiling Z^2 by a translation sublattice L (the deck group of the universal cover acts by
orientation-preserving fixed-point-free automorphisms of the {4,4} tiling, i.e. translations).
With 10 faces, [Z^2 : L] = 10.  Enumerate all 18 index-10 sublattices (Hermite normal form),
build the face lists, test the 'geometric' necessary conditions, and classify up to isomorphism.

Necessary conditions tested (all proved in PROOF-SKELETON.md, Lemma 1/2):
  (C1) each face has 4 distinct vertices;
  (C2) two distinct faces share at most 2 vertices, and if 2 then they form an edge of both;
  (C3) at every vertex v the 2d+1 = 9 vertices of its star (v, 4 neighbours, 4 opposite corners)
       are pairwise distinct.
"""
import itertools, json, sys
from canon import canon_form

N = 10
def sublattices(n):
    for a in range(1, n + 1):
        if n % a: continue
        d = n // a
        for b in range(a):
            # basis rows (a,0)?  use HNF columns: L = {(a*i + b*j, d*j)}
            yield ((a, 0), (b, d))

def quotient_map(L):
    (a, _), (b, d) = L
    # canonical representative of (x,y) mod L: y mod d, then x mod a after removing j*(b,d)
    def rep(x, y):
        j, y0 = divmod(y, d)
        x0 = (x - j * b) % a
        return (x0, y0)
    pts = sorted({rep(x, y) for x in range(a) for y in range(d)})
    idx = {p: i for i, p in enumerate(pts)}
    faces = []
    for (x, y) in pts:
        faces.append(tuple(idx[rep(*q)] for q in [(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)]))
    return faces

def check(faces):
    why = []
    for f in faces:
        if len(set(f)) < 4: why.append("C1"); break
    if why: return why
    edges = lambda f: {frozenset((f[i], f[(i + 1) % 4])) for i in range(4)}
    for f, g in itertools.combinations(faces, 2):
        s = set(f) & set(g)
        if len(s) > 2 or (len(s) == 2 and (frozenset(s) not in edges(f) or frozenset(s) not in edges(g))):
            why.append("C2"); break
    for v in range(N):
        star = set()
        cnt = 0
        for f in faces:
            if v in f:
                star |= set(f); cnt += 3
        if len(star) != 2 * (cnt // 3) + 1:
            why.append("C3"); break
    return why

def main():
  res = {}
  for L in sublattices(N):
      F = quotient_map(L)
      w = check(F)
      key = canon_form(F) if not w or w == ["C3"] else None
      print(L, "fails" if w else "PASSES", w, (key[:1] if key else ""))
      if not w:
          res.setdefault(canon_form(F), []).append(L)
  print("lattices:", len(list(sublattices(N))))
  print("isomorphism classes passing C1-C3:", len(res))
  for k, v in res.items():
      print("  class with", len(v), "lattices:", v, "faces:", json.dumps(quotient_map(v[0])))

if __name__=="__main__": main()

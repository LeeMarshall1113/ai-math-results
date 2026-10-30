# Proof outline

Definitions are those of GS 2009, p. 21; see README. Throughout, T is an acoptic polyhedral torus in R³ with 10
strictly convex quadrilateral faces.

## Lemma 1 (counts, intersections, stars)
1. 4F = 2E and V − E + F = 0 give E = 20 and V = 10.
2. **(C2)** Two faces meet in nothing, one common vertex or one common edge. The intersection of two convex faces
   is convex, so it has one component. Being acoptic, that component is a common vertex or a common edge.
3. **(C3)** Let v have degree d. Then v, its d neighbours and its d opposite corners are 2d+1 distinct points.
   * Two neighbours coincide: two edges with the same endpoints would be the same segment, contradicting the
     simplicity of the faces.
   * A neighbour w equals the opposite corner of a face Q at v: then vw is an edge, and also a diagonal of Q, so
     its relative interior lies in the relative interior of Q. The face carrying the edge vw then meets Q in a
     set that is not a common cell, contradicting (C2).
   * Two faces at v have the same opposite corner w: both contain the segment vw (by convexity) as a diagonal,
     which again contradicts (C2).
4. So 2d+1 ≤ V = 10, hence d ≤ 4. The degree sum is 2E = 40 = 4V, so every vertex has degree 4.

*Convexity is used in 2 and 3; the lemma is false for non-convex faces.*

## Lemma 2 (unique map)
1. Glue a unit square for each face. Every vertex has degree 4, so this gives a flat torus whose universal cover
   is the square tiling.
2. Deck transformations are fixed-point-free and orientation-preserving, hence translations. So T ≅ Z²/L with
   [Z² : L] = 10.
3. By (C3) the nine points of a 3×3 block are distinct mod L, so L ∩ {−2..2}² = {0}.
4. Write L = ⟨(a,0),(b,d)⟩ in Hermite normal form, with ad = 10 and 0 ≤ b < a.
   * a ∈ {1, 2}: (a,0) ∈ L lies in the box.
   * a = 5: some (b−5t, 2) lies in the box.
   * a = 10: b ∈ {0,1,2,8,9} puts a vector (·,1) in the box, and b ∈ {4,5,6} puts a vector (·,2) in the box.
     This leaves b ∈ {3, 7}, two mirror-image lattices.
5. With the labelling (x,y) ↦ x − 3y mod 10, the faces are Q_i = (i, i+1, i+8, i+7). This is {4,4}_(3,1).
6. Computer cross-checks:
   * `enum_lattice.py`: all 18 index-10 lattices give 1 class.
   * `enum_backtrack.py`: labelled backtracking from a fixed star, independent of steps 1–4, gives 2 labelled
     completions (a chiral pair), 1 class, and no non-orientable completion.

## Lemma 3 (AD)
* **Statement.** Let F and G be strictly convex faces with distinct planes and a common edge ab, and let c be any
  other vertex of F. Then the two vertices x, y of G off ab satisfy χ(a,b,c,x) = χ(a,b,c,y) ≠ 0.
* **Proof.**
  * G \ ab lies in one open half-plane of plane(G) bounded by the line ab.
  * plane(G) ∩ plane(F) = line ab, so that half-plane lies in one open half-space of plane(F).
  * χ(a,b,c,·) is the sign of the side of plane(F).

## Lemma 4 (GP)
For points a, b and four further points y₁..y₄, the real Plücker identity is
[ab y₁y₂][ab y₃y₄] − [ab y₁y₃][ab y₂y₄] + [ab y₁y₄][ab y₂y₃] = 0.
So the three sign products are either all zero or contain both signs.

## Lemma 5 (computer)
The CNF encoding GP (all 3150 relations), PL (χ(Q_i) = 0) and AD (all 40 ordered adjacent face pairs) for
{4,4}_(3,1) is unsatisfiable. See `certs/`.

## Conclusion
Lemmas 1–2 reduce T to {4,4}_(3,1). Any realisation's chirotope would satisfy Lemmas 3–4, contradicting Lemma 5.
Lemmas 3–5 use neither embeddedness nor the absence of intersections, which gives Theorem′.

## Remarks
* **Not a forced flattening.** With planarity alone (no AD), non-degenerate numerical solutions exist; with AD,
  none were found (`code/numtest_ad.py`).
* **8-face core.** Keeping all GP relations but only the PL and AD constraints of the 8 faces other than Q₀ and
  Q₃ (which share an edge) is already unsatisfiable. This is certificate-checked: `certs/m8_core.lrat` (from
  `encode.py`) and `certs/m8_indep.lrat` (from `encode2.py`) both give `s VERIFIED UNSAT` from cake_lpr.
* **Minimal core.** Within it, a group-minimal unsatisfiable core uses the 8 PL conditions, 16 AD conditions and 37 GP
  relations (`code/mus.log`). That minimality is solver-reported, not certificate-checked.

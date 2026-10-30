# No torus with ten convex quadrilateral faces

**Status: computer-assisted; LRAT certificates checked by cake_lpr (CakeML-verified) and a second LRAT checker; two independently written encoders agree; hostile AI referee: SURVIVES; no human peer review.**

## Statement and sources

**Definitions.** Grünbaum–Szilassi (GS), *Geometric realizations of special toroidal complexes*, Contrib.
Discrete Math. 4(1) (2009) 21–39, doi:10.55016/ojs/cdm.v4i1.61986, p. 21:
* Faces are planar simple polygons with no collinear adjacent edges. For convex faces this means *strictly*
  convex.
* The faces at each vertex form one circuit, and adjacent faces are not coplanar.
* The polyhedron is *acoptic*: every component of the intersection of two faces is a common vertex or edge.
  For convex faces, GS note this means the intersection is empty, one vertex or one edge.

**Theorem.** No acoptic polyhedral torus in R³ has exactly 10 faces that are all strictly convex quadrilaterals.

This is the F = 10 half of GS **Conjecture 3.2** ("There exist no quadrangulations, of the kind described in
Proposition 3.1, with F = 10 or 11").

**Theorem′ (stronger, about one map).** The quadrangulation {4,4}_(3,1) has no realisation in R³ with planar,
strictly convex faces and non-coplanar adjacent faces, *even if faces may intersect*.
* The faces are Q_i = (i, i+1, i+8, i+7) mod 10, and the graph is the circulant C₁₀(1,3).
* Already the 8 faces left after deleting two adjacent faces (Q₀ and Q₃) admit no such realisation. This is
  certificate-checked in the same way as the 10-face result (`m8_core`, `m8_indep` below).

**Corollary (conditional).** Assume two things:
* J. W. G. Seamons's 11-face torus is correct (*A Torus with Eleven Convex Quadrilateral Faces*, Zenodo,
  doi:10.5281/zenodo.23750303, 29 Oct 2026, not refereed). Our independent exact checker accepts it.
* GS Prop. 3.1 holds for F = 9 and F ≥ 12.

Then a torus of this kind with F convex quadrilateral faces exists **iff F = 9 or F ≥ 11**.

## Prior work

* **GS 2009**: Prop. 3.1 (F = 9 and F ≥ 12), Conjecture 3.2 (F = 10, 11), and Prop. 3.3 (overarching,
  non-convex examples for F = 10, 11).
* **B. Grünbaum**, *Acoptic polyhedra*, Contemp. Math. 223 (1999) 163–199.
  * Fig. 14 lists {4,4}_(3,1) among the maps for which it is not known whether they have acoptic realisations.
  * The map that Simutis (PhD thesis, UC Davis, 1977) proved has no convex realisation is the 9-face Ljubić
    torus (p. 13, Fig. 13). It is a different map.
* **J. W. G. Seamons**, Zenodo 23750303 (2026): settles F = 11 with an explicit witness.
* **J. Schwörbel** and **J. M. Wills**: flag-transitive toroidal maps have acoptic realisations with possibly
  non-convex faces. {4,4}_(3,1) is chiral (not flag-transitive), so this does not apply.
* **U. Brehm and W. Kühnel**, *Equivelar maps on the torus*, Europ. J. Combin. 29 (2008) 1843–1861.
* **Novelty search**, 2026-10-29: Zenodo API, arXiv API, GitHub, OSF, the openai/math CONTENTS.md and the web.
  It found no result for F = 10 (`NOVELTY.md`). Please tell us if we missed one.

## Method (details in `PROOF.md`)

1. **Counting.** 4F = 2E and V − E + F = 0 give V = 10 and E = 20.
2. **Star lemma** (uses convexity and acoptic intersections).
   * At a vertex of degree d, the 2d+1 vertices of its star are distinct: the vertex itself, its d neighbours
     and its d opposite corners.
   * So 2d+1 ≤ 10, hence d ≤ 4. The degree sum is 40, so every vertex has degree 4.
3. **Unique map.**
   * A degree-4 quadrangulation of the torus is Z²/L, the square grid modulo a lattice L of index 10.
   * The star lemma forces L ∩ {−2,…,2}² = {0}.
   * The Hermite normal form then leaves only L = ⟨(10,0),(3,1)⟩ and its mirror image, which is {4,4}_(3,1).
   * Checked twice by computer: a lattice enumeration, and an exhaustive labelled backtracking that does not use
     the lattice argument.
4. **Necessary sign conditions.** χ(p₁..p₄) is the sign of the 4×4 determinant [1 pᵢ]. Every realisation
   satisfies:
   * **GP**: the 3-term Grassmann–Plücker sign relations, including zeros (3150 relations);
   * **PL**: χ(Q_i) = 0 for each face;
   * **AD**: if faces F and G share the edge ab, the two other vertices x, y of G satisfy
     χ(a,b,c,x) = χ(a,b,c,y) ≠ 0. Reason: G is strictly convex and plane(G) ≠ plane(F).
5. **SAT.** GP + PL + AD is unsatisfiable: CaDiCaL, with an LRAT proof checked by cake_lpr.
   Self-intersection is not encoded, which is why Theorem′ holds.

## How to verify

**Requirements.**
* Python 3. The core checks use the standard library only; `numtest_ad.py` and `mus.py` need numpy, scipy and
  python-sat.
* CaDiCaL ≥ 1.9.
* Optionally, an LRAT checker such as cake_lpr.

**Command.**

```
CADICAL=/path/cadical LRAT_CHECK="/path/cake_lpr" PY=python3 bash code/check_all.sh
```

**What it does.**
* Runs both map enumerations.
* Runs `soundness_test.py`. This takes two genuine tori (a 9-face torus and Seamons's 11-face torus), checks
  them with the exact checker `geomcheck.py`, and requires their exact chirotopes to satisfy every generated
  clause of both encoders. It also requires their own CNFs to be SAT.
* Runs `mutation_test.py`, which corrupts tori and chirotopes and requires them to be rejected.
* Regenerates the five CNFs (three 10-face, two 8-face) and compares their SHA-256 hashes with `certs/CHECKSUMS.sha256`.
* Solves each CNF and checks the LRAT proofs.

**Recorded run.** 2026-10-29, one thread, about 1.5 min, `certs/replay_2026-10-29.log`. All three LRATs gave
`s VERIFIED UNSAT` from cake_lpr and passed a second LRAT checker. The two 8-face certificates were added
later the same day and checked the same way (cake_lpr: `s VERIFIED UNSAT`; second checker: pass).

**Certificates.**

| CNF | contents | LRAT |
|---|---|---|
| `m10_core` | GP + PL + AD; this is the proof | 92 KB |
| `m10_full` | adds convexity, intersection patterns and mirror symmetry breaking | 320 KB |
| `m10_indep` | separately written encoder `encode2.py`; no auxiliary variables, different AD triple, brute-force GP tables | 31 MB |
| `m8_core` | as `m10_core`, without the PL and AD constraints of Q₀ and Q₃ (`encode.py ... GP,PL,AD 0,3`); all GP relations kept | 89 KB |
| `m8_indep` | as `m10_indep`, without the PL and AD constraints of Q₀ and Q₃ (`encode2.py ... --drop=0,3`) | 50 KB |

## Limits

**What the check covers.**
* Steps 4–5, given the CNF generators.
* The generators are only partly trusted, and partly cross-checked:
  * They are short. Two independent implementations agree.
  * The true chirotopes of two genuine tori satisfy every clause.
  * Mutated inputs are rejected.

**What it does not cover.**
* There is no formal proof that each generator faithfully encodes GP, PL and AD.
* The star lemma and the grid-quotient step are hand proofs. They are cross-checked only by the combinatorial
  enumerations.
* There is no human-readable proof of the SAT step. A group-minimal unsatisfiable core has 37 GP relations,
  16 AD conditions and 8 faces (`code/mus.log`). Its minimality is solver-reported, not certificate-checked;
  the unsatisfiability of the 8-face formula is certificate-checked (`m8_core`, `m8_indep`).
* Residual literature risk: Simutis's 1977 thesis could not be obtained. Grünbaum knew it and still listed
  {4,4}_(3,1) as unknown in 1999, and the F = 10 case as open in 2009.
* The corollary depends on Seamons's unrefereed witness.
* No human has refereed this work.

## Files

* `PROOF.md`: the argument.
* `NOVELTY.md`: the searches.
* `code/`: generators, enumerations, checkers, tests and `check_all.sh`.
* `certs/`: the LRAT proofs, `CHECKSUMS.sha256` and the replay log.
* `zenodo_note/`: the LaTeX source of the technical note and its metadata.

## AI disclosure

This work was generated with AI (Claude, from Anthropic), directed by Lee Marshall. It was checked by two
independent encoders and a formally verified LRAT checker, and attacked by a hostile AI referee. It has not been
through human peer review. Please check the computations yourself with the checkers provided.

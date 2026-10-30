# Referee report: "θ̄(K14) = 4" (thickness-k14)

Referee: independent Claude session, 2026-10-29. Not public. Times US Eastern.

## Verdict

**SURVIVES.** I found no fatal or fixable flaw that bears on the claim.

- **My independent encoder (`renc.py`) gives UNSAT for K14 with 3 layers.** It was written from
  the definitions before I read their encoder. It has no lex reflection and no sweep lemmas, and
  it uses my own cube partition. The run gave 232 roots and 352 leaves, every leaf LRAT
  "s VERIFIED UNSAT" by cake_lpr. `cover_mine.py --regen` reports COVER OK with 0 hash
  mismatches. Cost: 5457 CPU-s of CaDiCaL and 21.0 GB of LRAT, checked and deleted. Wall time
  16:50–17:44 EDT on 2 threads.
- Their encoding, Euler bound, case split, sweep lemmas and lex reflection are all valid. Their
  cover_check (with --regen) passes on both trees, and 26 of their leaves re-checked end to end
  all pass.
- Upper bound θ̄(K14) ≤ 4 is confirmed, and I built and exactly checked an explicit 4-layer drawing.
- Novelty: I found no prior claim. DEH 2000 and Seamons 2026 both leave K14 open.

Cosmetic items (none affects the result):
(a) `results/tree.jsonl` keeps 1417 unreachable UNSAT records under the removed rc-3 splits.
(b) cover_check trusts the runner's `cake` flag and does not re-hash the LRAT that cake_lpr read.
(c) The LRATs are deleted, so a third party must regenerate them (about 1.5 CPU-h).
(d) "Hull order = label order" is used only to pick the sweep edges, and it is correct.

## 1. Independent encoder (`renc.py`, written before reading `gen.py`)

**Model.** For a straight-line drawing of K_n, no three vertices are collinear (the middle one would
lie on an edge), so the points are in general position. A small rotation makes the x-coordinates
distinct and changes no crossing. Label the points 1..n by x. X(a,b,c), a<b<c, means "ccw".

| block | constraint | why every real drawing satisfies it |
|---|---|---|
| S | for a<b<c<d, the sequence X(abc),X(abd),X(acd),X(bcd) changes sign at most once (8 clauses: no alternating 3-subsequence) | rank-3 signotope axiom (Felsner–Weil), necessary for x-sorted point sets. Also checked empirically below on 41 real sets *after* the normal-form relabelling. |
| K | for every 4-set and **all three** pairings pq\|rs: Y ↔ [O(pqr)≠O(pqs)] ∧ [O(rsp)≠O(rsq)]; Y → no layer contains both | exact segment-crossing test in general position. I keep the pairing ab\|cd too and do not rely on x-ranges. |
| L | each edge gets ≥1 (and ≤1) layer | a partition |
| N | X(1,a,b) = + for all 1<a<b | see §1.1 |
| P | layer value-precedence along the lex edge order | rename layers by first use |
| E | for each layer k: #{e : e∈L_k or e uncrossed} + #{hull vertices} ≤ 3n−3; plus its sum over k, (t−1)·#U + t·h ≤ t(3n−3) − C(n,2); H_v ↔ v lies in no triangle, U_e ↔ e crosses no edge (full equivalences), W_ek ≥ C_ek ∨ U_e; my own Sinz counter; units H_1, H_2, H_n | see §3. H_2 and H_n need no geometry: "2 inside abc" with 1∉{a,b,c} forces X(2ab),X(2ac),X(2bc),X(abc) = s,−s,s,s (two sign changes), and with a = 1 it contradicts N; the case of n is symmetric |

No lex/reflection breaking and no sweep lemmas. *Lesson from two aborted attempts:* with only
one-sided indicators (H, U forced true and never forced false) and only the per-layer bounds, the
cubes with ≥ 3 hull points did not finish within 300 s. The per-layer counters never see the sum,
and an "H true" literal does not make v a hull vertex. With full equivalences and the summed bound
they refute in under 0.1 s. Their h5 and ge6 cases are refuted by unit propagation alone. I traced
why, and the argument is valid: h ≥ 5 hull edges are uncrossed, and so are (1,2) and (n−1,n),
because x-sorted segments with disjoint x-ranges cannot cross. That gives 2U + 3h ≥ 2·6 + 15 > 26. Cubes: my own partition of the hull indicators
H_3..H_13 (a binary tree over H_3, H_4, … that stops at the third positive literal; 232 roots,
checked by brute force over all 2^11 patterns to be a partition), then adaptive 16-way splits on
orientation variables after a 300 s timeout. Each leaf is solved by CaDiCaL 3.0.1 with LRAT,
checked by cake_lpr and deleted; SHA-256 values are kept in `runs/k14t3/tree.jsonl`.
`cover_mine.py 14 3 runs/k14t3 --regen` re-checks the partition, the tree, the cake flags and every
leaf CNF hash.

### 1.1 Soundness of the symmetry breaking

Soundness needs only this: **every real drawing has at least one image that satisfies the CNF**.
The brief asks for the normal form to hold for abstract signotopes, but that is not necessary. The
direction used is "UNSAT over abstract ⇒ UNSAT over real". A symmetry-breaking constraint therefore
has to keep at least one representative of each *realisable* drawing, and nothing more. I checked
it at that level. I also checked that each map used sends real drawings to real drawings.

- **Normal form (N).** Let p be a hull vertex and q a hull neighbour. Take p′ just outside the hull,
  near p on the extension of qp, so that around p′ the point p comes first and the others keep their
  angular order around p. Map to infinity a line through p′ that misses the hull. Segments go to
  segments and crossings are preserved. Orientations all keep their sign or all flip; a mirror
  fixes that. The pencil at p′ becomes the vertical lines, so the x-order is p first and then the
  angular order around p, which gives X(1,a,b)=+. **Empirical check:** `controls.py real`
  relabels each real set this way (p = lowest point, the rest sorted ccw), uses the *original*
  orientations, and requires the full CNF (S included) to be satisfiable under those units. This
  passed on Seamons' K13 and on 40 random exact drawings (n = 9–12, h from 3 up).
- **Authors' lex reflection (gen.lex_reflection).** I derived the reflection map myself. Mirror the
  set, keep p1 and relabel i ↦ σ(i) = n+2−i for i ≥ 2. Then X′(a,b,c) = X(σc,σb,σa), so for n = 14 it
  exchanges their sequences A = X(a,a+1,a+2) for a = 6..2 and B for a = 8..12 exactly. The encoding
  is a one-sided lex A ≤ B; the auxiliary e can only add constraints, and it is satisfiable when
  A ≤ B. So one of each mirror pair survives. **Sound.**
- **Colour precedence.** Sound for any assignment, in their version and in mine.

## 2. Controls (referee encoder)

`controls.py table` (pysat CaDiCaL 1.9.5, full encoding, with and without E):
K4/1 SAT, K5/1 UNSAT, K8/2 SAT, K9/2 UNSAT, K12/3 SAT, K13/3 SAT (40 s without E, 245 s with E; `controls_table.log`).
This reproduces the DEH table for n ≤ 12 (1 / 2 / 3), Seamons' K13 = 3 and θ̄(K9) = 3.
`controls.py real`: Seamons' K13 witness (Zenodo 10.5281/zenodo.22239253) checks exactly: 0
crossings within layers, and it is accepted with and without E. 40 random exact drawings are
accepted. `upper.py`: the authors' four `witnesses/k14_t4*.json` are exact 4-layer drawings and
are accepted by my encoder at t = 4.

## 3. Audit of the authors' lemmas and the Euler bound

- **Per-class Euler bound: valid.** A crossing-free straight-line graph on a point set with h hull
  points extends to a triangulation, and every triangulation has 3n−3−h edges (isolated points
  allowed). An edge of K_n that crosses no edge (an "uncrossed" edge, all hull edges included) can
  be added to any layer. So |L_k ∪ U| + h ≤ 3n−3 for every k. Summing over the t layers gives
  |E| + (t−1)|U| + t·h ≤ t(3n−3), i.e. 2U + 3h ≤ 26 for K14/3. That is exactly `gen.euler` mode 3:
  `excess` = U repeated t−1 times, H repeated t times, bound t(3n−3) − C(n,2). They also add a
  per-layer bound |L_k| + h ≤ 3n−3, which is weaker and valid. Duplicate literals in a Sinz
  counter are counted with multiplicity, which is correct. H is a full equivalence there, which is
  fine.
- **Hull case split: complete.** With N, 1, 2 and n are hull vertices. S ⊆ {3..13}: |S| ≤ 2 is
  handled exactly (1+11+55) and |S| ≥ 3 by the 165 covers "T ⊆ S". That makes 232 cases.
- **Sweep lemmas: valid, and not load-bearing.** For a hull edge BC there is a linear functional
  d with d(B) < d(C) < d(others) (the distance from BC, slightly tilted). (i) An interior point of
  a triangle is below the maximum vertex. (ii) If pq crosses xy, then max d(p,q) > min d(x,y): the
  crossing point is interior to pq and d(p) ≠ d(q). In the exact-S cases, consecutive hull labels
  are hull edges, because the angular order around the hull vertex 1 equals the hull's cyclic order.
  The run without sweeps is also UNSAT, so these lemmas carry no weight.
- **Omitted ab\|cd pairing**: dropping constraints is sound, and x-sorted segments with disjoint
  x-ranges cannot cross anyway.
- **cover_check.py**: re-run here. `results`: 232 cases, 307 verified leaves, 0 missing,
  `--regen` 0 hash mismatches. `results_nosweep` with `TK_NOSWEEP=1 --regen`: the same, 0
  mismatches. Note that `results/tree.jsonl` holds 1724 UNSAT records but only 307 are reachable.
  The rest sit under the removed rc-3 splits, as the README says. The reachable leaf set is
  identical in both runs, and 142 of the CNF hashes differ (the h-cases with sweeps), as expected.
- **20 random leaves end to end** (`sample_leaves.py`, seed 20261029, 10 from each tree): the CNF
  was regenerated from gen.py with matching SHA-256 (20/20), CaDiCaL returned UNSAT (20/20) and
  cake_lpr reported VERIFIED UNSAT (20/20). A uniform sample is mostly ge6/h5 leaves (2 were h3), so
  I added a stratified sample of 6 h3/h4 leaves: 6/6 hash match, UNSAT and cake-verified (up to 54 s each).

## 4. Upper bound

θ̄(K14) ≤ ⌈14/4⌉ = 4 by DEH Thm 2.1. An explicit construction also follows from θ̄(K_{n+1}) ≤
θ̄(K_n)+1: take Seamons' K13 points plus the point (1000003, 777781) and put its 13 edges, a star
with no crossings, into layer 4. `upper.py` checks this exactly (0 collinear triples, 0 crossings
within layers) and writes it to `k14_t4_referee.json`.

## 5. Literature and novelty (2026-10-29)

- DEH (JGAA 4(3) 2000; arXiv math/9910185), Table 1: values 1 (n ≤ 4), 2 (5–8), 3 (9–12),
  [3,4] for n = 13–14, 4 for 15–16. Open problem 1 asks for K13 and K14. The K15 lower bound is
  proved by triangulation counting plus Lemmas 4.2/4.3 (the authors' sweeps generalise these).
- Seamons, Zenodo 10.5281/zenodo.22239253 (2026-09-02): θ̄(K13) = 3 by an integer witness. Its
  abstract says DEH "left the exact values for K13 and K14 open" and it does not settle K14.
- Zenodo API (`"geometric thickness"`, sorted by most recent): only Seamons and an unrelated 2017
  record. arXiv API (title/abstract thickness ∧ complete, math.CO/cs.CG/cs.DM): nothing on K14.
  Recent related work does not touch it (Saturated drawings of geometric thickness k, 2503.03577;
  ∃R-completeness for multigraphs, 2312.05010; point-thickness of complete multipartite graphs,
  2601.01734). openai/math (origin/HEAD 2026-10-26): no geometric-thickness entry. Web search:
  K14 is still reported open. I found no prior claim of θ̄(K14) = 4.

## 6. Residual trust base

The signotope axiom S and the projective normal form are the only geometric facts the CNF relies
on, and both are standard. They are also checked on real sets after relabelling. On top of these
sit the Euler/uncrossed-edge bound (a short proof, above), the CaDiCaL→cake_lpr chain, and the
small Python generators, of which there are now two independent ones. Neither run keeps LRAT
files; regenerating them costs about 1 CPU-hour.

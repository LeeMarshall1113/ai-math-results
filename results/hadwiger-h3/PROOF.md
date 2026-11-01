# Proof of H₃ ≤ 13

## 0. Statement

For sets A, B ⊂ ℝⁿ, C(A, B) is the least number of translates of B that cover A. A convex body is a compact
convex set with non-empty interior. Hadwiger's covering number is H_n = max_K C(K, int K) over convex bodies
K ⊂ ℝⁿ; by Boltyanski it equals the illumination number.

**Theorem.** Every convex body K ⊂ ℝ³ is covered by 13 translates of its interior. Hence (compactness, §9)
K is covered by 13 homothets x_i + rK with 0 < r < 1, and H₃ ≤ 13.

Earlier bounds: 20 (Lassak 1998), 16 (Papadoperakis 1999), 14 (Prymak 2023). The lower bound is 8 (the cube).

The proof follows Prymak's framework (Papadoperakis's reduction plus an exact LP certificate over a
6-parameter space). It changes three things:

1. **Lemma M** (§3): minimality of the circumscribed box forces a long chord of K in a face-diagonal direction.
   One of the 13 covering bodies uses it. This removes the obstruction that caps Prymak's method at 14 (§10).
2. **Symmetry domain** (§6): Prymak's region D is not a fundamental domain for the cube symmetries. We use a
   corrected domain D\* and certify 570 extra start boxes, so that every orbit of configurations is covered.
3. **Strict containment** (§5): Prymak's Lemma 2.4(ii) fails at some box corners. Every leaf is instead checked
   with strict inequalities, which gives open covering sets directly from Lemma 2.4(i).

Everything in §§1–6 and §9 is machine-checked in Lean 4 (`lean/`), with the certificate as an axiom (§8).

## 1. Normalisation

**(N) Minimal box.** A parallelotope Π(U, α, β) = {x : α_m ≤ ⟨u_m, x⟩ ≤ β_m, m = 1, 2, 3} has volume
∏(β_m − α_m)/|det U|. For fixed U the smallest one containing K has volume ∏ w_K(u_m)/|det U|, where
w_K(u) = h_K(u) + h_K(−u) is the width. On unit triples this is continuous where det U ≠ 0 and tends to ∞ as
det U → 0 (w_K ≥ 2r for an inscribed ball of radius r), so a minimum exists. Minimality is affine-invariant, so
after an affine map the minimal box is C = [0, 1]³. Then w_K(e_i) = 1 for each i and every facet
F_ij = {x ∈ C : x_i = j} meets K.

**(P) Paired contact points** (Papadoperakis, Lemma 3(a)). One can choose q_i0 ∈ K ∩ F_i0 with
q_i1 = q_i0 + e_i ∈ K ∩ F_i1, for i = 1, 2, 3.

*Proof (axis 3).* Let π be the projection to (x₁, x₂) and A_3j = K ∩ F_3j. Suppose π(A_30) ∩ π(A_31) = ∅. Strict
separation gives g = (a₁, a₂, 0) and δ > 0 with max_{A_31}⟨g,·⟩ − min_{A_30}⟨g,·⟩ < −δ. For
u_ε = e₃ + εg the parallelotope with normals e₁, e₂, u_ε has det 1 and volume w_K(u_ε). The first-variation
bound h_K(u + εg) ≤ h_K(u) + ε(h_{F(K,u)}(g) + o(1)) gives w_K(u_ε) ≤ 1 + ε(−δ + o(1)) < 1 for small ε > 0,
contradicting minimality. A common point (s₁, s₂) gives q_30 = (s₁, s₂, 0) and q_31 = (s₁, s₂, 1). ∎

So the contact points are the columns of

        ⎛ 0   1   p₃  p₃  p₅  p₅ ⎞
  A_p = ⎜ p₁  p₁  0   1   p₆  p₆ ⎟ ,   V_p = {A_p e₁, …, A_p e₆},   O_p = conv V_p ⊆ K,
        ⎝ p₂  p₂  p₄  p₄  0   1  ⎠

for some configuration p ∈ [0, 1]⁶. Also ±e_l ∈ K − K (P′). The admissible p form a product of three compact
convex sets; any admissible p may be used.

## 2. Covering reduction

Let E be the 1-skeleton (12 edges) of C.

**Lemma G.** Let V contain one point of K on each facet of C. For every x ∈ bd K the set
W_x = x − cone(K − x) meets E ∪ V.

*Proof.* W_x − x is a convex cone. If x ∈ bd C put z = x; otherwise follow the ray from x away from some
k ∈ K \ {x} to its exit point z ∈ bd C, which lies in W_x. If z ∈ E, done. Otherwise z is in the relative
interior of one facet F_ij; let q ∈ V be the point on it. If q = z, done. Otherwise the ray z + μ(z − q), μ ≥ 0,
stays in the plane of F_ij and in W_x (since z − x and x − q lie in −cone(K − x)), and leaves F_ij through
its relative boundary, which is in E. ∎

**Lemma R** (Papadoperakis, Lemma 4, in the form for K's own p). If E ∪ V ⊆ ⋃_{i ≤ m}(v_i + X_i) with every
X_i ⊆ int K, then C(K, int K) ≤ m.

*Proof.* If y ∈ W_x ∩ (v + X), write y = x − t(k − x) with k ∈ K, t ≥ 0, and k′ = y − v ∈ int K. Then
x − εv = (1 − ε(1 + t))x + εtk + εk′ ∈ int K for 0 < ε ≤ 1/(1 + t). By compactness of bd K one ε₀ works for
a finite subcover, and then also for interior points (write y ∈ int K on a ray from a fixed o ∈ int K). So K is
covered by the m translates int K + εv_i. ∎

Neither lemma needs the pairing in (P); the pairing is used for the 6-parameter form of V_p, for O_p ⊆ K, and
in Lemma M.

## 3. Lemma M (long chord)

**Lemma M.** Let [0, 1]³ be a minimal-volume parallelotope containing K, and D = K − K. For each axis i, with
{j, k} the other two axes, there is s ∈ {+1, −1} such that t(e_j + s e_k) ∈ D for some
t ≥ τ := 1/2 + (√2 − 1)/6 = 0.569035…. So K contains a segment [z, z + τ(e_j + s e_k)].

*Proof (i = 1).* By (P′), ±e_l ∈ D, and D ⊆ [−1, 1]³. With u₁ = e₁, u₂ = (a, 1, 1), u₃ = (b, 1, −1) we have
|det U| = 2, and the slab parallelotope has volume h_D(u₁)h_D(u₂)h_D(u₃)/2. Minimality gives
h_D(a, 1, 1) · h_D(b, 1, −1) ≥ 2 for all a, b, so one of A = inf_a h_D(a, 1, 1), B = inf_b h_D(b, 1, −1) is ≥ √2;
say A (s = +1). Let (0, t, t) be where the ray meets bd D, with outer normal n = (a′, 1 + d, 1 − d) scaled so that
h_D(n) = 2t. From ±e₂, ±e₃ ∈ D: 2t ≥ 1 + |d|. By sublinearity and D ⊆ [−1, 1]³:
√2 ≤ h_D(a′, 1, 1) ≤ h_D(n) + h_D(0, −d, d) ≤ 2t + 2|d|. Hence 2t ≥ min_d max(1 + |d|, √2 − 2|d|) = 1 + (√2 − 1)/3. ∎

The sign s depends on K and is not known in advance, so the certificate treats both signs. For the octahedron O_c
(c = (1/2)⁶) the corresponding value is exactly 1/2. We found Lemma M neither in the literature nor in
Papadoperakis's or Prymak's papers; it is elementary and may be known.

## 4. Inner body

Fix an axis i, a sign s, w = e_j + s e_k, λ = 7/10 and t = 569/1000 < τ. If O_p ⊆ K and [z, z + tw] ⊆ K, then

  λz + (1 − λ) int O_p + λt[0, w] ⊆ int K,

because each point is a convex combination with weight 1 − λ > 0 on an interior point of K. The chord position z
only shifts the set, so it does not enter the covering problem, which stays 6-dimensional.

## 5. Boxes, leaves, and strict containment

For a box P = ∏[a_i, a_i + δ_i] ⊂ [0, 1]⁶ let Q_P = ⋂_{v vertex of P} O_v. (It equals the intersection over 8
sign-pattern vertices; the checkers use either form.) Let R_P(f) = {A_p e_f : p ∈ P}, a rectangle in face f.

**Prymak's Lemma 2.4(i).** Q_P ⊆ O_p for every p ∈ P. (Proved in Lean, `QP_sub_O`.) Hence
int Q_P ⊆ int O_p ⊆ int K.

**Structures.** A 13-structure has one translate a(v) per cube vertex (8), one translate b_f per face f ≠ f₀ (5),
the rectangle of face f₀ assigned to a(v₀), and an injective assignment τ of the five faces f ≠ f₀ to edges or
to nothing. Edge e = (u₀, u₁) not in the image of τ has one breakpoint r_e: the points u₀ and u₀ + r_e·e lie in
a(u₀), and u₀ + r_e·e and u₁ in a(u₁). An edge assigned to face f has two breakpoints r_e, s_e, with the middle
piece in b_f. By convexity the 13 translates then cover E. Rectangles are covered through their 4 corners.

**Leaves.** A leaf is a box P with stored rational translation vectors (integers over 10⁹) such that:

* *plain leaf:* the 13 sets T + int Q_P cover E ∪ ⋃_f R_P(f) with the given structure;
* *chord leaf* (axis i): for **each** s = ±1 there is a structure in which a(v₀) is replaced by
  T + μλt w_s + (1 − λ) int Q_P, μ ∈ [0, 1] (one extra variable μ per point assigned to it), and the 13 sets
  cover E ∪ ⋃_f R_P(f).

Each membership is a linear inequality c₀ + ⟨n, y − T⟩ **> 0** (resp. (1 − λ)c₀ + ⟨n, y − T − μλt w⟩ > 0) for
every row (c₀, n) of the exact H-representation of Q_P. Strict inequalities on all rows place y − T in int Q_P.

**Leaf soundness** (Lean `leaf_sound`). If P is a leaf and K has configuration p ∈ P, then E ∪ V_p is covered by
13 translates of subsets of int K: plain translates by §5 above; the chord translate by §4, using the sign s that
Lemma M gives for K on axis i. By Lemma R, C(K, int K) ≤ 13.

**Why strict containment.** Prymak passes from closed covers by Q_P to open sets with his Lemma 2.4(ii): some
translate of Q_P lies in int O_p. This fails at some box corners. Example (exact, `code/corner_example.py`):

* box P = ∏[a_i, a_i + 1/160] with a = (0, 0, 63/160, 5/8, 0, 57/80), corner p = (1/160, 1/160, 63/160, 5/8,
  1/160, 57/80);
* six facets of O_p contain a point of Q_P (exact vertices of Q_P);
* the inner normals (0, 80, 9040), (18080, −160, 0), (−480, 0, −496), (3680, −4960, 0) of four of them satisfy
  (215760/113)·n₁ + (104377/113)·n₂ + 34800·n₃ + 1·n₄ = 0;
* so for every translation u some touched facet has ⟨n_i, u⟩ ≤ 0, and Q_P + u is not inside int O_p.

With strict leaf inequalities Lemma 2.4(ii) is not used anywhere.

## 6. Symmetry: D, D\*, and the orbit extension

The 48 symmetries of the cube act on configurations as signed coordinate permutations (p′_m = p_σ(m) or
1 − p_σ(m)). Put s_x = p₃ + p₅, s_y = p₁ + p₆, s_z = p₂ + p₄. The reflection in axis i sends s_i to 2 − s_i and
fixes the other two sums; it flips (p₁, p₆) for y, (p₂, p₄) for z, (p₃, p₅) for x. A coordinate permutation
permutes (s_x, s_y, s_z).

**Prymak's region** D = {p₁, p₂, p₃ ≤ 1/2, s_z ≤ s_y ≤ s_x} **is not a fundamental domain.** The condition
p₁, p₂, p₃ ≤ 1/2 uses up the three reflections, and permutations order |s_i − 1|, not s_i. Exact example:

  p = (459/500, 1481/5000, 1689/10000, 9529/10000, 9359/10000, 3079/10000)

has no image in D among its 48 images, and none in the closed union of the 39,380 start boxes of Prymak's
U^D_20 (`code/symmetry_gap.py`; Lean `DPrymak_not_fundamental`). For uniform random p, about 2.7% have no image
in that union.

**Corrected domain.** D\* = {p₁, p₂, p₃ ≤ 1/2, |s_z − 1| ≤ |s_y − 1| ≤ |s_x − 1|} meets every orbit: sort the
|s_i − 1| by a permutation, then reflect each of p₁, p₂, p₃ that exceeds 1/2; the reflections do not change any
|s_i − 1|. (Lean `exists_sym_mem_DStar`.)

**Orbit extension.** Every p lies in some closed 1/10 grid box B, and each symmetry maps 1/10 grid boxes to
1/10 grid boxes. So it suffices that every orbit of 1/10 boxes (10⁶ boxes) contains a certified start box. The
39,380 old start boxes miss 570 orbits. One representative of each (the lexicographically smallest box of the
orbit whose interior meets int D\*) was certified with the same pipeline (`fullcert_ext/`). The 570 were derived
twice, independently (`code/newmethod/dstar_domain.py`, `code/replay/orbitboxes.py`), with identical results;
the finite orbit check is also done in Lean by `native_decide` (`orbit_cover`).

The certificate itself need not be symmetric. Given K with configuration p, pick g with g(p) in a certified box
and apply the leaf to gK: gK has the same minimal box, its paired contacts give g(p), and Lemma M holds for gK
on every axis. Certifying all of D\* instead would need about 35,000 more start boxes.

## 7. The certificate

* Start boxes: the 39,380 boxes of side 1/10 meeting Prymak's U^D_20 (`certs/startboxes_main.json`), and the 570
  extension boxes (`certs/startboxes_ext.json`).
* Each start box is bisected on its longest side (first index on ties) down to side ≥ 1/640, until a 13-structure
  is found for the leaf.
* Leaves: 7,473,755 (main) + 87,828 (extension) = 7,561,583. Chord leaves: 999, all in the 8 start boxes around
  c = (1/2)⁶; all other leaves, including every extension leaf, are plain.
* Data: 197 + 57 gzipped JSON-lines files, 2.1 GB, one line per start box with every leaf's box, structure and
  rational solution (both signs for chord leaves). Hashes in `certs/SHA256SUMS.*`; files on Zenodo.
* Production: float LPs (HiGHS) find candidate solutions, rounded to multiples of 10⁻⁹. Several search policies
  were used over the run (`certs/POLICY.*.md`); policy affects which structure is found, never soundness, since
  every leaf is checked exactly.

## 8. Checking

| check | what it verifies | result |
|---|---|---|
| `newmethod/certcheck_fast.py` | every leaf, exact integer arithmetic, non-strict; tiling of every start box by a bisection tree; start-box set = U^D_20 | all 7,473,755 + 87,828 leaves valid; 1% re-checked with the rational reference `certcheck.py`, 0 disagreements |
| `newmethod/certcheck_strict.py` | the same with strict inequalities: exactly the Lean predicate `LeafOK` | all 7,473,755 + 87,828 leaves valid |
| `newmethod/certcheck_ext.py` | orbit coverage: every 1/10 box of [0,1]⁶ has a symmetry image among the verified start boxes; symmetry action derived from labelled contact points | 0 boxes without an image |
| `replay/hcheck.py --strict` | written separately from the specification; ignores the stored structure and breakpoints; computes each translate's exact interval on each edge and requires the 13 intervals to cover it; rectangles by corners; Q_P from 8 sign-pattern vertices (64-vertex cross-check); int64 arithmetic with stated bounds and rigorous inner intervals | 7,561,583 leaves, 0 failures; 39,950/39,950 start boxes tiled; 0 orbits without a certified box |
| `lean/` | §§1–6 and §9, with `all_leaves_pass` (the certificate) as the only project axiom, plus one `native_decide` grid check | 0 `sorry`; `#print axioms` lists only these two plus Lean's standard three |

`certcheck_fast.py`, `certcheck_strict.py` and `certcheck_ext.py` share `certcheck.constraints()` (structure →
constraints) and `hw.py`; their independence from the certifier is in the arithmetic and the Q_P construction.
`hcheck.py` shares no code with them.

## 9. Homothets

If K ⊆ ⋃_{i ≤ 13}(x_i + int K), then K ⊆ ⋃(x_i + o + r(K − o)) for some o ∈ int K and some r < 1 depending on
K: the sets x_i + o + r(int K − o) increase to x_i + int K as r ↑ 1 (they are nested for convex K), and K is
compact. Lean: `homothets_of_interior`. No uniform r is claimed.

## 10. Why 13 needs the chord

At c = (1/2)⁶, O_c is the ℓ₁ ball of radius 1/2 about (1/2, 1/2, 1/2). The 8 cube vertices and the 6 face centres
are 14 points at pairwise ℓ₁ distance ≥ 1, so no translate of int O_c contains two of them (Prymak, Remark 2.3).
Every box containing c inherits this, so covering by translates of int Q_P alone cannot go below 14. Pairs at
distance exactly 1 differ by e_i or by ½(e_j ± e_k); the chord body is longer than O_c in the direction
e_j + s e_k (t > 1/2), and one chord translate merges one such pair. Away from c the bound 14 does not bind: plain
13-structures exist everywhere except near c, which is why only 999 leaves use the chord.

## 11. The H₃ ≤ 14 proof

The two points in §§5–6 apply to Prymak's proof of H₃ ≤ 14 as well.

* Its certificate covers U^D_20, so it misses the same 570 orbits. Our extension boxes are certified by plain
  leaves only (13 translates of int Q_P, no chord), which is inside Prymak's framework; they cover those orbits.
* Its passage to open sets uses Lemma 2.4(ii). The strict-containment argument of §5 replaces it for any leaf
  whose stored solution satisfies the leaf inequalities strictly. We have not re-run Prymak's 4.6M leaves under
  the strict test. When we re-solved 300 random leaf boxes of his certificate and the 8 leaves containing c with
  our own LP (exact check, non-strict), the float margins were positive in every case.

So with these two changes his argument goes through, provided his leaves pass a strict check (or are re-solved).
H₃ ≤ 14 also follows from H₃ ≤ 13.

## References

See `README.md`.

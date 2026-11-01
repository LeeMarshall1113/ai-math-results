import Hadwiger.Basic

/-!
# Convex-geometry helpers, and Prymak's Lemma 2.4(i)

* `dot u x = ∑ⱼ uⱼ xⱼ`; every continuous linear functional on `E3` is `dot u`.
* separation in the form "`x ∉ S` (closed convex) ⇒ some `dot u` strictly separates".
* `O c` is compact and convex.
* `Q_sub_O` (**Prymak, Lemma 2.4(i)**): for a box `P` of configurations,
  `Q_P := ⋂_{v ∈ vert P} O_v ⊆ O_p` for every `p ∈ P`.  Proved via support functions.
-/

namespace Hadwiger

open Set

/-- Dot product on `E3`. -/
def dot (u x : E3) : ℝ := ∑ j, u j * x j

lemma dot_add (u x y : E3) : dot u (x + y) = dot u x + dot u y := by
  simp [dot, mul_add, Finset.sum_add_distrib]

lemma dot_sub (u x y : E3) : dot u (x - y) = dot u x - dot u y := by
  simp [dot, mul_sub, Finset.sum_sub_distrib]

lemma dot_smul (u x : E3) (a : ℝ) : dot u (a • x) = a * dot u x := by
  simp [dot, Finset.mul_sum]; congr 1; ext j; ring

lemma dot_neg_left (u x : E3) : dot (-u) x = - dot u x := by
  simp [dot, Finset.sum_neg_distrib]

/-- Every continuous linear functional on `E3` is a dot product. -/
lemma strongDual_eq_dot (f : StrongDual ℝ E3) (x : E3) :
    f x = dot (fun j => f (Pi.single j 1)) x := by
  have hx : x = ∑ j, x j • (Pi.single j (1 : ℝ) : E3) := by
    funext i; simp [Finset.sum_apply, Pi.single_apply]
  conv_lhs => rw [hx]
  simp [map_sum, map_smul, dot, mul_comm]

/-- Strict separation of a point from a closed convex set by a dot product. -/
lemma exists_dot_sep {S : Set E3} (hc : Convex ℝ S) (hcl : IsClosed S) {x : E3} (hx : x ∉ S) :
    ∃ u : E3, ∀ y ∈ S, dot u y < dot u x := by
  obtain ⟨f, r, hf, hr⟩ := geometric_hahn_banach_closed_point hc hcl hx
  refine ⟨fun j => f (Pi.single j 1), fun y hy => ?_⟩
  rw [← strongDual_eq_dot, ← strongDual_eq_dot]
  exact (hf y hy).trans hr

/-- If a point is not in the interior of a convex set with non-empty interior, a non-zero
dot product attains its supremum over the set at that point. -/
lemma exists_dot_support {S : Set E3} (hc : Convex ℝ S) (hint : (interior S).Nonempty) {x : E3}
    (hx : x ∉ interior S) : ∃ u : E3, u ≠ 0 ∧ ∀ y ∈ S, dot u y ≤ dot u x := by
  obtain ⟨f, hf⟩ := geometric_hahn_banach_open_point hc.interior isOpen_interior hx
  set u : E3 := fun j => f (Pi.single j 1)
  have hfu : ∀ z, f z = dot u z := fun z => strongDual_eq_dot f z
  refine ⟨u, ?_, ?_⟩
  · intro hu
    obtain ⟨a, ha⟩ := hint
    have := hf a ha
    rw [hfu, hfu, hu] at this
    simp [dot] at this
  · intro y hy
    have hcl : y ∈ closure (interior S) := by
      rw [hc.closure_interior_eq_closure_of_nonempty_interior hint]; exact subset_closure hy
    have hclosed : IsClosed {z : E3 | f z ≤ f x} := isClosed_le f.continuous continuous_const
    have hsub : interior S ⊆ {z : E3 | f z ≤ f x} := fun z hz => (hf z hz).le
    have := (closure_minimal hsub hclosed) hcl
    simpa [hfu] using this

/-- A dot product bounded on the generators is bounded on the convex hull. -/
lemma dot_le_of_mem_convexHull {s : Set E3} {u x : E3} {M : ℝ} (hs : ∀ y ∈ s, dot u y ≤ M)
    (hx : x ∈ convexHull ℝ s) : dot u x ≤ M := by
  have hconv : Convex ℝ {z : E3 | dot u z ≤ M} := by
    intro a ha b hb α β hα hβ hαβ
    simp only [mem_ofPred_eq, dot_add, dot_smul] at *
    have hM : α * M + β * M = M := by rw [← add_mul, hαβ, one_mul]
    nlinarith [mul_le_mul_of_nonneg_left ha hα, mul_le_mul_of_nonneg_left hb hβ]
  exact convexHull_min hs hconv hx

lemma dot_lt_of_mem_convexHull {s : Set E3} {u x : E3} {M : ℝ} (hs : ∀ y ∈ s, dot u y < M)
    (hx : x ∈ convexHull ℝ s) : dot u x < M := by
  have hconv : Convex ℝ {z : E3 | dot u z < M} := by
    intro a ha b hb α β hα hβ hαβ
    simp only [mem_ofPred_eq, dot_add, dot_smul] at *
    rcases hα.lt_or_eq with h | h
    · have hM : α * M + β * M = M := by rw [← add_mul, hαβ, one_mul]
      nlinarith [mul_lt_mul_of_pos_left ha h, mul_le_mul_of_nonneg_left hb.le hβ]
    · subst h; simp at hαβ; subst hαβ; simpa using hb
  exact convexHull_min hs hconv hx

lemma V_finite (c : Cfg) : (V c).Finite := Set.finite_range _

lemma O_convex (c : Cfg) : Convex ℝ (O c) := convex_convexHull ℝ _

lemma O_compact (c : Cfg) : IsCompact (O c) := (V_finite c).isCompact_convexHull ℝ

lemma O_closed (c : Cfg) : IsClosed (O c) := (O_compact c).isClosed

lemma V_sub_O (c : Cfg) : V c ⊆ O c := subset_convexHull ℝ _

/-! ### Boxes of configurations and `Q_P` -/

/-- A (real) box of configurations in Prymak's coordinates: `lo k ≤ p k ≤ hi k`. -/
structure RBox where
  lo : Fin 6 → ℝ
  hi : Fin 6 → ℝ

/-- Membership of a configuration in a box. -/
def RBox.Mem (P : RBox) (c : Cfg) : Prop := ∀ k, P.lo k ≤ toP c k ∧ toP c k ≤ P.hi k

/-- The vertex of `P` selected by `ε` (`true` = upper end). -/
def RBox.vert (P : RBox) (ε : Fin 6 → Bool) : Fin 6 → ℝ := fun k => if ε k then P.hi k else P.lo k

/-- Prymak's `Q_P = ⋂_{v ∈ vert P} O_v`. -/
def QP (P : RBox) : Set E3 := {x | ∀ ε : Fin 6 → Bool, x ∈ O (ofP (P.vert ε))}

lemma QP_convex (P : RBox) : Convex ℝ (QP P) := by
  have : QP P = ⋂ ε : Fin 6 → Bool, O (ofP (P.vert ε)) := by ext; simp [QP]
  rw [this]; exact convex_iInter fun ε => O_convex _

/-- Pointwise comparison behind Lemma 2.4(i). -/
lemma dot_pt_vert_le (P : RBox) (c : Cfg) (hc : P.Mem c) (u : E3) (ib : Fin 3 × Bool) :
    dot u (pt (ofP (P.vert (fun k => decide (u (kCoord k) ≤ 0)))) ib) ≤ dot u (pt c ib) := by
  unfold dot
  apply Finset.sum_le_sum
  intro j _
  simp only [pt]
  split_ifs with hj
  · exact le_rfl
  · have hne : ib.1 ≠ j := fun h => hj h.symm
    obtain ⟨h1, h2⟩ := kPair_kIdx ib.1 j hne
    simp only [ofP, hne, ite_false, RBox.vert, h2]
    have hm := hc (kIdx ib.1 j)
    simp only [toP, h1, h2] at hm
    by_cases hu : u j ≤ 0
    · simp only [hu, decide_true, ite_true]
      exact mul_le_mul_of_nonpos_left hm.2 hu
    · simp only [hu, decide_false]
      exact mul_le_mul_of_nonneg_left hm.1 (le_of_not_ge hu)

/-- **Prymak, Lemma 2.4(i).**  `Q_P ⊆ O_p` for every configuration `p` in the box `P`. -/
theorem QP_sub_O (P : RBox) (c : Cfg) (hc : P.Mem c) : QP P ⊆ O c := by
  intro x hx
  by_contra hxO
  obtain ⟨u, hu⟩ := exists_dot_sep (O_convex c) (O_closed c) hxO
  set ε : Fin 6 → Bool := fun k => decide (u (kCoord k) ≤ 0)
  have hxv := hx ε
  have : dot u x < dot u x := by
    refine dot_lt_of_mem_convexHull ?_ hxv
    rintro y ⟨ib, rfl⟩
    exact lt_of_le_of_lt (dot_pt_vert_le P c hc u ib) (hu _ (V_sub_O c ⟨ib, rfl⟩))
  exact lt_irrefl _ this

/-- Interior version: `int Q_P ⊆ int O_p`. -/
theorem interior_QP_sub (P : RBox) (c : Cfg) (hc : P.Mem c) : interior (QP P) ⊆ interior (O c) :=
  interior_mono (QP_sub_O P c hc)

end Hadwiger

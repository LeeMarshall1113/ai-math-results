import Hadwiger.LemmaM

/-!
# Papadoperakis' Lemma 3(a): touching chords `e_i ∈ K - K`

If `[0,1]³` is a minimal-volume parallelotope of the convex body `K`, then for each axis `i`
there is `x ∈ K` with `x + e_i ∈ K`.

Proof (written for this formalisation; Papadoperakis' original was not available).  Suppose
not.  Then `A = K ∩ {xᵢ = 0}` and `B = (K ∩ {xᵢ = 1}) - eᵢ` are disjoint compact convex subsets of
the plane `xᵢ = 0`, so a linear form `ℓ` not involving `xᵢ` separates them: `ℓ < α` on `A`,
`ℓ > α` on `K ∩ {xᵢ = 1}`.  The shear `yᵢ = xᵢ - ε (ℓ(x) - α)` (determinant `1`) maps `K`
strictly inside the slab `0 < yᵢ < 1` for small `ε > 0`, by compactness.  Rescaling the slab gives
an affine map sending `K` into the cube with determinant `> 1`, contradicting minimality.
-/

namespace Hadwiger

open Set

/-- Compactness helper: a positive continuous function on the part of a compact set where a
continuous condition holds has a positive lower bound there. -/
lemma exists_pos_lb {K : Set E3} (hK : IsCompact K) (ℓ φ : E3 → ℝ) (hℓ : Continuous ℓ)
    (hφ : Continuous φ) (α : ℝ) (hpos : ∀ x ∈ K, α ≤ ℓ x → 0 < φ x) :
    ∃ δ > 0, ∀ x ∈ K, α ≤ ℓ x → δ ≤ φ x := by
  set C := K ∩ {x | α ≤ ℓ x}
  have hC : IsCompact C := hK.inter_right (isClosed_le continuous_const hℓ)
  rcases C.eq_empty_or_nonempty with h | h
  · refine ⟨1, one_pos, fun x hx hαx => ?_⟩
    have : x ∈ C := ⟨hx, hαx⟩
    rw [h] at this; simp at this
  · obtain ⟨x0, hx0, hmin⟩ := hC.exists_isMinOn h hφ.continuousOn
    exact ⟨φ x0, hpos x0 hx0.1 hx0.2, fun x hx hαx => hmin ⟨hx, hαx⟩⟩

/-- **Papadoperakis, Lemma 3(a).** -/
theorem touching_chord {K : Set E3} (hK : IsConvexBody K) (hmin : MinimalBox K) (i : Fin 3) :
    ∃ x ∈ K, x + e i ∈ K := by
  by_contra hno
  push Not at hno
  obtain ⟨hconv, hcomp, hint⟩ := hK
  have hne : K.Nonempty := hint.mono interior_subset
  have hcube := hmin.1
  have hc : ∀ x ∈ K, 0 ≤ x i ∧ x i ≤ 1 := by
    intro x hx; have := hcube hx
    simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at this
    exact ⟨this.1 i, this.2 i⟩
  -- the two faces, both put in the plane xᵢ = 0
  set A := {x | x ∈ K ∧ x i = 0}
  set B := {y | y + e i ∈ K ∧ y i = 0}
  have hei : (e i) i = 1 := by simp [e]
  have hAc : IsCompact A := hcomp.inter_right (isClosed_eq (continuous_apply i) continuous_const)
  have hBc : IsCompact B := by
    have : B = (fun x => x - e i) '' {x | x ∈ K ∧ x i = 1} := by
      ext y; simp only [B, mem_setOf_eq, mem_image]
      constructor
      · rintro ⟨h1, h2⟩; exact ⟨y + e i, ⟨h1, by simp [h2, hei]⟩, by simp⟩
      · rintro ⟨x, ⟨h1, h2⟩, rfl⟩; exact ⟨by simpa using h1, by simp [h2, hei]⟩
    rw [this]
    exact (hcomp.inter_right (isClosed_eq (continuous_apply i) continuous_const)).image
      (by fun_prop)
  have hAconv : Convex ℝ A := by
    intro x hx y hy a b ha hb hab
    exact ⟨hconv hx.1 hy.1 ha hb hab, by simp [hx.2, hy.2]⟩
  have hBconv : Convex ℝ B := by
    intro x hx y hy a b ha hb hab
    refine ⟨?_, by simp [hx.2, hy.2]⟩
    have := hconv hx.1 hy.1 ha hb hab
    convert this using 1
    have hT : (a + b) • e i = e i := by rw [hab, one_smul]
    calc a • x + b • y + e i = a • x + b • y + (a + b) • e i := by rw [hT]
      _ = _ := by module
  have hdisj : Disjoint A B := by
    rw [Set.disjoint_left]
    rintro x ⟨hxK, _⟩ ⟨hxe, _⟩
    exact hno x hxK hxe
  -- a separating form `w` (with `w i = 0`) and level `α`
  obtain ⟨f, u, v, hfA, huv, hfB⟩ :=
    geometric_hahn_banach_compact_closed hAconv hAc hBconv hBc.isClosed hdisj
  set w0 : E3 := fun j => f (Pi.single j 1)
  have hf : ∀ z, f z = dot w0 z := fun z => strongDual_eq_dot f z
  set w : E3 := Function.update w0 i 0
  have hwi : w i = 0 := by simp [w]
  have hdw : ∀ z, dot w z = dot w0 z - w0 i * z i := by
    intro z
    fin_cases i <;> simp [dot, w, Fin.sum_univ_three, Function.update_apply] <;> ring
  set α := u
  have hA' : ∀ x ∈ K, x i = 0 → dot w x < α := by
    intro x hx hxi
    rw [hdw, hxi, mul_zero, sub_zero, ← hf]; exact hfA x ⟨hx, hxi⟩
  have hB' : ∀ x ∈ K, x i = 1 → α < dot w x := by
    intro x hx hxi
    have hy : x - e i ∈ B := ⟨by simpa using hx, by simp [hxi, hei]⟩
    have h1 := hfB _ hy
    rw [hf] at h1
    have : dot w x = dot w0 (x - e i) := by
      rw [hdw, dot_sub, dot_e, hxi]; ring
    rw [this]; linarith
  -- compactness bounds
  set ℓ : E3 → ℝ := fun x => dot w x
  have hℓc : Continuous ℓ := by simp only [ℓ, dot]; fun_prop
  obtain ⟨δ0, hδ0, hlb0⟩ := exists_pos_lb hcomp ℓ (fun x => x i) hℓc (continuous_apply i) α
    (fun x hx hαx => by
      rcases (hc x hx).1.lt_or_eq with h | h
      · exact h
      · exfalso; have := hA' x hx h.symm; simp only [ℓ] at hαx; linarith)
  obtain ⟨δ1, hδ1, hlb1⟩ := exists_pos_lb hcomp (fun x => - ℓ x) (fun x => 1 - x i) hℓc.neg
    (continuous_const.sub (continuous_apply i)) (-α)
    (fun x hx hαx => by
      rcases (hc x hx).2.lt_or_eq with h | h
      · linarith
      · exfalso; have := hB' x hx h; simp only [ℓ] at hαx; linarith)
  obtain ⟨xR, hxR, hRmax⟩ := hcomp.exists_isMaxOn hne
    (show ContinuousOn (fun x => |ℓ x - α|) K by fun_prop)
  set R := |ℓ xR - α| + 1
  have hR0 : 0 < R := by positivity
  have hRb : ∀ x ∈ K, |ℓ x - α| ≤ R := fun x hx => (hRmax hx).trans (by simp [R])
  set ε := min δ0 δ1 / (2 * R)
  have hε0 : 0 < ε := by positivity
  have hεR : ε * R = min δ0 δ1 / 2 := by simp only [ε]; field_simp
  set y : E3 → ℝ := fun x => x i - ε * (ℓ x - α)
  have hyc : Continuous y := by simp only [y]; fun_prop
  have hy : ∀ x ∈ K, 0 < y x ∧ y x < 1 := by
    intro x hx
    have hb := hRb x hx
    rw [abs_le] at hb
    have hmin0 := min_le_left δ0 δ1; have hmin1 := min_le_right δ0 δ1
    have hmpos : 0 < min δ0 δ1 := lt_min hδ0 hδ1
    constructor
    · by_cases h : α ≤ ℓ x
      · have := hlb0 x hx h
        have : ε * (ℓ x - α) ≤ ε * R := mul_le_mul_of_nonneg_left hb.2 hε0.le
        simp only [y]; linarith
      · have : 0 < ε * (α - ℓ x) := mul_pos hε0 (by linarith)
        simp only [y]; linarith [(hc x hx).1]
    · by_cases h : ℓ x ≤ α
      · have := hlb1 x hx (by linarith)
        have : ε * (α - ℓ x) ≤ ε * R := mul_le_mul_of_nonneg_left (by linarith) hε0.le
        simp only [y] at *; linarith
      · have : 0 < ε * (ℓ x - α) := mul_pos hε0 (by linarith)
        simp only [y]; linarith [(hc x hx).2]
  obtain ⟨xm, hxm, hymin⟩ := hcomp.exists_isMinOn hne hyc.continuousOn
  obtain ⟨xM, hxM, hymax⟩ := hcomp.exists_isMaxOn hne hyc.continuousOn
  set m := y xm
  set W := max (y xM - m) (1/2)
  have hW0 : 0 < W := lt_of_lt_of_le (by norm_num) (le_max_right _ _)
  have hW1 : W < 1 := by
    apply max_lt
    · linarith [(hy xM hxM).2, (hy xm hxm).1]
    · norm_num
  -- the affine map
  let Mt : Matrix (Fin 3) (Fin 3) ℝ := fun r c =>
    if r = i then ((if c = i then 1 else 0) - ε * w c) / W else (if c = r then 1 else 0)
  let b : E3 := fun r => if r = i then (ε * α - m) / W else 0
  have hmap : ∀ x ∈ K, Matrix.toLin' Mt x + b ∈ cube := by
    intro x hx
    have hxc := hcube hx
    simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at hxc ⊢
    have hrow : ∀ r, (Matrix.toLin' Mt x + b) r = if r = i then (y x - m) / W else x r := by
      intro r
      show (Matrix.mulVec Mt x) r + b r = _
      simp only [Matrix.mulVec, dotProduct, Mt, b, Fin.sum_univ_three, y, ℓ, dot]
      fin_cases i <;> fin_cases r <;> simp [w] <;> field_simp <;> ring
    have h1 : m ≤ y x := hymin hx
    have h2 : y x ≤ y xM := hymax hx
    constructor <;> intro r <;> rw [hrow r] <;> split_ifs
    · exact div_nonneg (by linarith) hW0.le
    · exact hxc.1 r
    · rw [div_le_one hW0]; linarith [le_max_left (y xM - m) (1/2)]
    · exact hxc.2 r
  have hdet := hmin.abs_det_le _ _ hmap
  rw [LinearMap.det_toLin'] at hdet
  have hdetM : Mt.det = 1 / W := by
    rw [Matrix.det_fin_three]
    fin_cases i <;> simp [Mt, w] <;> field_simp <;> ring
  rw [hdetM, abs_of_pos (by positivity), div_le_one hW0] at hdet
  linarith

end Hadwiger

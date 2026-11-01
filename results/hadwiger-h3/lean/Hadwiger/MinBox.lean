import Hadwiger.LemmaM

/-!
# Existence of a minimal-volume parallelotope

Every convex body `K ⊂ ℝ³` has an affine image whose minimal-volume circumscribed parallelotope
is `[0,1]³`.

Proof: among affine maps `x ↦ M x + b` sending `K` into the cube (a closed and, because `K`
contains a ball, bounded set of pairs `(M, b)`), choose one maximising `|det M|` (compactness).
Its image of `K` has minimal box `[0,1]³`: a parallelotope `L(cube) + b'` containing it with
`|det L| < 1` would give an admissible map with larger determinant, and a degenerate `L` is
excluded because the image of `K` has positive volume.
-/

namespace Hadwiger

open Set MeasureTheory

lemma affine_image_interior (A : E3 ≃ᵃ[ℝ] E3) (S : Set E3) :
    A '' interior S = interior (A '' S) :=
  (A.toContinuousAffineEquiv.toHomeomorph).image_interior S

lemma IsConvexBody.affine_image {K : Set E3} (hK : IsConvexBody K) (A : E3 ≃ᵃ[ℝ] E3) :
    IsConvexBody (A '' K) := by
  refine ⟨hK.1.affine_image (A : E3 →ᵃ[ℝ] E3), hK.2.1.image A.continuous_of_finiteDimensional, ?_⟩
  rw [← affine_image_interior]; exact hK.2.2.image _

/-- The affine map `x ↦ M x + b`, `M` given by its entries. -/
def mapAt (M : Fin 3 → Fin 3 → ℝ) (b x : E3) : E3 := fun i => ∑ j, M i j * x j + b i

/-- Affine maps sending `K` into the cube. -/
def Adm (K : Set E3) : Set ((Fin 3 → Fin 3 → ℝ) × E3) := {p | ∀ x ∈ K, mapAt p.1 p.2 x ∈ cube}

/-- Determinant of `M`. -/
noncomputable def detM (M : Fin 3 → Fin 3 → ℝ) : ℝ := (Matrix.of M).det

lemma mapAt_eq (M : Fin 3 → Fin 3 → ℝ) (b x : E3) :
    mapAt M b x = Matrix.toLin' (Matrix.of M) x + b := by
  funext i; simp [mapAt, Matrix.mulVec, dotProduct]

lemma Adm_closed (K : Set E3) : IsClosed (Adm K) := by
  have : Adm K = ⋂ x ∈ K, (fun p : (Fin 3 → Fin 3 → ℝ) × E3 => mapAt p.1 p.2 x) ⁻¹' cube := by
    ext p; simp [Adm]
  rw [this]
  refine isClosed_biInter fun x _ => IsClosed.preimage ?_ isClosed_Icc
  unfold mapAt; fun_prop

lemma mem_cube_coord {y : E3} (h : y ∈ cube) (i : Fin 3) : 0 ≤ y i ∧ y i ≤ 1 := by
  simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at h
  exact ⟨h.1 i, h.2 i⟩

lemma norm_e (j : Fin 3) : ‖e j‖ = 1 := by
  apply le_antisymm
  · refine (pi_norm_le_iff_of_nonneg zero_le_one).2 fun i => ?_
    by_cases h : i = j
    · subst h; simp [e]
    · simp [e, Pi.single_apply, h]
  · have := norm_le_pi_norm (e j) j
    simpa [e] using this

lemma Adm_bounded {K : Set E3} (hK : IsConvexBody K) : Bornology.IsBounded (Adm K) := by
  obtain ⟨z, hz⟩ := hK.2.2
  obtain ⟨r, hr, hball⟩ := Metric.mem_nhds_iff.1 (mem_interior_iff_mem_nhds.1 hz)
  have hzK : z ∈ K := interior_subset hz
  have hpm : ∀ (j : Fin 3) (σ : ℝ), |σ| ≤ 1 → z + (σ * (r/2)) • e j ∈ K := by
    intro j σ hσ
    apply hball
    rw [Metric.mem_ball, dist_eq_norm, add_sub_cancel_left, norm_smul, norm_e, mul_one,
      Real.norm_eq_abs, abs_mul, abs_of_pos (by positivity : 0 < r/2)]
    calc |σ| * (r/2) ≤ 1 * (r/2) := mul_le_mul_of_nonneg_right hσ (by positivity)
      _ < r := by linarith
  set C := max (1/r) (1 + 3 * (1/r) * ‖z‖)
  rw [isBounded_iff_forall_norm_le]
  refine ⟨C, fun p hp => ?_⟩
  obtain ⟨M, b⟩ := p
  have hMij : ∀ i j, |M i j| ≤ 1/r := by
    intro i j
    have h1 := mem_cube_coord (hp _ (hpm j 1 (by simp))) i
    have h2 := mem_cube_coord (hp _ (hpm j (-1) (by simp))) i
    simp only [mapAt, Pi.add_apply, Pi.smul_apply, smul_eq_mul, e, Pi.single_apply, mul_add,
      Finset.sum_add_distrib] at h1 h2
    simp only [mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true] at h1 h2
    have hrM : |M i j * r| ≤ 1 := by
      rw [abs_le]; constructor <;> nlinarith
    rw [abs_mul, abs_of_pos hr] at hrM
    rw [le_div_iff₀ hr]; exact hrM
  have hMn : ‖M‖ ≤ C := by
    refine (pi_norm_le_iff_of_nonneg (by positivity)).2 fun i =>
      (pi_norm_le_iff_of_nonneg (by positivity)).2 fun j => ?_
    rw [Real.norm_eq_abs]; exact (hMij i j).trans (le_max_left _ _)
  have hbn : ‖b‖ ≤ C := by
    refine (pi_norm_le_iff_of_nonneg (by positivity)).2 fun i => ?_
    rw [Real.norm_eq_abs]
    have h0 := mem_cube_coord (hp z hzK) i
    simp only [mapAt] at h0
    have hs : |∑ j, M i j * z j| ≤ 3 * (1/r) * ‖z‖ := by
      calc |∑ j, M i j * z j| ≤ ∑ j, |M i j * z j| := Finset.abs_sum_le_sum_abs _ _
        _ ≤ ∑ _j : Fin 3, (1/r) * ‖z‖ := by
          apply Finset.sum_le_sum; intro j _
          rw [abs_mul]
          exact mul_le_mul (hMij i j) (by rw [← Real.norm_eq_abs]; exact norm_le_pi_norm z j)
            (abs_nonneg _) (by positivity)
        _ = 3 * (1/r) * ‖z‖ := by simp; ring
    rw [abs_le] at hs ⊢
    constructor <;> linarith [le_max_right (1/r) (1 + 3 * (1/r) * ‖z‖)]
  exact max_le hMn hbn

lemma detM_continuous : Continuous (fun p : (Fin 3 → Fin 3 → ℝ) × E3 => |detM p.1|) := by
  have : (fun p : (Fin 3 → Fin 3 → ℝ) × E3 => detM p.1) = fun p =>
      p.1 0 0 * p.1 1 1 * p.1 2 2 - p.1 0 0 * p.1 1 2 * p.1 2 1 - p.1 0 1 * p.1 1 0 * p.1 2 2 +
      p.1 0 1 * p.1 1 2 * p.1 2 0 + p.1 0 2 * p.1 1 0 * p.1 2 1 - p.1 0 2 * p.1 1 1 * p.1 2 0 := by
    funext p; simp [detM, Matrix.det_fin_three]
  exact continuous_abs.comp (this ▸ by fun_prop)

lemma Adm_pos {K : Set E3} (hK : IsConvexBody K) : ∃ p ∈ Adm K, 0 < |detM p.1| := by
  obtain ⟨R, hR⟩ := hK.2.1.isBounded.exists_norm_le
  set R' := max R 1
  have hR' : 0 < R' := lt_of_lt_of_le one_pos (le_max_right _ _)
  refine ⟨(fun i j => if i = j then 1 / (2 * R') else 0, fun _ => 1/2), ?_, ?_⟩
  · intro x hx
    simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply, mapAt, ite_mul, zero_mul,
      Finset.sum_ite_eq, Finset.mem_univ, ite_true]
    have hxi : ∀ i, |x i| ≤ R' := fun i =>
      ((Real.norm_eq_abs _) ▸ norm_le_pi_norm x i).trans ((hR x hx).trans (le_max_left _ _))
    constructor <;> intro i <;> have := abs_le.1 (hxi i)
    · rw [show 1 / (2 * R') * x i = x i / (2 * R') by ring]
      have : -(1/2) ≤ x i / (2 * R') := by rw [le_div_iff₀ (by positivity)]; linarith
      linarith
    · rw [show 1 / (2 * R') * x i = x i / (2 * R') by ring]
      have : x i / (2 * R') ≤ 1/2 := by rw [div_le_iff₀ (by positivity)]; linarith
      linarith
  · simp [detM, Matrix.det_fin_three]; positivity

/-- **Existence of a minimal-volume parallelotope.** -/
theorem exists_minimalBox' (K : Set E3) (hK : IsConvexBody K) :
    ∃ A : E3 ≃ᵃ[ℝ] E3, MinimalBox (A '' K) := by
  have hcpt : IsCompact (Adm K) := Metric.isCompact_of_isClosed_isBounded (Adm_closed K) (Adm_bounded hK)
  obtain ⟨p0, hp0, hpos⟩ := Adm_pos hK
  obtain ⟨⟨M, b⟩, hpA, hmax⟩ := hcpt.exists_isMaxOn ⟨p0, hp0⟩ detM_continuous.continuousOn
  have hMpos : 0 < |detM M| := lt_of_lt_of_le hpos (hmax hp0)
  have hdet : LinearMap.det (Matrix.toLin' (Matrix.of M)) ≠ 0 := by
    rw [LinearMap.det_toLin']; exact abs_pos.1 hMpos
  set Le := LinearMap.equivOfDetNeZero _ hdet
  let A : E3 ≃ᵃ[ℝ] E3 := Le.toAffineEquiv.trans (AffineEquiv.constVAdd ℝ E3 b)
  have hA : ∀ x, A x = mapAt M b x := by
    intro x; rw [mapAt_eq]; simp [A, Le, add_comm]
  refine ⟨A, ?_, ?_⟩
  · rintro _ ⟨x, hx, rfl⟩; rw [hA]; exact hpA x hx
  · intro L b' hsub
    rw [volume_parallelotope, volume_cube]
    by_cases hL : LinearMap.det L = 0
    · -- degenerate: the image of K has positive volume
      exfalso
      have hvol := measure_mono (μ := volume) hsub
      rw [volume_parallelotope, hL, abs_zero, ENNReal.ofReal_zero] at hvol
      have hint : (interior (A '' K)).Nonempty := (hK.affine_image A).2.2
      have hpos' : 0 < volume (interior (A '' K)) := isOpen_interior.measure_pos volume hint
      have := (measure_mono (μ := volume) interior_subset).trans hvol
      exact absurd (lt_of_lt_of_le hpos' this) (lt_irrefl _)
    · set Lq := LinearMap.equivOfDetNeZero L hL
      -- the admissible map  x ↦ L⁻¹ (A x - b')
      set f : E3 →ₗ[ℝ] E3 := (Lq.symm : E3 →ₗ[ℝ] E3).comp (Matrix.toLin' (Matrix.of M))
      set N : Fin 3 → Fin 3 → ℝ := fun i j => LinearMap.toMatrix' f i j
      have hNf : Matrix.toLin' (Matrix.of N) = f := by
        have : Matrix.of N = LinearMap.toMatrix' f := rfl
        rw [this, Matrix.toLin'_toMatrix']
      have hadm : (N, Lq.symm (b - b')) ∈ Adm K := by
        intro x hx
        obtain ⟨y, hy, hyx⟩ := hsub ⟨x, hx, rfl⟩
        rw [mapAt_eq, hNf]
        have : f x + Lq.symm (b - b') = y := by
          have h1 : Lq y = L y := rfl
          have hAx : A x = Matrix.toLin' (Matrix.of M) x + b := by rw [hA, mapAt_eq]
          apply Lq.injective
          have h2 : Lq (f x) = Matrix.toLin' (Matrix.of M) x := by simp [f]
          rw [map_add, h2, LinearEquiv.apply_symm_apply, h1]
          rw [hAx] at hyx
          simp only at hyx
          rw [show L y = (L y + b') - b' by abel, hyx]; abel
        rw [this]; exact hy
      have hle := hmax hadm
      simp only [mem_setOf_eq] at hle
      have hdN : detM N = (LinearMap.det L)⁻¹ * detM M := by
        have : detM N = LinearMap.det f := by
          rw [← hNf, LinearMap.det_toLin']; rfl
        rw [this, LinearMap.det_comp, LinearMap.det_toLin']
        have hsymm : LinearMap.det (Lq.symm : E3 →ₗ[ℝ] E3) = (LinearMap.det L)⁻¹ := by
          have h1 := LinearMap.det_comp (Lq.symm : E3 →ₗ[ℝ] E3) (Lq : E3 →ₗ[ℝ] E3)
          have h2 : (Lq.symm : E3 →ₗ[ℝ] E3).comp (Lq : E3 →ₗ[ℝ] E3) = LinearMap.id := by
            ext x; simp
          rw [h2, LinearMap.det_id] at h1
          have h3 : LinearMap.det (Lq : E3 →ₗ[ℝ] E3) = LinearMap.det L := rfl
          rw [h3] at h1
          field_simp; linarith
        rw [hsymm]; rfl
      rw [hdN, abs_mul, abs_inv] at hle
      have hLpos : 0 < |LinearMap.det L| := abs_pos.2 hL
      have : 1 ≤ |LinearMap.det L| := by
        by_contra hlt; push Not at hlt
        have : |detM M| < (|LinearMap.det L|)⁻¹ * |detM M| := by
          have : 1 < (|LinearMap.det L|)⁻¹ := (one_lt_inv₀ hLpos).2 hlt
          nlinarith
        linarith
      exact ENNReal.one_le_ofReal.2 this

end Hadwiger

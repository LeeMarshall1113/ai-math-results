import Hadwiger.Geometry

/-!
# The long-chord lemma (Lemma M)

`MinimalBox K`: `K ⊆ [0,1]³` and the cube has minimal volume among all parallelotopes
(affine images of the cube) containing `K`.

**Lemma M.**  If `K` is a convex body with minimal box `[0,1]³` and with touching chords
`±e_l ∈ K - K` (Papadoperakis' Lemma 3(a)), then for every pair `j ≠ k` there are a sign
`s = ±1` and `x ∈ K` with `x + τ (e_j + s e_k) ∈ K`, where `τ = 1/2 + (√2 - 1)/6`.

The proof follows the paper, but in a form that needs neither attained infima nor
boundary normals: if both points `τ(e_j ± e_k)` were outside `D = K - K`, separating
hyperplanes produce two directions `v±` with `v±_j = 1`, `v±_k = ±1` and width `< √2`;
the slab parallelotope with rows `e_m, v₊, v₋` (determinant `±2`) would then have volume
`< 1`.
-/

namespace Hadwiger

open Set MeasureTheory

/-- `e i`: the `i`-th unit vector. -/
def e (i : Fin 3) : E3 := Pi.single i 1

/-- `[0,1]³` is a parallelotope of minimal volume containing `K`. -/
def MinimalBox (K : Set E3) : Prop :=
  K ⊆ cube ∧ ∀ (L : E3 →ₗ[ℝ] E3) (b : E3), K ⊆ (fun y => L y + b) '' cube →
    volume cube ≤ volume ((fun y => L y + b) '' cube)

lemma volume_cube : volume cube = 1 := by
  simp [cube, Real.volume_Icc_pi]

/-- Volume of a parallelotope: `vol(L(cube) + b) = |det L|`. -/
lemma volume_parallelotope (L : E3 →ₗ[ℝ] E3) (b : E3) :
    volume ((fun y => L y + b) '' cube) = ENNReal.ofReal |LinearMap.det L| := by
  have : (fun y => L y + b) '' cube = (fun x => x + b) '' (L '' cube) := by
    rw [image_image]
  rw [this, image_add_right, measure_preimage_add_right, Measure.addHaar_image_linearMap,
    volume_cube, mul_one]

/-- Minimality in determinant form: if an affine map `x ↦ M x + b` sends `K` into the cube,
then `|det M| ≤ 1`. -/
lemma MinimalBox.abs_det_le {K : Set E3} (h : MinimalBox K) (M : E3 →ₗ[ℝ] E3) (b : E3)
    (hK : ∀ x ∈ K, M x + b ∈ cube) : |LinearMap.det M| ≤ 1 := by
  by_cases hd : LinearMap.det M = 0
  · simp [hd]
  set Me := LinearMap.equivOfDetNeZero M hd
  have hsub : K ⊆ (fun y => (Me.symm : E3 →ₗ[ℝ] E3) y + (- Me.symm b)) '' cube := by
    intro x hx
    refine ⟨M x + b, hK x hx, ?_⟩
    simp only [LinearEquiv.coe_coe, map_add]
    have : Me.symm (M x) = x := by
      have : M x = Me x := rfl
      rw [this, LinearEquiv.symm_apply_apply]
    rw [this]; abel
  have hv := h.2 _ _ hsub
  rw [volume_parallelotope, volume_cube] at hv
  have hdet : LinearMap.det (Me.symm : E3 →ₗ[ℝ] E3) = (LinearMap.det M)⁻¹ := by
    have h1 := LinearMap.det_comp (Me.symm : E3 →ₗ[ℝ] E3) (Me : E3 →ₗ[ℝ] E3)
    have h2 : (Me.symm : E3 →ₗ[ℝ] E3).comp (Me : E3 →ₗ[ℝ] E3) = LinearMap.id := by
      ext x; simp
    rw [h2, LinearMap.det_id] at h1
    have h3 : LinearMap.det (Me : E3 →ₗ[ℝ] E3) = LinearMap.det M := rfl
    rw [h3] at h1
    field_simp
    linarith
  rw [hdet, abs_inv] at hv
  have hpos : 0 < |LinearMap.det M| := abs_pos.mpr hd
  have := (ENNReal.one_le_ofReal).mp hv
  rw [le_inv_comm₀ zero_lt_one hpos, inv_one] at this
  exact this

/-- The difference body `D = K - K`. -/
def diffBody (K : Set E3) : Set E3 := {d | ∃ x ∈ K, ∃ y ∈ K, d = y - x}

lemma diffBody_convex {K : Set E3} (hK : Convex ℝ K) : Convex ℝ (diffBody K) := by
  rintro _ ⟨x1, hx1, y1, hy1, rfl⟩ _ ⟨x2, hx2, y2, hy2, rfl⟩ a b ha hb hab
  refine ⟨a • x1 + b • x2, hK hx1 hx2 ha hb hab, a • y1 + b • y2, hK hy1 hy2 ha hb hab, ?_⟩
  simp only [smul_sub]; abel

lemma diffBody_compact {K : Set E3} (hK : IsCompact K) : IsCompact (diffBody K) := by
  have : diffBody K = (fun p : E3 × E3 => p.2 - p.1) '' (K ×ˢ K) := by
    ext d; constructor
    · rintro ⟨x, hx, y, hy, rfl⟩; exact ⟨(x, y), ⟨hx, hy⟩, rfl⟩
    · rintro ⟨⟨x, y⟩, ⟨hx, hy⟩, rfl⟩; exact ⟨x, hx, y, hy, rfl⟩
  rw [this]
  exact (hK.prod hK).image (by fun_prop)

lemma diffBody_coord_le {K : Set E3} (hK : K ⊆ cube) {d : E3} (hd : d ∈ diffBody K) (i : Fin 3) :
    |d i| ≤ 1 := by
  obtain ⟨x, hx, y, hy, rfl⟩ := hd
  have h1 := hK hx; have h2 := hK hy
  simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at h1 h2
  simp only [Pi.sub_apply, abs_le]
  constructor <;> linarith [h1.1 i, h1.2 i, h2.1 i, h2.2 i]

lemma dot_e (u : E3) (i : Fin 3) : dot u (e i) = u i := by
  simp [dot, e, Pi.single_apply]

/-- The Lemma-M constant `τ = 1/2 + (√2 - 1)/6 = 0.569035…`. -/
noncomputable def tauM : ℝ := 1/2 + (Real.sqrt 2 - 1)/6

lemma sqrt2_facts : 1 < Real.sqrt 2 ∧ Real.sqrt 2 * Real.sqrt 2 = 2 ∧ Real.sqrt 2 < 3/2 := by
  refine ⟨?_, Real.mul_self_sqrt (by norm_num), ?_⟩
  · rw [show (1:ℝ) = Real.sqrt 1 by simp]; exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  · rw [show (3/2:ℝ) = Real.sqrt (9/4) by
      rw [show (9/4:ℝ) = (3/2)^2 by norm_num, Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)

/-- Separation step: if `τ (e_j + s e_k) ∉ D`, there is a direction `v` with `v_j = 1`,
`v_k = s` and `dot v < √2` on `D`. -/
lemma short_direction {K : Set E3} (hconv : Convex ℝ K) (hcomp : IsCompact K) (hcube : K ⊆ cube)
    (htouch : ∀ i, e i ∈ diffBody K ∧ -e i ∈ diffBody K) {j k : Fin 3} (hjk : j ≠ k) {s : ℝ}
    (hs : s = 1 ∨ s = -1) (hout : tauM • (e j + s • e k) ∉ diffBody K) :
    ∃ v : E3, v j = 1 ∧ v k = s ∧ ∀ d ∈ diffBody K, dot v d < Real.sqrt 2 := by
  obtain ⟨u, hu⟩ := exists_dot_sep (diffBody_convex hconv) (diffBody_compact hcomp).isClosed hout
  obtain ⟨hs1, hs2, hs3⟩ := sqrt2_facts
  have hss : s * s = 1 := by rcases hs with rfl | rfl <;> norm_num
  have hkj : k ≠ j := fun h => hjk h.symm
  have hval : dot u (tauM • (e j + s • e k)) = tauM * (u j + s * u k) := by
    rw [dot_smul, dot_add, dot_smul, dot_e, dot_e]
  set N := u j + s * u k with hN
  have hej := hu _ (htouch j).1; have hej' := hu _ (htouch j).2
  have hek := hu _ (htouch k).1; have hek' := hu _ (htouch k).2
  rw [hval] at hej hej' hek hek'
  rw [dot_e] at hej hek
  rw [show -e j = (-1 : ℝ) • e j by simp, dot_smul, dot_e] at hej'
  rw [show -e k = (-1 : ℝ) • e k by simp, dot_smul, dot_e] at hek'
  have htau : 0 < tauM := by unfold tauM; linarith
  have hNpos : 0 < N := by
    by_contra hneg; push Not at hneg
    have : tauM * N ≤ 0 := mul_nonpos_of_nonneg_of_nonpos htau.le hneg
    linarith
  -- normalised normal n = (2/N) u, δ = n_j - 1
  set n : E3 := (2 / N) • u with hn
  set δ : ℝ := n j - 1 with hδ
  have hnj : n j = 2 * u j / N := by simp [hn] <;> ring
  have hnk : n k = 2 * u k / N := by simp [hn] <;> ring
  have hsum : n j + s * n k = 2 := by
    rw [hnj, hnk]; field_simp; rw [hN]
  have hnd : ∀ d ∈ diffBody K, dot n d < 2 * tauM := by
    intro d hd
    have := hu d hd
    rw [hval] at this
    rw [hn, show ((2 / N) • u) = (2 / N) • u from rfl]
    have : dot ((2 / N) • u) d = (2 / N) * dot u d := by
      simp [dot, Finset.mul_sum]; congr 1; ext i; ring
    rw [this]
    have h2N : 0 < 2 / N := by positivity
    calc 2 / N * dot u d < 2 / N * (tauM * N) := mul_lt_mul_of_pos_left (hu d hd |>.trans_eq hval) h2N
      _ = 2 * tauM := by field_simp
  -- |n_j| < 2τ and |n_k| < 2τ
  have hnj_lt : |n j| < 2 * tauM := by
    have h1 := hnd _ (htouch j).1; have h2 := hnd _ (htouch j).2
    rw [dot_e] at h1
    rw [show -e j = (-1 : ℝ) • e j by simp, dot_smul, dot_e] at h2
    rw [abs_lt]; constructor <;> linarith
  have hnk_lt : |n k| < 2 * tauM := by
    have h1 := hnd _ (htouch k).1; have h2 := hnd _ (htouch k).2
    rw [dot_e] at h1
    rw [show -e k = (-1 : ℝ) • e k by simp, dot_smul, dot_e] at h2
    rw [abs_lt]; constructor <;> linarith
  have hsnk : s * n k = 1 - δ := by rw [hδ]; linarith
  have habs_snk : |s * n k| = |n k| := by
    rw [abs_mul]; rcases hs with rfl | rfl <;> simp
  have hδlt : |δ| < 2 * tauM - 1 := by
    rw [abs_lt]; constructor
    · have := hsnk ▸ habs_snk ▸ hnk_lt
      rw [abs_lt] at this; linarith
    · have := hnj_lt; rw [abs_lt] at this; rw [hδ]; linarith
  -- v = n - δ e_j + s δ e_k
  refine ⟨n - δ • e j + (s * δ) • e k, ?_, ?_, ?_⟩
  · simp [e, Pi.single_apply, hkj, hδ]
  · have : n k = s * (1 - δ) := by rw [← hsnk, ← mul_assoc, hss, one_mul]
    simp [e, hjk, this]
    try ring
  · intro d hd
    have hdj := diffBody_coord_le hcube hd j
    have hdk := diffBody_coord_le hcube hd k
    have hexp : dot (n - δ • e j + (s * δ) • e k) d = dot n d - δ * d j + s * δ * d k := by
      simp [dot, e, Pi.single_apply, Finset.sum_add_distrib, Finset.sum_sub_distrib, sub_mul,
        add_mul, Finset.sum_ite_eq']
      try ring
    rw [hexp]
    have h1 : -(δ * d j) ≤ |δ| := by
      calc -(δ * d j) ≤ |δ * d j| := neg_le_abs _
        _ = |δ| * |d j| := abs_mul _ _
        _ ≤ |δ| * 1 := mul_le_mul_of_nonneg_left hdj (abs_nonneg _)
        _ = |δ| := mul_one _
    have h2 : s * δ * d k ≤ |δ| := by
      calc s * δ * d k ≤ |s * δ * d k| := le_abs_self _
        _ = |s| * |δ| * |d k| := by rw [abs_mul, abs_mul]
        _ ≤ 1 * |δ| * 1 := by
          rcases hs with rfl | rfl <;> simp <;> exact mul_le_of_le_one_right (abs_nonneg _) hdk
        _ = |δ| := by ring
    have := hnd d hd
    have htau2 : 2 * tauM + 2 * (2 * tauM - 1) = Real.sqrt 2 := by unfold tauM; ring
    linarith

/-- Volume contradiction: two short directions `v₊`, `v₋` contradict minimality. -/
lemma no_two_short {K : Set E3} (hK : IsConvexBody K) (hmin : MinimalBox K)
    (htouch : ∀ i, e i ∈ diffBody K ∧ -e i ∈ diffBody K)
    {j k m : Fin 3} (hjk : j ≠ k) (hjm : j ≠ m) (hkm : k ≠ m)
    {vp vm : E3} (hvp : vp j = 1 ∧ vp k = 1) (hvm : vm j = 1 ∧ vm k = -1)
    (hp : ∀ d ∈ diffBody K, dot vp d < Real.sqrt 2) (hm : ∀ d ∈ diffBody K, dot vm d < Real.sqrt 2) :
    False := by
  obtain ⟨hconv, hcomp, hint⟩ := hK
  have hne : K.Nonempty := hint.mono interior_subset
  obtain ⟨hs1, hs2, hs3⟩ := sqrt2_facts
  -- maxima of dot v on D and minima of dot v on K
  have hDc := diffBody_compact hcomp
  have hDne : (diffBody K).Nonempty := ⟨_, (htouch j).1⟩
  obtain ⟨dp, hdp, hdpmax⟩ := hDc.exists_isMaxOn hDne
    (show ContinuousOn (fun d => dot vp d) _ by unfold dot; fun_prop)
  obtain ⟨dm, hdm, hdmmax⟩ := hDc.exists_isMaxOn hDne
    (show ContinuousOn (fun d => dot vm d) _ by unfold dot; fun_prop)
  obtain ⟨xp, hxp, hxpmin⟩ := hcomp.exists_isMinOn hne
    (show ContinuousOn (fun d => dot vp d) _ by unfold dot; fun_prop)
  obtain ⟨xm, hxm, hxmmin⟩ := hcomp.exists_isMinOn hne
    (show ContinuousOn (fun d => dot vm d) _ by unfold dot; fun_prop)
  set Wp := dot vp dp; set Wm := dot vm dm
  have hWp1 : 1 ≤ Wp := by
    have := hdpmax (htouch j).1; simp only [mem_ofPred_eq] at this
    rw [dot_e, hvp.1] at this; exact this
  have hWm1 : 1 ≤ Wm := by
    have := hdmmax (htouch j).1; simp only [mem_ofPred_eq] at this
    rw [dot_e, hvm.1] at this; exact this
  have hWp2 : Wp < Real.sqrt 2 := hp dp hdp
  have hWm2 : Wm < Real.sqrt 2 := hm dm hdm
  -- the affine map
  let A : Matrix (Fin 3) (Fin 3) ℝ := fun r c =>
    if r = m then (if c = m then 1 else 0) else if r = j then vp c / Wp else vm c / Wm
  let b : E3 := fun r => if r = m then 0 else if r = j then - dot vp xp / Wp else - dot vm xm / Wm
  have hmap : ∀ x ∈ K, Matrix.toLin' A x + b ∈ cube := by
    intro x hx
    have hxc := hmin.1 hx
    simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at hxc ⊢
    have hrow : ∀ r, (Matrix.toLin' A x + b) r =
        if r = m then x m else if r = j then (dot vp x - dot vp xp) / Wp
        else (dot vm x - dot vm xm) / Wm := by
      intro r
      show (Matrix.mulVec A x) r + b r = _
      simp only [Matrix.mulVec, dotProduct, A, b]
      split_ifs with h1 h2
      · simp [Finset.sum_ite_eq']
      · simp only [dot]; rw [sub_div, Finset.sum_div]; congr 1
        · apply Finset.sum_congr rfl; intro c _; ring
        · ring
      · simp only [dot]; rw [sub_div, Finset.sum_div]; congr 1
        · apply Finset.sum_congr rfl; intro c _; ring
        · ring
    have hdp' : dot vp (x - xp) ≤ Wp := hdpmax ⟨xp, hxp, x, hx, rfl⟩
    have hdm' : dot vm (x - xm) ≤ Wm := hdmmax ⟨xm, hxm, x, hx, rfl⟩
    have hxp' : dot vp xp ≤ dot vp x := hxpmin hx
    have hxm' : dot vm xm ≤ dot vm x := hxmmin hx
    rw [dot_sub] at hdp' hdm'
    have hWp0 : 0 < Wp := by linarith
    have hWm0 : 0 < Wm := by linarith
    constructor <;> intro r <;> rw [hrow r] <;> split_ifs
    · exact hxc.1 m
    · exact div_nonneg (by linarith) hWp0.le
    · exact div_nonneg (by linarith) hWm0.le
    · exact hxc.2 m
    · rw [div_le_one hWp0]; linarith
    · rw [div_le_one hWm0]; linarith
  have hdet := hmin.abs_det_le _ _ hmap
  rw [LinearMap.det_toLin'] at hdet
  -- compute det A = ±2 / (Wp Wm)
  have hWp0 : 0 < Wp := by linarith
  have hWm0 : 0 < Wm := by linarith
  have hdetA : |A.det| = 2 / (Wp * Wm) := by
    rw [Matrix.det_fin_three]
    rw [abs_eq (by positivity)]
    fin_cases j <;> fin_cases k <;> fin_cases m <;> simp_all [A] <;>
      first | (left; field_simp; ring_nf; done) | (right; field_simp; ring_nf; done)
  rw [hdetA, div_le_one (by positivity)] at hdet
  have : Wp * Wm < 2 := by
    calc Wp * Wm < Real.sqrt 2 * Real.sqrt 2 :=
          mul_lt_mul'' hWp2 hWm2 (by linarith) (by linarith)
      _ = 2 := hs2
  linarith

/-- **Lemma M (long-chord lemma).** -/
theorem lemmaM {K : Set E3} (hK : IsConvexBody K) (hmin : MinimalBox K)
    (htouch : ∀ i, ∃ x ∈ K, x + e i ∈ K) (j k : Fin 3) (hjk : j ≠ k) :
    ∃ s : ℝ, (s = 1 ∨ s = -1) ∧ ∃ x ∈ K, x + tauM • (e j + s • e k) ∈ K := by
  have ht : ∀ i, e i ∈ diffBody K ∧ -e i ∈ diffBody K := by
    intro i; obtain ⟨x, hx, hx'⟩ := htouch i
    exact ⟨⟨x, hx, x + e i, hx', by abel⟩, ⟨x + e i, hx', x, hx, by abel⟩⟩
  have key : ∀ s : ℝ, (s = 1 ∨ s = -1) → tauM • (e j + s • e k) ∈ diffBody K →
      ∃ x ∈ K, x + tauM • (e j + s • e k) ∈ K := by
    rintro s - ⟨x, hx, y, hy, h⟩
    exact ⟨x, hx, by rw [h]; simpa using hy⟩
  -- the third axis
  obtain ⟨m, hjm, hkm⟩ : ∃ m : Fin 3, j ≠ m ∧ k ≠ m := by
    fin_cases j <;> fin_cases k <;> simp_all <;> first
      | exact ⟨0, by decide, by decide⟩ | exact ⟨1, by decide, by decide⟩
      | exact ⟨2, by decide, by decide⟩
  by_cases hp : tauM • (e j + (1:ℝ) • e k) ∈ diffBody K
  · exact ⟨1, Or.inl rfl, key 1 (Or.inl rfl) hp⟩
  by_cases hm : tauM • (e j + (-1:ℝ) • e k) ∈ diffBody K
  · exact ⟨-1, Or.inr rfl, key (-1) (Or.inr rfl) hm⟩
  exfalso
  obtain ⟨vp, hvp1, hvp2, hvp⟩ := short_direction hK.1 hK.2.1 hmin.1 ht hjk (Or.inl rfl) hp
  obtain ⟨vm, hvm1, hvm2, hvm⟩ := short_direction hK.1 hK.2.1 hmin.1 ht hjk (Or.inr rfl) hm
  exact no_two_short hK hmin ht hjk hjm hkm ⟨hvp1, hvp2⟩ ⟨hvm1, hvm2⟩ hvp hvm

end Hadwiger

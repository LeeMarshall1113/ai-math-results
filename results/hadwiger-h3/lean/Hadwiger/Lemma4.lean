import Hadwiger.MinBox
import Hadwiger.Leaf

/-!
# Papadoperakis' Lemma 4 (interior form): `C(K, int K) ≤ C(E ∪ V_p, int K)`

For a convex body `K ⊆ [0,1]³` containing a configuration `V_p` (one point on each face, the
points on opposite faces differing by `e_i`; in fact only "one point of `K` on each face" is
used) any cover of `E ∪ V_p` by `N` translates of `int K` yields a cover of `K` by `N`
translates of `int K`.  Minimality of the box is not needed.

Proof (written for this formalisation):
1. *Cone covering.*  Every `x ∈ [0,1]³` is a convex combination `θ a + (1-θ) o` with
   `a ∈ E ∪ V_p`, `θ > 0`, `o ∈ O_p`: follow the ray from `o₁ = q_{10}` through `x` to the
   boundary point `b` on a face `F`; then follow the ray in `F` from the touching point
   `q_F` through `b` to a point `a` on an edge of `F`.
2. *Illumination.*  If `a - tᵢ ∈ int K` then `x - θ tᵢ = θ (a - tᵢ) + (1-θ) o ∈ int K`, hence
   `x - ε tᵢ ∈ int K` for all `0 < ε ≤ θ`.
3. *Compactness* gives one `ε` for all `x ∈ K`, so `K ⊆ ⋃ᵢ (ε tᵢ + int K)`.
-/

namespace Hadwiger

open Set

/-- Exit point of a ray in the cube: starting at `o ∈ cube` in direction `d ≠ 0` with
`o + d ∈ cube`, there is `s ≥ 1` with `o + s d ∈ cube` lying on a face `x_l ∈ {0,1}` for some
coordinate `l` with `d_l ≠ 0`. -/
lemma ray_exit {o d : E3} (ho : o ∈ cube) (hod : o + d ∈ cube) (hd : d ≠ 0) :
    ∃ s : ℝ, 1 ≤ s ∧ o + s • d ∈ cube ∧ ∃ l, d l ≠ 0 ∧ OnBdry (o + s • d) l := by
  have hS : (Finset.univ.filter fun l : Fin 3 => d l ≠ 0).Nonempty := by
    by_contra h
    rw [Finset.not_nonempty_iff_eq_empty, Finset.filter_eq_empty_iff] at h
    exact hd (funext fun l => by simpa using h (Finset.mem_univ l))
  let sl : Fin 3 → ℝ := fun l => if 0 < d l then (1 - o l) / d l else (0 - o l) / d l
  obtain ⟨l0, hl0, hmin⟩ := Finset.exists_min_image _ sl hS
  have hl0d : d l0 ≠ 0 := (Finset.mem_filter.1 hl0).2
  have hoc := mem_cube_coord ho
  have hodc := mem_cube_coord hod
  -- every s_l ≥ 1
  have hsl1 : ∀ l, d l ≠ 0 → 1 ≤ sl l := by
    intro l hl
    have h1 := hoc l; have h2 := hodc l
    simp only [Pi.add_apply] at h2
    simp only [sl]
    split_ifs with hpos
    · rw [le_div_iff₀ hpos]; linarith
    · have hneg : d l < 0 := lt_of_le_of_ne (not_lt.1 hpos) hl
      rw [le_div_iff_of_neg hneg]; linarith
  refine ⟨sl l0, hsl1 l0 hl0d, ?_, l0, hl0d, ?_⟩
  · simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply, Pi.add_apply,
      Pi.smul_apply, smul_eq_mul]
    have hs0 : 0 ≤ sl l0 := le_trans zero_le_one (hsl1 l0 hl0d)
    constructor <;> intro l <;> have h1 := hoc l
    · by_cases hl : d l = 0
      · simp [hl, h1.1]
      · have hle := hmin l (Finset.mem_filter.2 ⟨Finset.mem_univ _, hl⟩)
        rcases lt_or_gt_of_ne hl with hneg | hpos
        · have : sl l = (0 - o l) / d l := by simp [sl, not_lt.2 hneg.le]
          rw [this, le_div_iff_of_neg hneg] at hle
          nlinarith
        · nlinarith [mul_nonneg hs0 hpos.le]
    · by_cases hl : d l = 0
      · simp [hl, h1.2]
      · have hle := hmin l (Finset.mem_filter.2 ⟨Finset.mem_univ _, hl⟩)
        rcases lt_or_gt_of_ne hl with hneg | hpos
        · nlinarith [mul_nonpos_of_nonneg_of_nonpos hs0 hneg.le]
        · have : sl l = (1 - o l) / d l := by simp [sl, hpos]
          rw [this, le_div_iff₀ hpos] at hle
          nlinarith
  · unfold OnBdry
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, sl]
    split_ifs with hpos
    · right; field_simp; ring
    · left; field_simp; ring

/-- Convex combination helper: `α u + β v = (α + β) w` with `w ∈ S`. -/
lemma combo_scale {S : Set E3} (hS : Convex ℝ S) {u v : E3} (hu : u ∈ S) (hv : v ∈ S)
    {α β : ℝ} (hα : 0 ≤ α) (hβ : 0 ≤ β) : ∃ w ∈ S, α • u + β • v = (α + β) • w := by
  rcases (add_nonneg hα hβ).lt_or_eq with h | h
  · refine ⟨(α / (α + β)) • u + (β / (α + β)) • v, hS hu hv (by positivity) (by positivity)
      (by field_simp), ?_⟩
    rw [smul_add, smul_smul, smul_smul]; congr 1 <;> congr 1 <;> field_simp
  · refine ⟨u, hu, ?_⟩
    have hα0 : α = 0 := by linarith
    have hβ0 : β = 0 := by linarith
    simp [hα0, hβ0]

/-- **Cone covering.**  Every point of the cube is `θ a + (1-θ) o` with `a ∈ E ∪ V_p`,
`0 < θ ≤ 1` and `o ∈ O_p`. -/
theorem cube_cone_cover (c : Cfg) (hc : c ∈ CfgSpace) {x : E3} (hx : x ∈ cube) :
    ∃ a ∈ skeleton ∪ V c, ∃ θ : ℝ, 0 < θ ∧ θ ≤ 1 ∧ ∃ o ∈ O c, x = θ • a + (1 - θ) • o := by
  have hptc : ∀ ib, pt c ib ∈ cube := by
    intro ⟨i, b⟩
    simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply, pt]
    constructor <;> intro j <;> split_ifs with h
    · cases b <;> simp [bv]
    · exact (hc i j (fun h' => h h'.symm)).1
    · cases b <;> simp [bv]
    · exact (hc i j (fun h' => h h'.symm)).2
  set o := pt c (0, false)
  have hoO : o ∈ O c := V_sub_O c ⟨_, rfl⟩
  by_cases hxo : x = o
  · exact ⟨o, Or.inr ⟨_, rfl⟩, 1, one_pos, le_rfl, o, hoO, by simp [hxo]⟩
  -- first ray: from o through x to a face point b
  obtain ⟨s, hs1, hbc, l, -, hbl⟩ := ray_exit (hptc _) (by convert hx using 1; abel) (sub_ne_zero.2 hxo)
  set b := o + s • (x - o)
  have hs0 : 0 < s := lt_of_lt_of_le one_pos hs1
  have hxb : x = (1/s) • b + (1 - 1/s) • o := by
    simp only [b, smul_add, smul_smul, sub_smul, one_smul]
    rw [show 1 / s * s = 1 by field_simp]; module
  have hinv0 : 0 < 1/s := by positivity
  have hinv1 : 1/s ≤ 1 := by rw [div_le_one hs0]; exact hs1
  obtain ⟨β, hβ⟩ := bv_of_onBdry hbl
  set q := pt c (l, β)
  have hqO : q ∈ O c := V_sub_O c ⟨_, rfl⟩
  have hql : q l = b l := by simp [q, pt, hβ]
  by_cases hbq : b = q
  · exact ⟨q, Or.inr ⟨_, rfl⟩, 1/s, hinv0, hinv1, o, hoO, by rw [hxb, hbq]⟩
  -- second ray: in the face, from q through b to an edge point a
  obtain ⟨r, hr1, hac, l', hl'd, hal'⟩ := ray_exit (hptc _) (by convert hbc using 1; abel)
    (sub_ne_zero.2 hbq)
  set a := q + r • (b - q)
  have hr0 : 0 < r := lt_of_lt_of_le one_pos hr1
  have hll' : l ≠ l' := by
    intro h; subst h; apply hl'd; simp [hql]
  have hal : OnBdry a l := by
    unfold OnBdry; simp only [a, Pi.add_apply, Pi.smul_apply, Pi.sub_apply, hql, sub_self,
      smul_eq_mul, mul_zero, add_zero]
    exact hbl
  have haE : a ∈ skeleton := ⟨hac, l, l', hll', hal, hal'⟩
  have hba : b = (1/r) • a + (1 - 1/r) • q := by
    simp only [a, smul_add, smul_smul, sub_smul, one_smul]
    rw [show 1 / r * r = 1 by field_simp]; module
  have hrinv0 : 0 < 1/r := by positivity
  have hrinv1 : 1/r ≤ 1 := by rw [div_le_one hr0]; exact hr1
  obtain ⟨w, hw, hwq⟩ := combo_scale (O_convex c) hqO hoO
    (show 0 ≤ 1/s * (1 - 1/r) by nlinarith) (show 0 ≤ 1 - 1/s by linarith)
  refine ⟨a, Or.inl haE, 1/s * (1/r), by positivity, by nlinarith, w, hw, ?_⟩
  rw [hxb, hba, smul_add, smul_smul, smul_smul, add_assoc, ← smul_smul (1/s) (1 - 1/r),
    smul_smul, hwq]
  congr 1; ring_nf

/-- Interior points are stable along the segment back to a point of `K`. -/
lemma interior_shift_mono {K : Set E3} (hK : Convex ℝ K) {x w : E3} (hx : x ∈ K) {θ ε : ℝ}
    (hθ : x - θ • w ∈ interior K) (hε0 : 0 < ε) (hεθ : ε ≤ θ) : x - ε • w ∈ interior K := by
  have hθ0 : 0 < θ := lt_of_lt_of_le hε0 hεθ
  have := hK.combo_interior_self_subset_interior (show 0 < ε / θ by positivity)
    (show 0 ≤ 1 - ε / θ by rw [sub_nonneg, div_le_one hθ0]; exact hεθ) (by ring)
    (Set.add_mem_add (Set.smul_mem_smul_set hθ) (Set.smul_mem_smul_set hx))
  convert this using 1
  simp only [smul_sub, smul_smul]
  rw [show ε / θ * θ = ε by field_simp]; module

/-- **Papadoperakis, Lemma 4** (interior form, own configuration). -/
theorem lemma4 {K : Set E3} (hK : IsConvexBody K) (hcube : K ⊆ cube) {c : Cfg}
    (hc : c ∈ CfgSpace) (hV : V c ⊆ K) {N : ℕ} {t : Fin N → E3}
    (hcov : skeleton ∪ V c ⊆ ⋃ i, tr (t i) (interior K)) : CoveredByInteriorTranslates N K := by
  have hO : O c ⊆ K := convexHull_min hV hK.1
  -- pointwise illumination
  have hpt : ∀ x ∈ K, ∃ i, ∃ θ : ℝ, 0 < θ ∧ x - θ • t i ∈ interior K := by
    intro x hx
    obtain ⟨a, ha, θ, hθ0, hθ1, o, ho, rfl⟩ := cube_cone_cover c hc (hcube hx)
    obtain ⟨i, hi⟩ := mem_iUnion.1 (hcov ha)
    refine ⟨i, θ, hθ0, ?_⟩
    have := hK.1.combo_interior_self_subset_interior hθ0 (show 0 ≤ 1 - θ by linarith) (by ring)
      (Set.add_mem_add (Set.smul_mem_smul_set (mem_tr.1 hi)) (Set.smul_mem_smul_set (hO ho)))
    convert this using 1
    simp only [smul_sub]; abel
  -- open sets U (i, n)
  let U : Fin N × ℕ → Set E3 := fun p => {y | y - (1 / ((p.2 : ℝ) + 1)) • t p.1 ∈ interior K}
  have hUo : ∀ p, IsOpen (U p) := fun p =>
    isOpen_interior.preimage (by fun_prop)
  have hKU : K ⊆ ⋃ p, U p := by
    intro x hx
    obtain ⟨i, θ, hθ, hxθ⟩ := hpt x hx
    obtain ⟨n, hn⟩ := exists_nat_one_div_lt hθ
    refine mem_iUnion.2 ⟨(i, n), interior_shift_mono hK.1 hx hxθ (by positivity) hn.le⟩
  obtain ⟨s, hs⟩ := hK.2.1.elim_finite_subcover U hUo hKU
  set n0 := s.sup (fun p => p.2)
  refine ⟨fun i => (1 / ((n0 : ℝ) + 1)) • t i, fun x hx => ?_⟩
  obtain ⟨p, hp, hxp⟩ := mem_iUnion₂.1 (hs hx)
  refine mem_iUnion.2 ⟨p.1, ?_⟩
  rw [mem_tr]
  apply interior_shift_mono hK.1 hx hxp (by positivity)
  have : (p.2 : ℝ) ≤ n0 := by exact_mod_cast Finset.le_sup (f := fun p : Fin N × ℕ => p.2) hp
  gcongr

end Hadwiger

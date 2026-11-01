import Hadwiger.Basic

/-!
# Symmetries of the cube and the reduction to `D*`

A symmetry `g = (σ, f)` (`σ` a permutation of the axes, `f j` = "reflect axis `j`") acts on
points by `(g x)_j = flip (f j) (x_{σ⁻¹ j})`, where `flip true t = 1 - t`.  It acts on
configurations by `(g c) k j = flip (f j) (c (σ⁻¹ k) (σ⁻¹ j))`, and `g '' V c = V (g c)`.

Main results:
* `exists_sym_mem_DStar`: every configuration in `[0,1]⁶` has an image in `D*`.
* `gAff g`: the action on points as an affine equivalence; it maps the cube and its
  1-skeleton onto themselves.
-/

namespace Hadwiger

open Set

/-- `flip true t = 1 - t`, `flip false t = t`. -/
def flip (b : Bool) (t : ℝ) : ℝ := if b then 1 - t else t

@[simp] lemma flip_false (t : ℝ) : flip false t = t := rfl
@[simp] lemma flip_true (t : ℝ) : flip true t = 1 - t := rfl

lemma flip_flip (b : Bool) (t : ℝ) : flip b (flip b t) = t := by cases b <;> simp

lemma flip_bv (b f : Bool) : flip f (bv b) = bv (xor b f) := by
  cases b <;> cases f <;> simp [bv]

lemma flip_mem_Icc {b : Bool} {t : ℝ} (h : 0 ≤ t ∧ t ≤ 1) : 0 ≤ flip b t ∧ flip b t ≤ 1 := by
  cases b <;> simp <;> constructor <;> linarith [h.1, h.2]

/-- A symmetry of the cube. -/
structure Sym where
  σ : Equiv.Perm (Fin 3)
  f : Fin 3 → Bool

/-- Action on points. -/
def gmap (g : Sym) (x : E3) : E3 := fun j => flip (g.f j) (x (g.σ.symm j))

/-- Action on configurations. -/
def gcfg (g : Sym) (c : Cfg) : Cfg := fun k j => flip (g.f j) (c (g.σ.symm k) (g.σ.symm j))

lemma gmap_pt (g : Sym) (c : Cfg) (i : Fin 3) (b : Bool) :
    gmap g (pt c (i, b)) = pt (gcfg g c) (g.σ i, xor b (g.f (g.σ i))) := by
  funext j
  simp only [gmap, pt, gcfg]
  by_cases hj : j = g.σ i
  · subst hj; simp [flip_bv]
  · have : g.σ.symm j ≠ i := by
      intro h; apply hj; rw [← h]; simp
    simp [hj, this]

lemma gmap_V (g : Sym) (c : Cfg) : gmap g '' V c = V (gcfg g c) := by
  ext y; constructor
  · rintro ⟨x, ⟨⟨i, b⟩, rfl⟩, rfl⟩
    exact ⟨_, (gmap_pt g c i b).symm⟩
  · rintro ⟨⟨k, b⟩, rfl⟩
    refine ⟨pt c (g.σ.symm k, xor b (g.f k)), ⟨_, rfl⟩, ?_⟩
    rw [gmap_pt]; simp

/-- The inverse action on points. -/
def gmapInv (g : Sym) (y : E3) : E3 := fun i => flip (g.f (g.σ i)) (y (g.σ i))

lemma gmapInv_gmap (g : Sym) (x : E3) : gmapInv g (gmap g x) = x := by
  funext i; simp [gmapInv, gmap, flip_flip]

lemma gmap_gmapInv (g : Sym) (y : E3) : gmap g (gmapInv g y) = y := by
  funext j; simp [gmapInv, gmap, flip_flip]

/-- Sign of a coordinate reflection. -/
def sgn (b : Bool) : ℝ := if b then -1 else 1

/-- Linear part of the action. -/
def gLin (g : Sym) : E3 ≃ₗ[ℝ] E3 where
  toFun x := fun j => sgn (g.f j) * x (g.σ.symm j)
  invFun y := fun i => sgn (g.f (g.σ i)) * y (g.σ i)
  map_add' x y := by funext j; simp [mul_add]
  map_smul' a x := by funext j; simp; ring
  left_inv x := by
    funext i; simp only
    cases g.f (g.σ i) <;> simp [sgn]
  right_inv y := by
    funext j; simp only [Equiv.apply_symm_apply]
    cases g.f j <;> simp [sgn]

/-- Translation part of the action. -/
def gOff (g : Sym) : E3 := fun j => bv (g.f j)

lemma gmap_eq (g : Sym) (x : E3) : gmap g x = gLin g x + gOff g := by
  funext j; simp only [gmap, gLin, gOff, LinearEquiv.coe_mk, LinearMap.coe_mk, AddHom.coe_mk,
    Pi.add_apply]
  cases g.f j <;> simp [sgn, bv] <;> ring

/-- The action on points as an affine equivalence. -/
def gAff (g : Sym) : E3 ≃ᵃ[ℝ] E3 where
  toFun := gmap g
  invFun := gmapInv g
  left_inv := gmapInv_gmap g
  right_inv := gmap_gmapInv g
  linear := gLin g
  map_vadd' p v := by
    show gmap g (v + p) = gLin g v + gmap g p
    rw [gmap_eq, gmap_eq, map_add]; abel

@[simp] lemma gAff_apply (g : Sym) (x : E3) : gAff g x = gmap g x := rfl

lemma gmap_mem_cube (g : Sym) {x : E3} (hx : x ∈ cube) : gmap g x ∈ cube := by
  simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at hx ⊢
  refine ⟨fun j => (flip_mem_Icc ⟨hx.1 _, hx.2 _⟩).1, fun j => (flip_mem_Icc ⟨hx.1 _, hx.2 _⟩).2⟩

lemma flip_onBdry {b : Bool} {t : ℝ} (h : t = 0 ∨ t = 1) : flip b t = 0 ∨ flip b t = 1 := by
  cases b <;> rcases h with h | h <;> simp [h]

lemma gmap_mem_skeleton (g : Sym) {x : E3} (hx : x ∈ skeleton) : gmap g x ∈ skeleton := by
  obtain ⟨hc, i, j, hij, hi, hj⟩ := hx
  refine ⟨gmap_mem_cube g hc, g.σ i, g.σ j, fun h => hij (g.σ.injective h), ?_, ?_⟩
  · unfold OnBdry gmap; simpa using flip_onBdry hi
  · unfold OnBdry gmap; simpa using flip_onBdry hj

lemma gmapInv_eq (g : Sym) (y : E3) : gmapInv g y = gmap ⟨g.σ.symm, fun i => g.f (g.σ i)⟩ y := by
  funext i; simp [gmapInv, gmap]

lemma gmap_image_cube (g : Sym) : gmap g '' cube = cube := by
  apply Subset.antisymm
  · rintro _ ⟨x, hx, rfl⟩; exact gmap_mem_cube g hx
  · intro y hy
    refine ⟨gmapInv g y, ?_, gmap_gmapInv g y⟩
    rw [gmapInv_eq]; exact gmap_mem_cube _ hy

lemma gmap_image_skeleton (g : Sym) : gmap g '' skeleton = skeleton := by
  apply Subset.antisymm
  · rintro _ ⟨x, hx, rfl⟩; exact gmap_mem_skeleton g hx
  · intro y hy
    refine ⟨gmapInv g y, ?_, gmap_gmapInv g y⟩
    rw [gmapInv_eq]; exact gmap_mem_skeleton _ hy

lemma gcfg_mem_CfgSpace (g : Sym) {c : Cfg} (hc : c ∈ CfgSpace) : gcfg g c ∈ CfgSpace := by
  intro k j hkj
  apply flip_mem_Icc
  exact hc _ _ (fun h => hkj (g.σ.symm.injective h))

/-! ### Reduction to `D*` -/

/-- Column sum `S_j = ∑_{i ≠ j} c i j`. -/
def colSum (c : Cfg) (j : Fin 3) : ℝ := ∑ i, if i = j then 0 else c i j

/-- `refl true s = 2 - s`, `refl false s = s` (effect of a reflection on a column sum). -/
def refl (b : Bool) (s : ℝ) : ℝ := if b then 2 - s else s

lemma colSum_gcfg (g : Sym) (c : Cfg) (j : Fin 3) :
    colSum (gcfg g c) j = refl (g.f j) (colSum c (g.σ.symm j)) := by
  have key : ∀ (b : Bool) (u v : ℝ), flip b u + flip b v = refl b (u + v) := by
    intro b u v; cases b <;> simp [refl] <;> ring
  have h1 : colSum c (g.σ.symm j) =
      ∑ k, if k = j then 0 else c (g.σ.symm k) (g.σ.symm j) := by
    unfold colSum
    rw [← Equiv.sum_comp g.σ.symm]
    simp
  rw [h1]
  unfold colSum gcfg
  fin_cases j <;> simp [Fin.sum_univ_three, key] <;> ring_nf

lemma abs_refl_sub (b : Bool) (t : ℝ) : |refl b t - 1| = |t - 1| := by
  cases b
  · simp [refl]
  · simp only [refl, ite_true]; rw [show 2 - t - 1 = -(t - 1) by ring, abs_neg]

/-- Every triple of reals can be sorted by a permutation. -/
lemma exists_perm_sorted (a : Fin 3 → ℝ) :
    ∃ τ : Equiv.Perm (Fin 3), a (τ 2) ≤ a (τ 1) ∧ a (τ 1) ≤ a (τ 0) := by
  rcases le_total (a 0) (a 1) with h01 | h01 <;> rcases le_total (a 1) (a 2) with h12 | h12 <;>
    rcases le_total (a 0) (a 2) with h02 | h02
  · exact ⟨Equiv.swap 0 2, by simp [Equiv.swap_apply_of_ne_of_ne]; exact ⟨h01, h12⟩⟩
  · exact ⟨Equiv.swap 0 2, by simp [Equiv.swap_apply_of_ne_of_ne]; exact ⟨h01, h12⟩⟩
  · -- a0 ≤ a1, a2 ≤ a1, a0 ≤ a2 : order 1,2,0
    refine ⟨Equiv.swap 0 1 * Equiv.swap 1 2, ?_⟩
    simp [Equiv.Perm.mul_apply, Equiv.swap_apply_of_ne_of_ne]; exact ⟨h02, h12⟩
  · -- a0 ≤ a1, a2 ≤ a1, a2 ≤ a0 : order 1,0,2
    refine ⟨Equiv.swap 0 1, ?_⟩
    simp [Equiv.swap_apply_of_ne_of_ne]; exact ⟨h02, h01⟩
  · -- a1 ≤ a0, a1 ≤ a2, a0 ≤ a2 : order 2,0,1
    refine ⟨Equiv.swap 0 2 * Equiv.swap 1 2, ?_⟩
    simp [Equiv.Perm.mul_apply, Equiv.swap_apply_of_ne_of_ne]; exact ⟨h01, h02⟩
  · -- a1 ≤ a0, a1 ≤ a2, a2 ≤ a0 : order 0,2,1
    refine ⟨Equiv.swap 1 2, ?_⟩
    simp [Equiv.swap_apply_of_ne_of_ne]; exact ⟨h12, h02⟩
  · exact ⟨1, by simp; exact ⟨h12, h01⟩⟩
  · exact ⟨1, by simp; exact ⟨h12, h01⟩⟩

lemma flip_decide_le (t : ℝ) : flip (decide (1/2 < t)) t ≤ 1/2 := by
  by_cases h : 1/2 < t
  · rw [decide_eq_true h]; simp; linarith
  · rw [decide_eq_false h]; simp; linarith

/-- **Symmetry reduction.**  Every configuration in `[0,1]⁶` is mapped into `D*` by some
symmetry of the cube. -/
theorem exists_sym_mem_DStar (c : Cfg) (hc : c ∈ CfgSpace) :
    ∃ g : Sym, toP (gcfg g c) ∈ DStar := by
  obtain ⟨τ, h21, h10⟩ := exists_perm_sorted (fun j => |colSum c j - 1|)
  let f : Fin 3 → Bool := fun j =>
    if j = 0 then decide (1/2 < c (τ 1) (τ 0))
    else if j = 1 then decide (1/2 < c (τ 0) (τ 1))
    else decide (1/2 < c (τ 0) (τ 2))
  refine ⟨⟨τ.symm, f⟩, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro k
    have := gcfg_mem_CfgSpace ⟨τ.symm, f⟩ hc (kPair k) (kCoord k) (kPair_ne_kCoord k)
    exact this
  · show flip (f 1) (c (τ 0) (τ 1)) ≤ 1/2
    exact flip_decide_le _
  · show flip (f 2) (c (τ 0) (τ 2)) ≤ 1/2
    exact flip_decide_le _
  · show flip (f 0) (c (τ 1) (τ 0)) ≤ 1/2
    exact flip_decide_le _
  · have e2 : toP (gcfg ⟨τ.symm, f⟩ c) 1 + toP (gcfg ⟨τ.symm, f⟩ c) 3 =
        colSum (gcfg ⟨τ.symm, f⟩ c) 2 := by
      simp [toP, colSum, Fin.sum_univ_three, kPair, kCoord]
    have e1 : toP (gcfg ⟨τ.symm, f⟩ c) 0 + toP (gcfg ⟨τ.symm, f⟩ c) 5 =
        colSum (gcfg ⟨τ.symm, f⟩ c) 1 := by
      simp [toP, colSum, Fin.sum_univ_three, kPair, kCoord]
    rw [e2, e1, colSum_gcfg, colSum_gcfg, abs_refl_sub, abs_refl_sub]
    simpa using h21
  · have e1 : toP (gcfg ⟨τ.symm, f⟩ c) 0 + toP (gcfg ⟨τ.symm, f⟩ c) 5 =
        colSum (gcfg ⟨τ.symm, f⟩ c) 1 := by
      simp [toP, colSum, Fin.sum_univ_three, kPair, kCoord]
    have e0 : toP (gcfg ⟨τ.symm, f⟩ c) 2 + toP (gcfg ⟨τ.symm, f⟩ c) 4 =
        colSum (gcfg ⟨τ.symm, f⟩ c) 0 := by
      simp [toP, colSum, Fin.sum_univ_three, kPair, kCoord]
    rw [e1, e0, colSum_gcfg, colSum_gcfg, abs_refl_sub, abs_refl_sub]
    simpa using h10

end Hadwiger

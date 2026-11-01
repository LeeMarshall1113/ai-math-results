import Hadwiger.Leaf
import Hadwiger.Touching
import Hadwiger.MinBox
import Hadwiger.Lemma4
import Hadwiger.Symmetry

/-!
# From the tiling of the parameter space to `H₃ ≤ 13`

* `QBox`: rational boxes `∏ [a_k, a_k + d_k]` in Prymak's coordinates (the certificate format).
* `Certified B`: `B` is tiled by a bisection tree (halving one coordinate at a time) whose
  leaves all pass `LeafOK`.  `certcheck.py` checks the longest-side tree, a special case.
* `oldStartIdx`: the 39,380 start boxes of side `1/10` meeting Prymak's `U^D_20`.

No axioms here: `H3_le_13_of_symcover` takes the certified box family as a parameter.  The
certificate axiom and the unconditional theorem are in `Final.lean`.
-/

namespace Hadwiger

open Set

/-! ### Rational boxes and the bisection tree -/

/-- A rational box `∏ [a_k, a_k + d_k]` (Prymak's coordinates, 0-based). -/
structure QBox where
  a : Fin 6 → ℚ
  d : Fin 6 → ℚ

/-- The real box of a rational box. -/
def QBox.toR (B : QBox) : RBox := ⟨fun k => (B.a k : ℝ), fun k => (B.a k : ℝ) + (B.d k : ℝ)⟩

/-- Membership of Prymak's `p`. -/
def QBox.MemP (B : QBox) (p : Fin 6 → ℝ) : Prop := ∀ k, (B.a k : ℝ) ≤ p k ∧ p k ≤ B.a k + B.d k

/-- Lower half of `B` in coordinate `k`. -/
def QBox.lower (B : QBox) (k : Fin 6) : QBox := ⟨B.a, Function.update B.d k (B.d k / 2)⟩
/-- Upper half of `B` in coordinate `k`. -/
def QBox.upper (B : QBox) (k : Fin 6) : QBox :=
  ⟨Function.update B.a k (B.a k + B.d k / 2), Function.update B.d k (B.d k / 2)⟩

/-- `B` is tiled by a bisection tree all of whose leaves pass the leaf check. -/
inductive Certified : QBox → Prop
  | leaf {B : QBox} : LeafOK B.toR → Certified B
  | split {B : QBox} (k : Fin 6) : Certified (B.lower k) → Certified (B.upper k) → Certified B

/-- Every point of a certified box lies in a passing leaf. -/
theorem Certified.exists_leaf {B : QBox} (h : Certified B) {p : Fin 6 → ℝ} (hp : B.MemP p) :
    ∃ L : QBox, LeafOK L.toR ∧ L.MemP p := by
  induction h with
  | leaf hl => exact ⟨_, hl, hp⟩
  | @split B k _ _ ih1 ih2 =>
    by_cases hk : p k ≤ (B.a k : ℝ) + (B.d k : ℝ) / 2
    · apply ih1
      intro k'
      by_cases h : k' = k
      · subst h; simp only [QBox.lower, Function.update_self]; push_cast
        exact ⟨(hp k').1, hk⟩
      · simp only [QBox.lower, Function.update_of_ne h]; exact hp k'
    · apply ih2
      intro k'
      by_cases h : k' = k
      · subst h; simp only [QBox.upper, Function.update_self]; push_cast
        constructor
        · linarith
        · linarith [(hp k').2]
      · simp only [QBox.upper, Function.update_of_ne h]; exact hp k'

lemma QBox.memP_toR {B : QBox} {c : Cfg} (h : B.MemP (toP c)) : B.toR.Mem c := h

/-! ### The certified start boxes -/

/-- Prymak's index condition for `U^D_{2M}` (`certcheck.in_UD`), 0-based indices. -/
def inUD (k : Fin 6 → ℕ) : Prop := k 1 + k 3 ≤ k 0 + k 5 + 1 ∧ k 0 + k 5 + 1 ≤ k 2 + k 4 + 2

instance : DecidablePred inUD := fun k => by unfold inUD; infer_instance

/-- The start-box index set of the existing run: `1/10`-boxes (`k₀,k₁,k₂ < 5`,
`k₃,k₄,k₅ < 10`) containing a `1/20`-box of `U^D_20`. -/
def oldStartIdx : Set (Fin 6 → ℕ) :=
  {k | (∀ i : Fin 6, k i < (if i.val < 3 then 5 else 10)) ∧
    ∃ ε : Fin 6 → Bool, inUD (fun i => 2 * k i + (if ε i then 1 else 0))}

/-- Start box with index `k`: `∏ [k_i/10, (k_i+1)/10]`. -/
def startBox (k : Fin 6 → ℕ) : QBox := ⟨fun i => (k i : ℚ) / 10, fun _ => 1 / 10⟩

/-! ### Affine transport -/

/-- Pulling a cover back along an affine equivalence. -/
lemma cover_pullback (A : E3 ≃ᵃ[ℝ] E3) {X K : Set E3} {N : ℕ} {t : Fin N → E3}
    (h : A '' X ⊆ ⋃ i, tr (t i) (interior (A '' K))) :
    X ⊆ ⋃ i, tr (A.symm.linear (t i)) (interior K) := by
  intro x hx
  obtain ⟨i, hi⟩ := mem_iUnion.1 (h ⟨x, hx, rfl⟩)
  rw [mem_tr, ← affine_image_interior] at hi
  obtain ⟨y, hy, hyx⟩ := hi
  refine mem_iUnion.2 ⟨i, ?_⟩
  rw [mem_tr]
  have : x = y + A.symm.linear (t i) := by
    apply A.injective
    have h1 : A (y + A.symm.linear (t i)) = A y + A.linear (A.symm.linear (t i)) := by
      rw [show y + A.symm.linear (t i) = A.symm.linear (t i) +ᵥ y from add_comm _ _,
        A.map_vadd]; simp [add_comm]
    rw [h1, hyx]; simp
  rw [this]; simpa using hy

lemma chord_image (A : E3 ≃ᵃ[ℝ] E3) {K : Set E3} {w : E3} {τ : ℝ} (h : HasChordDir K w τ) :
    HasChordDir (A '' K) (A.linear w) τ := by
  obtain ⟨x, hx, hxw⟩ := h
  refine ⟨A x, ⟨x, hx, rfl⟩, ⟨x + τ • w, hxw, ?_⟩⟩
  rw [show x + τ • w = (τ • w) +ᵥ x from add_comm _ _, A.map_vadd]
  simp [add_comm]

lemma chord_neg {K : Set E3} {w : E3} {τ : ℝ} (h : HasChordDir K (-w) τ) : HasChordDir K w τ := by
  obtain ⟨x, hx, hxw⟩ := h
  exact ⟨x + τ • -w, hxw, by simpa using hx⟩

lemma gLin_e (g : Sym) (i : Fin 3) : gLin g (e i) = sgn (g.f (g.σ i)) • e (g.σ i) := by
  funext j
  simp only [gLin, LinearEquiv.coe_mk, LinearMap.coe_mk, AddHom.coe_mk, e, Pi.smul_apply,
    smul_eq_mul, Pi.single_apply]
  by_cases h : j = g.σ i
  · subst h; simp
  · have : g.σ.symm j ≠ i := fun h' => h (by rw [← h']; simp)
    simp [h, this]

lemma sgn_sq (b : Bool) : sgn b * sgn b = 1 := by cases b <;> simp [sgn]

/-- Chords survive a cube symmetry (with possibly changed sign). -/
lemma chord_sym (g : Sym) {K : Set E3} {j k : Fin 3}
    (h : ∃ s : ℝ, (s = 1 ∨ s = -1) ∧
      HasChordDir K (e (g.σ.symm j) + s • e (g.σ.symm k)) tauM) :
    ∃ s : ℝ, (s = 1 ∨ s = -1) ∧ HasChordDir (gAff g '' K) (e j + s • e k) tauM := by
  obtain ⟨s, hs, hch⟩ := h
  have h2 := chord_image (gAff g) hch
  have hlin : (gAff g).linear (e (g.σ.symm j) + s • e (g.σ.symm k)) =
      sgn (g.f j) • (e j + (s * sgn (g.f j) * sgn (g.f k)) • e k) := by
    show gLin g (e (g.σ.symm j) + s • e (g.σ.symm k)) = _
    rw [map_add, map_smul, gLin_e, gLin_e]
    simp only [Equiv.apply_symm_apply, smul_add, smul_smul]
    congr 1
    congr 1
    calc s * sgn (g.f k) = s * sgn (g.f k) * (sgn (g.f j) * sgn (g.f j)) := by rw [sgn_sq, mul_one]
      _ = _ := by ring
  rw [hlin] at h2
  refine ⟨s * sgn (g.f j) * sgn (g.f k), ?_, ?_⟩
  · rcases hs with rfl | rfl <;> cases g.f j <;> cases g.f k <;> simp [sgn]
  · cases hfj : g.f j
    · have hpos : sgn (g.f j) • (e j + (s * sgn (g.f j) * sgn (g.f k)) • e k) =
          e j + (s * sgn false * sgn (g.f k)) • e k := by rw [hfj]; simp [sgn]
      rw [hpos] at h2; exact h2
    · apply chord_neg
      have hneg : sgn (g.f j) • (e j + (s * sgn (g.f j) * sgn (g.f k)) • e k) =
          -(e j + (s * sgn true * sgn (g.f k)) • e k) := by
        rw [hfj]; cases g.f k <;> simp [sgn] <;> abel
      rw [hneg] at h2; exact h2

/-! ### Main theorem -/

/-- The certificate reduction, with the box family as a parameter: if a family of certified
boxes covers `D*`, then every convex body in `ℝ³` is covered by 13 translates of its interior.
Uses no axioms beyond Lean's standard ones: the minimal box (`exists_minimalBox'`),
Papadoperakis' Lemma 3(a) (`touching_chord`) and Lemma 4 (`lemma4`) are proved. -/
theorem H3_le_13_of_symcover (𝓑 : Set QBox) (hcert : ∀ B ∈ 𝓑, Certified B)
    (hsc : ∀ c ∈ CfgSpace, ∃ g : Sym, ∃ B ∈ 𝓑, B.MemP (toP (gcfg g c)))
    (K : Set E3) (hK : IsConvexBody K) : CoveredByInteriorTranslates 13 K := by
  -- normalisation
  obtain ⟨A, hmin⟩ := exists_minimalBox' K hK
  set K1 := A '' K with hK1
  have hK1b : IsConvexBody K1 := hK.affine_image A
  choose x hx hxe using touching_chord hK1b hmin
  let c : Cfg := fun i j => x i j
  have hcube := hmin.1
  have hc : c ∈ CfgSpace := by
    intro i j _
    have := hcube (hx i)
    simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at this
    exact ⟨this.1 j, this.2 j⟩
  have hxi0 : ∀ i, x i i = 0 := by
    intro i
    have h1 := hcube (hx i); have h2 := hcube (hxe i)
    simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at h1 h2
    have := h2.2 i
    simp [e] at this
    linarith [h1.1 i]
  have hptc : ∀ i b, pt c (i, b) = if b then x i + e i else x i := by
    intro i b; funext j
    by_cases hj : j = i
    · subst hj; cases b <;> simp [pt, bv, e, hxi0]
    · cases b <;> simp [pt, hj, e, c]
  have hV : V c ⊆ K1 := by
    rintro _ ⟨⟨i, b⟩, rfl⟩
    rw [hptc]; cases b
    · exact hx i
    · exact hxe i
  -- Lemma M on K1
  have hM := lemmaM hK1b hmin (fun i => ⟨x i, hx i, hxe i⟩)
  -- symmetry into a certified box
  obtain ⟨g, B, hB, hpB⟩ := hsc c hc
  set K2 := gAff g '' K1
  set c' := gcfg g c
  have hK2b : IsConvexBody K2 := hK1b.affine_image (gAff g)
  have hV2 : V c' ⊆ K2 := by
    rw [← gmap_V]; exact image_mono hV
  have hO2 : O c' ⊆ K2 := convexHull_min hV2 hK2b.1
  have hch2 : ∀ j k : Fin 3, j ≠ k → ∃ s : ℝ, (s = 1 ∨ s = -1) ∧
      HasChordDir K2 (e j + s • e k) tauM := by
    intro j k hjk
    apply chord_sym g
    obtain ⟨s, hs, y, hy, hys⟩ := hM (g.σ.symm j) (g.σ.symm k)
      (fun h => hjk (g.σ.symm.injective h))
    exact ⟨s, hs, y, hy, hys⟩
  -- the leaf containing c'
  obtain ⟨L, hL, hpL⟩ := (hcert B hB).exists_leaf hpB
  obtain ⟨t, ht⟩ := leaf_sound hL (QBox.memP_toR hpL) hK2b hO2 hch2
  -- pull back to K1
  have ht1 : gAff g '' (skeleton ∪ V c) ⊆ ⋃ i, tr (t i) (interior (gAff g '' K1)) := by
    rw [image_union]
    have e1 : gAff g '' skeleton = skeleton := by
      rw [show (gAff g : E3 → E3) = gmap g from rfl]; exact gmap_image_skeleton g
    have e2 : gAff g '' V c = V c' := by
      rw [show (gAff g : E3 → E3) = gmap g from rfl]; exact gmap_V g c
    rw [e1, e2]; exact ht
  have hcov1 := cover_pullback (gAff g) ht1
  obtain ⟨t1, ht1'⟩ := lemma4 hK1b hmin.1 hc hV hcov1
  -- pull back to K
  exact ⟨_, cover_pullback A ht1'⟩

/-- The `D*` form: a certified family covering `D*` suffices (no axioms). -/
theorem H3_le_13_of_cover (𝓑 : Set QBox) (hcert : ∀ B ∈ 𝓑, Certified B)
    (hcov : ∀ p ∈ DStar, ∃ B ∈ 𝓑, B.MemP p) (K : Set E3) (hK : IsConvexBody K) :
    CoveredByInteriorTranslates 13 K := by
  refine H3_le_13_of_symcover 𝓑 hcert (fun c hc => ?_) K hK
  obtain ⟨g, hg⟩ := exists_sym_mem_DStar c hc
  obtain ⟨B, hB, hpB⟩ := hcov _ hg
  exact ⟨g, B, hB, hpB⟩

/-- Composition of symmetries. -/
def Sym.comp (g2 g1 : Sym) : Sym := ⟨g1.σ.trans g2.σ, fun j => xor (g2.f j) (g1.f (g2.σ.symm j))⟩

lemma flip_flip_xor (a b : Bool) (t : ℝ) : flip a (flip b t) = flip (xor a b) t := by
  cases a <;> cases b <;> simp

lemma gcfg_comp (g2 g1 : Sym) (c : Cfg) : gcfg (g2.comp g1) c = gcfg g2 (gcfg g1 c) := by
  funext k j
  simp [gcfg, Sym.comp, flip_flip_xor, Bool.xor_comm]

/-- Covering by translates of the interior gives covering by smaller homothets. -/
theorem homothets_of_interior {K : Set E3} (hK : IsConvexBody K) {N : ℕ}
    (h : CoveredByInteriorTranslates N K) :
    ∃ r : ℝ, 0 < r ∧ r < 1 ∧ ∃ x : Fin N → E3, K ⊆ ⋃ i, (fun y => x i + r • y) '' K := by
  obtain ⟨t, ht⟩ := h
  obtain ⟨z, hz⟩ := hK.2.2
  -- open sets: y ∈ t i + z + r_n (int K - z), r_n = 1 - 1/(n+2)
  let rn : ℕ → ℝ := fun n => 1 - 1 / ((n : ℝ) + 2)
  have hrn0 : ∀ n, 0 < rn n := fun n => by
    simp only [rn]; rw [sub_pos, div_lt_one (by positivity)]; linarith [n.cast_nonneg (α := ℝ)]
  have hrn1 : ∀ n, rn n < 1 := fun n => by simp only [rn]; linarith [show 0 < 1 / ((n:ℝ) + 2) by positivity]
  let U : Fin N × ℕ → Set E3 := fun p => {y | z + (rn p.2)⁻¹ • (y - t p.1 - z) ∈ interior K}
  have hUo : ∀ p, IsOpen (U p) := fun p => isOpen_interior.preimage (by fun_prop)
  have hcov : K ⊆ ⋃ p, U p := by
    intro y hy
    obtain ⟨i, hi⟩ := mem_iUnion.1 (ht hy)
    rw [mem_tr] at hi
    set w := y - t i
    -- the map r ↦ z + r⁻¹ (w - z) is continuous at r = 1 and equals w ∈ int K there
    have hcont : ContinuousAt (fun r : ℝ => z + r⁻¹ • (w - z)) 1 := by
      apply ContinuousAt.add continuousAt_const
      exact (continuousAt_inv₀ one_ne_zero).smul continuousAt_const
    have hmem : (fun r : ℝ => z + r⁻¹ • (w - z)) 1 ∈ interior K := by simpa using hi
    have hev := hcont.preimage_mem_nhds (isOpen_interior.mem_nhds hmem)
    obtain ⟨ε, hε, hball⟩ := Metric.mem_nhds_iff.1 hev
    obtain ⟨n, hn⟩ := exists_nat_one_div_lt hε
    refine mem_iUnion.2 ⟨(i, n), ?_⟩
    have : rn n ∈ Metric.ball (1:ℝ) ε := by
      rw [Metric.mem_ball, Real.dist_eq]
      simp only [rn]
      rw [show 1 - 1 / ((n:ℝ) + 2) - 1 = -(1 / ((n:ℝ) + 2)) by ring, abs_neg,
        abs_of_pos (by positivity)]
      calc 1 / ((n:ℝ) + 2) ≤ 1 / ((n:ℝ) + 1) := by gcongr; linarith
        _ < ε := hn
    have := hball this
    simpa [U, w, sub_sub] using this
  obtain ⟨s, hs⟩ := hK.2.1.elim_finite_subcover U hUo hcov
  set n0 := s.sup (fun p => p.2)
  refine ⟨rn n0, hrn0 n0, hrn1 n0, fun i => t i + (1 - rn n0) • z, fun y hy => ?_⟩
  obtain ⟨p, hp, hyp⟩ := mem_iUnion₂.1 (hs hy)
  refine mem_iUnion.2 ⟨p.1, ?_⟩
  -- monotonicity in n: rn p.2 ≤ rn n0
  have hle : rn p.2 ≤ rn n0 := by
    have : (p.2 : ℝ) ≤ n0 := by exact_mod_cast Finset.le_sup (f := fun p : Fin N × ℕ => p.2) hp
    simp only [rn]; gcongr
  -- the point k := z + (rn n0)⁻¹ (y - t - z) lies in K (segment between z and the U-point)
  set u := z + (rn p.2)⁻¹ • (y - t p.1 - z)
  have huK : u ∈ K := interior_subset hyp
  have hzK : z ∈ K := interior_subset hz
  set θ := rn p.2 / rn n0
  have hθ0 : 0 ≤ θ := div_nonneg (hrn0 _).le (hrn0 _).le
  have hθ1 : θ ≤ 1 := by rw [div_le_one (hrn0 _)]; exact hle
  have hk : (1 - θ) • z + θ • u ∈ K := hK.1 hzK huK (by linarith) hθ0 (by ring)
  refine ⟨(1 - θ) • z + θ • u, hk, ?_⟩
  simp only [u, θ, smul_add, smul_smul]
  have h1 : rn n0 * (rn p.2 / rn n0 * (rn p.2)⁻¹) = 1 := by
    field_simp [(hrn0 p.2).ne', (hrn0 n0).ne']

  rw [h1]
  module

end Hadwiger

import Hadwiger.Symmetry

/-!
# Prymak's region `D` is not a fundamental domain

The configuration `p = (459/500, 1481/5000, 1689/10000, 9529/10000, 9359/10000, 3079/10000)`
(found by the hostile referee) has no image under the 48 symmetries of the cube in Prymak's
region `D = {p ∈ [0,1/2]³ × [0,1]³ : p₂+p₄ ≤ p₁+p₆ ≤ p₃+p₅}`.  So Prymak's equation (4) (the
reduction of `max_p` over `[0,1]⁶` to `max_p` over `D`) does not follow from the symmetries, and
the certificate must cover `D*` (`exists_sym_mem_DStar`) instead.

The check is a finite computation over `ℤ` (values scaled by `10⁴`, `decide`), transported
to `ℝ` by casting.
-/

namespace Hadwiger

/-- Integer mirror of `flip` (values scaled by `10⁴`). -/
def flipZ (b : Bool) (t : ℤ) : ℤ := if b then 10000 - t else t

/-- The counterexample, in Prymak's (0-based) coordinates, times `10⁴`. -/
def pZ : Fin 6 → ℤ := ![9180, 2962, 1689, 9529, 9359, 3079]

/-- Integer mirror of `10⁴ · toP (gcfg g (ofP p))`, with `τ = σ⁻¹`. -/
def imgZ (τ : Fin 3 → Fin 3) (f : Fin 3 → Bool) (p : Fin 6 → ℤ) : Fin 6 → ℤ := fun k =>
  flipZ (f (kCoord k))
    (if τ (kPair k) = τ (kCoord k) then 0 else p (kIdx (τ (kPair k)) (τ (kCoord k))))

/-- Integer mirror of membership in `D`. -/
def inDZ (q : Fin 6 → ℤ) : Prop :=
  (∀ k, 0 ≤ q k ∧ q k ≤ 10000) ∧ q 0 ≤ 5000 ∧ q 1 ≤ 5000 ∧ q 2 ≤ 5000 ∧
    q 1 + q 3 ≤ q 0 + q 5 ∧ q 0 + q 5 ≤ q 2 + q 4

instance (q : Fin 6 → ℤ) : Decidable (inDZ q) := by unfold inDZ; infer_instance

lemma no_image_Z : ∀ (τ : Fin 3 → Fin 3) (f : Fin 3 → Bool), Function.Injective τ →
    ¬ inDZ (imgZ τ f pZ) := by
  decide

/-- The configuration as a real configuration. -/
noncomputable def pR : Fin 6 → ℝ := fun k => (pZ k : ℝ) / 10000

lemma img_cast (g : Sym) (k : Fin 6) :
    toP (gcfg g (ofP pR)) k = (imgZ g.σ.symm g.f pZ k : ℝ) / 10000 := by
  simp only [toP, gcfg, ofP, imgZ, pR]
  split_ifs <;> cases g.f (kCoord k) <;> simp [flipZ] <;> ring

/-- **Prymak's `D` is not a fundamental domain**: this configuration of `[0,1]⁶` has no
symmetric image in `D`. -/
theorem DPrymak_not_fundamental :
    ofP pR ∈ CfgSpace ∧ ∀ g : Sym, toP (gcfg g (ofP pR)) ∉ DPrymak := by
  constructor
  · intro i j hij
    fin_cases i <;> fin_cases j <;> simp_all [ofP, kIdx, pR, pZ] <;> norm_num
  · intro g hg
    apply no_image_Z g.σ.symm g.f g.σ.symm.injective
    obtain ⟨h0, h1, h2, h3, h4, h5⟩ := hg
    simp only [img_cast] at h0 h1 h2 h3 h4 h5
    refine ⟨fun k => ⟨?_, ?_⟩, ?_, ?_, ?_, ?_, ?_⟩
    · have := (h0 k).1; exact_mod_cast (by linarith : (0:ℝ) ≤ (imgZ g.σ.symm g.f pZ k : ℝ))
    · have := (h0 k).2; exact_mod_cast (by linarith : (imgZ g.σ.symm g.f pZ k : ℝ) ≤ 10000)
    · exact_mod_cast (by linarith : (imgZ g.σ.symm g.f pZ 0 : ℝ) ≤ 5000)
    · exact_mod_cast (by linarith : (imgZ g.σ.symm g.f pZ 1 : ℝ) ≤ 5000)
    · exact_mod_cast (by linarith : (imgZ g.σ.symm g.f pZ 2 : ℝ) ≤ 5000)
    · exact_mod_cast (by linarith : (imgZ g.σ.symm g.f pZ 1 : ℝ) + imgZ g.σ.symm g.f pZ 3 ≤
        imgZ g.σ.symm g.f pZ 0 + imgZ g.σ.symm g.f pZ 5)
    · exact_mod_cast (by linarith : (imgZ g.σ.symm g.f pZ 0 : ℝ) + imgZ g.σ.symm g.f pZ 5 ≤
        imgZ g.σ.symm g.f pZ 2 + imgZ g.σ.symm g.f pZ 4)

end Hadwiger

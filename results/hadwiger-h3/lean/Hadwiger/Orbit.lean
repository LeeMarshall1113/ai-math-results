import Hadwiger.Main
import Hadwiger.ExtraBoxes

/-!
# Orbit coverage of the `1/10`-grid

The cube symmetries act on the `1/10`-boxes of `[0,1]⁶` as signed coordinate permutations
(`gridAct`).  The certified start boxes are the old ones (`oldStartIdx`) and the 570 extension
boxes (`extraList`).  Coverage of all configurations up to symmetry is reduced to the finite
statement `gridCheck = true`: every box with `k₀,k₁,k₂ ≤ 4` has a symmetric image among the
certified start boxes (any configuration is first moved to `p₀,p₁,p₂ ≤ 1/2`).

`gridCheck = true` is established by `native_decide`, i.e. by running compiled code; this adds the
trust assumption `gridCheck_true._native.native_decide.ax_1_1` (an auxiliary axiom asserting what
the compiled evaluation of `gridCheck` returned) to the axioms of the main theorem.  Kernel `decide` on `1.25 · 10⁵` boxes × 48 symmetries is out of reach.
The same fact was checked independently by `replay/orbitboxes.py`.
-/

namespace Hadwiger

open Set

/-- Source coordinate of the signed permutation induced by `g` on Prymak's coordinates. -/
def symSrc (g : Sym) (m : Fin 6) : Fin 6 := kIdx (g.σ.symm (kPair m)) (g.σ.symm (kCoord m))

lemma toP_gcfg (g : Sym) (c : Cfg) (m : Fin 6) :
    toP (gcfg g c) m = flip (g.f (kCoord m)) (toP c (symSrc g m)) := by
  have hne : g.σ.symm (kPair m) ≠ g.σ.symm (kCoord m) :=
    fun h => kPair_ne_kCoord m (g.σ.symm.injective h)
  obtain ⟨h1, h2⟩ := kPair_kIdx _ _ hne
  simp only [toP, gcfg, symSrc, h1, h2]

/-- Action of `g` on grid indices. -/
def gridAct (g : Sym) (k : Fin 6 → ℕ) : Fin 6 → ℕ := fun m =>
  if g.f (kCoord m) then 9 - k (symSrc g m) else k (symSrc g m)

/-- `p` lies in the grid box `k`. -/
def InGrid (k : Fin 6 → ℕ) (p : Fin 6 → ℝ) : Prop :=
  ∀ m, k m ≤ 9 ∧ (k m : ℝ) / 10 ≤ p m ∧ p m ≤ ((k m : ℝ) + 1) / 10

lemma inGrid_gcfg (g : Sym) {k : Fin 6 → ℕ} {c : Cfg} (h : InGrid k (toP c)) :
    InGrid (gridAct g k) (toP (gcfg g c)) := by
  intro m
  obtain ⟨h9, hlo, hhi⟩ := h (symSrc g m)
  rw [toP_gcfg]
  simp only [gridAct]
  cases g.f (kCoord m)
  · simp only [flip_false, Bool.false_eq_true, ite_false]; exact ⟨h9, hlo, hhi⟩
  · simp only [flip_true, ite_true]
    refine ⟨by omega, ?_, ?_⟩ <;> rw [Nat.cast_sub h9] <;> push_cast <;> linarith

lemma startBox_memP {k : Fin 6 → ℕ} {p : Fin 6 → ℝ} (h : InGrid k p) : (startBox k).MemP p := by
  intro m
  obtain ⟨-, h1, h2⟩ := h m
  simp only [startBox]; push_cast
  constructor <;> linarith

/-- Every configuration with `p₀,p₁,p₂ ≤ 1/2` lies in a grid box with `k₀,k₁,k₂ ≤ 4`. -/
lemma exists_grid {c : Cfg} (hc : c ∈ CfgSpace) (h0 : toP c 0 ≤ 1/2) (h1 : toP c 1 ≤ 1/2)
    (h2 : toP c 2 ≤ 1/2) : ∃ k : Fin 6 → ℕ, (∀ m, k m < (if m.val < 3 then 5 else 10)) ∧
      InGrid k (toP c) := by
  have hp : ∀ m, 0 ≤ toP c m ∧ toP c m ≤ 1 := fun m => hc _ _ (kPair_ne_kCoord m)
  set B : Fin 6 → ℕ := fun m => if m.val < 3 then 5 else 10 with hBdef
  refine ⟨fun m => min (B m - 1) ⌊10 * toP c m⌋₊, fun m => ?_, fun m => ?_⟩
  · have : B m - 1 < B m := by simp only [B]; split_ifs <;> omega
    exact lt_of_le_of_lt (min_le_left _ _) this
  · dsimp only
    have h0' := (hp m).1; have h1' := (hp m).2
    have hfl := Nat.floor_le (show 0 ≤ 10 * toP c m by linarith)
    have hlt := Nat.lt_floor_add_one (10 * toP c m)
    have hB9 : B m - 1 ≤ 9 := by simp only [B]; split_ifs <;> omega
    refine ⟨(min_le_left _ _).trans hB9, ?_, ?_⟩
    · have : ((min (B m - 1) ⌊10 * toP c m⌋₊ : ℕ) : ℝ) ≤ ⌊10 * toP c m⌋₊ := by
        exact_mod_cast min_le_right _ _
      linarith
    · rcases le_total ⌊10 * toP c m⌋₊ (B m - 1) with hle | hle
      · rw [min_eq_right hle]; linarith
      · rw [min_eq_left hle]
        by_cases hm : m.val < 3
        · have hpm : toP c m ≤ 1/2 := by
            have : m = 0 ∨ m = 1 ∨ m = 2 := by omega
            rcases this with rfl | rfl | rfl <;> assumption
          simp only [B, hm, ite_true]; norm_num; linarith
        · simp only [B, hm, ite_false]; norm_num; linarith

/-! ### The finite check -/

/-- Sufficient test for membership in `oldStartIdx` (closed form of the `∃ ε` condition). -/
def oldFast (k : Fin 6 → ℕ) : Bool :=
  k 0 < 5 && k 1 < 5 && k 2 < 5 && k 3 < 10 && k 4 < 10 && k 5 < 10 &&
    ((2 * (k 1 + k 3) ≤ 2 * (k 0 + k 5) + 1 && 2 * (k 0 + k 5) + 1 ≤ 2 * (k 2 + k 4) + 4) ||
     (2 * (k 1 + k 3) ≤ 2 * (k 0 + k 5) + 2 && 2 * (k 0 + k 5) + 2 ≤ 2 * (k 2 + k 4) + 4) ||
     (2 * (k 1 + k 3) ≤ 2 * (k 0 + k 5) + 3 && 2 * (k 0 + k 5) + 3 ≤ 2 * (k 2 + k 4) + 4))

lemma oldFast_sound {k : Fin 6 → ℕ} (h : oldFast k = true) : k ∈ oldStartIdx := by
  simp only [oldFast, Bool.and_eq_true, Bool.or_eq_true, decide_eq_true_eq] at h
  refine ⟨fun i => ?_, ?_⟩
  · fin_cases i <;> simp <;> omega
  · by_cases hA : 2 * (k 1 + k 3) ≤ 2 * (k 0 + k 5) + 1 ∧ 2 * (k 0 + k 5) + 1 ≤ 2 * (k 2 + k 4) + 4
    · refine ⟨![false, false, true, false, true, false], ?_⟩
      simp only [inUD]; simp; omega
    by_cases hB : 2 * (k 1 + k 3) ≤ 2 * (k 0 + k 5) + 2 ∧ 2 * (k 0 + k 5) + 2 ≤ 2 * (k 2 + k 4) + 4
    · refine ⟨![true, false, true, false, true, false], ?_⟩
      simp only [inUD]; simp; omega
    · refine ⟨![true, false, true, false, true, true], ?_⟩
      simp only [inUD]; simp; omega

/-- The six permutations of three axes (only soundness matters: any list of symmetries). -/
def permList : List (Equiv.Perm (Fin 3)) :=
  [1, Equiv.swap 0 1, Equiv.swap 0 2, Equiv.swap 1 2, Equiv.swap 0 1 * Equiv.swap 1 2,
    Equiv.swap 1 2 * Equiv.swap 0 1]

/-- The 48 symmetries. -/
def symList : List Sym :=
  permList.flatMap fun σ =>
    [false, true].flatMap fun a => [false, true].flatMap fun b => [false, true].map fun c =>
      ⟨σ, ![a, b, c]⟩

/-- The grid index `![a,b,c,d,e,f]`. -/
def mkIdx (a b c d e f : ℕ) : Fin 6 → ℕ := ![a, b, c, d, e, f]

/-- The finite orbit-coverage check over all boxes with `k₀,k₁,k₂ ≤ 4`. -/
def gridCheck : Bool :=
  (List.range 5).all fun a => (List.range 5).all fun b => (List.range 5).all fun c =>
  (List.range 10).all fun d => (List.range 10).all fun e => (List.range 10).all fun f =>
    symList.any (fun g => oldFast (gridAct g (mkIdx a b c d e f))) ||
    symList.any (fun g => extraList.contains (List.ofFn (gridAct g (mkIdx a b c d e f))))

/-- **Trust assumption** (`native_decide`, adds the auxiliary axiom `gridCheck_true._native.native_decide.ax_1_1`): the finite check passes.
Independently confirmed by `code/replay/orbitboxes.py`. -/
theorem gridCheck_true : gridCheck = true := by native_decide

/-- Certified start-box indices: old start boxes and the 570 extension boxes. -/
def CertIdx (k : Fin 6 → ℕ) : Prop := k ∈ oldStartIdx ∨ List.ofFn k ∈ extraList

lemma grid_orbit_cover {k : Fin 6 → ℕ} (hk : ∀ m, k m < (if m.val < 3 then 5 else 10)) :
    ∃ g : Sym, CertIdx (gridAct g k) := by
  have h := gridCheck_true
  simp only [gridCheck, List.all_eq_true, List.mem_range] at h
  have hk' : k = mkIdx (k 0) (k 1) (k 2) (k 3) (k 4) (k 5) := by
    funext m; fin_cases m <;> rfl
  have := h (k 0) (by simpa using hk 0) (k 1) (by simpa using hk 1) (k 2) (by simpa using hk 2)
    (k 3) (by simpa using hk 3) (k 4) (by simpa using hk 4) (k 5) (by simpa using hk 5)
  rw [← hk'] at this
  simp only [Bool.or_eq_true, List.any_eq_true] at this
  rcases this with ⟨g, -, hg⟩ | ⟨g, -, hg⟩
  · exact ⟨g, Or.inl (oldFast_sound hg)⟩
  · exact ⟨g, Or.inr (by simpa using hg)⟩

/-- **Orbit coverage.**  Every configuration is moved by two cube symmetries into a certified
start box. -/
theorem orbit_cover (c : Cfg) (hc : c ∈ CfgSpace) :
    ∃ g1 g2 : Sym, ∃ k, CertIdx k ∧ (startBox k).MemP (toP (gcfg g2 (gcfg g1 c))) := by
  obtain ⟨g1, hg1⟩ := exists_sym_mem_DStar c hc
  obtain ⟨-, h0, h1, h2, -⟩ := hg1
  obtain ⟨k, hkb, hk⟩ := exists_grid (gcfg_mem_CfgSpace g1 hc) h0 h1 h2
  obtain ⟨g2, hg2⟩ := grid_orbit_cover hkb
  exact ⟨g1, g2, _, hg2, startBox_memP (inGrid_gcfg g2 hk)⟩

end Hadwiger

import Mathlib

/-!
# Basic definitions for the formalisation of `H₃ ≤ 13`

* `E3 = Fin 3 → ℝ` is three-space (its product topology is the Euclidean one; only
  topology, convexity, translations and Lebesgue volume are used).
* `tr t S` is the translate `t + S`.
* `cube = [0,1]³`, `skeleton` its 1-skeleton `E`.
* A *configuration* is `c : Cfg = Fin 3 → Fin 3 → ℝ`; `c i j` (for `i ≠ j`) is the `j`-th
  coordinate of the touching point on the face `xᵢ = 0`; the diagonal is ignored.
  `pt c (i, b)` is the touching point on the face `xᵢ = b`.  In Prymak's notation
  `p = (p₁,…,p₆) = (c 0 1, c 0 2, c 1 0, c 1 2, c 2 0, c 2 1)`, so that `pt c` runs over the
  columns of Prymak's matrix `A_p`, i.e. over `V_p`.
* `O c = conv V_p`.
-/

namespace Hadwiger

open Set

/-- Three-space. -/
abbrev E3 := Fin 3 → ℝ

/-- The translate `t + S`. -/
def tr (t : E3) (S : Set E3) : Set E3 := {y | y - t ∈ S}

@[simp] lemma mem_tr {t y : E3} {S : Set E3} : y ∈ tr t S ↔ y - t ∈ S := Iff.rfl

lemma tr_mono {t : E3} {S S' : Set E3} (h : S ⊆ S') : tr t S ⊆ tr t S' := fun _ hy => h hy

lemma tr_tr (t u : E3) (S : Set E3) : tr t (tr u S) = tr (t + u) S := by
  ext y; simp [sub_sub]

/-- A convex body: compact convex set with non-empty interior. -/
def IsConvexBody (K : Set E3) : Prop := Convex ℝ K ∧ IsCompact K ∧ (interior K).Nonempty

/-- `K` is covered by `N` translates of its interior, i.e. `C(K, int K) ≤ N`. -/
def CoveredByInteriorTranslates (N : ℕ) (K : Set E3) : Prop :=
  ∃ t : Fin N → E3, K ⊆ ⋃ i, tr (t i) (interior K)

/-- The unit cube `[0,1]³`. -/
def cube : Set E3 := Icc 0 1

/-- `x` lies on a face of the cube in coordinate `i`. -/
def OnBdry (x : E3) (i : Fin 3) : Prop := x i = 0 ∨ x i = 1

/-- The 1-skeleton `E` of the cube: points of the cube with two coordinates in `{0,1}`. -/
def skeleton : Set E3 := {x | x ∈ cube ∧ ∃ i j, i ≠ j ∧ OnBdry x i ∧ OnBdry x j}

/-- Configurations (the diagonal entries are ignored). -/
abbrev Cfg := Fin 3 → Fin 3 → ℝ

/-- `0/1` as a real number. -/
def bv (b : Bool) : ℝ := if b then 1 else 0

/-- The touching point of configuration `c` on the face `xᵢ = b`. -/
def pt (c : Cfg) (ib : Fin 3 × Bool) : E3 := fun j => if j = ib.1 then bv ib.2 else c ib.1 j

/-- `V_p`: the six touching points. -/
def V (c : Cfg) : Set E3 := range (pt c)

/-- `O_p = conv V_p`. -/
def O (c : Cfg) : Set E3 := convexHull ℝ (V c)

/-- The configuration space `[0,1]⁶`. -/
def CfgSpace : Set Cfg := {c | ∀ i j, i ≠ j → 0 ≤ c i j ∧ c i j ≤ 1}

/-! ### Prymak's coordinates `p ∈ ℝ⁶` -/

/-- Pair (face index) of Prymak's parameter `p_{k+1}`. -/
def kPair : Fin 6 → Fin 3 := ![0, 0, 1, 1, 2, 2]
/-- Coordinate of Prymak's parameter `p_{k+1}`. -/
def kCoord : Fin 6 → Fin 3 := ![1, 2, 0, 2, 0, 1]

lemma kPair_ne_kCoord (k : Fin 6) : kPair k ≠ kCoord k := by
  fin_cases k <;> decide

/-- Prymak's `p` of a configuration. -/
def toP (c : Cfg) : Fin 6 → ℝ := fun k => c (kPair k) (kCoord k)

/-- Index of the parameter `c i j` (`i ≠ j`) in Prymak's `p`. -/
def kIdx (i j : Fin 3) : Fin 6 :=
  if i = 0 then (if j = 1 then 0 else 1) else if i = 1 then (if j = 0 then 2 else 3)
  else (if j = 0 then 4 else 5)

lemma kPair_kIdx (i j : Fin 3) (h : i ≠ j) : kPair (kIdx i j) = i ∧ kCoord (kIdx i j) = j := by
  fin_cases i <;> fin_cases j <;> simp_all [kIdx, kPair, kCoord] <;> decide

/-- The configuration of Prymak's `p` (diagonal set to `0`). -/
def ofP (p : Fin 6 → ℝ) : Cfg := fun i j => if i = j then 0 else p (kIdx i j)

lemma toP_ofP (p : Fin 6 → ℝ) : toP (ofP p) = p := by
  funext k; fin_cases k <;> simp [toP, ofP, kPair, kCoord, kIdx]

lemma pt_ofP_toP (c : Cfg) : pt (ofP (toP c)) = pt c := by
  funext ib j
  simp only [pt]
  split_ifs with h
  · rfl
  · have hne : ib.1 ≠ j := fun h' => h h'.symm
    simp only [ofP, hne, if_false, toP]
    rw [(kPair_kIdx _ _ hne).1, (kPair_kIdx _ _ hne).2]

lemma V_ofP_toP (c : Cfg) : V (ofP (toP c)) = V c := by simp [V, pt_ofP_toP]

lemma O_ofP_toP (c : Cfg) : O (ofP (toP c)) = O c := by simp [O, V_ofP_toP]

/-- The symmetry-reduced region `D*` (in Prymak's 1-based notation:
`p₁,p₂,p₃ ≤ 1/2`, `|s_z - 1| ≤ |s_y - 1| ≤ |s_x - 1|` with `s_z = p₂+p₄`, `s_y = p₁+p₆`,
`s_x = p₃+p₅`), inside `[0,1]⁶`.  Here `p` is 0-based. -/
def DStar : Set (Fin 6 → ℝ) :=
  {p | (∀ k, 0 ≤ p k ∧ p k ≤ 1) ∧ p 0 ≤ 1/2 ∧ p 1 ≤ 1/2 ∧ p 2 ≤ 1/2 ∧
    |p 1 + p 3 - 1| ≤ |p 0 + p 5 - 1| ∧ |p 0 + p 5 - 1| ≤ |p 2 + p 4 - 1|}

/-- Prymak's region `D = {p ∈ [0,1/2]³ × [0,1]³ : p₂+p₄ ≤ p₁+p₆ ≤ p₃+p₅}` (0-based here).
It is *not* a fundamental domain of the cube symmetries (see `Symmetry.lean`). -/
def DPrymak : Set (Fin 6 → ℝ) :=
  {p | (∀ k, 0 ≤ p k ∧ p k ≤ 1) ∧ p 0 ≤ 1/2 ∧ p 1 ≤ 1/2 ∧ p 2 ≤ 1/2 ∧
    p 1 + p 3 ≤ p 0 + p 5 ∧ p 0 + p 5 ≤ p 2 + p 4}

end Hadwiger

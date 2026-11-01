import Hadwiger.LemmaM

/-!
# The leaf check and its soundness

A *leaf* is a box `P` of configurations together with a certificate.  The certificate data is
exactly what `certcheck.py` / `certcheck_fast.py` check (`constraints()` there):

* `A v` (`v : Fin 8`, the cube vertices in `itertools.product((0,1), repeat=3)` order):
  the translate attached to vertex `v`;
* `B f` (`f : Fin 6`, faces `x=0, x=1, y=0, y=1, z=0, z=1`): the translate of face `f`;
  face `f0`'s rectangle is given to vertex `v0` instead, so `8 + 5 = 13` translates are used;
* for each of the 12 edges (`EDGES` order: pairs `a < b` of vertex indices at distance one)
  breakpoints `t e`, `s e`, and optionally a face `mid e` whose translate covers the middle.

The checked conditions say that each vertex, each corner of each rectangle `R_P ∩ F_f`, and each
edge breakpoint lies in the relevant translate.  **The Lean predicate requires the points to lie
in the *interior* of the translates of `Q_P`** (equivalently: every facet inequality holds
*strictly*).  `certcheck_fast.py` tests `≥ 0`; see the README for why the strict test is needed
(Prymak's Lemma 2.4(ii) fails at some box corners) and for the sampled strict re-check.

* *plain leaf*: 13 translates of `int Q_P`;
* *chord leaf* with axis `i`: for **both** signs `s = ±1`, the translate of vertex `v0` is
  replaced by the merged body `(1-λ) int Q_P + λ t [0, w_s]`, `w_s = e_j + s e_k`
  (`j < k` the other two axes), `λ = 7/10`, `t = 569/1000`.
-/

namespace Hadwiger

open Set

/-! ### Combinatorics of the cube, in `certcheck.py`'s order -/

/-- `CUBE = list(itertools.product((0,1), repeat=3))`. -/
def cubeBit (v : Fin 8) (i : Fin 3) : Bool := (v.val / 2 ^ (2 - i.val)) % 2 = 1

/-- Cube vertex `v`. -/
def cubeVert (v : Fin 8) : E3 := fun i => bv (cubeBit v i)

/-- `EDGES`: pairs `a < b` of vertex indices at L1-distance one, in lexicographic order. -/
def edgeList : Fin 12 → Fin 8 × Fin 8 :=
  ![(0,1), (0,2), (0,4), (1,3), (1,5), (2,3), (2,6), (3,7), (4,5), (4,6), (5,7), (6,7)]

/-- Face `f` (`0..5` = `x=0, x=1, y=0, y=1, z=0, z=1`) as (axis, side). -/
def faceIdx (f : Fin 6) : Fin 3 × Bool := (⟨f.val / 2, by omega⟩, f.val % 2 = 1)

/-- The edge point `a + r (b - a)`. -/
def edgePt (e : Fin 12) (r : ℝ) : E3 :=
  cubeVert (edgeList e).1 + r • (cubeVert (edgeList e).2 - cubeVert (edgeList e).1)

/-- Every edge runs in a coordinate direction `m`, from bit `0` to bit `1`. -/
lemma edgeList_spec : ∀ e : Fin 12, ∃ m : Fin 3, cubeBit (edgeList e).1 m = false ∧
    cubeBit (edgeList e).2 m = true ∧ ∀ i, i ≠ m → cubeBit (edgeList e).1 i = cubeBit (edgeList e).2 i := by
  decide

/-- Every (direction, other bits) pair is an edge of the list. -/
lemma edgeList_surj : ∀ (m : Fin 3) (β : Fin 3 → Bool), ∃ e : Fin 12,
    cubeBit (edgeList e).1 m = false ∧ cubeBit (edgeList e).2 m = true ∧
    ∀ i, i ≠ m → cubeBit (edgeList e).1 i = β i ∧ cubeBit (edgeList e).2 i = β i := by
  decide

lemma bv_of_onBdry {t : ℝ} (h : t = 0 ∨ t = 1) : ∃ b, t = bv b := by
  rcases h with h | h
  · exact ⟨false, by simp [h, bv]⟩
  · exact ⟨true, by simp [h, bv]⟩

/-- Every point of the 1-skeleton lies on one of the 12 listed edges. -/
lemma skeleton_sub_edges {x : E3} (hx : x ∈ skeleton) :
    ∃ e : Fin 12, ∃ r : ℝ, 0 ≤ r ∧ r ≤ 1 ∧ x = edgePt e r := by
  obtain ⟨hc, i, j, hij, hi, hj⟩ := hx
  -- the third coordinate
  obtain ⟨m, him, hjm, hall⟩ : ∃ m : Fin 3, i ≠ m ∧ j ≠ m ∧ ∀ l : Fin 3, l = i ∨ l = j ∨ l = m := by
    fin_cases i <;> fin_cases j <;> simp_all <;> first
      | exact ⟨0, by decide, by decide, by decide⟩ | exact ⟨1, by decide, by decide, by decide⟩
      | exact ⟨2, by decide, by decide, by decide⟩
  obtain ⟨bi, hbi⟩ := bv_of_onBdry hi
  obtain ⟨bj, hbj⟩ := bv_of_onBdry hj
  let β : Fin 3 → Bool := fun l => if l = i then bi else bj
  obtain ⟨e, he1, he2, he3⟩ := edgeList_surj m β
  simp only [cube, mem_Icc, Pi.le_def, Pi.zero_apply, Pi.one_apply] at hc
  refine ⟨e, x m, hc.1 m, hc.2 m, ?_⟩
  funext l
  simp only [edgePt, cubeVert, Pi.add_apply, Pi.smul_apply, Pi.sub_apply, smul_eq_mul]
  rcases hall l with rfl | rfl | rfl
  · rw [(he3 l him).1, (he3 l him).2]; simp [β, hbi]
  · have hli : l ≠ i := fun h => hij h.symm
    rw [(he3 l hjm).1, (he3 l hjm).2]; simp [β, hli, hbj]
  · rw [he1, he2]; simp [bv]

/-! ### Certificate data and the leaf predicate -/

/-- Certificate data for one structure (one sign). -/
structure CoverData where
  A : Fin 8 → E3
  B : Fin 6 → E3
  f0 : Fin 6
  v0 : Fin 8
  mid : Fin 12 → Option (Fin 6)
  t : Fin 12 → ℝ
  s : Fin 12 → ℝ

/-- The covering conditions, for covering sets `S v T` (vertex translates; `T` the
translation vector) and `Sb T` (face translates). -/
def CoverData.Holds (D : CoverData) (P : RBox) (S : Fin 8 → E3 → Set E3) (Sb : E3 → Set E3) :
    Prop :=
  (∀ v, cubeVert v ∈ S v (D.A v)) ∧
  (∀ f ε, pt (ofP (P.vert ε)) (faceIdx f) ∈
      (if f = D.f0 then S D.v0 (D.A D.v0) else Sb (D.B f))) ∧
  (∀ e, match D.mid e with
    | none => edgePt e (D.t e) ∈ S (edgeList e).1 (D.A (edgeList e).1) ∧
              edgePt e (D.t e) ∈ S (edgeList e).2 (D.A (edgeList e).2)
    | some f => f ≠ D.f0 ∧ edgePt e (D.t e) ∈ S (edgeList e).1 (D.A (edgeList e).1) ∧
              edgePt e (D.s e) ∈ S (edgeList e).2 (D.A (edgeList e).2) ∧
              edgePt e (D.t e) ∈ Sb (D.B f) ∧ edgePt e (D.s e) ∈ Sb (D.B f))

/-- `λ = 7/10`. -/
noncomputable def lamC : ℝ := 7 / 10
/-- `t = 569/1000`. -/
noncomputable def tC : ℝ := 569 / 1000

/-- Chord direction for axis `i` and sign `s`: `e_j + s e_k`, `j < k` the other two axes
(`certcheck.chord_dir`). -/
def chordDir (i : Fin 3) (s : ℝ) : E3 :=
  if i = 0 then e 1 + s • e 2 else if i = 1 then e 0 + s • e 2 else e 0 + s • e 1

/-- The other two axes of `i`, in increasing order. -/
def chordAxes (i : Fin 3) : Fin 3 × Fin 3 :=
  if i = 0 then (1, 2) else if i = 1 then (0, 2) else (0, 1)

lemma chordDir_eq (i : Fin 3) (s : ℝ) : chordDir i s = e (chordAxes i).1 + s • e (chordAxes i).2 := by
  fin_cases i <;> simp [chordDir, chordAxes]

lemma chordAxes_ne (i : Fin 3) : (chordAxes i).1 ≠ (chordAxes i).2 := by
  fin_cases i <;> decide

/-- Translate of `int Q_P` by `T`. -/
def intQ (P : RBox) (T : E3) : Set E3 := tr T (interior (QP P))

/-- Translate of the merged body `(1-λ) int Q_P + λ t [0, w]` by `T`. -/
def merged (P : RBox) (w : E3) (T : E3) : Set E3 :=
  {y | ∃ μ : ℝ, 0 ≤ μ ∧ μ ≤ 1 ∧ ∃ q ∈ interior (QP P), y = T + (μ * (lamC * tC)) • w + (1 - lamC) • q}

/-- A plain leaf: 13 translates of `int Q_P` cover `E ∪ R_P` with the given structure. -/
def PlainOK (P : RBox) : Prop :=
  ∃ D : CoverData, D.Holds P (fun _ T => intQ P T) (intQ P)

/-- A chord leaf with axis `i`: for both signs, the structure covers with the vertex-`v0`
translate replaced by the merged body. -/
def ChordOK (P : RBox) (i : Fin 3) : Prop :=
  ∀ s : ℝ, (s = 1 ∨ s = -1) → ∃ D : CoverData,
    D.Holds P (fun v T => if v = D.v0 then merged P (chordDir i s) T else intQ P T) (intQ P)

/-- **The leaf predicate.** -/
def LeafOK (P : RBox) : Prop := PlainOK P ∨ ∃ i, ChordOK P i

/-! ### Soundness -/

lemma intQ_convex (P : RBox) (T : E3) : Convex ℝ (intQ P T) := by
  intro x hx y hy a b ha hb hab
  have := (QP_convex P).interior hx hy ha hb hab
  simp only [intQ, mem_tr] at *
  convert this using 1
  have hT : (a + b) • T = T := by rw [hab, one_smul]
  calc a • x + b • y - T = a • x + b • y - (a + b) • T := by rw [hT]
    _ = a • (x - T) + b • (y - T) := by module

lemma merged_convex (P : RBox) (w T : E3) : Convex ℝ (merged P w T) := by
  rintro _ ⟨μ1, h1, h1', q1, hq1, rfl⟩ _ ⟨μ2, h2, h2', q2, hq2, rfl⟩ a b ha hb hab
  refine ⟨a * μ1 + b * μ2, by positivity, ?_, a • q1 + b • q2,
    (QP_convex P).interior hq1 hq2 ha hb hab, ?_⟩
  · nlinarith
  · have hT : T = a • T + b • T := by rw [← add_smul, hab, one_smul]
    conv_rhs => rw [hT]
    simp only [smul_add, smul_smul, add_smul]
    module

/-- 1-D covering of `[0,1]` by the chain `[0,t] ∪ [t,1]`. -/
lemma unit_cover2 (r t : ℝ) (h0 : 0 ≤ r) (h1 : r ≤ 1) : r ∈ segment ℝ 0 t ∨ r ∈ segment ℝ t 1 := by
  rw [segment_eq_uIcc, segment_eq_uIcc]
  simp only [Set.mem_uIcc]
  by_cases h : r ≤ t
  · left; left; exact ⟨h0, h⟩
  · right; left; exact ⟨(lt_of_not_ge h).le, h1⟩

/-- 1-D covering of `[0,1]` by the chain `[0,t] ∪ [t,s] ∪ [s,1]`. -/
lemma unit_cover3 (r t s : ℝ) (h0 : 0 ≤ r) (h1 : r ≤ 1) :
    r ∈ segment ℝ 0 t ∨ r ∈ segment ℝ t s ∨ r ∈ segment ℝ s 1 := by
  rw [segment_eq_uIcc, segment_eq_uIcc, segment_eq_uIcc]
  simp only [Set.mem_uIcc]
  by_cases h : r ≤ t
  · left; left; exact ⟨h0, h⟩
  · by_cases h' : r ≤ s
    · right; left; left; exact ⟨(lt_of_not_ge h).le, h'⟩
    · right; right; left; exact ⟨(lt_of_not_ge h').le, h1⟩

lemma edgePt_mem_of_segment {S : Set E3} (hS : Convex ℝ S) (e : Fin 12) {r a b : ℝ}
    (ha : edgePt e a ∈ S) (hb : edgePt e b ∈ S) (hr : r ∈ segment ℝ a b) : edgePt e r ∈ S := by
  obtain ⟨α, β, hα, hβ, hαβ, rfl⟩ := hr
  have : α • edgePt e a + β • edgePt e b = edgePt e (α • a + β • b) := by
    simp only [edgePt, smul_eq_mul]
    rw [show ∀ v u : E3, α • (v + a • u) + β • (v + b • u) = (α + β) • v + (α * a + β * b) • u
      from fun v u => by module, hαβ, one_smul]
  rw [← this]; exact hS ha hb hα hβ hαβ

lemma edgePt_zero (e : Fin 12) : edgePt e 0 = cubeVert (edgeList e).1 := by simp [edgePt]
lemma edgePt_one (e : Fin 12) : edgePt e 1 = cubeVert (edgeList e).2 := by simp [edgePt]

/-- A convex set containing all `64` vertex-images of a configuration-affine point map
contains the image of the whole box. -/
lemma pt_mem_of_vertices (P : RBox) (ib : Fin 3 × Bool) {S : Set E3} (hS : Convex ℝ S)
    (hv : ∀ ε, pt (ofP (P.vert ε)) ib ∈ S) (c : Cfg) (hc : P.Mem c) : pt c ib ∈ S := by
  -- affinity of `p ↦ pt (ofP p) ib` in each coordinate
  have haff : ∀ (p : Fin 6 → ℝ) (k : Fin 6) (x y α β : ℝ), α + β = 1 →
      pt (ofP (Function.update p k (α * x + β * y))) ib =
        α • pt (ofP (Function.update p k x)) ib + β • pt (ofP (Function.update p k y)) ib := by
    intro p k x y α β hab
    funext j
    simp only [pt, ofP, Pi.add_apply, Pi.smul_apply, smul_eq_mul, Function.update_apply]
    split_ifs <;> first | (rw [← add_mul, hab, one_mul]) | ring
  -- induction on the number of free coordinates
  have key : ∀ n : ℕ, n ≤ 6 → ∀ p : Fin 6 → ℝ, (∀ k, P.lo k ≤ p k ∧ p k ≤ P.hi k) →
      (∀ k : Fin 6, n ≤ k.val → (p k = P.lo k ∨ p k = P.hi k)) → pt (ofP p) ib ∈ S := by
    intro n
    induction n with
    | zero =>
      intro _ p _ hp
      have : p = P.vert (fun k => decide (p k = P.hi k)) := by
        funext k
        simp only [RBox.vert]
        by_cases h : p k = P.hi k
        · simp [h]
        · simp [h]; exact (hp k (Nat.zero_le _)).resolve_right h
      rw [this]; exact hv _
    | succ n ih =>
      intro hn p hbox hp
      set k : Fin 6 := ⟨n, by omega⟩
      have hseg : p k ∈ segment ℝ (P.lo k) (P.hi k) := by
        rw [segment_eq_Icc ((hbox k).1.trans (hbox k).2)]; exact hbox k
      obtain ⟨α, β, hα, hβ, hαβ, hpk⟩ := hseg
      have hp' : p = Function.update p k (α * P.lo k + β * P.hi k) := by
        rw [← smul_eq_mul, ← smul_eq_mul, hpk]; simp
      rw [hp', haff p k _ _ α β hαβ]
      have hlo : pt (ofP (Function.update p k (P.lo k))) ib ∈ S := by
        apply ih (by omega)
        · intro k'; by_cases h : k' = k
          · subst h; simp [(hbox k).1.trans (hbox k).2]
          · simp [Function.update_apply, h, hbox k']
        · intro k' hk'
          by_cases h : k' = k
          · subst h; simp
          · simp only [Function.update_apply, h, ite_false]
            apply hp k'
            have : k'.val ≠ n := fun h' => h (Fin.ext h')
            omega
      have hhi : pt (ofP (Function.update p k (P.hi k))) ib ∈ S := by
        apply ih (by omega)
        · intro k'; by_cases h : k' = k
          · subst h; simp [(hbox k).1.trans (hbox k).2]
          · simp [Function.update_apply, h, hbox k']
        · intro k' hk'
          by_cases h : k' = k
          · subst h; simp
          · simp only [Function.update_apply, h, ite_false]
            apply hp k'
            have : k'.val ≠ n := fun h' => h (Fin.ext h')
            omega
      exact hS hlo hhi hα hβ hαβ
  have := key 6 le_rfl (toP c) hc (fun k hk => absurd k.isLt (by omega))
  rwa [pt_ofP_toP] at this

/-- A chord of `K` of direction `w` and length factor `τ`. -/
def HasChordDir (K : Set E3) (w : E3) (τ : ℝ) : Prop := ∃ x ∈ K, x + τ • w ∈ K

lemma tC_le_tauM : tC ≤ tauM := by
  unfold tC tauM
  have : (1.414 : ℝ) ≤ Real.sqrt 2 := by
    rw [show (1.414 : ℝ) = Real.sqrt (1.414 ^ 2) by rw [Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_le_sqrt (by norm_num)
  linarith

/-- The merged body sits inside a translate of `int K`. -/
lemma merged_sub {P : RBox} {c : Cfg} (hc : P.Mem c) {K : Set E3} (hK : IsConvexBody K)
    (hO : O c ⊆ K) {w : E3} (hch : HasChordDir K w tauM) (T : E3) :
    ∃ T' : E3, merged P w T ⊆ tr T' (interior K) := by
  obtain ⟨x, hx, hxw⟩ := hch
  refine ⟨T - lamC • x, ?_⟩
  rintro y ⟨μ, hμ0, hμ1, q, hq, rfl⟩
  simp only [mem_tr]
  have hqK : q ∈ interior K := interior_mono hO (interior_QP_sub P c hc hq)
  -- the chord point x + (μ t) w ∈ K
  have hseg : x + (μ * tC) • w ∈ K := by
    have hτ : 0 < tauM := lt_of_lt_of_le (by norm_num [tC]) tC_le_tauM
    have hθ0 : 0 ≤ μ * tC / tauM := by unfold tC at *; positivity
    have hθ1 : μ * tC / tauM ≤ 1 := by
      rw [div_le_one hτ]
      calc μ * tC ≤ 1 * tC := mul_le_mul_of_nonneg_right hμ1 (by norm_num [tC])
        _ = tC := one_mul _
        _ ≤ tauM := tC_le_tauM
    have := hK.1 hx hxw (by linarith : 0 ≤ 1 - μ * tC / tauM) hθ0 (by ring)
    convert this using 1
    rw [smul_add, smul_smul, ← add_assoc, ← add_smul]
    congr 1
    · simp
    · congr 1; field_simp
  have hcombo := hK.1.combo_interior_self_subset_interior (show (0:ℝ) < 1 - lamC by norm_num [lamC])
    (show (0:ℝ) ≤ lamC by norm_num [lamC]) (by ring)
    (Set.add_mem_add (Set.smul_mem_smul_set hqK) (Set.smul_mem_smul_set hseg))
  convert hcombo using 1
  simp only [smul_add, smul_smul]
  module

lemma intQ_sub {P : RBox} {c : Cfg} (hc : P.Mem c) {K : Set E3} (hO : O c ⊆ K) (T : E3) :
    intQ P T ⊆ tr T (interior K) :=
  tr_mono (interior_mono hO |>.trans' (interior_QP_sub P c hc))

/-- Generic soundness of the covering conditions: if each covering set is convex and lies in a
translate of `int K`, then `E ∪ V_c` is covered by `13` translates of `int K`. -/
lemma CoverData.cover (D : CoverData) (P : RBox) (S : Fin 8 → E3 → Set E3) (Sb : E3 → Set E3)
    (hD : D.Holds P S Sb) (hSc : ∀ v T, Convex ℝ (S v T)) (hSbc : ∀ T, Convex ℝ (Sb T))
    {K : Set E3} (hS : ∀ v T, ∃ T', S v T ⊆ tr T' (interior K))
    (hSb : ∀ T, ∃ T', Sb T ⊆ tr T' (interior K)) (c : Cfg) (hc : P.Mem c) :
    ∃ t : Fin 13 → E3, skeleton ∪ V c ⊆ ⋃ i, tr (t i) (interior K) := by
  obtain ⟨hV, hR, hE⟩ := hD
  choose TA hTA using fun v => hS v (D.A v)
  choose TB hTB using fun f => hSb (D.B f)
  -- index the 13 translates by `Fin 8 ⊕ {f // f ≠ f0}`
  let ι := Fin 8 ⊕ {f : Fin 6 // f ≠ D.f0}
  have hcard : Fintype.card ι = 13 := by
    simp only [ι, Fintype.card_sum, Fintype.card_fin, Fintype.card_subtype_compl,
      Fintype.card_unique]
  let eq : ι ≃ Fin 13 := Fintype.equivFinOfCardEq hcard
  let tt : ι → E3 := fun i => match i with
    | Sum.inl v => TA v
    | Sum.inr f => TB f.1
  refine ⟨fun i => tt (eq.symm i), ?_⟩
  have hA : ∀ v, ∀ y ∈ S v (D.A v), y ∈ ⋃ i, tr (tt (eq.symm i)) (interior K) := by
    intro v y hy
    exact mem_iUnion.2 ⟨eq (Sum.inl v), by simpa using hTA v hy⟩
  have hB : ∀ f, f ≠ D.f0 → ∀ y ∈ Sb (D.B f), y ∈ ⋃ i, tr (tt (eq.symm i)) (interior K) := by
    intro f hf y hy
    exact mem_iUnion.2 ⟨eq (Sum.inr ⟨f, hf⟩), by simpa using hTB f hy⟩
  rintro x (hx | hx)
  · -- edges
    obtain ⟨e, r, hr0, hr1, rfl⟩ := skeleton_sub_edges hx
    have hE' := hE e
    have h0 : edgePt e 0 ∈ S (edgeList e).1 (D.A (edgeList e).1) := by
      rw [edgePt_zero]; exact hV _
    have h1 : edgePt e 1 ∈ S (edgeList e).2 (D.A (edgeList e).2) := by
      rw [edgePt_one]; exact hV _
    cases hm : D.mid e with
    | none =>
      rw [hm] at hE'
      obtain ⟨ha, hb⟩ := hE'
      rcases unit_cover2 r (D.t e) hr0 hr1 with h | h
      · exact hA _ _ (edgePt_mem_of_segment (hSc _ _) e h0 ha h)
      · exact hA _ _ (edgePt_mem_of_segment (hSc _ _) e hb h1 h)
    | some f =>
      rw [hm] at hE'
      obtain ⟨hf, ha, hb, hta, hsb⟩ := hE'
      rcases unit_cover3 r (D.t e) (D.s e) hr0 hr1 with h | h | h
      · exact hA _ _ (edgePt_mem_of_segment (hSc _ _) e h0 ha h)
      · exact hB f hf _ (edgePt_mem_of_segment (hSbc _) e hta hsb h)
      · exact hA _ _ (edgePt_mem_of_segment (hSc _ _) e hb h1 h)
  · -- the touching points: rectangles
    obtain ⟨ib, rfl⟩ := hx
    obtain ⟨i, b⟩ := ib
    let f : Fin 6 := ⟨2 * i.val + (if b then 1 else 0), by have := i.isLt; split_ifs <;> omega⟩
    have hf : faceIdx f = (i, b) := by
      simp only [faceIdx, f, Prod.mk.injEq]
      constructor
      · ext; simp; split_ifs <;> omega
      · cases b <;> simp <;> omega
    have hR' := hR f
    rw [hf] at hR'
    by_cases h0 : f = D.f0
    · simp only [h0, ite_true] at hR'
      exact hA _ _ (pt_mem_of_vertices P (i, b) (hSc _ _) hR' c hc)
    · simp only [h0, ite_false] at hR'
      exact hB f h0 _ (pt_mem_of_vertices P (i, b) (hSbc _) hR' c hc)

/-- **Soundness of a leaf.**  If `P` passes the leaf check and `p ∈ P`, then for every convex
body `K ⊇ O_p` that has chords `τ (e_j ± e_k)` for every pair of axes (Lemma M), `E ∪ V_p`
is covered by `13` translates of `int K`. -/
theorem leaf_sound {P : RBox} (hP : LeafOK P) {c : Cfg} (hc : P.Mem c) {K : Set E3}
    (hK : IsConvexBody K) (hO : O c ⊆ K)
    (hch : ∀ j k : Fin 3, j ≠ k → ∃ s : ℝ, (s = 1 ∨ s = -1) ∧ HasChordDir K (e j + s • e k) tauM) :
    ∃ t : Fin 13 → E3, skeleton ∪ V c ⊆ ⋃ i, tr (t i) (interior K) := by
  rcases hP with ⟨D, hD⟩ | ⟨i, hi⟩
  · exact D.cover P _ _ hD (fun _ T => intQ_convex P T) (intQ_convex P)
      (fun _ T => ⟨T, intQ_sub hc hO T⟩) (fun T => ⟨T, intQ_sub hc hO T⟩) c hc
  · obtain ⟨s, hs, hchs⟩ := hch (chordAxes i).1 (chordAxes i).2 (chordAxes_ne i)
    obtain ⟨D, hD⟩ := hi s hs
    refine D.cover P _ _ hD ?_ (intQ_convex P) ?_ (fun T => ⟨T, intQ_sub hc hO T⟩) c hc
    · intro v T; split_ifs
      · exact merged_convex P _ T
      · exact intQ_convex P T
    · intro v T; split_ifs
      · rw [chordDir_eq]; exact merged_sub hc hK hO hchs T
      · exact ⟨T, intQ_sub hc hO T⟩

end Hadwiger

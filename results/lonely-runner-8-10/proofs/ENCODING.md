# The gate CNF and why UNSAT proves a gate

Shared by every package in `proofs/k*/`. The generator is `common/gen_cnf.py`, which is
about 160 lines (half comments) of plain Python with no dependencies. Tags: [known] means cited,
[proved] means proved here, [computed] means produced by a script in this folder.

## 1. What is proved per prime

Fix k, a prime-power level L = q^e (in practice L = k+1), and a prime P with P ∤ L. Set
N = L·P and dist_N(r) = min(r mod N, N − r mod N).

Take a k-tuple x of residues mod N with P ∤ x_i for all i. It is **proper** when either
* (a) gcd(x_j : j ≠ i, N) > 1 for some i, or
* (b) some t ∈ Z/N has (k+1)·dist_N(t·x_i) ≥ N for every i, i.e. ‖t x_i / N‖ ≥ 1/(k+1).

I(k, L, P) is the set of improper tuples.

This is Trakulthongchai's Definition 3 (2511.22427, "(k,ℓ,p)-proper") with ℓ = L. With L = k+1 it
is the hypothesis of Rosenfeld's Lemma 6 (2509.14111): his condition (i), "gcd(S ∪ {(k+1)p}) = 1
for every (k−1)-subset S", is the negation of (a). [known]

**Gate lemma** [known: Rosenfeld 2509.14111 Lemma 6 for ℓ = k+1; Trakulthongchai 2511.22427
Lemma 5 for any ℓ]. Let k ≥ 3, assume LRC(k−1), and let P be a prime with I(k, ℓ, P) = ∅. Then
P divides v_1⋯v_k for every set of k distinct positive integers without a lonely time.

## 2. Theorem: if G(k, L, P) is UNSAT, then I(k, L, P) = ∅  [proved]

The theorem is used in one direction only. The proof is a model construction: from an
improper tuple, build an assignment and check it against every clause block. The comments in
`gen_cnf.py` restate the same argument next to each block of code.

### 2.1 Normal form

Let x be improper.
* **Two coordinates are prime to q.** The prime divisors of N are q and P, and P divides no
  coordinate. So (a) fails iff for every i some j ≠ i has q ∤ x_j, which says that at least
  two coordinates are prime to q.
* **Normalise one of them to 1.** Pick such a coordinate x_i. It is a unit mod N. Let
  u = x_i^{-1} mod N and set x'_j = fold_N(u·x_j) ∈ {1..⌊N/2⌋}, where fold_N(r) = min(r mod N,
  N − r mod N).
* **x' is still improper.**
  * (a) is unchanged, because u is a unit and gcd(−r, N) = gcd(r, N).
  * (b) is unchanged because t ↦ t·u^{-1} permutes Z/N, (t u^{-1})(u x_j) = t x_j, and
    dist_N(−r) = dist_N(r).
  * P ∤ x'_j still holds.
* **Define S.** Now x'_i = 1, and every x'_j is either 1 or lies in D = {2..⌊N/2⌋} ∖ PZ. Let
  S = {x'_j} ∖ {1} ⊆ D. Then |S| ≤ k − 1.
* **The full case.** If |S| = k − 1, then x' is (1, s_1, …, s_{k−1}) up to order, with
  distinct s_j. Condition (a) fails at the position of the 1, so gcd(s_1..s_{k−1}, N) = 1, and
  some s ∈ S is prime to q.
* **Bad coordinates.** Call value v *bad* at time t when (k+1)·dist_N(t v) < N. Since x' is
  improper, every t has a bad coordinate. Badness depends only on the value. The value 1 is
  good exactly at the t with (k+1)·dist_N(t) ≥ N.

### 2.2 The assignment

* x_v := [v ∈ S] for v ∈ D.
* y_c := [some s ∈ S has fold_P(s) = c] for c = 2..(P−1)/2.
* Every counter auxiliary takes its intended meaning: "at least j of the first i inputs are
  true".

### 2.3 Every clause holds

**A1 (level-1 time clauses).** These are the clauses ⋁_{c bad at s} y_c, for s ≤ (P−1)/2 with
(k+1)s ≥ P, where c is bad at s when (k+1)·dist_P(s c) < P.
* Take the time t = L s. Then dist_N(L s v) = L·dist_P(s v). So v is bad at t iff its class is
  bad at s.
* Values of class 1, including the value 1, are good at t, because (k+1)s ≥ P.
* So the bad coordinate at t lies in S, with class c ≥ 2, and y_c is true.

**A2 (at most k−1 of the y_c).** S meets at most |S| ≤ k−1 classes.

**B1 (time clauses).** These are the clauses ⋁_{v ∈ D bad at t} x_v, for ⌈N/(k+1)⌉ ≤ t ≤ ⌊N/2⌋.
* At these t the value 1 is good, because dist_N(t) = t.
* So the bad coordinate is some s ∈ S, and x_s is true.
* Times t in (N/2, N) give the same clauses as N − t. Times t < N/(k+1) are omitted, and
  omitting clauses is always sound.

**B2 and B3 (links).** These are x_v → y_{cls v} and y_c → ⋁_{cls v = c} x_v. Both hold by the
definitions of x and y.

**B4 (|S| ≤ k−1).**
* For each class c ≥ 2 there is an at-most-(k−1) counter over its members. Its outputs
  "≥ j members", for j ≥ 2, are the *extras*.
* A final at-most-(k−1) counter runs over Z = (y_c)_c ++ (x_v)_{cls v = 1} ++ extras.
* Let m_c = |S ∩ class c|. Under the assignment the true literals of Z number
  Σ_{c≥2}[m_c ≥ 1] + m_1 + Σ_{c≥2}(m_c − 1)_+ = Σ_c m_c = |S| ≤ k−1.
* A sequential counter whose auxiliaries carry their intended meaning satisfies all of its
  clauses whenever the bound holds. See `Builder.at_most`: each of its four clause shapes is
  checked against that meaning there.

**B5 (full support needs a unit).** The clause is "Z-count ≥ k−1 → ⋁_{q ∤ v} x_v". By B4
the count equals |S|, and by 2.1, |S| = k−1 forces some s ∈ S prime to q.

So the assignment is a model. Contrapositively, UNSAT means there is no improper tuple. ∎

Only the existence of the model matters. Clauses that are merely implied, such as A1 and A2,
or that only break symmetry, cannot make the argument unsound, as long as the constructed
model satisfies them.

## 3. Evidence that the code implements the text [computed]

`common/test_encoding.py` runs 37 small cases in 12 s. It covers L ∈ {2, 3, 4, 5, 7, 8, 9},
both L = k+1 and L ≠ k+1, k = 2..8, and P = 3..19.

**FULL cases (17).**
* It enumerates **all** of B(L,P)^k from the definition, with no normal form. There are up
  to 810,000 tuples per case.
* For every improper tuple it builds the assignment of §2.2 by fixing x and y, and confirms
  with a SAT call that the CNF has a model extending it.
* It also enumerates every model projected onto x, and checks that each decodes to an
  improper tuple, so the CNF is *exact*, not just sound.
* The two support sets coincide.
* Improper tuples occur in 5 of these cases: for example (k,L,P) = (4,5,7) with 13,824
  improper tuples, and (3,9,5) with 864.

**NF cases (20).**
* It enumerates every candidate support S ⊆ D with |S| ≤ k−1, using an independent
  bitmask coverage test, and compares against the model projection in the same two
  directions.
* This includes (7,8,11) with 526 improper supports and (7,8,13) with 72. These are the
  real level of the k = 7 proof, at primes where the gate fails.

**Large cases.**
* The SAT/UNSAT verdicts reproduce the published gate tables exactly: k = 7 at p = 11..167, and
  k = 8 at p = 11..241. Details are in the per-k PROOF.md.
* Every SAT answer is decoded and re-checked against the definition, by `prove.py` and
  `improper_def`.

## 4. What a Lean proof of the encoding would need [speculation: design only]

The goal is `gate k L P` ← `¬ Satisfiable (G k L P)`, plus a verified LRAT checker. One
option is to reflect cake_lpr's guarantee. Another is a Lean LRAT checker: Lean 4 core has
`Std.Tactic.BVDecide.LRAT` with a soundness proof, and `lean-sat`-style checkers exist.
The pieces:

1. **Definitions.**
   * `improper k L P x` over `Fin k → ZMod (L*P)`.
   * The CNF as a Lean function `G : ℕ → ℕ → ℕ → CNF`, which is a transcription of
     `gen_cnf.py`, about 80 lines of code.
   * The bridge between `G` and the DIMACS file is checked by evaluation: `decide` or
     `native_decide` on the hash, or the Lean checker reads `G` directly and never touches the file.
2. **Normal-form lemma (2.1).** Units in `ZMod N` and invariance of properness under unit
   scaling and negation. This needs `ZMod.unitOfCoprime`, gcd facts for N = q^e·P, and
   `Finset` cardinality for |S| ≤ k−1. It is moderate: about 200–400 lines.
3. **Model lemma (2.2–2.3).** For each block, the constructed assignment satisfies it.
   * A1 needs dist_N(L s v) = L·dist_P(s v). This is elementary `Int.emod` arithmetic.
   * B4 needs a generic proof that the sequential counter with intended auxiliaries satisfies
     its clauses. That is reusable and about 150 lines.
   * The counting identity Σ = |S| is a `Finset.sum_fiberwise` argument.
4. **Glue.** UNSAT of `G` gives `I = ∅`, and the published gate lemma then gives `P ∣ ∏ v`.
   The gate lemma, MSS Theorem A, and Rosenfeld's Lemma 4 would remain cited axioms unless
   they are formalised too. They are the larger part of the work.

The bottleneck is not the encoding. It is checking several GB of LRAT inside Lean: Lean's
`bv_decide` checker is used on proofs of about 10^6 steps. The cheaper route is to trust
cake_lpr, which is verified in HOL4, and to prove items 1–3 only.

## 5. Levels with several prime factors (k = 9, L = 10): `common/gen_cnf_composite.py` [proved; not yet used in a package]

When L has two or more prime factors, an improper tuple need not contain a unit, so the
anchor "x_1 = 1" of §2.1 is not available. The changes:

* **When (a) fails.** It fails iff, for **every** prime q | L, at least two coordinates are
  prime to q.
* **Anchor.**
  * Let Q = q*^e be the largest prime-power part of L, and M = L/Q.
  * Pick a coordinate x_i prime to q*. It is a unit mod QP.
  * Use the unit u with u ≡ x_i^{-1} (mod QP) and u ≡ 1 (mod M). It exists by CRT and is
    a unit mod N.
  * After multiplying by u and folding, the anchor lies in A = fold_N{1 + jQP : 0 ≤ j < M}.
    For L = 10 that is {1, 5P−1}.
* **Support.** S is the set of folded values of *all* coordinates, so |S| ≤ k and S ∩ A ≠ ∅.
* **gcd clauses.**
  * For each q | L, some s ∈ S is prime to q.
  * If |S| = k, then all values are distinct, so two elements of S are prime to q.
  * "Two of" uses prefix chains with *downward-sound* definitions: o_i → o_{i−1} ∨ c_i,
    t_i → t_{i−1} ∨ c_i, t_i → t_{i−1} ∨ o_{i−1}. A model cannot set them true spuriously,
    and their intended values satisfy them.
* **Times.** The time clauses cover every t in 1..⌊N/2⌋, because the anchor 5P−1 is not bad
  where 1 is.
* **Unchanged.** Block A is unchanged, because the anchor has class 1. Each model check of §2.3
  goes through as before.

The CNF is claimed sound only, not exact. Its gcd condition is weaker than "some padding
makes (a) fail". It was tested by full brute force over B(L,P)^K for 10 cases with L ∈ {6, 10},
including (5,6,7) with 552,960 improper tuples (`test_encoding.py --composite`). Every
improper tuple's support extends to a model, and in these cases no spurious model appeared.

## 6. Split mode: a level-1 case split plus one LRAT proof per leaf (k ≥ 8) [proved]

For k = 8 the single whole-gate proof is too large: at p ≥ 199 CaDiCaL had written more than
17 GB of LRAT per prime after 52 minutes, and the run was stopped. So the gate is split.

* **The split.** `common/level1.py`, about 60 lines, is trusted. Its docstring states and proves
  the theorem. From (k, P) alone it computes a list of leaves (I_j, E_j), which are disjoint
  class sets.
* **What it guarantees.** The y-part of every model of G satisfies block A, so it is a
  *cover*: at most k−1 classes whose level-1 bad sets meet every s in T1. Every cover
  contains some I_j and avoids the matching E_j.
* **The search.** It is a depth-first search that branches on the uncovered time with the
  fewest candidates, with two dead-end rules:
  * |I| = k−1 with a time still uncovered;
  * the counting bound r·max|cov ∩ U| < |U|.
* **Leaf CNFs.** The j-th leaf CNF is G with the unit clauses y_c (c ∈ I_j) and ¬y_c
  (c ∈ E_j) appended. It is plain concatenation, done in verify.py/split.py.
* **The checks.** CaDiCaL refutes each leaf CNF. Both checkers must accept every leaf proof
  against the regenerated leaf CNF.
* **Conclusion.** If every leaf CNF is UNSAT then G is UNSAT. Any model of G would extend some
  leaf, because its y-part is a cover. By §2 or §5, I(k, L, P) = ∅.

Only the leaf list's *completeness* is trusted, and that rests on the docstring proof. A
wrong leaf list could only lose covers. It cannot make a refutation accept.

`checkers/level1.cpp` is a C++ port. It is about 20× faster and is used for k = 9, where the
Python search would take about 30 h.
* On the cases compared it prints byte-identical leaf lists: k = 7 p = 163; k = 8 p = 53, 89,
  151, 241; k = 9 p = 137, 139, 181.
* For k = 9 the port is part of the trusted base.

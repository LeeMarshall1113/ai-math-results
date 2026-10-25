# Hostile referee: π(K₄⁻) ≤ 14337150771693633/5·10¹⁶ ≈ 0.2867430154

Referee pass on 2026-10-24 by an agent that had not seen the work before. My brief was to break the claim.

**Verdict: SURVIVES.** I found no error in the mathematics, the encoding, the checker or the arithmetic. Two independent evaluations agree to the exact rational:
- the original checker, which I reran;
- my own re-implementation, which uses a different algorithm.

The weakest points are listed at the end.

## 1. Degree regularity: the mathematics

**What the certificate actually uses.** cert_C has no `linear_stat` rows. Its non-SOS terms are of two kinds.

*Interval rows* use 24 five-vertex 1-flags g:

  p_g([[d·g]]₁ − t₁[[g]]₁) − q_g([[d·g]]₁ − t₂[[g]]₁), with p, q ≥ 0.

Here d is the rooted edge through the root.

*Product rows* use 3 four-vertex graphs g:

  pp_g(e·g − t₁g) − pq_g(e·g − t₂g), with pp, pq ≥ 0.

The bound holds if every one of these terms has φ ≥ 0 at an extremal limit φ with φ(e) = π ∈ [t₁, t₂]. The SOS blocks need no hypothesis.

**The lemma.** Some extremal limit has degree d_v = π for almost every root v. This is true, and there are two standard routes to it.

- **Deletion.** Start from extremal n-vertex graphs and repeatedly delete a vertex of degree below (π−ε)·C(n−1,2). Because ex(n)/C(n,3) decreases to π, only o(n) vertices are ever deleted. Take ε → 0 diagonally. The limit has φ(e) = π and minimum degree at least π almost surely. The degrees average to φ(e) = π, so d_v = π almost surely.
- **Razborov's differential method.** arXiv 2609.27495 (the tetrahedron paper), Prop. 2.6, cites Razborov, Cor. 4.6(a), and Baber §5.1.1. It shows that every maximiser of φ(ρ) satisfies φ(S) = 0, with S = 3(ρ² − [[d·d]]₁). S is minus three times the variance of the root degree, so d_v is constant almost surely. That constant is E d_v = π. The hypothesis is that the theory is vertex-uniform, which holds for K₄⁻-free 3-graphs, just as it does for K₄-free ones.

Either way, φ^v(g)·d_v = π·φ^v(g) almost surely, so φ([[d·g]]₁) = π·φ([[g]]₁) for every 1-flag g.

The tetrahedron paper encodes the same fact as the single linear equality S = 0. This work instead uses the consequence for every g, which is bilinear in π, and linearises it on an interval. The two approaches are logically equivalent ways of using the same property. The linearised one needs π ∈ [t₁, t₂].

**The encoding.** I checked this against check_cert.py and against my own code.
- cG_g/105 counts (v, 4-set S) with (H[v+S], v) ≅ g. There are 7·C(6,4) = 105 such pairs, so this is [[g]]₁ in 7-vertex form.
- cD_g adds the condition that v plus the remaining pair is an edge, which is [[d·g]]₁: the product of a 3-vertex and a 5-vertex 1-flag has 3+5−1 = 7 vertices.
- In the limit, φ(cD − t cG)/105 = (π − t)·φ([[g]]₁). This is ≥ 0 for t = t₁ and ≤ 0 for t = t₂, so the signs p ≥ 0 on the t₁ term and −q on the t₂ term are the right ones.
- The product rows rely on φ(e·g) = φ(e)φ(g), which holds for every limit. With the edge on V∖U, the 35 sets U give exactly the 7-vertex expansion of e·g.

**The interval logic.** Let π = φ(e) at a regular extremal limit.
- If π > t₂ = 0.28689: impossible, because the plain certificate (no regularity used) gives π ≤ 0.2868885176 < t₂.
- If π ∈ [t₁, t₂]: then π ≤ max_H val_C(H) = 0.28674301543.
- If π < t₁ = 0.2867: the claim holds trivially.

So π ≤ max(t₁, v) = v. This is exactly the logic of check_bound.py: the cover must reach the plain bound, the worst case is max(t₀, v_i), and an interval with v < t₁ is treated as excluded.

**The SOS blocks.** Each block computes c_b(H)/((7)_s·C(7−s, k)). This is the 7-vertex expansion of [[F_a F_b]]_σ summed over the Gram matrix G = Σ v vᵀ, with integer v, so G is PSD by construction. The checker requires each injective θ to induce the labelled type exactly, so it uses the type correctly. Leaving out some types (only 9 of the 5-vertex types appear) is allowed: any subset of SOS terms is valid.

## 2. The exact certificate: independent reproduction

**The original checker.** I reran it on 1 thread:

    check_bound.py 0.2867433 certs/cert_plain.json.gz certs/cert_C.json.gz

Result: PASS, "BOUND 2867433/10000000 ESTABLISHED". The log is in `notes/referee_rerun_check_bound.out`.

**My own code.** This ran; the log is in `notes/referee_indep.out`.

`notes/indep_catalogue.py` checks the catalogue:
- It brute-forces the labelled K₄⁻-free counts: 11, 253 and 33652 for n = 4, 5, 6.
- It counts n = 7 as pairs (labelled 6-vertex G, link L) by a direct numpy test, without the superset-sum transform. The count is 35333735.
- For the 8157 classes, it confirms each is K₄⁻-free, that no two are isomorphic (min over all 7! relabellings) and that the orbit sum is 35333735. **Match.**

`notes/indep_check.py` evaluates val(H) by a different algorithm:
- It sums over all 5040 vertex orderings and divides by the split stabiliser (k!², 48 or 144).
- It identifies flags through a lookup table of every labelled copy, not a canonical form.
- It computes the Gram matrices exactly from split int64 products.
- It evaluates all 8157 graphs exactly, with Fractions.

| certificate | max val | attained at |
|---|---|---|
| plain | 3012329434941157/10500000000000000 | H = 4132800 |
| C | 14337150771693633/50000000000000000 | H = 0 |

Both maxima, and the graphs where they are attained, are identical to check_cert.py.

**Hunting for checker bugs.**
- Negative p, q, pp or pq is rejected (asserted).
- μ is free, which is correct for equality rows, and cert_C uses none.
- Flags are checked to extend the type, to be K₄⁻-free, canonical and distinct.
- A float factor would be truncated silently by `int()`. That is harmless: the truncated G is still PSD and is the one used for the evaluation.
- The checker does not assert t₁ ≤ t₂. With t₁ > t₂ a certificate could lower its values, but check_bound then requires the cover to reach the plain bound and reports max(t₀, ·), so no false bound results. It holds here anyway (t₁ < t₂). Hardening: add `assert t1 < t2`.

## 3. The final number

The value 14337150771693633/5·10¹⁶ = 0.28674301543387… comes straight out of the exact checkers above.

- It is below 0.2867433, the bound check_bound was asked to establish.
- It is below t₂, and above t₁, so the interval is not excluded and the value is what counts.
- The plain value 0.28688851761… is below t₂, so the interval cover is complete.
- The improvement over 0.286889 is 1.46·10⁻⁴, as stated.

## 4. Novelty (quick check only)

I ran web searches for 2013–2026 bounds on π(K₄⁻). I found nothing below 0.2868 (0.286889, Falgas-Ravry–Vaughan). The newest relevant item, arXiv 2609.27495, concerns K₄, not K₄⁻, and its text does not mention K₄⁻.

This was not an exhaustive literature search. Before going public, it is worth checking specifically whether anyone has already published "FRV plus degree regularity" for K₄⁻, e.g.:
- Baber's thesis, §5.1, which used degree regularity for K₄;
- Flagmatic example files.

## Weakest points (none fatal)

1. **REPORT.md states the regularity lemma without proof.** It says "min-degree deletion lemma" and nothing more. A write-up should give the deletion argument or cite Razborov's Cor. 4.6, as arXiv 2609.27495 does. A sketch is in §1 above.
2. **The two checkers share a pipeline.** Both read the same certificate format and conventions, and the certificate came from one author's pipeline. My evaluation code is independent, but it read the same JSON. There is no formal (Lean) check.
3. **The novelty check was quick**, and specifically did not cover prior use of degree regularity for K₄⁻.
4. **Small checker hardening:** add `t1 < t2`, and reject non-integer factors rather than truncating them.

If the regularity step were rejected for some reason, the plain certificate alone still gives the exact bound 0.2868885176, which matches the old record.

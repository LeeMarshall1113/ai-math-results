# LRC(8), nine runners: certificate-checked replication

**Status [computed, 2026-10-21]: VERIFIED.** All 35 gates were checked, and every one of their
38,734 leaf certificates was accepted by **both** cake_lpr and gate_lrat_check. The exact size
inequality holds.

    python proofs/common/verify.py proofs/k8

This replicates Trakulthongchai, arXiv 2511.22427 (gate set S_8), and Rosenfeld,
arXiv 2512.01912, using independent code.

## 1. The reduction [known; cited]

k = 8 speeds (9 runners). Assume a counterexample v_1..v_8: distinct positive integers with
gcd 1 and no lonely time.

1. **LRC(7)** is the package `../k7` (certificate-checked here). It also rests on LRC(6),
   Barajas–Serra 2008.
2. **Size bound.** MSS Theorem A (Forum Math. Sigma 13 (2025) e164) plus AM–GM, as in
   Rosenfeld 2509.14111 Cor. 3 = 2512.01912 Cor. 2 and Trakulthongchai 2511.22427 footnote 1:
   ∏ v_i ≤ B_8 = (36^7/8)^8, with ln B_8 = 184.0415.
3. **Forced divisors.** lcm(2..9) = 2520 divides ∏ v_i (Rosenfeld 2509.14111 Lemma 4).
4. **Gates.** Trakulthongchai 2511.22427 Lemma 5 with ℓ = 9 = k+1 (= Rosenfeld 2509.14111
   Lemma 6 = 2512.01912 Lemma 3), under LRC(7) and p prime: I(8, 9, p) = ∅ ⇒ p | ∏ v_i.
5. **Conclusion.** 2520·∏ gates ≤ ∏ v_i ≤ B_8 contradicts the exact test
   2520·∏gates·8^8 > 36^56 (`common/mass.py`).

| | value |
|---|---|
| ln B_8 − ln 2520 | 184.0415 − 7.8320 = 176.2095 needed |
| our 35 gates | **178.7309** |
| the full S_8, 47..241 (39 primes) | 189.4834 |

**Gate set.** 35 primes, chosen to minimise the work subject to the exact bound.
* 27 of them are from S_8: 47, 53, 59, 61, 71, 79, 107, 131, 139, 149, 151, 157, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241.
* The other 8, 251, 257, 263, 269, 271, 277, 281, 283, lie beyond Trakulthongchai's list. They
  are certified here as gates.
* Left out: the primes in 67..163 with the largest level-1 splits. For example 73 has 27,693
  leaves and 109 has 20,636, while the large primes need only 8 to 128 leaves each.

## 2. Certificates [computed]

**Why split mode.** A single whole-gate LRAT proof does not scale to k = 8.
* At p = 199..241, CaDiCaL had written 14–17 GB of LRAT per prime after 52 minutes. It was
  stopped before finishing, to protect the shared disk.
* So each gate uses the level-1 split of `../ENCODING.md` §6. The trusted `common/level1.py`
  lists leaves (I, E). Each leaf CNF G(8,9,p) + units gets its own CaDiCaL LRAT proof, and
  both cake_lpr and gate_lrat_check must accept every one.

| | |
|---|---|
| leaves (= LRAT proofs) | 38,734: 28,018 at the six gates 47..79, and 8 to 2,000 per gate at p ≥ 107 |
| proof volume | **21.1 GB** of text LRAT; ≥ 83.4M additions (step counts were not logged for the six run-1 gates) |
| level-1 split (level1.py, Python) | 4,252 s in total |
| solve | 3,536 CPU-s |
| gate_lrat_check | 1,799 CPU-s |
| cake_lpr | 2,163 s wall (parallel in the final run; serial over /mnt/c in run 1) |
| total wall | 7,121 s, summed over gates |

**Runs.** `verify_gates.jsonl` keeps one line per gate.
* Gates 47, 53, 59, 61, 71 and 79 were checked in run 1, on Windows, with cake_lpr reading over
  /mnt/c. See `verify_run1_partial.log`.
* The other 29 were checked in the final run, inside WSL. Its summary line is in
  `verify_log.jsonl`, and the full log is `verify_final.log`.
* Each gate was accepted by both checkers within a single run.
* The first pass `split_records.jsonl` (p = 11..59) is a separate earlier production run. It
  also certifies the non-gates.

**Non-gates 11..43 reproduced.** Every prime 11 ≤ p ≤ 43 has a SAT leaf. The decoded tuple
was re-checked as improper from the definition, e.g. p = 43 in split_records.jsonl. This agrees
with Trakulthongchai: I(8,9,p) ≠ ∅ exactly for 11 ≤ p ≤ 43.

## 3. Trusted base

The k7 trusted base (cited lemmas, ENCODING.md, gen_cnf.py, mass.py, cake_lpr /
gate_lrat_check), plus:
* `common/level1.py` (~60 lines, with the completeness proof in its docstring);
* the leaf-CNF concatenation in `split.py` (`leaf_cnf_text`, 4 lines).

Not trusted: CaDiCaL, `level1.cpp`, lr.py, and the hashes.

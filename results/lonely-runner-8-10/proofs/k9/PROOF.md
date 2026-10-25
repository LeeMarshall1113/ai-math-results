# LRC(9), ten runners: certificate-checked replication

**Status [computed, 2026-10-21]: VERIFIED.** All 45 gates were checked, and every one of their
31,593 leaf certificates was accepted by **both** cake_lpr and gate_lrat_check. The exact size
inequality holds. LRC(9) therefore follows from the cited lemmas, together with LRC(8)
(`../k8`) and LRC(7) (`../k7`).

    LEVEL1_EXE=proofs/common/checkers/bin/level1 python proofs/common/verify.py proofs/k9

This replicates Trakulthongchai, arXiv 2511.22427 (k = 9, gate set S_9 = primes 137..401),
using independent code.

## 1. The reduction [known; cited]

1. **LRC(8)** is the package `../k8` (certificate-checked here).
2. **Size bound.** MSS Theorem A plus AM–GM, as in Trakulthongchai, footnote 1:
   ∏ v_i ≤ B_9 = (45^8/9)^9, with ln B_9 = 254.3047.
3. **Forced divisors.** lcm(2..10) = 2520 divides ∏ v_i (Rosenfeld 2509.14111 Lemma 4,
   stated for general k).
4. **Gates.** Trakulthongchai 2511.22427 Lemma 5, which holds for *any* level ℓ. We use
   ℓ = 10 = k+1, under LRC(8) and p prime: I(9, 10, p) = ∅ ⇒ p | ∏ v_i.
   * 10 is not a prime power, so the CNF is `common/gen_cnf_composite.py`
     (`../ENCODING.md` §5). Its gcd condition uses both primes 2 and 5.
   * Trakulthongchai computes I(9,10,p) inside the fibres of the chains 1→2→10 and 1→5→10.
     We certify emptiness at level 10 directly.
5. **Conclusion.** The exact test is 2520·∏gates·9^9 > 45^72 (`common/mass.py`).

| | value |
|---|---|
| ln B_9 − ln 2520 | 246.4727 needed |
| our 45 gates | **250.9207** |
| S_9 = 137..401 (47 primes) | 259.9379 |

**Gate set.** The 45 primes 139..409 except 151 and 163.
* 137, 151 and 163 were left out because their level-1 splits are the largest (8,280,
  33,282 and 9,045 leaves).
* 409 lies beyond Trakulthongchai's list. It was included because its split has only 9 leaves.

The level-1 split survey for 131..409 is summarised in §2.

## 2. Certificates [computed]

* **Split mode,** as for k = 8 (`../ENCODING.md` §6).
  * The leaves come from `checkers/level1.cpp`. Its output is byte-identical to `level1.py` on
    every case compared. The Python version would take about 30 h at k = 9.
  * Each leaf CNF is G(9,10,p) from `gen_cnf_composite.py`, with the leaf's unit clauses appended.
  * CaDiCaL 3.0.1 (Linux build, same commit c607304) writes the LRAT, and both checkers must
    accept it.

| | |
|---|---|
| leaves (= LRAT proofs) | 31,593, from 9 (p = 409) to 4,683 (p = 139) per gate |
| level-1 search | 1.6M nodes (p = 139) to 220M (p = 401); 7,688 s in total |
| proof volume | **50.9 GB** of text LRAT, ≥ 290.8M additions (step counts were not logged for the four run-1 gates) |
| solve | 3,748 CPU-s |
| gate_lrat_check | 2,080 CPU-s |
| cake_lpr | 5,881 s wall, as six parallel checkers; about 4 CPU-h at about 3.4 MB/s |
| total wall | 14,496 s, summed over gates |

**Runs.** `verify_gates.jsonl` keeps one line per gate.
* Gates 139, 149, 157 and 167 were checked in run 1. That run used Windows CaDiCaL and cake_lpr
  reading over WSL /mnt/c. It was stopped afterwards because the 9p I/O was slow. See
  `verify_run1_partial.log`.
* The remaining 41 gates were checked in the final run, entirely inside WSL (`common/run_wsl.sh`, with
  `--resume`). Its summary is in `verify_log.jsonl`.
* Each gate was accepted by both checkers in a single run. No gate's verdict was pieced
  together from partial runs.

**Not done.** The non-gates below 137 were not certified for k = 9. The only small-prime
evidence is the level-1 survey, which gives the leaf counts for 131..409.

## 3. Trusted base

The k8 base, with `gen_cnf_composite.py` in place of `gen_cnf.py` and `level1.cpp` in place
of `level1.py`.

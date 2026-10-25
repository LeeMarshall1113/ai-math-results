# LRC(7), eight runners: certificate-checked replication

**Status [computed, 2026-10-21]: VERIFIED.** All 28 gate certificates are accepted by both
checkers: cake_lpr (CakeML-verified) and gate_lrat_check (ours). The exact size inequality holds.

One command re-checks everything:

    python proofs/common/verify.py proofs/k7

On this box it is run as `CAKE_LPR=<path to cake_lpr> python ...`.

This replicates Rosenfeld, arXiv 2509.14111 (v2, Oct 2025), using independent code. The
gate CNF and the checker were written from the definitions, and no code from the paper was used.

## 1. The reduction [known; cited]

Notation: k = 7 speeds (8 runners), and LRC(k) is the statement for k speeds. A
*counterexample* is a set of k distinct positive integers v_1..v_k with no t satisfying
‖t v_i‖ ≥ 1/(k+1) for all i. Scaling by the gcd preserves this, so we may assume gcd = 1.

1. **Induction hypothesis LRC(6).** This is Barajas & Serra, *The lonely runner with seven
   runners*, Electron. J. Combin. 15 (2008) R48. It is not re-checked here.
2. **Size bound.**
   * Malikiosis–Santos–Schymura, Forum Math. Sigma 13 (2025) e164 (arXiv 2411.06903),
     Theorem A, under LRC(k−1): a counterexample with gcd = 1 has
     Σ_{S⊆[k]} gcd(v_S) ≤ C(k+1,2)^{k−1}.
   * By AM–GM, as in Rosenfeld 2509.14111 Corollary 3:
     ∏ v_i ≤ B_k := (C(k+1,2)^{k−1}/k)^k.
   * For k = 7, ln B_7 = 126.3312.
3. **Forced divisors** (Rosenfeld 2509.14111 Lemma 4). lcm(2..k+1) divides ∏ v_i. For k = 7,
   lcm(2..8) = 840.
4. **Prime gates** (Rosenfeld 2509.14111 Lemma 6, under k ≥ 3 and LRC(k−1); equivalently
   Trakulthongchai 2511.22427 Lemma 5 with ℓ = k+1).
   * If every 7-tuple of residues mod 8p, none divisible by p, with gcd(S ∪ {8p}) = 1 for all
     6-subsets S, has a time t with ‖t v_i/(8p)‖ ≥ 1/8 for all i, then p | ∏ v_i.
   * That hypothesis says exactly I(7, 8, p) = ∅ (ENCODING.md §1).
5. **Conclusion.** The 28 gates below are distinct primes > 8, so they are coprime to 840 and
   to each other. By items 3 and 4, 840·∏ gates divides ∏ v_i. Then
   ∏ v_i ≥ 840·∏ gates > B_7 ≥ ∏ v_i, which is a contradiction. So LRC(7) holds.

The inequality 840·∏ gates > B_7 is checked in exact integer arithmetic by `common/mass.py`.
It compares 840·∏gates·7^7 with 28^42.

| | value |
|---|---|
| ln B_7 | 126.3312 |
| ln 840 | 6.7334 |
| needed from gates | > 119.5978 |
| ln ∏ gates (28 primes) | **125.6289** |
| Rosenfeld's 27 primes alone (31..163, no 41) | 120.5109, which also suffices |

**Strictness.** MSS Theorem A is strict (">"), and Rosenfeld's Corollary 3 is non-strict.
We prove a strict ">" against ∏ ≤ B_k, which covers both readings.

**Gate set.** The 28 gates are {31, 37, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167}. This is every prime in
31..167 except 41.

## 2. The certificates [computed]

Per gate p there is one CNF, G(7, 8, p), from `common/gen_cnf.py`, and one LRAT refutation of
it. The theorem "G UNSAT ⇒ I(7, 8, p) = ∅" is proved in `../ENCODING.md` §2.

| | |
|---|---|
| solver (untrusted) | CaDiCaL 3.0.1 (git c607304), `-q --lrat --no-binary`; deterministic, so proofs reproduce byte for byte |
| CNF sizes | 1,195 vars / 2,382 clauses (p=31) to 6,907 / 13,857 (p=167) |
| total proof | 14.5M LRAT additions, **5.48 GB** text (716 MB as .xz) |
| largest proof | p=167: 1.80M steps, 872 MB |
| solve time | 481 CPU-s total |
| check time (end-to-end verify, 2026-10-21) | wall **355 s**; cake_lpr 346 s; gate_lrat_check 42 s |

`manifest.json` lists every gate with its CNF SHA-256, LRAT SHA-256, step count and the
solver options used. The proofs are not in git, because each is between 13 and 872 MB.
* The .xz files live in `proofs/k7/artifacts/` on Lee's machine (gitignored).
* When a proof file is absent, verify.py regenerates it with CaDiCaL. A regenerated proof is
  checked the same way, and its hash is compared with the manifest for reproducibility.

**Non-gates.** These are also reproduced, and agree with Rosenfeld and with lr.py. For each
of the 7 primes 11, 13, 17, 19, 23, 29 and 41 the CNF is SAT. The decoded tuple was
re-checked as improper from the definition. The tuples are in `manifest.json` under
"non_gates", e.g. p=41: (1, 56, 66, 94, 122, 131, 150) mod 328.

## 3. Trusted base

| component | size | why trusted |
|---|---|---|
| cited lemmas: MSS Thm A, Rosenfeld Cor 3, Lemma 4, Lemma 6; LRC(6) (Barajas–Serra) | | peer-reviewed / published, not re-checked |
| `../ENCODING.md` §2 | 1 page | the model-construction argument, by hand |
| `common/gen_cnf.py` | 163 lines of Python (half comments) | generates the CNF the checkers read; validated by `common/test_encoding.py` (37 brute-force cases) |
| `common/mass.py` | 42 lines | exact integer inequality |
| **cake_lpr** | CakeML-verified LRAT checker | soundness proved in HOL4 down to the x64 machine code. Trusted: `basis_ffi.c` (the C FFI wrapper), gcc/the assembler, and the OS. Built from upstream commit a36874a; `cake_lpr.S` SHA-256 2f3af32d… matches upstream `cake_lpr.sha256`. Upstream's listed hash for basis_ffi.c is stale (the file changed on 22 Jul 2026) |
| gate_lrat_check.cpp | 140 lines of C++ | ours, an independent second check; either checker alone suffices, given its own trust |
| Python `lzma`, the OS | | they only move bytes, and a corrupted proof is rejected |

Not trusted: CaDiCaL, lr.py, the manifest, the hashes, and prove.py.

## 4. Reproduce from scratch (Linux)

    sh proofs/common/checkers/build_unix.sh      # builds gate_lrat_check + cake_lpr (pinned, hash-checked)
    python3 proofs/common/test_encoding.py       # optional: encoding vs brute force, 12 s (needs python-sat)
    CADICAL=/path/to/cadical python3 proofs/common/verify.py proofs/k7

Without stored proofs, add about 8 CPU-min of CaDiCaL to regenerate them.

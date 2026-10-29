# Lonely runner: the tight instances for six and seven speeds

**Status: computer-assisted; two independent implementations agree; hostile AI referee: SURVIVES; no human peer review.**

## Statement and sources

**Notation.**
* For positive integers $v_1,\dots,v_k$ (the speeds), let
  $\mathrm{ML}(v)=\max_t\min_i\|tv_i\|$, where $\|x\|$ is the distance to the nearest
  integer.
* The Lonely Runner Conjecture for $k$ speeds ($k+1$ runners), LRC($k$), says
  $\mathrm{ML}(v)\ge 1/(k+1)$.
* A speed set is *tight* if equality holds.
* Problem 1 of Perarnau–Serra's survey (*The Lonely Runner Conjecture turns 60*,
  Comput. Sci. Rev. 58 (2025) 100798, arXiv:2409.20160) asks for a characterisation of
  the tight instances.

**Theorem B1 (six speeds, seven runners).** If $\mathrm{ML}(v_1,\dots,v_6)\le 1/7$, then
$v$ is a dilate of a permutation of $(1,2,3,4,5,6)$.

**Theorem B2 (seven speeds, eight runners).** If $\mathrm{ML}(v_1,\dots,v_7)\le 1/8$, then
$v$ is a dilate of a permutation of one of
* $(1,2,3,4,5,6,7)$;
* $(1,2,3,4,5,7,12)$ (Wills);
* $(1,4,5,6,7,11,13)$ (Goddyn–Wong).

**Consequences.**
* **The tight instances.** For 7 and 8 runners, the tight instances are exactly these
  dilates. This answers Perarnau–Serra Problem 1 in those cases.
* **New proofs of LRC(6) and LRC(7).** B1 and LRC(5) give LRC(6). B2 and LRC(6) give
  LRC(7). In each case: if some $v$ had ML below the bound, the theorem would make it a
  dilate of a tuple whose ML equals the bound exactly.
* **Kravitz's Conjecture 8.3 for 7 and 8 speeds.** This is the "one very fast runner"
  form of his loneliness spectrum conjecture (N. Kravitz, *Barely lonely runners and very
  lonely runners*, Comb. Theory 1 (2021) #17, arXiv:1912.06034).
  * His Theorem 10.3 reduces the conjecture for $n$ speeds to two inputs: LRC($n-1$) and
    the full list of tight $(n-1)$-sets.
  * Each primitive tight set must then satisfy his condition (II.b).
  * For $n=7$, B1 supplies the list and LRC(6), and his Prop. 10.4 covers $(1,\dots,6)$.
  * For $n=8$, B2 supplies the list and LRC(7). Prop. 10.4 covers $(1,\dots,7)$, and an
    exact computation (`code/original/kravitz_IIb.py`) checks (II.b) for the other two
    sets.

## Prior work

* **L. Goddyn and E. B. Wong**, *Tight instances of the lonely runner*, Integers 6 (2006)
  A38. Constructs families of tight instances, including the two non-trivial septuples
  above.
* **Y. G. Chen** (J. Number Theory 37 (1991)) and **Bohman–Holzman–Kleitman** (EJC 8(2)
  (2001) R3) list the tight instances for 4 and 5 speeds.
* **N. Kravitz** (arXiv:1912.06034). The spectrum conjecture. Conjecture 8.3 was proved
  for 4 and 6 speeds (his Cor. 10.5, 10.7), with Thm 10.3 reducing it to the tight sets.
* **Y. Zhang**, arXiv:2608.13599 (v3, 18 Oct 2026). Classifies the one-entry
  modifications $([n-1]\setminus\{r\})\cup\{w\}$, and says that for at most five speeds
  the tight sets are known completely.
* **F. Cordella**, arXiv:2609.03444 (v2, 8 Sep 2026; code Zenodo 22345079). The
  six-speed spectrum up to finitely many exceptions. States that the seven-speed case is
  blocked by the missing classification of the tight sextuples. His search up to speed
  130 found only $(1,\dots,6)$.
* **G. Perarnau and O. Serra** (arXiv:2409.20160). The survey and Problem 1.
* **R. D. Malikiosis, F. Santos and M. Schymura** (MSS), *Linearly-exponential checking
  is enough for the lonely runner conjecture and some of its variants*, Forum Math. Sigma
  13 (2025) e164, arXiv:2411.06903. Bounds the size of counterexamples. Lemma T below is a
  strict form of their argument.
* **V. Bonfioli**, Zenodo 20534835 and 23445368. Work on the 4-speed spectrum and an
  8-runner LRC proof. Lists known tight vectors, and states that the list is not a
  classification.
* **Other proofs of LRC(6) and LRC(7).**
  * Barajas–Serra, EJC 15 (2008) R48 (LRC(6)), with H. Chen's repair of its §6.6,
    Zenodo 22852990.
  * Rosenfeld, arXiv:2509.14111 (LRC(7)).

A fresh novelty search on 2026-10-29 (`NOVELTY.md`: Zenodo, arXiv, GitHub, OSF, web)
found no classification of the tight instances for six or seven speeds. Please tell us if
we missed one.

## Method (full details in `PROOF.md`)

1. **Lemma T, a strict MSS bound.**
   * Rerunning the MSS proof with the scale $\lambda=\frac{n-2}{n}+\frac1m$ in place of
     $\frac{n-1}{n+1}$ gives: if LRC($n-1$) holds, $v$ is primitive and
     $\sum_{J\ne\emptyset}\gcd(v_J)>m^{n-1}$, then $\mathrm{ML}(v)\ge\frac1n-\frac1{2m}$.
   * $n=6$, $m=22$ gives $\frac{19}{132}>\frac17$. $n=7$, $m=29$ gives
     $\frac{51}{406}>\frac18$.
   * So a primitive $v$ with ML at most $\frac1{k+1}$ has $\sum v_i\le22^5$ for $k=6$ and
     $\sum v_i\le29^6$ for $k=7$.
   * `PROOF.md` §3 writes out the whole argument: the segment length $m(\lambda-\mu)=1$,
     and why the scale enters nowhere else.
2. **Prime gates.**
   * A prime $p$ is a *gate* if every residue tuple mod $(k+1)p$ that is covered at every
     time with the closed threshold $\le\frac1{k+1}$, and not excluded by an elementary
     gcd lemma (Lemma G), reduces mod $p$ to a unit/sign/permutation image of a known
     tight tuple.
   * Exhaustive enumeration shows the following.

     | k | primes | non-gates | gates |
     |---|---|---|---|
     | 6 | $11\le p\le997$ | 17, 19, 23, 29, 31, 37, 47, 53, 61, 79 | 154 |
     | 7 | $11\le p\le1153$ | 11, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 67 | 175 |

   * Level-one search is reduced by a proved private-time normalisation (Lemma N in
     `PROOF.md`).
3. **Invariant divisibility.**
   * At a gate not dividing $\prod v$, the vector $v$ matches some tight $T$ mod $p$.
   * Then $p$ divides every $I^T_j(v)=e_j(v^2)E_1(T)^j-e_1(v^2)^jE_j(T)$.
   * These invariants all vanish only when $v$ is a permutation of $T$.
   * Bounding $\prod v$ and $|I^T_j|$ by the size bound from step 1, the product of the
     gates exceeds what a non-tight $v$ allows. This is checked exactly in integers.

     | k | needed (nats) | gate mass (nats) | margin (nats) |
     |---|---|---|---|
     | 6 | 283.80 | 915.27 (primes ≤ 997) | **631.47** |
     | 7 | 1050.94 | 1070.19 (primes ≤ 1153) | **19.24** |

   * With primes only up to 1129, the $k=7$ margin would be 5.14 nats, and every gate
     from 173 to 1129 would be individually load-bearing.
   * B2 uses primes up to 1153. Then any two gate verdicts could be wrong without
     affecting the theorem.

## How to verify

**Requirements.**
* Python 3 with `numpy`, `numba` and `mpmath`. `sympy` is not needed.
* Pure CPU, one thread, peak memory well under 1 GB.

**Commands.**

```
cd code
bash original/check.sh        # implementation 1: sieve k=6 and k=7, exact bookkeeping, compare with gates/*.json
bash referee/check.sh         # implementation 2 (independent code), same checks
bash run_all.sh               # both, plus all brute-force validations; transcripts in code/logs/
```

Set `PY=/path/to/python` if `python` is not the right interpreter.

**Expected output.** `original/check.sh` ends with:

```
k=6: m=22, S0=m^5=5153632, needed gate mass log(B) = 283.799574 nats
  primes 11..997: 154 gates, mass 915.272049 nats, margin 631.472475 nats, exact inequality holds: True
k=7: m=29, S0=m^6=594823321, needed gate mass log(B) = 1050.943096 nats
  primes 11..1129: 173 gates, mass 1056.086713 nats, margin 5.143617 nats, exact inequality holds: True
  primes 11..1153: 175 gates, mass 1070.185222 nats, margin 19.242126 nats, exact inequality holds: True
  k=6: 164 published primes, verdict/orbit-count mismatches: [], unpublished primes: []
  k=7: 187 published primes, verdict/orbit-count mismatches: [], unpublished primes: []
ORIGINAL CHECK OK
```

`referee/check.sh` ends with the same margins and `REFEREE CHECK OK`.

**Runtimes** (recorded run, 2026-10-29, one thread):
* `original/check.sh`: about 55 min (k=6 sieve 1.8 min, k=7 sieve 53.6 min).
* `referee/check.sh`: about 58 min (k=6 sieve 1.6 min, k=7 sieve 56.6 min).
* `run_all.sh`, including validations: 2 h 3 min (63.6 + 59.5 min), with peak commit 0.64 GB.

The transcripts of that run are in `code/logs/`. `code/CHECKSUMS.sha256` lists every file
in `code/`.

**What the validations check** (`--validate`, run by `run_all.sh`):
* `tight6.py validate`/`validate2` and `validate7.py`: brute force over every multiset,
  every lift and every time at small primes.
* `ref_validate.py`: the pruned level-one search against unpruned enumeration, at all
  primes up to 97 ($k=6$) and up to 61 ($k=7$), and the bitset lift against a naive count.
* `ref_ml.py`: exact ML of the tight tuples and of controls.
* `kravitz_IIb.py`: Kravitz's condition (II.b).

**Files.**
* `code/gates/gates_k6.json`, `code/gates/gates_k7.json`: for each prime, the verdict,
  the level-one orbit count, the non-exempt orbit count and the survivor counts from both
  implementations.
* `code/make_gates.py`: rebuilds these lists from the two runs and refuses to write them
  if the runs disagree.

## Limits and dependencies

**Dependencies.**
* **B1** uses LRC(5) (Bohman–Holzman–Kleitman 2001; also Renault 2004), in Lemma T and in
  Lemma G.
* **B2** uses LRC(6), in Lemma T and in Lemma G. Either published proof works:
  * Barajas–Serra 2008 with H. Chen's repair (Zenodo 22852990);
  * Corollary 1 here (B1 + LRC(5)). B2 then also depends on the B1 computation.
* **MSS** (arXiv:2411.06903). Lemma T reruns their §3 argument with a different scale.
  The full proof is in `PROOF.md` §3. It has not been refereed by a journal.
* **The Kravitz corollary** uses his Theorem 10.3 and Prop. 10.4. The numbering follows
  arXiv:1912.06034 v1; check it against the published version. For 8 speeds it needs
  LRC(7), from B2 + LRC(6) or from Rosenfeld.

**Limits.**
* The theorems are computer-assisted. They rest on the finite gate verdicts, computed by
  two independent exhaustive enumerations that agree at every prime, plus brute-force
  cross-checks at small primes.
* There is no formal proof certificate. SAT/LRAT certificates were produced and checked
  only for a few small six-speed gates ($p=13,41,43,59$). Whole-gate certificates for
  all gates would run to terabytes.
* The statements are for integer speeds.
* No human has refereed this work.

## Referee summary

A hostile AI referee re-derived Lemma T from the MSS source and found it correct. It
re-proved the normalisation and checked Lemma G, the gate reduction and the invariant
bookkeeping. It also wrote the second implementation from the definitions; that code
imports nothing from the first.

Both implementations agree on every verdict and every level-one orbit count:
* $k=6$: 164 primes;
* $k=7$: 187 primes, up to 1153.

The referee's required changes are all applied here:
* the full proof of Lemma T;
* B2 stated for ML ≤ 1/8, so that it supplies LRC(7) for Kravitz's Theorem 10.3;
* primes up to 1153;
* the dependencies made explicit;
* a written proof of the normalisation;
* the gate lists and checkers.

**Verdicts:** B1 SURVIVES; B2 SURVIVES; Kravitz 8.3 for 7 and 8 speeds SURVIVES.

## Files

* `PROOF.md`: the full argument.
* `NOVELTY.md`: the novelty search, with queries and timestamps.
* `code/original/`: implementation 1 (`lr.py`, `tight6.py`, `tight7.py`, checkers,
  validations, Kravitz scripts).
* `code/referee/`: implementation 2 (`ref_level1.py`, `ref_lift.py`, `ref_gates.py`,
  `ref_bounds.py`, `ref_ml.py`, `ref_validate.py`, checker).
* `code/gates/`: the gate lists. `code/logs/`: the recorded run. `code/run_all.sh`,
  `code/CHECKSUMS.sha256`.
* `zenodo_note/`: the LaTeX source of the short preprint and the Zenodo metadata.

## AI disclosure

This work was generated with AI (Claude, from Anthropic), directed by Lee Marshall. It
was checked by two independent implementations and attacked by a hostile AI referee. It
has not been through human peer review. Please check the computations yourself with the
checkers provided.

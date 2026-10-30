# The geometric thickness of K14 is 4

**Status: computer-assisted. Two independently written SAT encoders, each with its own case tree, give UNSAT. Every leaf LRAT proof was checked by cake_lpr (CakeML-verified). Hostile AI referee: SURVIVES. No human peer review.**

## Statement and sources

**Definition.** The geometric thickness θ̄(G) is the least t such that G has a straight-line drawing (vertices at
distinct points in the plane, edges as segments) whose edges can be split into t classes, each crossing-free.
Dillencourt–Eppstein–Hirschberg (DEH), *Geometric thickness of complete graphs*, J. Graph Algorithms Appl. 4(3)
(2000) 5–17, arXiv:math/9910185.

**Theorem.** θ̄(K14) = 4.
* **Lower bound (new):** no straight-line drawing of K14 has its 91 edges split into 3 crossing-free layers.
* **Upper bound (known):** θ̄(K_n) ≤ ⌈n/4⌉ (DEH Thm 2.1), so θ̄(K14) ≤ 4.
  * Explicit exact witness: `certs/witnesses/k14_t4_referee.json`. It is Seamons's 3-layer K13 drawing plus the point
    (1000003, 777781), whose 13 edges (a star, hence crossing-free) form layer 4. Integer coordinates, checked exactly.
  * Four further exact 4-layer K14 drawings, realised from SAT models with hull sizes 3–6:
    `certs/witnesses/k14_t4*.json`.

**Consequence.** With DEH (n ≤ 12, and n = 15, 16) and Seamons (n = 13), θ̄(K_n) is known for every n ≤ 16:
1 (n ≤ 4), 2 (5–8), 3 (9–13), 4 (14–16).

## Prior work

* **DEH 2000**: θ̄(K_n) ≤ ⌈n/4⌉ (Thm 2.1); a counting lower bound (Thm 3.1); Table 1 gives θ̄(K13), θ̄(K14) ∈ {3, 4}
  and θ̄(K15) = θ̄(K16) = 4. Their open problem 1 asks for K13 and K14. The K15 lower bound uses triangulation counting
  plus Lemmas 4.2/4.3.
* **J. W. G. Seamons**, Zenodo record doi:10.5281/zenodo.22239253, 2 Sep 2026: θ̄(K13) = 3 by an integer witness. The record says DEH left K13 and K14 open;
  it does not settle K14. Its `k13.json` is included as `certs/witnesses/k13_seamons.json` (SHA-256 `8e36330a…1ef50`,
  equal to the record's `SHA256SUMS.txt`).
* **Method sources**: signotope encoding of point sets, the p1 normal form and the reflection predicate follow
  Heule–Scheucher (TACAS 2024). Signotope axioms: Felsner–Weil.
* **Novelty search**, 2026-10-29: Zenodo API ("geometric thickness", newest first), arXiv API (thickness ∧ complete in
  math.CO, cs.CG, cs.DM), GitHub, OSF, the openai/math CONTENTS.md and the web. Recent related work does not touch K14
  (arXiv 2503.03577, 2312.05010, 2601.01734). No earlier claim of θ̄(K14) = 4 was found. Please tell us if we missed one.

## Method

1. **Points to signotopes.**
   * A straight-line K_n drawing has no 3 collinear vertices (the middle one would lie on an edge). A small rotation
     makes x-coordinates distinct and changes no crossing.
   * Label points 1..n by x. Variable X(a,b,c), a<b<c, means "counter-clockwise".
   * Rank-3 signotope axioms on every 4-set. UNSAT over signotopes implies UNSAT over point sets.
2. **Layers and crossings.**
   * C(e,k): edge e in layer k, exactly one layer per edge.
   * Y(pq,rs) ↔ [o(pqr) ≠ o(pqs)] ∧ [o(rsp) ≠ o(rsq)], the exact crossing test in general position.
     Y → pq and rs not in the same layer.
3. **Symmetry breaking** (each keeps at least one image of every real drawing).
   * p1 normal form: X(1,a,b) = + for all 1<a<b (p1 a hull vertex, others in angular order; projective map).
   * Reflection lex predicate (authors' encoder only).
   * Layer value-precedence.
4. **Euler bound.**
   * H_v ↔ v lies in no triangle (v is a hull vertex). U_e ↔ e crosses no edge.
   * A crossing-free straight-line graph on n points with h hull points has ≤ 3n−3−h edges, and an uncrossed edge
     can be added to any layer. So |L_k ∪ U| + h ≤ 3n−3 for each k.
   * Summed over the t layers: (t−1)·|U| + t·h ≤ t(3n−3) − C(n,2), i.e. 2|U| + 3h ≤ 26 for K14 with t = 3.
5. **Hull-size case split.**
   * The normal form makes 1, 2 and n hull vertices. S = other hull vertices among 3..13.
   * Authors' tree: |S| = 0, 1, 2 exactly (1 + 11 + 55 cases) and |S| ≥ 3 by the 165 cases "T ⊆ S", |T| = 3. 232 cases.
   * Referee's tree: a binary tree over H_3, H_4, … that stops at the third positive literal. 232 roots, checked by
     brute force over all 2^11 patterns to be a partition.
6. **Sweep lemmas** (authors' encoder, |S| ≤ 2 cases only; redundant).
   * For each hull edge BC, a linear functional d with d(B) < d(C) < d(others).
   * Interior point of a triangle lies below its max vertex; if pq crosses xy then max d(p,q) > min d(x,y) and
     vice versa. Generalises DEH Lemmas 4.2/4.3.
   * **A second run without them** (`TK_NOSWEEP=1`) is also UNSAT on all 232 cases, so they carry no weight.
7. **Cube-and-conquer.** Each node solved by CaDiCaL 3.0.1 with LRAT, 300 s limit. On timeout, split on the next 4
   free orientation variables (16 children). Each UNSAT leaf's LRAT checked by cake_lpr, hashed, then deleted.

## Verification

| run | encoder | symmetry breaking | sweeps | leaves | CaDiCaL CPU-s (leaves) | LRAT total |
|---|---|---|---|---|---|---|
| `tree_authors_sweep` | `gen.py` | p1, lex, precedence | yes | 307 | 2093 | 16.6 GB |
| `tree_authors_nosweep` | `gen.py` | p1, lex, precedence | no | 307 | 2429 | 18.1 GB |
| `tree_referee` | `renc.py` | p1, precedence | no | 352 | 3056 | 21.0 GB |

* **Two independent encoders.** `renc.py` was written by the referee from the definitions before reading `gen.py`.
  It has no lex predicate and no sweeps, keeps all 3 pairings per 4-set, uses its own counters and its own cube
  partition (`cnc.py`).
* **Every leaf LRAT** of all three runs: `s VERIFIED UNSAT` from cake_lpr. 26 authors' leaves were also re-checked
  end to end by the referee.
* **LRATs are not shipped**: about 35 GB for the two authors' runs and 21 GB for the referee's. `certs/tree_*.jsonl`
  record, per leaf, the SHA-256 of the CNF and of the LRAT that was checked, its size and the solve time. `FULL=1`
  regenerates them.
* **Controls.**
  * Referee encoder: K4/1 SAT, K5/1 UNSAT, K8/2 SAT, K9/2 UNSAT, K12/3 SAT, K13/3 SAT (DEH table, Seamons).
  * Seamons's K13 and 40 random exact drawings (n = 9–12) satisfy the referee CNF after normalisation.
  * Seamons's K13, the five K14 4-layer drawings and three exact K10 drawings (hull 4, 5, 6) satisfy the authors'
    production CNF of their own hull case.
* **Tree hygiene.** The authors' original tree file had 1417 unreachable UNSAT records and 2 unreachable split records
  under 135 spurious splits (CaDiCaL 3.0.1 rc-3 aborts, since removed). The published `tree_authors_sweep.jsonl` keeps
  only the 312 reachable records (5 splits, 307 leaves).

## How to verify

**Requirements.**
* Python 3 with python-sat (`pip install python-sat`) for the controls. The cover checks and `replay.py` use the
  standard library only.
* For `FULL=1`: CaDiCaL (tested: 3.0.1) and an LRAT checker (cake_lpr).

**Default check** (about 1 min, under 1 GB RAM, no solver binary needed):

```
cd code && PY=python3 bash check_all.sh
```

* Certificate hashes (`certs/CHECKSUMS.sha256`).
* Upper bound: rebuilds and exactly checks the 4-layer drawing; checks the other four.
* Soundness controls for both encoders (the K13/3 SAT control adds about 5 min: `CONTROLS_K13=1`).
* Completeness of all three trees: every case/root present, every split has all 16 children, every leaf UNSAT with
  a verified flag; the referee's root partition by brute force; **every leaf CNF regenerated** and compared with
  its published SHA-256.

**Full re-certification:**

```
cd code && FULL=1 CADICAL=cadical LRAT_CHECK="cake_lpr --CML_HEAP_SIZE=2500 --CML_STACK_SIZE=1000" bash check_all.sh
```

* For each of the 966 reachable leaves: regenerate the CNF, compare its hash, solve with
  `cadical --lrat=true --binary=false --congruence=false`, check the LRAT, re-hash the LRAT the checker read, delete.
* Then the cover check is re-run against the replay manifest (`work/replay_*.jsonl`): every reachable leaf must have
  a hash-matched, verified replay record.
* The LRAT hash is compared with the published one and the agreement is reported. A mismatch is not a failure:
  LRAT output can differ across CaDiCaL builds and options. Observed on 2026-10-29 (same 3.0.1 build): byte-identical
  for the sampled `nosweep` and `referee` leaves; different for the sampled `sweep` leaf (65 MB published, 74 MB now).
* `CADICAL_OPTS` defaults to `--congruence=false`; CaDiCaL 3.0.1 aborts with rc 3 in congruence closure on some
  leaves without it.
* Resumable: rerunning skips leaves already verified. `TREES=referee` selects trees; `SAMPLE=k SEED=s` replays k
  random leaves per tree.
* **Cost** (one thread, estimate): about 3.5 hours for all three trees.
  * CaDiCaL: 2.1 CPU-h (sum of the recorded leaf times: 2093 + 2429 + 3056 s).
  * cake_lpr: about 1.2 h for 56 GB of LRAT, at the 10–16 MB/s measured in the smoke test.
* **Disk**: one leaf at a time; largest LRAT 1.7 GB, CNF ≤ 3 MB. Keep about 4 GB free in `work/` if your checker
  wrapper copies the files. **RAM**: cake_lpr with a 2.5 GB heap.

**Recorded runs** (2026-10-29, one thread, from a clean copy):
* Default `check_all.sh`: 56 s, all checks pass.
* `FULL=1 SAMPLE=1 SEED=1` (one leaf per tree, 74–256 MB LRATs): 2.5 min; 3/3 CNF hashes match, UNSAT, `s VERIFIED UNSAT`.

## Limits

**What the check covers.**
* Given the two CNF generators: every case of both case trees is UNSAT with a verified LRAT, and the trees are
  complete partitions.
* The leaf CNFs that were certified are exactly the ones the shipped encoders generate (hash comparison).

**What it does not cover.**
* The soundness of each encoding rests on hand-proved lemmas and the referee's audit, not on a formal proof:
  * signotope axioms for x-sorted point sets; the crossing test; the p1 normal form (projective map) and, for
    `gen.py`, the reflection predicate;
  * the Euler/uncrossed-edge bound and the hull indicators;
  * completeness of the hull case split; the sweep lemmas (not load-bearing: the no-sweep run is also UNSAT).
* These are cross-checked empirically only: real drawings satisfy both encodings, and the small-n table is reproduced.
* There is no formal proof that `gen.py` or `renc.py` encodes the definitions faithfully. Two independent
  implementations agree.
* There is no human-readable proof of the SAT step.
* No human has refereed this work.

## Files

* `code/gen.py`, `code/run_all.py`: authors' encoder, hull cases and cube-and-conquer driver.
* `code/renc.py`, `code/cnc.py`: referee's independent encoder and driver.
* `code/cover_check.py`, `code/cover_mine.py`: completeness checks (with `--regen`, `--replay`).
* `code/replay.py`: full re-certification of a tree. `code/tools.py`: solver and checker calls.
* `code/controls.py`, `code/check_known_prod.py`: soundness controls. `code/upper.py`: upper-bound witness and its check.
* `code/check_all.sh`: one-command check.
* `certs/tree_*.jsonl`: case trees with CNF and LRAT hashes. `certs/witnesses/`: exact drawings.
  `certs/CHECKSUMS.sha256`.
* `zenodo_note/`: LaTeX source of the technical note and its metadata.

## AI disclosure

This work was generated with AI (Claude, from Anthropic), directed by Lee Marshall. It was checked by two independent
encoders and a formally verified LRAT checker, and attacked by a hostile AI referee. It has not been through human peer
review. Please check the computations yourself with the checkers provided.

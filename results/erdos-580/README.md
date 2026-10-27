# Erdős problem #580 (vertex form) holds for n <= 21

## Claim

**Statement (vertex form).** Let G be a graph on n vertices in which at least ⌈n/2⌉ vertices have degree at least ⌈n/2⌉. Then G contains every tree with at most ⌊n/2⌋ vertices as a subgraph.

**Result.** This holds for all n <= 21. Our own checks cover n <= 18, which reproduces the earlier work of Zeraoulia. The cases **n = 20 and n = 21 are new**. n = 19 is covered by n = 18 plus the standard step below, and n = 21 by n = 20.

The case n = 20 (m = 10) is settled as follows. There are 106 trees on 10 vertices.
- 75 are covered by cited theorems.
- 31 are refuted as counterexample candidates by SAT, with LRAT certificates checked by cake_lpr (a formally verified LRAT checker) and, as a second check, `gate_lrat_check`.

A hostile referee pass (`notes/REFEREE-n20.md`) concluded SURVIVES, with only documentation gaps and none affecting soundness.

## Prior work

- Loebl–Komlós–Sós conjecture (the edge form of the problem). Zhao proved it for large n. The vertex form is the version stated on erdosproblems.com as #580.
- Zeraoulia, Zenodo record 21348157 and a forum post of 2026-07-16, claims n <= 19. Nothing beyond n = 19 was found in the forum thread, on the web or on GitHub (as of 2026-10-26).
- Piguet–Stein (arXiv 0712.3382): the LKS conjecture for trees of diameter <= 5 (Thm 2) and for the family T(k,l,c) (Thm 3). Bazgan–Li–Woźniak: paths and a path with a star.

## Method

1. **Edge-minimal host reduction.** A minimal counterexample on n = 2m vertices can be taken with L (the m high-degree vertices) of exact degree m and S (the other m vertices) independent. n odd follows from n - 1 by deleting a vertex. A tree with fewer than m vertices embeds in an m-vertex tree.
2. **Leaf-extension (Hall) lemma.** It lets us require only that the tree's core (the tree minus its leaves) embed with the leaf-parents mapped into L. This makes the encoding sufficient but not necessary, which is fine for refutation.
3. **Tree classification** (`classify.py`, output in `classification_m*.json`). The trees on m vertices are enumerated and sorted by the cited results. At m = 10 the split is D 58 (diameter <= 5, Piguet–Stein Thm 2), P 11 (path or path plus star, Bazgan–Li–Woźniak), C 6 (T(k,l,c) with c even, Piguet–Stein Thm 3) and E 31 (needs SAT). 58 + 11 + 6 = 75 trees are cited.
4. **SAT + LRAT + cake_lpr.** For each E tree, `lks.py` writes a CNF whose models are hosts avoiding the tree. CaDiCaL proves it UNSAT and emits an LRAT proof. `gate_lrat_check` and cake_lpr both check the proof. The hardest cases are the spiders S(3,3,3) and S(2,3,4). S(2,3,4) has 14.8M clauses and a 650 MB LRAT proof. An independent audit (`audit.py`) checks every generated CNF clause by clause against the intended encoding.

Trusted base: cake_lpr, `lks.py`, `classify.py` with the cited theorems, and the two lemmas above. See `REPORT.md` for details, including a harmless misquotation of Piguet–Stein in Zeraoulia's argument.

## Files

- `REPORT.md`: working report (the paths to local job logs were removed).
- `notes/REFEREE-n20.md`: hostile referee report for n = 20.
- `lks.py`, `classify.py`, `audit.py`, `audit_cnf.py`, `run_certs.py`: encoding, classification, audit, certificate pipeline. `gen_eager.py`, `eager.py`, `cegar.py`, `hybrid.py`, `find_variant.py`, `plan_variants.py`: alternative encodings and variant generators used for the hard trees.
- `classification_m7.json` to `classification_m10.json`: tree classification.
- `results/m4.jsonl` to `results/m10.jsonl`: one row per certificate run (hashes, solver and checker return codes, `verified`). At m = 10 the 31 E trees each have a final row with `verified: true`. Earlier rows for the same tree are failed or alternative attempts (for example cake_lpr heap exhaustion).
- `manifest.json`: CNF and LRAT hashes and sizes for the certificates.
- `diff_zenodo.py`, `check_zbase.py`, `check_zenodo.sh`: audit of Zeraoulia's certificates. They need his Zenodo files in `zenodo/`, which are not included.
- `cake_check.sh`, `d73_tree45.log`: the cake_lpr wrapper, and the log of the final S(2,3,4) check (`s VERIFIED UNSAT`).

The CNF and LRAT files are **not included**. They are multi-GB and cheap to regenerate.

## How to verify

Needs Python 3, CaDiCaL, `gate_lrat_check` and `cake_lpr`. Tool locations are set through the environment variables `CADICAL`, `GATE_LRAT_CHECK` and `CAKE_LPR`. cake_lpr runs in WSL or Linux, and `~/lrchk/bin/` is the default location of the checkers. If you are not on Windows with WSL, adapt the `wsl` calls in `run_certs.py`.

```
python verify.py --quick   # re-classifies m = 4..10, then regenerates and re-checks every certificate that solves in < 30 s
python verify.py           # all certificates; about 6 CPU-hours, tens of GB of disk, up to 10 GB RAM for the largest cake_lpr runs
```

To regenerate one certificate, for example m = 10, tree 45: `python run_certs.py 10:45`. Set `KEEP_LRAT=1` to keep the proof, and `CAKE_HEAP_MB` and `CAKE_STACK_MB` for the large cases. The Zenodo audit is skipped when `zenodo/` is absent. `python classify.py 10` alone re-derives the 75 / 31 split in seconds.

## Status and level

Verified by SAT with formally checked certificates, plus cited theorems, and survived one AI hostile referee. No human peer review. The earlier n <= 19 claim is unrefereed, and we reproduced n <= 18 independently.

**Calibrated level: 2.5 (range 2 to 3).** It is a finite extension of an existing method, but it is a first for n = 20 and 21 and uses a verified checker.

No post or contact about this result has been made on erdosproblems.com or elsewhere.

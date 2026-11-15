# No set of type (2,5) with respect to lines in PG(3,9)

**Status: computer-assisted proof, not peer reviewed. Reviewed adversarially by an AI referee. An independent human check is invited.**

## Statement

**Theorem.** There is no set of points of PG(3,9) that meets every line in exactly 2 or exactly 5 points. Equivalently there is no set of type (2,5) and none of type (5,8) in PG(3,9).

**Corollary (depends on cited theorems).** For every prime power q < 25 and every r >= 3, a set of points of PG(r,q) that meets every line in exactly a or b points (a < b, both occurring) is a point, a hyperplane, or the complement of one. This uses Calderbank-Kantor, Bull. London Math. Soc. 18 (1986), Theorem 12.5 (attributed there to Tallini Scafati) and the type (0,n) result as cited by Napolitano (Contrib. Discrete Math. 3 (2008), Cor. 4.4). We have read Calderbank-Kantor and Napolitano; we have not read Tallini Scafati's own papers.

## Source of the problem and prior work

The question is stated as open in Gavrilyuk-Matkin (arXiv:1803.10442, 2018: first open case PG(3,9), (a,b) = (2,5)), Kurz-Landjev-Pavese-Rousseva (arXiv:2211.16793, 2022: no such set known for r > 2) and Ihringer (arXiv:2411.16288, Section 5.2: first open case, open for more than five decades). The planar sections of a hypothetical set are the (2,5)-sets of PG(2,9), classified by computer by Penttila-Royle (Des. Codes Cryptogr. 6 (1995) 229-245). We did not have access to their paper and have **not** compared our planar classification (3 classes of size 26, 7 of size 35) with theirs. Searches are recorded in `LITERATURE-CHECK.md`; they are repeated on the day of release.

## Method

A hypothetical set has 260 or 287 points. Fixing a point of it, the directions of the 2-secants through the point form a (2,5)-set of PG(2,9); by a certified classification there are 10 such sets up to collineation. Normalising by a collineation and an elation gives 10 SAT cases. In each, the section on one plane is one of 16 or 18 sections up to homothety (a certified completeness proof), giving 174 cases, each refuted by SAT with an LRAT proof checked by the CakeML-verified checker `cake_lpr`; 7 of them could not be closed directly and are split once more; two more, C2_c00 and C1_c04, were split although they had closed directly (9 splits, 1,224 children in all). Details and the status of every step (hand proof, certified, script check, cited) are in `PROOF-CHAIN.md`; the semantics of the encoder are proved in `ENCODER-SEMANTICS.md`; the note is `note/pg39_note.pdf`.

## How to check it

Everything under `replay/` uses only the Python standard library. From the folder `pg39/`:

```
python replay/rebuild_all.py --selftest   # rebuilds all 1,424 formulas, compares SHA-256 with certs/manifest.jsonl (about 40 s)
python replay/cover_asserts.py --lrat     # cover logic (about 40 s); needs the 19 certs/cover* formulas and proofs from Zenodo in certs/
python replay/planar_validity.py          # planar solution lists (about 6 s)
python replay/planar_count.py
python replay/planar_aut.py
python replay/lemma_checks.py
python replay/level0_arith.py
python replay/order3_check.py
python replay/encoder_exhaustive.py       # about 2.5 min (--quick for less)
```

Together under an hour on a laptop (one thread, under 4 GB). To re-check the stored certificates (19 completeness and 4 planar proofs; download from Zenodo; Linux binary of `cake_lpr`): `replay/cake_recheck_local.sh` (about 4 minutes, heap 2 GB; log `logs/cake_recheck_local.log`: 22 of 22 `s VERIFIED UNSAT` on 11 November 2026, before the 19th completeness proof existed; `logs/cake_recheck_c1c04cover.log`: the 19th, `cover2_C1_c04`, `s VERIFIED UNSAT`).

The 167 cube proofs and 1,084 child proofs of the cloud run (571 GB) and the 140 child proofs of the local split of C1_c04 (33.3 GB) were **not retained**. They are regenerated deterministically: `replay/rebuild_all.py` rebuilds every formula in memory and compares its hash with the manifest (it does not write the formula); `code/a_slice_driver.py` writes the formulas from a job list (`results/full_rerun_stage{1,2,3}.json`, `results/c1c04_jobs.json`), runs `cadical -t 7200 --lrat=true --binary=true --congruence=false`, then `cake_lpr`, and compares with the recorded hashes. Total about 125 CPU-hours.

All of them have been regenerated locally except one: the full local regeneration finished on 2026-11-13 (`logs/full_rerun_chain.log`, records `results/full_rerun_dry.jsonl`, comparison `code/full_vs_local.py` -> `results/full_vs_local.json`). **1,250 of the 1,251 cloud certificates (166 of 167 cubes, 1,084 of 1,084 children) were re-run locally and are bit-identical to the cloud records** (status, conflicts, proof size, proof and formula SHA-256), none timed out or failed the checker, 77.5 CPU-hours (`certs/PROOF-INDEX.tsv`, column `local_identical`). The one exception is the monolithic cube C1_c04: its 4.96 GB proof needs a 6500 MB checker heap, more than the local machine allows. That cube is closed locally instead through a second split (`certs/split_C1_c04.json`: plane 10, 140 children, completeness proof `certs/cover2_C1_c04.{cnf,lrat}`); all 140 children and the completeness proof were solved and checked locally (`results/c1c04_manifest_local.jsonl`, 140 of 140 `s VERIFIED UNSAT`, 4.2 CPU-hours). Every one of the 174 cubes therefore has a refutation verified on the local machine.

## Files

| path | contents |
|---|---|
| `note/pg39_note.pdf`, `.tex` | technical note |
| `PROOF-CHAIN.md` | the argument, with the status of each step |
| `ENCODER-SEMANTICS.md` | proof of the encoding |
| `LITERATURE-CHECK.md` | literature, dependencies, order-3 cross-check |
| `replay/` | independent stdlib scripts |
| `code/` | generators and drivers used to produce the certificates |
| `certs/` | `manifest.jsonl`, `PROOF-INDEX.tsv`, `cubes_*.json`, `split_*.json` (cube and child lists) |
| `results/` | planar classes and frame solutions, ledger, parameter table |
| `logs/` | run logs |
| Zenodo | formulas (CNF) and stored LRAT proofs: https://doi.org/10.5281/zenodo.24259948 (4 planar proofs, 19 cover proofs with formulas, the small formulas; put them in `certs/`) |

## Hashes

| item | sha256 |
|---|---|
| `certs/manifest.jsonl` | 3329644a1db61827d40933d493018b01aab2a937793c7b672ac9d9bc2e41f1ae (recomputed 2026-11-13 after the C1_c04 split; recompute on release) |
| `certs/PROOF-INDEX.tsv` | c284c39e13cded5005b0f2b14699a13f1f7914d2d8f8319a55be04ff5dc2e32e |
| `results/planar_classes.json` | b258b2ac933009df57fb4e18d7ae312f35e2a9788c3436079ee1a515fb9bd14f |
| CaDiCaL 3.0.1 binary | f30ecc0fc740ea7ab8bb928dc4c0b3e065faa8cfcfad07eb1da59c4e08d502a5 |
| `cake_lpr` binary | 272a75de5b8add99371f7098a1af076e7c419c9eba2e7f2915e67a530671acaf |

## Limitations

- `cake_lpr` and the CakeML/HOL4 proof of the checker are cited, not checked. CaDiCaL is not trusted.
- The refutation proofs were deleted after checking. What a reader trusts, unless they re-run them, are the drivers' records (formula and proof hashes, the checker's last line). For 1,250 of the 1,251 cloud certificates the local run agrees bit for bit with the cloud run; the remaining cube is closed by a split that was verified locally only (no cloud record of its 140 children).
- The hand lemmas are not formalised; they are proved in the text and machine-checked on the actual geometry.
- Encoder exactness at full size rests on a proof plus exhaustive small cases and sampling.
- The planar classification has not been compared with Penttila-Royle.
- The corollary depends on Tallini Scafati's theorem as stated in Calderbank-Kantor.

## How to cite

Lee Marshall, *No set of type (2,5) with respect to lines in PG(3,9)*, 2026, folder `results/pg39-type25` of https://github.com/LeeMarshall1113/ai-math-results; certificates: https://doi.org/10.5281/zenodo.24259948 (v1).

## AI statement

The search, encodings, certificates and drafts were produced by AI assistants (Claude, Anthropic) under the direction of Lee Marshall, with an adversarial review by a second AI system. None of it has been through human peer review. Please check the certificates yourself, and please report errors or earlier work we missed (open an issue).

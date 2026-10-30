# CHECKLIST: torus-quad10 (copied from math-comms/PIPELINE.md, 2026-10-29)

## 0. Before compute
- [x] Significance check. Kept in the private scout notes, not here, because ratings stay private.
- [x] Novelty pre-check on Zenodo API, arXiv, web, GitHub, OSF and openai/math CONTENTS.md (`NOVELTY.md`).
- [x] Selection score. Scout file, item 1.
- [ ] Row in GATES.md. The coordinator handles this.
- [x] Compute plan: negligible (under 2 min, under 1 GB), so no lease was needed.
- [x] Compute is conservative: no paid compute.
- [n/a] Open hand-off.

## 1. The result
- [x] Precise statement, with the GS 2009 p. 21 definitions quoted and Conjecture 3.2 checked against the source.
- [x] PROOF.md, with every computer step named.
- [x] Certificates with SHA-256 checksums (`certs/CHECKSUMS.sha256`).

## 2. Verification
- [x] Independent checker: `encode2.py` (written separately) plus cake_lpr (verified) and a second LRAT checker.
  `geomcheck.py` is an exact geometric checker that fails closed.
- [x] Coverage assertion:
  - the map is unique by two enumerations;
  - all 3150 GP relations are included;
  - all 40 ordered adjacent pairs are included.
- [x] Mutation test (`code/mutation_test.py`).
- [x] End-to-end replay with one command (`code/check_all.sh`, about 1.5 min). Log: `certs/replay_2026-10-29.log`.
  Peak memory under 1 GB, not separately recorded.
- [x] Hostile referee: SURVIVES (erdos-triage/torus-quad10/referee/REFEREE.md). Fixes applied.
- [x] Verified checker: cake_lpr on all three 10-face LRATs and both 8-face LRATs. There is no Lean proof of the generators or of Lemmas 1–2.

## 3. Novelty and context
- [ ] Fresh novelty re-check on the day of publishing.
- [ ] Simutis 1977 thesis (UC Davis) obtained and checked. Not found online on 2026-10-29.
- [ ] Level rating rechecked (private).

## 4. Ship (Lee decides; nothing pushed)
- [ ] Zenodo fast-stake note. Lee fills the `% LEE` sections of `zenodo_note/note.tex` and uploads it himself.
- [ ] Public folder in ai-math-results. Staged here only. Before publishing, remove `code/tools/cake_wsl.sh` (machine-specific) and strip local paths from `certs/replay_2026-10-29.log`.
- [ ] Zenodo DOI via the GitHub integration.
- [ ] LEDGER.md, BOARD.md and GATES.md updated.
- [ ] Same-day courtesy email draft for Lee. Lee's rule: no contact with Seamons until Lee decides.

## 5. After
- [ ] Follow-up date in LEDGER.md (about 2 weeks).
- [ ] Paper only at level >= 3, ideally with an expert co-author.

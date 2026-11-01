# CHECKLIST: hadwiger-h3 (from math-comms/PIPELINE.md, 2026-11-01)

## 0. Before compute
- [x] Significance check. Kept in the private scout notes, not here, because ratings stay private.
- [x] Novelty pre-check (`NOVELTY.md`).
- [x] Compute plan: certificate run and checks under broker leases, preemptible, no paid compute.
- [n/a] Open hand-off.

## 1. The result
- [x] Precise statement (`README.md`, `PROOF.md` §0).
- [x] PROOF.md with every computer step named, including D\*, the 570-box orbit extension and strict containment.
- [x] Certificate hashes in `certs/SHA256SUMS.*`; data archived separately (Zenodo dataset record).

## 2. Verification
- [x] Exact checker `certcheck_fast.py` on all leaves, 1% cross-checked by the rational reference `certcheck.py`.
- [x] Strict checker `certcheck_strict.py` (the Lean leaf predicate) on all 7,473,755 + 87,828 leaves; the
  extension re-run on 2026-11-01 with the published version of the script.
- [x] Independent replay `replay/hcheck.py --strict`, written from the specification: 7,561,583 leaves, 0 failures.
- [x] Coverage: `certcheck_ext.py` (orbit coverage), `replay/orbitboxes.py`, Lean `orbit_cover`.
- [x] Counterexamples checked exactly: `symmetry_gap.py`, `corner_example.py`, Lean `DPrymak_not_fundamental`.
- [x] Lean: 0 `sorry`; axioms = `all_leaves_pass` + one `native_decide` + standard three. Public copy rebuilt
  and audited on 2026-11-01 (`certs/lean_audit_2026-11-01.txt`).
- [x] Quick check `code/check_sample.sh`: passes (about 4 min).
- [ ] Full `code/check_all.sh` from the archives: not run end to end in this layout (about 45 CPU-hours). Tested
  2026-11-01: archive and chunk hashes, extraction, and the sharded strict check of `fullcert_ext` (57 files,
  87,828 leaves, 0 bad) pass. The recorded verdicts come from the same scripts in the working repository.
- [x] Hostile referee, two rounds: first round found D and the stale draft; second round found the Lean coverage
  hypothesis mismatch (G1) and the partial strict run (G2). Both closed. Reports are in the private working
  repository.

## 3. Novelty and context
- [ ] Fresh novelty re-check on the day of publishing (arXiv, Zenodo, GitHub, OSF, web).
- [ ] Papadoperakis 1999 obtained (ask Prymak, or interlibrary loan) and attribution checked.
- [ ] Level rating rechecked (private).

## 4. Ship (Lee decides; nothing pushed)
- [ ] Zenodo dataset record (certificate archives) and note record: drafts via the private draft script; Lee
  reviews and publishes.
- [ ] Replace the two `zenodo.XXXXXXXX` placeholders in `README.md` and the note.
- [ ] Push ai-math-results; new GitHub release, so the repository concept DOI gets a new version.
- [ ] Same-day courtesy email to A. Prymak and A. Arman (skeleton in the private papers folder).
- [ ] LEDGER.md, BOARD.md and GATES.md updated.

## 5. After
- [ ] Follow-up date in LEDGER.md (about 2 weeks).
- [ ] Paper with full proofs; ideally with an expert co-author.

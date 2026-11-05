# Heptagon packing: every packing of congruent regular heptagons has density at most the Kuperberg double-lattice density

**Status: Computer-assisted proof, AI-generated, not peer reviewed, not formally verified.** The geometric lemmas (L1, L2, L3, L5 below) are written up as a working draft (`PROOFS-DRAFT.md`, not paper prose) that has been attacked by two hostile AI audits but not read by a human referee. The interval computation (L4) is a rerunnable program, not a stored per-box certificate. Six sweep slices that had resumed from checkpoints were replayed from scratch (run `rp1`): all finished with 0 stuck.

## Statement and sources

**Theorem (claimed).** Every packing of congruent regular heptagons in the plane has density at most

delta_DL = (2/97)(-111 + 492u - 356u^2), where u = cos(pi/7), so delta_DL = 0.892690686126508...

This is the density of the Kuperberg double-lattice packing of the regular heptagon. The packing was constructed by G. Kuperberg and W. Kuperberg (*Double-lattice packings of convex bodies in the plane*, Discrete Comput. Geom. 5 (1990) 389-397) and is conjectured to be optimal among all packings (see Kallus and Kusner). Its local optimality (among nearby packings) was proved by Kallus and Kusner ([arXiv:1509.02241](https://arxiv.org/abs/1509.02241)). The method here follows Hales and Kusner's proof of the analogous statement for the regular pentagon ([arXiv:1602.07220](https://arxiv.org/abs/1602.07220)). Related background: [arXiv:1305.0289](https://arxiv.org/abs/1305.0289) and [arXiv:2207.08959](https://arxiv.org/abs/2207.08959). The papers are not copied here; please read them at the links.

## Proof outline

Details are in `PROOFS-DRAFT.md` (proofs of the structural lemmas and of the interface to the computation) and `NOTES.md` (lab notebook: definitions, design decisions, run history).

- **L1, saturation and density.** A packing can be assumed saturated (no heptagon can be added). Then Delaunay circumradii of the centres are at most 2 (heptagon circumradius 1), and the density is an average of local quantities over Delaunay triangles.
- **L2, Delaunay facts.** The facts about Delaunay triangles that the bookkeeping uses (edge lengths between 2u and 4, circumradius at most 2, and so on). The draft proves these in finite windows, with a separate finite triangulation per window.
- **L3, transfer.** Across a long edge of an obtuse triangle, area is transferred to the neighbouring triangle by an amount bounded by a table, D5 (`runs/D5_nodes.npy`, `runs/D5_bins.json`, built by `src/build_D5.py`). This lets obtuse triangles be handled by a separate sweep.
- **L4, computer sweeps.** A rigorous interval branch-and-bound (`src/bnb2.py`, on top of the outward-rounded interval primitives in `src/rig.py`) shows that, for every admissible triangle configuration, the adjusted local area is at least the critical value A/(2 delta_DL) = 1.5326754446782212... (A = area of the regular heptagon of circumradius 1), tested against the target 1.5326754446782214. There are two sweeps: the **main sweep** (non-obtuse triangles) and the **obtuse-mode sweep** (using D5). Each is cut into 45 slices of one parameter. The main sweep splits slice 39, next to the double lattice, into ten sub-slices 39a-39j, and slice 38 was finished in six shards.
- **Certified tube near the double lattice.** At the double-lattice configuration the inequality is an equality, so branch-and-bound cannot close its neighbourhood. A tube of **1,530 charts** (`runs/charts_w030c.json`, made by `src/local_tube_arb.py`) covers a neighbourhood of the double-lattice contact curve. Inside the charts the inequality is proved by a local argument in arb ball arithmetic, and the sweeps exclude the charts (`src/chart_excl.py`). The charts are glued and rechecked by `src/tube_glue_arb.py` and `src/tube_recheck_arb.py`.
- **L5, averaging.** Summing the local bound over the Delaunay triangles of a window, with the L3 transfers cancelling in pairs, gives the density bound.

## What was checked

Numbers are taken from the logs in `runs/` and from the two audit reports and NOTES.

- **Interval branch-and-bound, run `a1`.** Main sweep: 54 jobs (44 slices, slice 39 as ten sub-slices), all `done: true, stuck: 0, open: 0`; slice 38 timed out twice and was finished from its dumped open set in six shards. Obtuse sweep: 45 of 45 slices `done: true, stuck: 0, open: 0`. Box counts: main 91,520,252, obtuse 34,047,698, together **125,567,950 boxes** (68.9 container-hours; for slices that were resumed from a checkpoint, only the work after the resume is counted). Every log begins with a provenance block (file hashes, numpy 2.5.4, environment, exact command line); the runs used `bnb2.py` sha256 `9f4fb628...`, `rig.py` `1ad0261b...`.
- **Local spot replays** (`runs/replay_a1/`, made by `src/local_replay.py`): obtuse slices 13 (1,066,185 boxes) and 30 (452,191) have box counts identical to the cloud runs; main slice 39c (2,083,897 boxes) is identical to the cloud run; main slice 5 (466,659 boxes) finished with 0 stuck, but its cloud run was resumed from a checkpoint, so only the verdict is comparable.
- **Tube.** `runs/tube_w030_fix.log` (1,531 charts computed, 1,530 certified; the walk stops at the first failure), `runs/tube_glue_w030c.log` (1,529 pairs, 0 failures for the gluing and window-overlap conditions) and `runs/tube_recheck_w030c.log` (all 1,530 charts rechecked after widening every exported number by 2 ulp; conditions (A)-(D), the inflated cube K', E', G, Ov and the window check W; 1,529 pairs, 0 failures).
- **Exact double-lattice check.** `runs/dl_exact.log` (from `src/dl_exact.py`, exact algebra): the double-lattice configuration satisfies the contact conditions and area(s2) * 2 * delta_DL = A exactly; three negative controls fail as they should.
- **Two independent hostile AI audits:**
  - [`AUDIT-2026-11-04.md`](AUDIT-2026-11-04.md) found no flaw that breaks the argument, but found three small rigour flaws (outward half-widths, parameter slivers, a chart offset) and some unclear items. The fixes were made and **both sweeps were rerun** (run `a1`); the logs here are the post-fix runs.
  - [`AUDIT2-2026-11-05.md`](AUDIT2-2026-11-05.md) audited the rerun. Verdict "FIX NEEDED (small)", with nothing found that breaks the argument. It closed the first audit's flaws, found one new formal flaw (FLAW-4: exported chart data rounded to nearest rather than outward, about 1e-18 against a slack of about 3e-21) and left the six resumed slices as unclear. FLAW-4 was afterwards closed without a rerun by `src/tube_recheck_arb.py` (log above; argument in `PROOFS-DRAFT.md`, Lemma 6.0). That later step was not itself re-audited.

**Caveats.**
- The branch-and-bound writes progress logs, **not** a per-box certificate, so "replay" means rerunning the program, not checking a stored proof object.
- Six slices that resumed from checkpoints (main 2, 39f, 39g; obtuse 7, 26, 44) were replayed from scratch, with checkpointing off and the same `bnb2.py` (`runs/sweep_main_rp1/`, `runs/sweep_obt_rp1/`). All are `done: true, stuck: 0`: main k2 499,793 boxes, 39f 2,367,945, 39g 8,430,765; obtuse k7 1,267,429, k26 505,137, k44 855,327. Every slice of both sweeps is therefore backed by a from-scratch run or a checked resume chain (slice 38: two dump resumes plus six shards, multiset-checked).
- `PROOFS-DRAFT.md` is a draft with internal markers; L1, L2, L3, L5 and the tube lemma are checked by audit only, not by a machine.
- Nothing is formalised in Lean or another proof assistant.

## How to replay one slice

Requirements: Python 3 and numpy (the cloud runs used CPython 3.13 and numpy 2.5.4 on Linux; the local replays used CPython 3.14 and numpy 2.5.4 on Windows). The tube and exact-check scripts also need python-flint, mpmath and sympy. Slices take minutes to hours on one thread. Example, obtuse slice 13 (1,066,185 boxes; about 25 minutes on one thread in the local replay), run from the folder that contains `src/` and `runs/`:

```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 HEPT_DUAL2=1 \
  HEPT_MU='[-1.1684855e-12, -1.1684853e-12]' \
  HEPT_OBT4=runs/D5_nodes.npy,runs/D5_bins.json
python src/bnb2.py runs/f_w030.npy 3 2 1.5326754446782214 40000 \
  '{"lo": [2.548, 0, 0.8, 0.26, -1e-12, -1e-12], "hi": [4.0, 2.0, 3.5, 0.28, 0.8975979010256552, 0.8975979010256552]}'
```

This is taken from the header of `runs/replay_a1/obt_k13.log`. The last line of the output should show `"done": true, "stuck": 0, "open": 0` and `"processed": 1066185`. For a main-mode slice, set `HEPT_CHARTS=runs/charts_w030c.json` instead of `HEPT_OBT4` and take the region from the header of the matching `runs/sweep_main_a1/k*.log`. `python src/local_replay.py obt 13` (or `main 5`, ...) builds the same command and writes into `runs/replay_a1/`. Compare the file hashes in your log header with the provenance blocks of the original logs.

Other checks: `python src/tube_glue_arb.py runs/charts_w030c.json`, `python src/tube_recheck_arb.py runs/charts_w030c.json` and `python src/dl_exact.py`.

## Files

- `PROOFS-DRAFT.md`: working proofs of L1, L2, L3, L5, the interface to L4 and the tube lemma. `NOTES.md`: lab notebook. `AUDIT-2026-11-04.md`, `AUDIT2-2026-11-05.md`: the two hostile audit reports. These four files are unedited except that local scratch-directory paths and one internal budget remark were removed; they still mention a shared compute broker and a cloud provider, which are operational notes.
- `src/`: all tracked source: interval primitives `rig.py`, branch-and-bound `bnb2.py`, tube code `local_tube_arb.py`, `chart_excl.py`, `tube_glue_arb.py`, `tube_recheck_arb.py`, table builder `build_D5.py`, exact check `dl_exact.py`, sweep drivers `modal_sweep.py`, `modal_obt.py`, `local_replay.py`, and earlier exploratory scripts.
- `runs/`: inputs of the certified runs (`f_w030.npy`, `dl_ref.json`, `charts_w030c.json`, `D5_nodes.npy`, `D5_bins.json`); logs of the final sweeps (`sweep_main_a1/`, `sweep_main_a1s/`, `sweep_main_a1r/`, `sweep_obt_a1/`, their driver logs, `replay_a1/`); the tube, gluing and recheck logs; `dl_exact.log`; the slice-38 shard script `k38_shards.sh`. The logs and scripts are unedited and contain the author's original machine paths.
- Not included: third-party PDFs (see the arXiv links above), older superseded runs and tables that predate the audit fixes, and the unfinished from-scratch replays.

## AI disclosure

This work was generated with AI (Claude, from Anthropic), directed by Lee Marshall. It was attacked by two hostile AI audits, whose reports are included, and cross-checked by local spot replays. It has not been through human peer review and is not formally verified. Please check the computations yourself.

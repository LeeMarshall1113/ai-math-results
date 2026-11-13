# Karp's evasiveness conjecture for ten vertices: eight of nine types excluded, type 5 open

**Status: computer-assisted partial result; hand proofs conditional on the fixed-point theorems of Smith and Oliver (cited, not re-proved); one AI referee read the proofs and source but did not run the replay; not independently replayed by a third party; no human peer review. It does NOT settle the n = 10 case.**

## Statement

A graph property on n vertices is *evasive* if every decision tree for it, in the worst case, asks about all n(n−1)/2 possible edges. Karp's conjecture (Aanderaa–Karp–Rosenberg) says that every nontrivial monotone graph property is evasive. It is known for prime powers n (Kahn–Saks–Sturtevant) and n = 6; n = 10 is the smallest open case.

Angel and Borja showed that a non-evasive monotone graph property on ten vertices has one of nine types (determined by which of ten special graph classes it contains) and excluded types 1, 3, 7 and 9, noting that they had found no groups for types 2, 4, 5, 6 and 8.

**What is proved here.** A non-evasive monotone graph property on ten vertices:
1. is of type 5: **eight of the nine types are excluded** (types 1–4 and 6–9; the exclusions of 1, 3, 7, 9 re-derive Angel–Borja with replayable certificates, those of 2, 4, 6, 8 are new);
2. contains every disconnected graph and no graph with disconnected complement;
3. contains nine explicit graphs (C₁₀, the Petersen graph, K₅ × K₂ and six others) and none of their complements;
4. has a trace on 32 graph classes that is one of eight explicit patterns, and is never self-dual.

**What is not proved.** Whether a type-5 property exists. Type 5 is open.

**Conditional on:** the cited fixed-point theorems (Smith 1941; Oliver 1975, as used by Angel–Borja and Engström) and Kahn–Saks–Sturtevant. Every other step is a finite computation with explicit permutation groups, replayed by the program below.

## Prior work

Angel and Borja (the nine types; types 1, 3, 7, 9), Borja's thesis, Engström (Oliver groups on ten vertices; his Theorem 14), Miller's survey (arXiv:1306.0110). Full references are in the note. Angel–Borja is arXiv:1603.04412 (J. Graph Theory 2019); Engström is arXiv:math/0512421. A web search on 2026-11-12 (arXiv, Google Scholar, GitHub, Zenodo, OSF) found no later work on n = 10 and no extension of the nine types; please tell us if we missed something.

## How to check it

Standard library only; no third-party packages, no network. Tested with CPython 3.14.0 on Windows 11. Do not run with `-O`.

```
cd results/karp10-partial
python replay/replay_all.py
```

Expected: a line `manifest: 20 certificates, 17 scripts, 23 rows; OK`, then the manifest digest, one line per row, and finally

```
ALL VALID (23 rows: 20 certificates + 3 further checks)
```

with exit code 0 (about 2 minutes on one thread, well under 1 GB of memory). The script first checks that `certs/` and `replay/` contain exactly the files listed in `replay/RELEASE-MANIFEST.json`, each with its pinned SHA-256 (CRLF normalised to LF), and then runs each of the 23 rows in a fresh isolated interpreter. It has no switch to skip a check.

**Manifest digest.** `sha256(replay/RELEASE-MANIFEST.json)` = `875ca163b77c0761b3be8e2b69fccfccc08da47451d39908fc796d472a1cb349`. The replay prints it on its second line; compare it with this README.

## Verification status (please read)

* The replay was written and run by the authors (Lee Marshall, with Claude). It passes on a clean copy of this folder.
* An AI referee read the proofs and the replay source and judged the method sound; it did **not** execute the replay.
* No independent third party has replayed it yet. The authors' own mutation controls (1,116 deliberately corrupted inputs, all rejected) and an independently written second implementation of the graph and homology computations (`replay/mutation_controls.py`, `replay/indep_graph.py`, `replay/indep_homology.py`, `replay/crosscheck_independent.py`) are included, but these are the authors' checks too.
* The group lists in the certificates were produced with GAP 4.12.1 from the conjugacy classes of subgroups of S₁₀; that generation script is not included. The replay does not rely on it: it re-verifies every group (order, normality, generation) from the generator lists in the certificates.
* The type-5 searches (solver searches over necessary conditions from about 1,500 subgroup classes, and a second bounded attempt) were inconclusive: they found neither an exclusion nor a counterexample. Their logs are not part of this release. A further 7.4-hour run on 11-12 November 2026 (Smith conditions of the eight 2-groups, one branch) also ended without a verdict.

## Files

| path | contents |
|---|---|
| `note/non_evasive_ten_vertices.pdf`, `.tex` | the technical note (11 pages) |
| `PARTIAL-THEOREM.md`, `TYPE5-ATTEMPT-2.md` | working documents exactly as pinned in the manifest on 11 November 2026 (statement and status of each part; the second bounded type-5 attempt). They are internal notes, kept unedited so that `python replay/make_release_manifest.py --verify` passes; where they say nothing is published, that was true when they were written |
| `note/check_explicit.py`, `note/z5_table.py`, `note/explicit/`, `note/z5_patterns.json`, `note/z5_table.tex`, `note/run_checks.sh`, `note/crosscheck_z3.py` | the scripts and data the note refers to (explicit groups of types 2, 4, 6, 8 and the pattern table); `crosscheck_z3.py` is optional and needs the `z3-solver` package, everything else is standard library |
| `certs/` | the 20 certificates (JSON) |
| `replay/` | `replay_all.py`, the checkers, the pinned manifest `RELEASE-MANIFEST.json`, the authors' mutation controls and independent cross-check, and `DEPENDENCIES.md` (audited list of every file, import and stdlib module the replay touches) |

The note was compiled from the `.tex`; the date line inside the note is a draft date.

## Corrections

If you find an error, a gap, or earlier work we missed, please open an issue.

## AI disclosure

The computations and the drafting were done by AI (Claude, from Anthropic) under Lee Marshall's direction. The proofs have not been through human peer review.

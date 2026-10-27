# Erdős #580: verification of n ≤ 18, and progress on n = 20

Written 2026-10-22 by the coordinator from the agent's report (last commit d94303d). To re-check: `python verify.py --quick`, or `python verify.py` (about 6 CPU-h; run it under a lease).

## Problem
In an n-vertex graph where at least ⌈n/2⌉ vertices have degree ≥ ⌈n/2⌉, is every tree with at most ⌊n/2⌋ vertices a subgraph? This is the vertex form; the edge form (Loebl–Komlós–Sós) is proved only for large n, by Zhao. Zeraoulia claims n ≤ 19 (Zenodo 21348157, forum post of 2026-07-16).

## Reduction: sound, and every case is covered
- **Edge-minimal host.** We may assume |L| = |S| = m, every L-vertex has degree exactly m, and S is independent. The argument is correct for every n = 2m, so we reused it at m = 10.
- **Leaf-extension (Hall) lemma:** correct.
- **Tree classification:** our own classifier reproduces his split exactly (32 / 7 / 3 / 5, with the same five exceptional trees).
- **One misquote, harmless.** He states Piguet–Stein Thm 3's odd-c condition as ℓ + c ≥ ⌊n/2⌋; the published condition is ⌊n/2⌋ + 1. All three of his trees in that class use the c-even branch, and the odd branch cannot apply at k = 8.
- **n ≤ 17, and n = 19 from n = 18:** correct.
- **A limitation of the method, not of his result.** The rooted-core lemma is only sufficient: for 25 trees at m = 6, 8, 9, 10 its formula is satisfiable even though the host contains the tree. His five trees are unaffected. At m = 10 we used partial leaf deletion.

## His certificates
- **Encoding.** We re-derived it. His avoidance clauses are identical to ours under the edge-variable map. His totalizer base matches "degree = 9" exactly (all 2^17 neighbourhood assignments per vertex checked).
- **Checking.** We converted all 4 DRUP traces to LRAT with drat-trim -L. **cake_lpr: `s VERIFIED UNSAT` on all 4**, and gate_lrat_check accepts all 4.
- **Independent re-proof.** With our own encoding, CaDiCaL 3.0.1 LRAT proofs and cake_lpr, we certified all 5 exceptional n = 18 trees, every non-star tree for m ≤ 8, and 45 of the 46 at m = 9. The exception is P9, which still relies on Bazgan–Li–Woźniak.

## n = 20 (which also gives n = 21)
- **Covered by published theorems:** 75 of the 106 trees (diameter ≤ 5, path-plus-star, even-c caterpillars).
- **Need SAT:** the other 31. **29 are certified** (CaDiCaL + LRAT + cake_lpr, plus gate_lrat_check). That includes the spider S(3,3,3): 3.07M clauses, 2.7 CPU-h, 3.0 GB LRAT.
- **Open:** the spiders S(2,2,5) (7.4M clauses, ~7 CPU-h) and S(2,3,4) (14.8M clauses, 13+ CPU-h; cake_lpr may hit the 6 GB WSL cap). Decision D69: the coordinator is proceeding locally with leased, preemptible jobs, S(2,2,5) first.

## Novelty
The forum and web searches show nothing beyond n = 19. This is the first verified-checker confirmation of his n = 18 certificates. Once the two spiders are done, n = 20 and 21 would be new.

## Trusted base
cake_lpr; `lks.py` (tree enumeration, base encoding, clause validator); `classify.py` plus the cited theorems; and the two lemmas.

## Pending (n = 20): m = 10 trees 30 and 45
Tree 30 first attempt: CaDiCaL UNSAT, gate_lrat_check ok, cake_lpr failed (heap 2500 MB exhausted). Re-running with CAKE_HEAP_MB=4500 (full regenerate + re-solve).
Logs (launched 2026-10-22 23:37 EDT): tree 30 (local job log, not included) ; tree 45 S(2,3,4) (local job log, not included) . Results append to results/m10.jsonl; n = 20/21 not settled until both show verified:true. If cake_lpr cannot fit 6 GB for tree 45, split by cube-and-conquer.

## Update 2026-10-23 (tree 45 still pending)
The tree-45 LRAT had been deleted (KEEP_LRAT was unset), so cake_lpr cannot be retried without re-solving. Re-run queued as lease erdos-3461 (KEEP_LRAT=1, CAKE_HEAP_MB=5300, CAKE_STACK_MB=400; cake_check.sh now honours CAKE_STACK_MB), log (local job log, not included) Not yet granted at 16:35 EDT (queue blocked behind other leases). Novelty sweep so far: erdosproblems.com thread 580 (only Zeraoulia n<=19, 16 Jul 2026, and Tang), web search, `gh search repos/code`: nothing for n>=20. n = 20/21 remain NOT settled.

## Update 2026-10-26: tree 45 VERIFIED; n = 20 and 21 settled (pending a final referee)
cake_lpr reports `s VERIFIED UNSAT` (rc 0) for tree 45, the spider S(2,3,4), on the kept LRAT from the 2026-10-23 re-solve. gate_lrat_check had already passed. This ran with WSL temporarily at 10 GB under a 10 GB lease (Lee's D73); log: d73_tree45.log. All 31 SAT trees at m = 10 are now certified, and the other 75 trees are covered by cited theorems. So the vertex form of #580 holds for n = 20, and by the existing step for n = 21. Novelty: nothing beyond n = 19 found (forum thread 580, web, GitHub). Status: CANDIDATE until one final hostile referee pass on the m = 10 classification and coverage is done.

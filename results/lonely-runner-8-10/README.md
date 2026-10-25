# Lonely runner conjecture, 8, 9 and 10 runners: formally verified certificates

**Status: machine-verified certificates (LRAT proofs checked by the CakeML-verified cake_lpr and by our own gate_lrat_check). The cases were already proved in unrefereed preprints; what is new is the verified certificates, from independent code.**

## Claim
The lonely runner conjecture for k = 7, 8, 9 speeds (8, 9, 10 runners) follows from certified computations plus cited theory:

| speeds k (runners) | gates | proof size | verified |
|---|---|---|---|
| 7 (8) | 28 primes (31..167) | 14.5M LRAT steps, 5.48 GB | 2026-10-21, both checkers |
| 8 (9) | 35 primes (47..283) | 38,734 leaf proofs, 21.1 GB | 2026-10-21, both checkers |
| 9 (10) | 45 primes (139..409) | 31,593 leaf proofs, 50.9 GB | 2026-10-21, both checkers |

## Prior work
- 8 runners: Rosenfeld, arXiv 2509.14111. 9 and 10 runners: Trakulthongchai, arXiv 2511.22427, and Rosenfeld, arXiv 2512.01912 (preprints; not refereed to our knowledge). 7 runners: Barajas-Serra 2008, used as an input.
- The reduction uses Malikiosis-Santos-Schymura (Forum Math. Sigma 13 (2025) e164).
- Ours: the CNF generator and checkers were written from the definitions, with no code from the papers; each gate is certified UNSAT and the checkers accept every proof.

## Method
Reduction to prime gates (cited); each gate is encoded as a CNF whose unsatisfiability is proved by CaDiCaL with an LRAT proof; the exact size inequality (`proofs/common/mass.py`) then gives a contradiction. For k = 8, 9 each gate is split into leaves (level-1 split) with one proof per leaf. The argument that UNSAT of the CNF implies a gate is `proofs/ENCODING.md`; the reduction for each k is `proofs/kN/PROOF.md`.

## How to verify
Not re-run in preparing this repository (the LRAT artifacts are tens of GB). The recorded verification logs and per-gate results are in `proofs/k7`, `proofs/k8`, `proofs/k9` (`verify_log*.jsonl`, `verify_gates.jsonl`, `verify_final.log`, `manifest.json` with gate lists and hashes).

To regenerate and re-check, from this folder:
```
bash proofs/common/checkers/build_unix.sh        # builds gate_lrat_check, fetches and hash-checks cake_lpr (Linux/WSL/macOS)
export CADICAL=/path/to/cadical GATE_LRAT_CHECK=... CAKE_LPR=...
python proofs/common/verify.py proofs/k7         # the one command; missing proofs are regenerated with CaDiCaL
bash proofs/common/run_wsl.sh k8                 # k8/k9, inside Linux/WSL (needs tens of GB of disk)
```
`run_wsl.sh` expects the binaries at the paths named in the script (edit them). For k = 9 the level-1 split uses a C++ program built from `proofs/common/checkers/level1.cpp` (set `LEVEL1_EXE`). Expect hours of CPU for k = 8, 9.

## Not included
The LRAT/CNF proof artifacts (`artifacts/`, 5 to 50 GB per k), compiled `.exe`/`.obj` files, and the k = 10 and later work (in progress).

## Caveats
- Trusted base: the CNF generator and exact inequality (`gen_cnf.py` / `gen_cnf_composite.py`, `mass.py`), the checkers, and the cited theorems (LRC(k-1), the MSS size bound, the gate lemma). Not formalised end to end in a proof assistant.
- Not a first proof of these cases.
- Some gates for k = 8, 9 were checked in a first run and the rest in a second; every gate was accepted by both checkers within a single run (see the PROOF.md files). Non-gate primes below 137 were not certified for k = 9.
- Log files had the machine-specific job-scheduler header lines removed.

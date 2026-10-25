# Certificate-checked proof packages (Q1)

| k (runners) | status | gates | proof size | verify time |
|---|---|---|---|---|
| 7 (8) | **VERIFIED** with cake_lpr and gate_lrat_check, 2026-10-21 | 28 (31..167) | 14.5M steps, 5.48 GB | 355 s wall |
| 8 (9) | **VERIFIED**, both checkers, 2026-10-21 | 35 (47..283) | 38,734 leaf proofs, 21.1 GB | 7,121 s summed over gates |
| 9 (10) | **VERIFIED**, both checkers, 2026-10-21 | 45 (139..409) | 31,593 leaf proofs, 50.9 GB | 14,496 s summed over gates |

Each package `kN/` contains:
* PROOF.md: the reduction with exact citations, the gate set and its exact mass check, the
  trusted base, and the timings.
* manifest.json: the gate list, plus hashes and solver options per gate.

The proofs themselves live in `kN/artifacts/`. They are gitignored because of their size. If
a proof is missing, verify.py regenerates it with CaDiCaL.

    python proofs/common/verify.py proofs/k7          # the one command
    bash proofs/common/run_wsl.sh k8                  # k8/k9 on this box (inside WSL; add --resume to reuse gates)

`common/`:
* gen_cnf.py: the trusted CNF generator.
* mass.py: the trusted exact inequality.
* verify.py: the driver.
* prove.py: the untrusted producer.
* test_encoding.py: brute-force validation of the generator.
* checkers/: the gate_lrat_check source, plus build scripts that also fetch and hash-check
  cake_lpr.

The argument that UNSAT of the CNF implies the gate is in ENCODING.md.

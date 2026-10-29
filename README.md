# ai-math-results

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23748158.svg)](https://doi.org/10.5281/zenodo.23748158)

Computer-assisted mathematical results by Lee Marshall, produced with AI assistance and checked by independent verifiers. Every result comes with its certificates and a checker you can run yourself.

Each folder under `results/` has its own README with:
- the precise statement and its sources;
- prior work, including any earlier or overlapping claims;
- the method;
- how to check it, and what the check does and does not cover.

## Results

| Result | Verification | How to check | Folder |
|---|---|---|---|
| **Collatz cycles:** no nontrivial positive cycle has 111 or fewer local minima. A further extension to 127 is computed but not formally verified | Lean 4, kernel-checked, conditional on two published inputs stated in the README | `lake build Collatz && lake env lean scripts/Statement.lean` in [collatz-m-cycles](https://github.com/LeeMarshall1113/collatz-m-cycles) (about 9–10 h) | [results/collatz-m-cycles](results/collatz-m-cycles) |
| **Agaev–Chebotarev polygon conjecture, n = 4:** Λ₄ = S₄. Together with n ≤ 3 (2005) and Chebotarev's disproof for n ≥ 5 (2026), this settles the conjecture. An earlier claim by D. Z. Liu is credited in the folder | Exact certificates for 836 face inequalities; three independent checkers, one using only the Python standard library | `cd results/polygon-n4/proof && bash fullcheck.sh` (about 10 min), or `cd indep3 && python -I check3.py` (about 5 min) | [results/polygon-n4](results/polygon-n4) |
| **Lonely runner, tight instances for 6 and 7 speeds:** every 6-speed set with maximal loneliness ≤ 1/7 is a dilate of (1,…,6); every 7-speed set with ≤ 1/8 is a dilate of (1..7), (1,2,3,4,5,7,12) or (1,4,5,6,7,11,13). Answers Perarnau–Serra's Problem 1 for 7 and 8 runners and gives Kravitz's Conjecture 8.3 for 7 and 8 speeds | Strict Malikiosis–Santos–Schymura bound (proved in full) plus prime gates; two independent implementations agree exactly | `cd results/lonely-runner-tight-6-7/code && bash run_all.sh` (about 2 h, one thread) | [results/lonely-runner-tight-6-7](results/lonely-runner-tight-6-7) |
| **No torus with 10 convex quadrilateral faces:** no acoptic polyhedral torus in R³ has exactly 10 strictly convex quadrilateral faces. This is the F = 10 case of Grünbaum–Szilassi Conjecture 3.2; with Seamons's 11-face torus (2026) it settles which face counts occur | LRAT certificates checked by cake_lpr (a formally verified checker); two independently written encoders agree | `cd results/torus-quad10/code && LRAT_CHECK=cake_lpr bash check_all.sh` (about 2 min) | [results/torus-quad10](results/torus-quad10) |
| **Erdős problem #580 (vertex form)** holds for n ≤ 21; n = 20 and 21 are new | SAT with LRAT proofs checked by cake_lpr (a formally verified checker), plus cited theorems | `cd results/erdos-580 && python verify.py --quick` (full run about 6 CPU-hours) | [results/erdos-580](results/erdos-580) |
| **Lonely runner conjecture, 8, 9 and 10 runners:** formally verified certificates for cases already proved in the literature | LRAT certificates checked by cake_lpr | `python proofs/common/verify.py proofs/k7` in the folder | [results/lonely-runner-8-10](results/lonely-runner-8-10) |
| **Ulam's packing conjecture** holds for origin-symmetric convex bodies within Banach–Mazur distance 1 + 6·10⁻⁸ of the ball | Interval arithmetic | `cd results/ulam-explicit-radius && python certify.py` (about 1 s) | [results/ulam-explicit-radius](results/ulam-explicit-radius) |
| **Turán density of K₄⁻:** π(K₄⁻) < 0.28674302. Since superseded by Lidický–Volec (0.2865, arXiv 2610.33850); kept as an exact, checkable certificate | Exact rational flag-algebra certificate | `cd results/turan-k4minus && python check_bound.py 0.2867433 certs/cert_plain.json.gz certs/cert_C.json.gz` (about 25 min) | [results/turan-k4minus](results/turan-k4minus) |

## Citing

Use the Zenodo DOI above, or the "Cite this repository" button (from `CITATION.cff`), and name the result folder. Each GitHub release is archived on Zenodo with its own version DOI.

## Corrections

If you find an error, a gap, or earlier work we missed, please open an issue. Earlier claims we learn about are credited in the relevant folder.

## AI disclosure

The work in this repository was generated with AI (Claude, from Anthropic), directed by Lee Marshall. Each result was checked by independent verifiers and attacked by hostile AI referees, whose reports are included where they exist. None of it has been through human peer review unless a result's README says so. Please check the certificates yourself where a check is provided.

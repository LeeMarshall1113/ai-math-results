# Agaev-Chebotarev polygon conjecture, n = 4: Lambda_4 = S_4

**Status: complete, computer-assisted, exact arithmetic. Passed one hostile AI referee ("SURVIVES"). No human peer review. No formal (Lean) check.**

## Statement and sources

Conjecture 1 of R. Agaev and P. Chebotarev, *On the spectra of nonsymmetric Laplacian matrices*, Linear Algebra Appl. 399 (2005) 157-168 (arXiv:math/0508176): for every n, Lambda_n = S_n. Here a standardized Laplacian is a real n x n matrix with zero row sums and off-diagonal entries in [-1/n, 0], Lambda_n is the set of all eigenvalues of all such matrices, and S_n is the convex hull of lambda_k = (1/n)(k - e - ... - e^k), e = exp(-2 pi i/n), and their conjugates.

**Result.** For n = 4, S_4 is the hexagon with vertices 0, (1 +- i)/4, (3 +- i)/4, 1, and every eigenvalue of every standardized 4 x 4 Laplacian satisfies |Im z| <= 1/4. With the published inclusions (2005 Theorems 5 and 6) this gives **Lambda_4 = S_4**. Equivalently, with M = 4L and weights w in [0,1]^12, every eigenvalue of M(w) has |Im| <= 1.

## History

- n <= 3: known (2005).
- n >= 5: the conjecture is false. P. Chebotarev, arXiv:2610.14100 (v1 8 Oct 2026, v2 20 Oct 2026), gives counterexamples for every n >= 5. It lists n = 4 as open (Question 1).
- n = 4: proved here, with open certificates that anyone can replay. This completes the conjecture.
- **Prior claim.** Our novelty search missed an earlier claim of the n = 4 case: D. Z. Liu, *The Agaev–Chebotarev Polygon Conjecture Holds Exactly for n ≤ 4*, Zenodo record [23630989](https://zenodo.org/records/23630989), created 2026-10-25 (01:43 UTC), before this work began. It uses a related approach: an LP reduction to at most six fractional weights, then exact Bernstein expansions on 9583 classes. As of 2026-10-28 the record contains only the PDF (its code is announced but not included). Liu's claim has priority; what this repository adds is an independent proof whose certificates are open and replay cleanly.

## Method (full details in `PROOF.md`)

1. **Identity.** For the cubic p(t) = t^3 - a t^2 + b t - c left after factoring out the eigenvalue 0, a polynomial Phi(a,b,c) = -Disc(p) - 4(a^2 - 3b + 4)^2 factors, with a real root r and complex pair x +- iy, as Phi = 4(y^2 - 1)(d^2 + (y+2)^2)(d^2 + (y-2)^2), d = r - x. By Bendixson, |Im| <= sqrt(3/2) < 2, so the last two factors are positive and Phi <= 0 rules out |Im| > 1.
2. **LP-vertex face reduction.** A hypothetical counterexample can be moved, by a linear program in the weights and the eigenvalue with the eigenvector fixed, to a vertex with at least 6 weights at 0 or 1 and at most 2 fractional weights per row. So it suffices to prove -Phi >= 0 on faces of [0,1]^12 with 6 free weights (38016 faces).
3. **Symmetry.** Relabelling vertices and the complement map M -> 4I - J - M reduce these to 836 canonical faces.
4. **Bernstein.** On 690 faces, exact Bernstein coefficients of the face polynomial are all nonnegative after bisection.
5. **Proof trees (136 faces) and tube certificates (10 faces).** The remaining faces are tight (max |Im| = 1 exactly). They are handled by exact recursive proof trees (monotone reduction, factorisation, splits, Duffy simplex maps, vertex blow-ups) and, for the last 10, by k-d "tube" certificates that compare N with its restriction to the zero-set components.

## How to verify

Requirements: Python 3, `sympy`, `python-flint` (pip install python-flint), numpy. Pure CPU, under 1 GB of memory.

```
cd proof
bash fullcheck.sh
```

This runs `verify.py` (re-derives identity (*), re-enumerates the 38016 faces and 836 representatives, **asserts that the Bernstein faces, trees and tubes validate as inputs and cover exactly the 836 canonical faces, disjointly**, then replays every Bernstein face and all 136 trees) and then `verify_tube.py` on the five tube files (replays all 10 tube certificates). It exits nonzero unless all 836 classes pass; the expected last line is `FULL CHECK PASSED: 826 + 10 = 836 canonical faces certified`. `proof/fullcheck.log` is the recorded run. Runtime: about 10 minutes on one thread, peak about 0.5 GB. `python verify.py` alone skips the tube certificates, so use `fullcheck.sh` or `indep3/check3.py all` for a full replay. `proof/CHECKSUMS.sha256` lists the certificate files (`cd proof && sha256sum -c CHECKSUMS.sha256`; the data files are stored byte-exact).

Certificates are included in full (all files under 5 MB, about 12 MB total); the search code that generated them is not. Checkers import only sympy, flint and each other.


### Third independent checker (pure Python, no CAS libraries)
`proof/indep3/` re-checks everything using **only the Python standard library** (Fractions and integers). It uses no python-flint, sympy or numpy, and none of the other checkers' code, so the two libraries the first two checkers rely on are no longer a single point of failure. It derives Phi, re-enumerates the 38016 faces and 836 orbits, rebuilds every face polynomial, replays every certificate with its own Bernstein/de Casteljau code, and checks exact, disjoint coverage by orbit.

```
cd proof/indep3
python -I check3.py           # about 4-5 min, under 0.1 GB; ends with ALL OK
python -I mutation_tests.py   # 47 deliberately corrupted inputs, all rejected
```
Limits are stated in `indep3/README.md`: it trusts CPython integer arithmetic, and the hand arguments (Bendixson, the LP-vertex reduction, why symmetry preserves |Im|) are proved in PROOF.md, not machine-checked.

## Referee summary

`notes/REFEREE.md` (hostile AI referee, 2026-10-28): verdict SURVIVES. The referee re-derived the identity and Bendixson bound (correcting the constant to sqrt(3/2); this folder uses that), checked the LP-vertex reduction and symmetry, independently recounted 38016 and 836, wrote its own checker (`notes/referee/refcheck.py`, own Bernstein and replay, no project code) and validated samples of all three certificate types plus mutation tests, and ran 400 numerical ascents finding no violation. Its minor points (no coverage assertion, unvalidated certificate inputs, stale text) are fixed here. Remaining caveats: both checkers rely on python-flint; the proof relies on the published 2005 Theorems 5 and 6 (Dmitriev-Dynkin bound).

## Files

- `PROOF.md`: the full argument. `notes/`: referee report and scripts, novelty check.
- `proof/`: checkers (`verify.py`, `vsome.py`, `verify_tube.py`, `fullcheck.sh`), `faces_canon.txt`, certificate `*.jsonl` files, checksums. Referee scripts take the proof directory as an argument, e.g. `python notes/referee/refcheck.py proof`.

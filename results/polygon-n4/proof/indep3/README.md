# indep3: a third, independent replay of the n = 4 proof

This directory holds a checker for `../../PROOF.md` that shares **no code and no
libraries** with the two earlier checkers (`../verify.py` + `../verify_tube.py`, and
`../../notes/referee/refcheck.py`). Both of those use python-flint, and the first also
uses sympy. P. Chebotarev asked for an independent replay in a second system.

## System used

- **CPython 3.14.0** (`C:\Python314\python.exe`, Windows) and **its standard library
  only**: `fractions`, `math` (`comb`, `gcd`, `lcm`, `isqrt`), `heapq`, `itertools`,
  `json`, `random`, `operator`.
- **Not used:** python-flint, sympy, numpy, gmpy2, or any other third-party package.
  Nothing from the project is imported (no `verify.py`, no search code).
- **No CAS.** I looked on 2026-10-28 and found none installed: no SageMath, Mathematica,
  PARI/GP, Maxima, Julia or Macaulay2, on Windows or in WSL (Ubuntu-24.04,
  docker-desktop). So all polynomial algebra is new code, written for this checker:
  - `poly.py`: sparse integer polynomials, ~340 lines;
  - `bernstein.py`: the Bernstein test, ~140 lines;
  - `check3.py`: the replay, ~690 lines.
- **Exact arithmetic throughout.** Python integers do all the polynomial work. Fractions
  are used only for box corners and for the rational sample points.

## What it shares with the other checkers

- **The certificate files (read-only), in `../`:** `result_0_836.jsonl`; the tree files
  `hard3_*`, `hard5_*`, `hard6_*`, `iso*_*`, `kd_*` (listed explicitly in `check3.py`);
  and the tube files `tube_t158`, `tube_tubeA`..`tube_tubeD`.
- **The data conventions of those files:**
  - face string position k is the matrix entry OFF[k], with
    OFF = [(i,j) for i in 0..3 for j in 0..3, i != j];
  - '1' means w_ij = 1, and the '*' entries are t0..t5 in string order;
  - the node semantics are those described in PROOF.md, Step 4 / 4b;
  - factor strings are written in python-flint's `str()` format. `poly.parse` reads them
    as plain text.
- **For information only:** `faces_canon.txt`. The checker prints whether the certified
  faces coincide with it, but nothing depends on that file.
- **The mathematics of PROOF.md Steps 1-3.** The hand arguments are not machine-checked
  here: Bendixson, the LP-vertex reduction to faces, and why symmetry preserves |Im|.
  What *is* machine-checked is listed below.

## What it checks (`python check3.py all`)

1. **Algebra.**
   - Disc(p) for p = t^3 - a t^2 + b t - c is computed as -Res(p, p'), using its own
     Sylvester determinant over Z[a,b,c]. It matches the textbook formula.
   - Phi = -Disc - 4(a^2-3b+4)^2 matches the expansion printed in PROOF.md.
   - Phi = Res_t(p(t+i), p(t-i)) / (8i), computed over Z[i][a,b,c].
   - Identity (*) holds as a polynomial identity in r, x, y.
2. **Faces.**
   - Faces are built row by row (at most 2 free entries per row, 6 in all, the rest in
     {0,1}): 38016 distinct faces.
   - Orbits under S4 x {id, complement}: 836.
   - The five Duffy tetrahedra used by `simplex` tile the cube (volume 1, plus a lattice
     test for gaps and overlaps).
3. **Face polynomials.**
   - N_F = -Phi(tr M, sum of 2x2 principal minors, sum of 3x3 principal minors). The
     minors are computed by the Leibniz formula, with det M = 0 asserted.
   - Cross-check: for every one of the 836 certified faces, N_F agrees exactly with
     Faddeev-LeVerrier (Fractions) at 2 random rational points.
   - The symmetry reduction is checked exactly: N_{gF}, composed with the
     relabelling/complement map of variables, equals N_F. This is done for 150 random
     (F, g) pairs, and for 20 faces with all 48 group elements.
4. **Certificates.** Every node type is replayed with this code's own composition
   routines and its own Bernstein test:
   - node types: bern, mono, split, splitat, blowup, simplex, factor, zero, const;
   - in the tube certificates: k-d split, bern, tree, reduce, and the TG
     sub-certificates (trees or nested tubes).

   The Bernstein test differs from both earlier checkers. It builds the tensor of
   Bernstein coefficients once, with integer weights lcm_j C(d,j)/C(d,j), and then
   subdivides with the **de Casteljau triangle** on that tensor. It never re-expands a
   composed polynomial. A box with a negative corner value fails at once.

   `factor` is checked **without factorisation software**:
   - each listed odd factor is parsed and divided out exactly, as often as it divides;
   - the remaining quotient must be (positive or negative constant) x (exact square). The
     square root is found term by term, and S*S == R is checked;
   - each factor of odd multiplicity has its sign proof replayed;
   - the product of the signs must be positive.

   A listed factor that does not divide is a rejection.

   Tube walks also check that the k-d leaves have total volume exactly 1.
5. **Coverage.**
   - The 690 Bernstein faces, 136 tree faces and 10 tube faces are pairwise disjoint.
   - Each one is a genuine face.
   - Each of the 836 orbits contains exactly one certified face.

   This is an orbit-based check. It does not rely on the project's choice of canonical
   representative.

Faces that have more than one certificate in the files (2 tree faces, and tube faces
with several `ok` records) have **every** certificate replayed.

## Results (2026-10-28)

See `fullcheck3.log` and `mutation_tests.log`.

Results of `python check3.py all`, run on one CPU thread. The run took 4 min 32 s, and
peak commit was 0.09 GB.

- identity (*), the Phi expansion, the resultant form, and the discriminant: **OK**;
- faces **38016**, orbits **836**;
- N_F vs Faddeev-LeVerrier on all 836 faces: **OK**;
- the symmetry identities: **OK**;
- Bernstein faces: **690/690**;
- proof trees: **138/138**. These cover the 136 tree faces; two faces have two
  certificates each in the files, and both of each pair pass;
- tube certificates: **10/10**, including every TG sub-certificate and nested tube, with
  k-d leaf volume exactly 1;
- coverage: 690 + 136 + 10 = 836, exactly one certified face per orbit, disjoint:
  **OK**;
- **ALL OK**. Nothing disagreed with the earlier checkers.

Nodes replayed:

| node    | count |
|---------|------:|
| bern    | 4925  |
| blowup  | 1359  |
| factor  | 1276  |
| mono    | 946   |
| splitat | 377   |
| simplex | 93    |

In all there were 7614 Bernstein tests on 9399 boxes.

The first run went through `lease-run.sh` (lease hunt-4094) and is kept in
`fullcheck3_lease_run1.log`; its results are identical. After that run, `check3.py`
changed by a pure refactor: identity (*) was moved into `star_identity()` for the
mutation tests. `fullcheck3.log` is the rerun with the committed code.

## Mutation tests (`python mutation_tests.py`)

Each test first replays the original input, which must be accepted, and then a
corrupted copy, which must be rejected. See `mutation_tests.log`.

**47/47 pass.** Each mutant below was rejected, and its original was accepted:

- **Bernstein routine:** (2t-1)^2 - 1/100, and (t0-t1)^2 - 10^-4.
- **Algebra:**
  - a coefficient of Phi changed: identity (*) fails;
  - Phi's constant term negated: identity (*) fails.
- **Face polynomial:** N_F with one coefficient off by 1 is caught by Faddeev-LeVerrier.
- **Shifted polynomials:**
  - a Bernstein face, shifted by -10^6;
  - a tight face, by plain Bernstein;
  - a valid tree replayed against N - 10^-6;
  - a valid tree replayed against N - 10^-12.
- **mono:** the direction flipped (3 trees); the variable changed (2 trees).
- **factor:**
  - the sign s of an odd factor flipped (3);
  - a factor string altered (3);
  - an odd factor dropped from the list (3).
- **blowup:** a corner bit flipped (3); a cone dropped; a variable removed from S (2).
- **splitat:** r outside (0,1); the cut point moved from 1/2 to 1/64 (2).
- **simplex:** a tetrahedron dropped; a corner bit flipped.
- **Tube certificates:**
  - replayed against N - 10^-6;
  - a reduce pin flipped;
  - a TG entry deleted;
  - all reduce leaves replaced by plain bern;
  - a k-d box dropped;
  - a nested TG tube weakened;
  - the tree leaves of a tube replaced by plain bern.
- **Coverage:**
  - a Bernstein face removed;
  - a tube face removed;
  - a tree face also listed as Bernstein;
  - a second face from an already-certified orbit added;
  - a non-face inserted.
- **Face rule:** a face with 3 free weights in one row is not among the 38016 faces.

One test design error was found and fixed before the final run. The first version used
(t0-t1)^2 itself as an "original that must be accepted". Plain bisection cannot certify
it, because it vanishes on the whole diagonal; that is why the `blowup` move exists. The
original is now (t0-t1)^2 + 1/100.

## Limits

- This checker trusts CPython's integer arithmetic, and the JSON files as the statement
  of the certificates.
- It does not machine-check the hand arguments of PROOF.md: the Bendixson bound, the
  LP-vertex reduction to 6-free faces, and that relabelling and complement preserve |Im|.
  Their algebraic consequence for N_F is checked exactly, but on a sample (above).

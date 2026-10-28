# Referee report: the computer-assisted proof of Lambda_4 = S_4

Hostile referee pass, 2026-10-28 (about 10:30 am EDT). I had not seen this work
before. My aim was to break it. Scripts and logs are in `notes/referee/`.

## Verdict: SURVIVES

I found no mathematical error. Every step either checks by hand, or was
re-verified by my own code, which imports nothing from the project. The weak
points are listed at the end. None of them is a hole in the argument.

## 1. Statement (checked against arXiv math/0508176 and arXiv 2610.14100v2)

- **Definition.** A standardized Laplacian has zero row sums and
  -1/n <= l_ij <= 0 off the diagonal (2005, eq. (5); 2026, Section 2).
  Lambda_n is the union of the spectra. This matches PROOF.md.
- **S_n.** In 2005, Section 4, S_n is the convex hull of
  lambda_k = (1/n)(k - mu - ... - mu^k), with mu = e^{-2 pi i/n}, and of their
  conjugates. For n = 4 this gives 0, (1+i)/4, (3+i)/4, 1 and their conjugates.
  That hexagon equals {|Im| <= Re, |Im| <= 1 - Re, |Im| <= 1/4}. Correct.
- **2005 Theorem 5.** The spectrum lies in two disks, two angles and a band. The
  angle at 0 has half-angle pi/2 - pi/n = pi/4, which gives |Im| <= Re. The angle
  at 1 is bounded by the rays through e^{+-2 pi i/n} = +-i, which gives
  |Im| <= 1 - Re. This is exactly what PROOF.md uses.
- **2005 Theorem 6.** S_n is contained in Lambda_n. The proof uses convex
  combinations of 0, L_k and L_{k+1}. These are polynomials in Q, and they are
  standardized Laplacians. Valid.
- **The equivalence** Lambda_4 = S_4 <=> sup |Im| = 1/4 is correct. The 2026
  paper states it explicitly in Section 6, and it lists n = 4 as open
  (Question 1).
- **Scaling.** M = 4L, with w in [0,1]^12. Correct.

## 2. Algebraic criterion

- **Factor t.** Since M1 = 0, det M = 0. So det(tI - M) = t(t^3 - a t^2 + b t - c),
  where a, b, c are the elementary symmetric functions of the eigenvalues, i.e.
  the trace and the sums of the 2x2 and 3x3 principal minors. A non-real
  eigenvalue is therefore a root of the cubic p.
- **Phi.** I re-derived Disc(p) = 18abc - 4a^3c + a^2b^2 - 4b^3 - 27c^2 by hand.
  The expansion of Phi in PROOF.md matches -Disc - 4(a^2 - 3b + 4)^2. With sympy
  I checked Res(p(t+i), p(t-i)) / (8i) = Phi and the factorisation (*).
- **Bendixson.** It gives |Im| <= ||(M - M^T)/2||_2. The paper bounds this by
  sqrt 3. In fact the bound is sqrt(3/2) (Frobenius^2 = 2(s1^2 + s2^2) <= 3).
  Either way it is below 2, so y > 1 implies Phi > 0.
- **Real roots.** If all roots are real, Phi = -Disc - (square) <= 0, so the
  criterion is an equivalence. Only the direction "Phi <= 0 implies |Im| <= 1"
  is needed.

## 3. LP-vertex reduction, symmetry, complement

- **The LP.** Fix the complex eigenvector z. The LP in (w', Re mu', Im mu') has
  8 real equality constraints. These are linear because z is fixed, so complex
  mu is handled correctly. The LP is bounded, because mu' is determined by a row
  with z_i != 0. It is pointed, because dmu * z = 0 forces dmu = 0. So the
  maximum is attained at a vertex, with Im mu* >= Im mu > 1.
- **At most 6 fractional weights.** A vertex needs 14 independent active
  constraints. The equalities give at most 8, so at least 6 of the w are at a
  bound.
- **At most 2 fractional per row.** Suppose row i had three fractional weights.
  Two real equations in three unknowns leave a null direction supported on
  row i, with mu fixed. That gives a feasible segment through w*, which
  contradicts it being a vertex.
- **On Proposition 17.** The global LP is needed for the "6 total" claim.
  Chebotarev's Proposition 17 works row by row with lambda fixed, and gives only
  the per-row bound. PROOF.md's LP gives both bounds. Correct.
- **Padding to exactly 6 free.** This is possible, since 4 rows x 2 = 8 >= 6.
- **Count.** The number of faces is 594 * 64 = 38016. Here 594 is the
  coefficient of x^6 in (1 + 3x + 3x^2)^4. I recomputed this independently.
- **Symmetry.** Relabelling is a permutation similarity.
- **Complement.** M^c = 4I - J - M. A left eigenvector u for mu != 0 has
  u^T 1 = 0, so mu is mapped to 4 - mu and |Im| is unchanged. The map sends
  faces to faces: fixed bits are flipped and t -> 1 - t.
- **Orbits.** My union-find over S_4 x Z_2 gives 836 orbits. Every orbit
  contains exactly one face of `faces_canon.txt`.
- **Coverage.** The certified sets are 690 Bernstein faces, 136 trees and 10
  tubes. They are pairwise disjoint, and their union is exactly the 836
  canonical faces (`referee/coverage_scan.py`).

## 4. Certificates: independent re-implementation (`referee/refcheck.py`)

What my checker does differently from the project's:

- **Its own N_F.** The trace and the minors are expanded explicitly with flint,
  and the hand-expanded Phi from PROOF.md is applied. The sympy route of
  verify.py is not used. I cross-checked N_F against numpy's characteristic
  polynomial on 60 faces x 5 random points: relative difference 7e-13.
- **Its own Bernstein test.** It checks the coefficients of
  (1 + x)^d p(x / (1 + x)). This is a different formula from verify.py's.
- **Its own replay** of every node type.
- **Its own k-d walk** for the tubes. This walk also checks that the leaf
  volumes sum to exactly 1, so the boxes cover each face with no gaps.

Results:

- **Bernstein faces:** 20 of 20 random faces OK.
- **Proof trees:** 10 of 10 OK. The sample includes the 4 slowest:
  `**0**0**1110`, `**0**0*011*0`, `**0*0*0*0*01` and `**0**0*0*110`.
- **Tube certificates:** 10 of 10 OK, and so are all their TG sub-certificates.
  Partition volume is exactly 1 at every level.
  - The `reduce` step is sound. With lo_l = 0, a pin is checked by dN/ds_l >= 0
    on the box; with hi_l = 1, by dN/ds_l <= 0. The pins are substituted in
    order, and N|_{pins} >= 0 is proved on all of [0,1]^6. Together these give
    N(t) >= N(pi_S(t)) >= 0.
  - The crossing moves (vertex blow-up and rational cut) live inside `tree`
    leaves, and my replay verifies them.
- **Mutation tests** (my checker must reject false claims). It rejects:
  - N on a tight face, by plain Bernstein;
  - N - 10^-6 on that face;
  - a valid tree replayed against N - 10^-6;
  - a tree with one `mono` direction flipped.
- **Numerical attack.** I ran 400 projected random-ascent runs of max |Im mu|
  over [0,1]^12. The best was 1.0000000000000029, which is the extremal value.
  No violation was found.

## 5. Replay of the project's own checker

I ran `vsome.py 500` and `verify_tube.py` on the 5 tube files under a lease
(log: `referee/fullcheck_ref.log`). Results:

- identity OK;
- 38016 faces, 836 canonical;
- Bernstein faces 690/690;
- trees 136/136;
- tubes 10/10.

On independence: verify.py, verify_tube.py and vsome.py import only sympy,
flint and each other. They import no search code (prover*, tube, iso*, kd,
phiface), and only read the JSON output of the search.

## Weakest points (none breaks the proof)

1. **Shared trusted base.** Both checkers rely on python-flint for mpoly
   compose, subs, derivative and factor. My checker does not use sympy for N_F,
   so that part is independent, but a flint bug in `compose` would hit both.
   For a journal version, a third replay in a different CAS would remove this,
   or a formal proof (e.g. Lean with exact Bernstein).
2. **fullcheck does not assert coverage.** `fullcheck.sh` runs `vsome.py`,
   which does not check that easy ∪ trees ∪ tubes = the 836 canonical faces.
   `verify.py main()` does, but it is not run. I checked coverage separately:
   it holds. **Fix:** add the assert.
3. **Unvalidated inputs.** verify.py does not check its inputs:
   - `mono.at` is not checked to be in {0,1};
   - `reduce` pins v are not checked to be in {0,1};
   - variable indices are not range-checked.

   A value outside {0,1} would be silently mis-handled. All values in the data
   are 0 or 1 (scanned), so no certificate is affected. **Fix:** add the asserts.
4. **Stale text in PROOF.md.**
   - Section 2 "Results" still says "103 of them are done; 43 remain". The
     Status section says all 836 are done.
   - "Step 4b" sits after Section 3.

   Clean these up before publication.
5. **Unproved extremiser claims.** The "Extremisers" list in Section 1 rests on
   numerics. It is not part of the proof and should be labelled that way.
6. **Proof dependence.** Λ4 ⊆ S4 depends on 2005 Thm 5 (angles) and on
   Thm 6 (S4 ⊆ Λ4), both published. Thm 5's angle part rests on the
   Dmitriev–Dynkin bound via Proposition 1. This is standard, but it is a
   dependency, not something re-proved here.

## Run note

The lease (hunt-4071, 1.5 GB) logged "lease lost (reaped)" mid-run, though the
job completed. Peak commit was 0.51 GB.

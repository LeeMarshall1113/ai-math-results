# The 2005 polygon conjecture at n = 4

**STATUS: complete, computer-assisted.** All 836 face inequalities have exact certificates, and an independent checker replays every one of them (see Section 4).

## 1. Statement and source

**Original source.** R. Agaev, P. Chebotarev, *On the spectra of nonsymmetric Laplacian
matrices*, Linear Algebra Appl. 399 (2005) 157-168 (arXiv:math/0508176), Conjecture 1.

**Current status.** P. Chebotarev, *Digraph Laplacian spectra: The 2005 polygon conjecture
fails for all orders above four*, arXiv:2610.14100 (v1 8 Oct 2026, v2 20 Oct 2026). It
proves that the conjecture fails for every n >= 5; it was already known for n <= 3. For
n = 4 it says "the question remains open" (Question 1, Section 7.3). Searches on
2026-10-25 (arXiv, the web) found no later work on n = 4.

**Definitions** (2005 paper, eq. (5); 2026 paper, Section 2).
- A *standardized Laplacian* of order n is a real n x n matrix L = (l_ij) with zero row
  sums and -1/n <= l_ij <= 0 for i != j.
- Lambda_n is the set of all eigenvalues of all such matrices.
- With eps = e^{-2 pi i/n}, set lambda_k = (1/n)(k - eps - ... - eps^k), k = 0..n-1.
  S_n is the convex hull of the lambda_k and their conjugates.
- Conjecture 1: Lambda_n = S_n.

**n = 4.**
- S_4 is the hexagon with vertices 0, (1 +- i)/4, (3 +- i)/4, 1. Equivalently,
  S_4 = { |Im z| <= Re z, |Im z| <= 1 - Re z, |Im z| <= 1/4 }.
- S_4 is contained in Lambda_4 (2005, Theorem 6).
- Lambda_4 lies in the two angles |Im z| <= Re z and |Im z| <= 1 - Re z (2005, Theorem 5,
  the Dmitriev-Dynkin bound).
- Hence (2026 paper, Section 6) **Lambda_4 = S_4 if and only if every eigenvalue of every
  standardized 4 x 4 Laplacian has |Im| <= 1/4.**

**Scaled form used below.** Put M = 4L. Then M = D - W, where W = (w_ij) has zero
diagonal and off-diagonal entries w_ij in [0,1]: 12 free parameters, the cube [0,1]^12.
D = diag(row sums of W). In this form the claim is:

> **Theorem.** For every w in [0,1]^12, every eigenvalue mu of M(w) satisfies |Im mu| <= 1.

**Extremisers (numerical observations, not part of the proof).** |Im mu| = 1 is attained on a large set:
- the directed 4-cycle (mu = 1 + i);
- the 4-cycle with all four chords i -> i+2 of equal weight t (mu = 1 + 2t + i);
- more generally the 4-cycle with chord weights satisfying t1 + t3 = t2 + t4. On this
  face the quantity in Section 2 factors as -(t1 - t2 + t3 - t4)^2 * (positive);
- the duplicated directed 3-cycle, with arbitrary weights between the twins (mu = 2 + i);
- in the eigenvector picture, a continuum of further configurations (numerics in
  `scratch/`).

So the maximum is reached on a positive-dimensional set, not at isolated points.

## 2. Proof

### Step 1: an exact algebraic test for |Im| <= 1

Since M has zero row sums, det(tI - M) = t * p(t), where

    p(t) = t^3 - a t^2 + b t - c,
    a = tr M,  b = sum of the 2x2 principal minors,  c = sum of the 3x3 principal minors.

The roots of p are the eigenvalues of M other than one copy of 0, so any non-real
eigenvalue is a root of p. Define

    Phi(a,b,c) = -Disc(p) - 4 (a^2 - 3b + 4)^2
               = -4a^4 + 4a^3c - a^2b^2 + 24a^2b - 32a^2 - 18abc + 4b^3 - 36b^2 + 96b + 27c^2 - 64.

(Phi equals Res_t(p(t+i), p(t-i)) / (8i).) If p has roots r and x +- iy, then

    Phi = 4 (y^2 - 1) (d^2 + (y+2)^2) (d^2 + (y-2)^2),   where d = r - x.      (*)

This is checked symbolically in `proof/verify.py`.

By Bendixson's theorem, |Im mu| <= ||(M - M^T)/2||_2 <= sqrt(3/2) < 2. (Here (M - M^T)/2 is real
antisymmetric with 12 entries in [-1/2, 1/2], so its Frobenius norm squared is at most 3.
Its eigenvalues are +-i s1, +-i s2, so F^2 = 2(s1^2 + s2^2) >= 2 s1^2, and the spectral
norm s1 is at most sqrt(3/2).) So the last factor in (*) is positive, and **y > 1 implies Phi > 0.** In
other words, Phi(w) <= 0 at a point w rules out an eigenvalue with |Im| > 1 there. (If
all roots of p are real there is nothing to prove.)

### Step 2: reduction to faces with at most 6 free weights

Suppose some w has an eigenvalue mu with Im mu > 1 and eigenvector z != 0. Fix z, and
solve the linear program in (w', mu') in [0,1]^12 x C:

    sum_{j != i} w'_ij (z_i - z_j) = mu' z_i   (i = 1..4),     maximise Im mu'.

- The LP is feasible (it contains (w, mu)).
- It is bounded: take an i with z_i != 0.
- Its feasible region is pointed, since a direction (0, d mu) needs d mu * z = 0.

So the maximum is attained at a vertex (w*, mu*), with Im mu* >= Im mu > 1, and M(w*) has
eigenvalue mu*. Two facts about this vertex:

- **At most 6 weights are fractional.** A vertex needs 14 active constraints of full
  rank. The equalities give at most 8, so at least 6 box constraints are active.
- **At most 2 weights per row are fractional.** If row i had all three fractional, then
  2 real equations in 3 unknowns, together with the other rows and mu* held fixed, leave
  a feasible line through w*, which contradicts its being a vertex. This is the argument
  of Chebotarev's Proposition 17.

Free some extra coordinates if needed (the per-row capacity is 8 >= 6). Then w* lies in
a *face*: exactly 6 free coordinates, at most 2 per row, the other 6 fixed in {0,1}.
There are 38016 such faces.

### Step 3: symmetry

Relabelling the vertices (a permutation similarity) preserves the spectrum.

Complementation w -> 1 - w gives M^c = 4I - J - M. If mu != 0 is an eigenvalue of M with
left eigenvector u, then u^T 1 = 0 and u^T M^c = (4 - mu) u^T, so |Im| is preserved.

Hence it suffices to treat one face per orbit of S_4 x {id, complement}. There are **836
canonical faces** (`proof/faces_canon.txt`).

### Step 4: Phi <= 0 on every canonical face (computer-assisted, exact)

For each canonical face F, N_F := -Phi restricted to F is a polynomial with integer
coefficients in 6 variables, of degree <= 4 in each. We prove N_F >= 0 on [0,1]^6 with
exact rational arithmetic only. Each proof is a tree made of the following steps, each
of which is a sufficient condition:

- `bern`: every Bernstein coefficient is >= 0, possibly after exact bisections.
- `factor`: N = c * prod f_i^{e_i} (checked by expanding); each odd-multiplicity factor
  gets a sign proof.
- `mono`: dN/dt_l >= 0 (or <= 0) on the box is proved by Bernstein, so
  N >= N|_{t_l=0} (or N|_{t_l=1}), and the restriction is proved recursively.
- `blowup`: the box is covered by the cones {x_m = max_{j in S} x_j}, where
  x_j = |t_j - c_j|. Each cone is the exact image of the box under x_j = x_m s_j, and
  N composed with that map is proved recursively. This resolves the zeros of N that sit
  on lower-dimensional faces of the box, and, for S = {a,b}, those on the diagonals
  t_a = t_b and t_a = 1 - t_b.
- `splitat`: cut t_l at a rational point.
- `simplex`: in local corner coordinates x on 3 variables, the cube is the simplex
  {sum x <= 1} together with 4 tetrahedra coned from (1,1,1). Each piece is the image of
  the box under a polynomial (Duffy) map, so a plane zero set sum x = 1 becomes a box
  face.

**Results.**
- 690 faces are settled by Bernstein alone, with no tight point.
- The other 146 faces all reach max |Im mu| = 1 exactly (numerically), so N_F has zeros
  in the box. These need the blow-up steps. 136 of them are settled by proof trees and the last 10 by
  tube certificates (Step 4b below); see Status.

### Step 4b: tube certificates for the last 10 tight faces (`proof/tube.py`)

On each of the last 10 faces, the exact zero set Z of N is a union of components:
box faces on which N vanishes identically, and squared linear factors (diagonals and
planes) inside box faces. `proof/comps.py` computes Z.

The face box [0,1]^6 is cut by a k-d partition into boxes. Each box gets one of three
certificates:
- `bern`: the Bernstein coefficients of N on the box are all >= 0.
- `reduce` S (the transversal argument):
  - S is a list of pins t_l = v shared by the components that touch the box.
  - In order, and with the earlier pins already substituted, dN/dt_l has the sign that
    makes N decrease towards t_l = v on the box. This is certified by Bernstein.
  - Hence N(t) >= N(pi_S(t)), where pi_S sets those coordinates to their pins.
  - N restricted to the face G_S = {pins} is proved >= 0 on the whole box by a separate
    certificate TG[S]. That is an ordinary proof tree, or recursively a tube certificate.
- `tree`: an ordinary proof tree for N on the box. This is used near crossings.
  Crossings need two extra moves: a vertex blow-up over the support of a mixed-sign
  squared linear form in the local quadratic part, and a cut at the rational root of a
  one-variable squared factor. Example: face `**0**0*01*10` at q = (1,1,0,1,0,1), with
  t3 = 1. The local form is 9[(x4+x5-x1)^2 + 4 x0 x1] + 36 x1 (x2+x5)^2 + ... ; blow up
  over {x1, x4, x5}, then cut at t4 = 1/2.

`proof/verify_tube.py` replays all of these using `verify.py`'s independent Bernstein
code and its tree replay.

Steps 1-4b give the Theorem, and hence Lambda_4 = S_4.

## 3. How to check

    cd proof
    bash fullcheck.sh       # the full replay: all 836 faces; nonzero exit unless every class passes

or, with the independent third checker, `python indep3/check3.py all`.

`fullcheck.sh` runs `verify.py` (Bernstein faces and trees) and then `verify_tube.py` (the 10 tube
certificates). `verify.py` on its own skips the tube certificates' content, so it is not a full replay.

`verify.py` does the following, without importing any of the search code:
1. It re-proves identity (*).
2. It re-enumerates the 38016 faces and their 836 orbit representatives.
3. It recomputes each N_F with sympy (trace and principal minors), independently of the
   flint code that searched.
4. It replays every proof tree, with its own Bernstein implementation (integer, scaled by
   d!) and its own bisection.

Requirements: Python 3, sympy, python-flint, numpy.

## 4. Status (2026-10-28)

**All 836 canonical faces are certified:**
- 690 by Bernstein;
- 136 tight faces by proof trees;
- 10 tight faces by tube certificates.

**End-to-end independent replay** (`proof/fullcheck.sh` -> `proof/fullcheck.log`):
- identity (*) OK;
- 38016 faces, 836 canonical;
- Bernstein faces 690/690;
- tight-face trees 136/136;
- tube certificates 10/10.

**Conclusion:** every eigenvalue of every standardized 4x4 Laplacian satisfies
|Im lambda| <= 1/4. Hence Lambda_4 = S_4, and the 2005 polygon conjecture holds for n = 4.

The earlier notes below (on SOS and on blow-up attempts) are kept for the record.

# Hostile referee report: PROOF.md (explicit Kallus neighbourhood, d_BM <= 1+6e-8)

Date: 2026-10-24. The referee had not seen this work before. Source for Kallus: the arXiv
source in `refs/k1/pessimum-v2.tex` (arXiv:1212.2551v2, "final submitted manuscript").

## Verdict: SURVIVES

I could not break the argument. Every lemma checks line by line, and every inequality points
the right way. An independent recomputation reproduces the margin. It uses different code, a
different cutoff, closed-form Legendre values and mpmath intervals.

| quantity | PROOF.md | referee (`notes/referee_check.py`) |
|---|---|---|
| a (Lemma 4) | 5.7774e-8 | 5.77736e-8 (condition holds) |
| eps | <= 1.7555e-7 | 1.755471e-7 |
| L | <= 3.46526e-4 | 3.465258e-4 |
| cbar | 8.4881e-4 | 8.488116e-4 |
| kappa | <= 4.37450e-4 | 4.374497e-4 |
| S (finite + tail) | 14097.31 + 0.196 | 14097.175 + 1.020 (cutoff l=4000); 14097.31 + 0.196 (cutoff 20000) |
| N_eff | <= 180.733 | 180.736 (cutoff 4000); 180.733 (cutoff 20000) |
| **12 kappa N_eff** | **0.94874** | **0.94875 / 0.94874**, a margin of about 5.1% |
| delta = 6.6e-8 | 0.9965, passes | 0.99653, passes |
| delta = 6.7e-8 | fails | 1.00429, fails |

I also checked the following numerically:

- Lemma 2: traces are 1/2, g1 = sqrt(3/8), G2 = sqrt 6 and ||Q^(k)||_op = 1.
- Every non-minimal vector of Lambda_0 has norm^2 >= 2.
- c_l > 0 for every even 4 <= l <= 20000. The minimum is c_10 = 65/65536, then c_4 = 7/32 and c_22 = 0.319.
- A Monte Carlo on SO(3) confirms the key identity A*A = 1.5 c_l on H_l. For l = 4, 10, 22 it gives
  0.3305 vs 0.3281, 0.001485 vs 0.001488, and 0.4781 vs 0.4785.
- The same run confirms the sup bound |AG| <= (1/4)||G||_2 sqrt(24(2l+1)c_l).

## The lemma chain as reconstructed, with attacks

1. **L1 (Lipschitz bound for log r_K).** Take a subgradient xi of the gauge. Then |xi| <= 1/s because
   sB is in K, and the tangential part divided by g(u) is at most sqrt(r^2/s^2 - 1). The proof uses
   one-sided derivatives along great circles and absolute continuity. No smoothness is assumed. OK.
2. **L2 (D3 linear algebra).** Recomputed independently by matrix inversion. OK.
3. **L3 (explicit Kallus Thm 4(ii)).**
   - Brouwer on the simplex D: I checked Psi(D) in D. Ybar = 12 g1 (e^eps - 1) >= 2 g1 s_rho because
     s_rho <= 6(e^eps - 1), and T = cbar s_rho/(1-cbar).
   - Admissibility of the 12 minimal images: <p_i,x_i> = r_i + tau_i and
     gamma_i <= tan gamma_i <= e^eps |Q'_perp x_i|. L1 applied along the arc, together with
     tau_i* = Psi_i(tau*), gives |p_i| >= r_K(p_i/|p_i|).
   - Non-minimal vectors: (1 - ||Q'||)sqrt2 >= e^eps.
   - Determinant: e2 <= 3m^2, e3 <= m^3, and the trace comes from the eutaxy coefficients 1/4.
     I rederived each kappa term: 1/4 eps e^eps; e^eps cbar/(2(1-cbar)); 72 eps e^{2eps};
     288 eps^2 e^{3eps}.
   - The rotation reduction r_{U^{-1}K}(x) = r_K(Ux) is correct, and the hypotheses are rotation invariant. OK.
4. **L4 (normalisation, explicit Kallus Lemma 2).**
   - Verified lambda_{(I+A)K} = lambda_0 o v_A - log|Mu|.
   - e1 is bounded by convexity of (1+alpha)^{-2} - 1 + 2 alpha, whose maximum is at alpha = -a.
     e2 is bounded by the half-spread of M.
   - Recomputed ||pi_2||_{inf->inf} = 5 int_0^1 |P_2| = 10/(3 sqrt 3).
   - The fixed point gives pi_2 lambda = 0. This is pi_2 of log r, not of r as in Kallus; it is used
     consistently in L5 and L6.
   - Rederived the oscillation bound and the Lipschitz bound (composition, plus the gradient of log|Mu|). OK.
5. **L5 (volume).** Jensen on e^{3 lambda}, using mean(lambda) = 0. OK.
6. **L6 (SO(3) duality, the new step).**
   - Measures: sigma and Haar are both probability measures, so int G lambda dU = 6||lambda||_1 and Ux
     is uniform on S^2.
   - The adjoint via the disintegration (F7) is correct. The pull-back of omega(Ux_1) gives
     (1/4)(2 + 8M_{1/2} + 2M_0) = Phi/2: the multiplicities are 1,4,2,4,1, and evenness folds -t onto t.
   - (I+nu)(I-nu+nu^2) = I+nu^3, so the remainder has multiplier (c_l-1)^3 on even l. The correction
     G_l = (1-c_l)^3 g_l/(1.5 c_l) cancels it for l >= 4. What is left lies in H_0 + H_2, which is
     orthogonal to lambda in Z_0.
   - The norm bound uses the reproducing kernel and Cauchy-Schwarz over l, with
     sum ||g_l||^2 <= ||sign lambda||_2^2 <= 1. The series converges absolutely in sup norm.
   - g = sign lambda is only L^inf. M_t, nu and A* are still defined a.e. by Fubini, and the pairing
     identities hold, so this is a gap in presentation only, not in the mathematics. **Valid.**
7. **Final averaging.** A lambda has Haar mean 3 m(lambda) = 0. Averaging the L3 bound over
   Omega = {A lambda < 0} gives -(1/2)||A lambda||_1 + 6 kappa ||lambda||_1, and L6 then needs
   12 kappa N < 1. When lambda != 0, L6 gives |Omega| > 0. The inequality is strict unless lambda = 0,
   that is, unless K is an ellipsoid. OK.

**Consistency.**

- "Lattice packing" means delta_L = phi = vol/(8 d_K) throughout, by (F1) for o-symmetric bodies.
- d_BM uses the multiplicative convention. GL_3 invariance reduces it to B <= K <= (1+delta)B, and
  the single certified run at delta* covers every smaller delta by inclusion.
- The conclusion holds for ALL o-symmetric K in the neighbourhood, not just a perturbation class.
  The proof uses no smoothness, no finite harmonic truncation, and no restriction on the shape of lambda.

## Citation check against Kallus (arXiv v2)

The content is quoted correctly:

- eutaxy coefficients n/|S| = 1/4 (Thm 1(iv));
- the radial-function construction, with (1-eps)B <= K <= (1+eps)B and beta <= 2 sqrt(eps);
- Phi^{-1} bounded on Z via ||K_l||_1 <= sqrt(2l+1), with the sum never evaluated;
- Remark 3, which speculates that the radius is governed by min |c_l|, with c_10 ~ 1e-3;
- the non-symmetric case is left open.

The exact value 65/65536 is PROOF.md's own computation (Kallus gives only ~1e-3), and I verified it.

**The numbering is wrong**, though, for the arXiv source (amsart, separate counters):

| PROOF.md says | Kallus arXiv v2 actually |
|---|---|
| Thm 3(ii) (radial construction) | Theorem 4(ii) (`cns1`) |
| Lemma 4 (normalisation) | Lemma 2 (`quadlem`) |
| Lemma 6 (Phi^{-1} bounded) | Lemma 4 (`philem`) |
| Thm 4 (main theorem) | Theorem 5 (`mainthm`) |
| Thm 1(iv), Remark 3 | correct |

The journal version may number differently. Either way, cite one stated version.

The paper's title is "The 3-ball is a local pessimum for packing". PROOF.md has it right. The title
in the review brief, "Least efficient packing shapes of the 3-dimensional Euclidean space", is wrong.

## Weakest points (none fatal)

1. **Overclaim in "Consequence for Ulam's conjecture".** The claim that no counterexample to Ulam's
   packing conjecture lies within d_BM 1+6e-8 of the ball is proved only for **origin-symmetric**
   bodies. Ulam's conjecture covers all convex bodies, and the non-symmetric case near the ball is
   open (Kallus says so in Remark 1). The claim must be restricted.
2. **Factual slip.** PROOF.md says "the l = 10 term alone is 14097.0". It is 14031.65
   (21 (1-c_10)^6 / (1.5 c_10)); 14097.3 is the whole finite sum. The totals are unaffected.
3. **"Hausdorff form" arithmetic.** (1+3e-8)/(1-3e-8) = 1 + 6e-8 + 1.8e-15 > 1 + 6e-8, so
   "<= 1+6e-8" is false as written. It does no harm, because the certificate passes at 6.6e-8,
   but the sentence should cite that. Also, the sandwich condition is not the Hausdorff metric.
4. Kallus citation numbering (see the table above).
5. **Thin margin.** The margin is about 5%. It rests on the c_l table up to l = 20000 plus the
   Bernstein tail. I reproduced it with an independent cutoff (4000), so it is robust, but if a
   later fix loosens kappa (for example, a missed term), there is little room.
6. L6 should say a word on measurability: sign lambda is only in L^inf, and A* is defined a.e.
7. The novelty sweep relied on OpenAlex and arXiv. Google Scholar and MathSciNet were not checked.

## What I tried that did not break it

- the inequality directions in the det/volume ratio;
- whether vol K >= vol B uses the right normalisation (it does: the mean of log r);
- whether taking pi_2 of log r instead of r matters (it is consistent throughout);
- whether admissibility of the non-minimal lattice vectors was missed (it is checked);
- whether d_BM is conflated with a support-function distance (it is not: everything is radial, and
  B <= K <= (1+d)B is exactly the BM sandwich);
- whether Lemma 6 secretly needs a smooth lambda (it does not);
- whether the averaging over Omega needs |Omega| > 0 (it does, and it holds);
- the equality case (both directions are fine).

To reproduce: `python notes/referee_check.py [delta] [l_cutoff]`. It takes about 2 s and needs mpmath and numpy.

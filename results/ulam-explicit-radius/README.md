# An explicit neighbourhood where Ulam's packing conjecture holds

**Status: computed with rigorous interval arithmetic, refereed (hostile AI referee, independent recomputation). No Lean check. The radius is tiny.**

## Claim
If K is an origin-symmetric convex body in R^3 with Banach-Mazur distance d_BM(K, ball) <= 1 + 6e-8, then its lattice packing density is at least pi/sqrt(18), the density of the ball, with equality only for ellipsoids. So no origin-symmetric counterexample to Ulam's packing conjecture lies within that distance of the ball. Non-symmetric bodies are not covered.
The radius 6e-8 is tiny: the certificate passes at 6.6e-8 and fails at 6.7e-8. It is of interest as the first explicit value, not as a practical range.

## Prior work
Kallus, "The 3-ball is a local pessimum for packing" (Adv. Math. 264 (2014), arXiv 1212.2551) proves this qualitatively ("there exists epsilon > 0"). We found no earlier explicit radius (citation sweep of 2026-10-24 in `PROOF.md`; Google Scholar was not accessible, so coverage of non-arXiv journals relies on OpenAlex).

## Method
Make Kallus's argument quantitative: Lipschitz bound for the log radial function, an explicit Brouwer-ball radius, spherical-harmonic estimates, and a final inequality 12 * kappa * N_eff < 1. The numerical steps are certified in arb interval arithmetic by `certify.py`. Details, and the obstacle to larger radii, are in `PROOF.md`.

## How to verify
```
pip install python-flint
python certify.py        # about 1 s; expect "ALL CHECKS PASSED"
```
`sanity.py` (numpy, scipy) is a Monte Carlo sanity check, not part of the proof. `notes/referee_check.py` (mpmath) is the referee's independent recomputation; the report is `notes/REFEREE.md`.

## Caveats
- Proof is by hand plus certified numerics; no formal verification.
- Relies on Kallus's framework and on Hales's proof for the ball.
- Numbering of Kallus's results in PROOF.md follows an earlier reading; the mapping is in a note near the top of PROOF.md.

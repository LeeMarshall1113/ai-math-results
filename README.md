# ai-math-results

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23748157.svg)](https://doi.org/10.5281/zenodo.23748157)

**Thirteen computer-assisted mathematical results, each shipped with its certificates and a checker you can run yourself.**

By Lee Marshall, produced with AI assistance (Claude, Anthropic). The rule of this repository: you should not have to trust us. Every folder states exactly what is proved, what earlier work exists, and what the check does and does not cover.

**Not a mathematician?** Try the [interactive page](https://leemarshall1113.github.io/ai-math-results/) (four of the problems as things you can play with), or read [In plain words](PLAIN-LANGUAGE.md): what each result says, why it mattered, and why a computer-checked proof can be trusted.

**Status at a glance:** 2 results reviewed by a domain expert · 2 with Lean 4 components · 6 with proofs checked by the formally verified checker `cake_lpr` · none peer reviewed by a journal yet.

## Start here: check one in a minute

```
git clone https://github.com/LeeMarshall1113/ai-math-results
cd ai-math-results/results/ulam-explicit-radius && python certify.py          # about 1 s
cd ../thickness-k14/code && bash check_all.sh                                  # about 1 min
cd ../../polygon-n4/proof/indep3 && python -I check3.py                        # about 5 min, standard library only
```

## Headline results

### Hadwiger covering in ℝ³: 13 copies suffice
Every convex body in ℝ³ is covered by 13 smaller homothetic copies, so H₃ ≤ 13. The previous published bound was 14 (Prymak, 2023). Also corrects two points in the H₃ ≤ 14 argument.

- **Evidence:** exact LP certificate (7.56 million leaves) checked by two separately written exact checkers; the reduction is formalised in Lean 4 with the certificate as an axiom.
- **Check:** `cd results/hadwiger-h3/code && bash check_sample.sh` (about 4 min; the full check is about 45 CPU-hours with the Zenodo data).
- **Links:** [folder](results/hadwiger-h3) · [standalone repository](https://github.com/LeeMarshall1113/hadwiger-h3) · [DOI 10.5281/zenodo.24451161](https://doi.org/10.5281/zenodo.24451161)

### Regular heptagons: the double lattice is the densest packing
Every packing of congruent regular heptagons in the plane has density at most the Kuperberg double-lattice density (2/97)(−111 + 492u − 356u²), u = cos(π/7), about 0.8927. This was the conjectured optimum; local optimality is due to Kallus–Kusner, and the method follows Hales–Kusner for the pentagon.

- **Evidence:** rigorous interval branch-and-bound (about 125.6 million boxes, none stuck) plus a certified tube of 1,530 charts; two hostile AI audits. Not formally verified. The certificates are rerunnable logs, not a stored per-box proof; the folder's caveats say exactly which slices were rerun from scratch.
- **Check:** `cd results/heptagon-packing`, then see the README (one slice takes minutes to hours).
- **Links:** [folder](results/heptagon-packing) · [standalone repository](https://github.com/LeeMarshall1113/heptagon-packing) · [DOI 10.5281/zenodo.24451167](https://doi.org/10.5281/zenodo.24451167)

### No two-intersection set of type (2,5) in PG(3,9)
No set of points of PG(3,9) meets every line in exactly 2 or exactly 5 points (so no set of type (2,5) or (5,8) with respect to lines). Stated as the first open case by Gavrilyuk–Matkin (2018) and Ihringer (2025).

- **Evidence:** SAT with cube splitting; every refutation has an LRAT proof checked by `cake_lpr`; full local regeneration; adversarial review by a second AI system. Not compared with Penttila–Royle's planar classification.
- **Check:** `cd results/pg39-type25`, then the replay scripts in the README (about 11 min, one thread).
- **Links:** [folder](results/pg39-type25) · [standalone repository](https://github.com/LeeMarshall1113/pg39-type25) · [DOI 10.5281/zenodo.24451165](https://doi.org/10.5281/zenodo.24451165)

### The Agaev–Chebotarev polygon conjecture, n = 4
Λ₄ = S₄. Together with n ≤ 3 (2005) and Chebotarev's disproof for n ≥ 5 (2026), this settles the conjecture. An earlier computer-assisted proof by D. Z. Liu is credited in the folder.

- **Evidence:** reviewed by a domain expert. Exact certificates for 836 face inequalities; two independent full replays (one using only the Python standard library) plus a partial third check.
- **Check:** `cd results/polygon-n4/proof && bash fullcheck.sh`, or `cd indep3 && python -I check3.py` (about 5 min).
- **Links:** [folder](results/polygon-n4) · [DOI 10.5281/zenodo.23748158](https://doi.org/10.5281/zenodo.23748158)

### Collatz: no cycle with 111 or fewer local minima
No nontrivial positive Collatz cycle has 111 or fewer local minima. An extension to 133 is computed (two independent codes per table value, two chain implementations per m) but not formally verified.

- **Evidence:** reviewed by a domain expert. Lean 4, kernel-checked, conditional on one published input (Barina's verification to 2^71), stated in the folder.
- **Check:** `lake build Collatz && lake env lean scripts/Statement.lean` in [collatz-m-cycles](https://github.com/LeeMarshall1113/collatz-m-cycles) (about 9–10 h).
- **Links:** [folder](results/collatz-m-cycles) · [DOI 10.5281/zenodo.23748153](https://doi.org/10.5281/zenodo.23748153)

## More results

| Result | Evidence | Check | Folder |
|---|---|---|---|
| **Karp's evasiveness conjecture, n = 10 (partial).** A non-evasive monotone graph property on ten vertices must be of Angel–Borja type 5; eight of the nine types are excluded, four of them new. Type 5 is open, so n = 10 is not settled | Exact finite checks on explicit permutation groups, conditional on cited theorems of Smith and Oliver; replayed by the authors only | about 2 min, standard library | [karp10-partial](results/karp10-partial) |
| **No torus with 10 convex quadrilateral faces.** The F = 10 case of Grünbaum–Szilassi Conjecture 3.2; with Seamons's 11-face torus (2026) it settles which face counts occur | LRAT proofs checked by `cake_lpr`; two independently written encoders agree | about 2 min | [torus-quad10](results/torus-quad10) |
| **Geometric thickness of K₁₄ is 4.** With the Dillencourt–Eppstein–Hirschberg upper bound (2000), this completes θ̄(Kₙ) for n ≤ 16 | LRAT proofs checked by `cake_lpr`; two encoders with separate case trees agree | about 1 min (full: 3.5 h) | [thickness-k14](results/thickness-k14) |
| **Lonely runner, tight instances for 6 and 7 speeds.** Classifies every speed set attaining the bound; answers Perarnau–Serra's Problem 1 for 7 and 8 runners and gives Kravitz's Conjecture 8.3 for 7 and 8 speeds | Strict Malikiosis–Santos–Schymura bound (proved in full) plus prime gates; two independent implementations agree exactly | about 2 h | [lonely-runner-tight-6-7](results/lonely-runner-tight-6-7) |
| **Erdős problem #580 (vertex form)** for n ≤ 21; n = 20 and 21 are new | LRAT proofs checked by `cake_lpr`, plus cited theorems | quick mode; full run about 6 CPU-hours | [erdos-580](results/erdos-580) |
| **Ulam's packing conjecture near the ball.** Holds for origin-symmetric convex bodies within Banach–Mazur distance 1 + 6·10⁻⁸ of the ball | Interval arithmetic | about 1 s | [ulam-explicit-radius](results/ulam-explicit-radius) |
| **Lonely runner, 8, 9 and 10 runners.** Formally checked certificates for cases already proved in the literature (not new theorems) | LRAT proofs checked by `cake_lpr` | see folder | [lonely-runner-8-10](results/lonely-runner-8-10) |
| **Turán density of K₄⁻:** π(K₄⁻) < 0.28674302. Since superseded by Lidický–Volec (0.2865, arXiv 2610.33850); kept as an exact, checkable certificate | Exact rational flag-algebra certificate | about 25 min | [turan-k4minus](results/turan-k4minus) |

Exact commands are in each folder's README.

## How a result gets here

1. **Novelty search** of the literature, arXiv, Zenodo and GitHub. Earlier or overlapping claims are credited in the folder.
2. **Certificates, not trust.** The search code produces a certificate (an LRAT proof, an exact LP solution, an interval log, a Lean proof). A separate checker verifies it; where possible a second checker is written independently.
3. **Replay from scratch** on a clean checkout.
4. **Hostile review** by a second AI system instructed to break the argument. Reports are included where they exist.
5. **Plain statement of limits.** Each README says what is conditional, what was not checked by machine, and what a reader must still take on trust.

## What each folder contains

- the precise statement and its sources;
- prior work, including any earlier or overlapping claims;
- the method;
- how to check it, and what the check does and does not cover.

## Citing

Use the Zenodo DOI above, or the "Cite this repository" button (from `CITATION.cff`), and name the result folder. Each GitHub release is archived on Zenodo with its own version DOI.

## Corrections

If you find an error, a gap, or earlier work we missed, please open an issue. Earlier claims we learn about are credited in the relevant folder.

## AI disclosure

The work in this repository was generated with AI (Claude, from Anthropic), directed by Lee Marshall. Each result was checked by independent verifiers and attacked by hostile AI referees, whose reports are included where they exist. None of it has been through human peer review unless a result's README says so. Please check the certificates yourself where a check is provided.

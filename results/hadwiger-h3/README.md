# Every convex body in ℝ³ is covered by 13 smaller homothetic copies

**Status: computer-assisted. Exact LP certificate (7,561,583 leaves) checked by two separately written exact
checkers; the reduction from the certificate to the theorem is formalised in Lean 4 with no `sorry`; hostile AI
referee, two rounds; no human peer review.**

## Statement

For a convex body K ⊂ ℝⁿ, C(K, int K) is the least number of translates of int K that cover K. Hadwiger's
covering number H_n is the maximum of C(K, int K) over convex bodies in ℝⁿ; by Boltyanski it equals the
illumination number. The Levi–Hadwiger–Gohberg–Markus conjecture is H_n = 2ⁿ.

**Theorem.** Every convex body K ⊂ ℝ³ is covered by 13 translates of its interior, and hence by 13 homothets
x_i + rK with 0 < r < 1 (r depending on K). So H₃ ≤ 13.

The previous bound was H₃ ≤ 14 (A. Prymak, SIAM J. Discrete Math. 37 (2023) 17–24, arXiv:2112.10698).
The lower bound is 8.

## Credit and prior work

The proof is built on Prymak's framework and uses his reduction, his bodies Q_P and his exact-LP certificate
method throughout. That framework in turn rests on I. Papadoperakis (Geom. Dedicata 75 (1999) 275–285,
H₃ ≤ 16): the minimal circumscribed parallelotope, the paired contact points (his Lemma 3(a)) and the reduction to
covering the cube's edges and contact points (his Lemma 4). We have not seen Papadoperakis's paper; we use it as
described by Prymak and by Prymak–Shepelska, and re-prove the two lemmas we use (`PROOF.md` §§1–2, and in Lean).

Other work in this line: M. Lassak (H₃ ≤ 20, 1998); A. Prymak and V. Shepelska (Papadoperakis's method in
dimensions 4–6, J. Geom. 2020); F. Diao (dimensions 5–6); A. Arman, A. Bondarenko and A. Prymak (small
dimensions, Canad. Math. Bull. 2025). Prymak's data and scripts are at github.com/andriyprm/h3atmost14.

## What is new

1. **Lemma M (long chord).** If [0,1]³ is a minimal-volume parallelotope containing K, then for each axis K
   contains a segment of vector τ(e_j ± e_k), for one of the two signs, with τ = 1/2 + (√2 − 1)/6 = 0.5690….
   One of the 13 covering bodies is a Minkowski combination of Prymak's Q_P and this chord. This removes the
   configuration p = (1/2)⁶ where Prymak's covering bodies need 14 (his Remark 2.3).
2. **A corrected symmetry domain and 570 extra start boxes** (see below).
3. **Strict containment** in every leaf, in place of Prymak's Lemma 2.4(ii) (see below).

## Two points in the H₃ ≤ 14 argument

Both came up while checking this work. We state them because the same fixes are needed here.

* **The region D is not a fundamental domain.** Prymak reduces by the 48 cube symmetries to
  D = {p₁, p₂, p₃ ≤ 1/2, p₂ + p₄ ≤ p₁ + p₆ ≤ p₃ + p₅}. The configuration
  p = (459/500, 1481/5000, 1689/10000, 9529/10000, 9359/10000, 3079/10000) has none of its 48 images in D, nor in
  the closed union of the start boxes of his certified region U^D_20. About 2.7% of random configurations behave
  this way. A correct domain is D\* = {p₁, p₂, p₃ ≤ 1/2, |s_z − 1| ≤ |s_y − 1| ≤ |s_x − 1|}, with s_x = p₃ + p₅,
  s_y = p₁ + p₆, s_z = p₂ + p₄. We do not use D\*
  directly: we certify one 1/10 box from each of the 570 grid orbits that U^D_20 misses, which with the old
  39,380 boxes meets every orbit.
* **Lemma 2.4(ii) fails at some box corners.** It states that some translate of Q_P lies in int O_p for p in
  the box P. For P = ∏[a_i, a_i + 1/160], a = (0, 0, 63/160, 5/8, 0, 57/80), and its corner
  p = a + (1, 1, 0, 0, 1, 0)/160, the body Q_P touches facets of O_p whose inner normals have 0 in their convex
  hull, so no translate of Q_P lies in int O_p. We check every leaf with strict inequalities instead; then
  int Q_P ⊆ int O_p (from Lemma 2.4(i)) is enough.

These fixes also repair the H₃ ≤ 14 argument: the 570 extension boxes are certified with plain leaves (13
translates of int Q_P, no chord), inside Prymak's framework, and strict containment replaces Lemma 2.4(ii) for any
leaf that satisfies its inequalities strictly. We have not re-run Prymak's own certificate under a strict test.
Details in `PROOF.md` §11.

## How to verify

**Requirements.** Python 3 with numpy and scipy. Optional: elan, for the Lean project.

**Quick check, no download** (about 4 min, one thread, under 1 GB):

```
cd results/hadwiger-h3/code && bash check_sample.sh
```

It checks both counterexamples above exactly, derives D\* and the 570 extension boxes, checks that every orbit of
1/10 boxes contains a certified start box, and runs both leaf checkers on `certs/sample/`: the start box at the
centre (1,771 leaves, 118 of them chord leaves) and six other start boxes, copied verbatim from the certificate.

**Full check** (about 45 CPU-hours; resumable; `WORKERS` processes of about 1.5 GB each):

1. Download `hadwiger-h3-fullcert.tar` (2.2 GB) and `hadwiger-h3-fullcert_ext.tar` (26 MB) from the Zenodo
   dataset record, doi:10.5281/zenodo.XXXXXXXX.
2. Run

   ```
   cd results/hadwiger-h3/code && DATA=/path/to/archives WORKERS=8 bash check_all.sh
   ```

It checks the archive and chunk hashes against `certs/SHA256SUMS.*`, then runs `certcheck_fast.py` (exact,
non-strict, with a 1% cross-check by the rational reference `certcheck.py`), `certcheck_ext.py` (orbit coverage),
`certcheck_strict.py` (the Lean leaf predicate), the independent replay `replay/hcheck.py --strict`, and the quick
check. `LEAN=1` also builds the Lean project.

**Lean.** `cd lean && lake exe cache get && lake build && lake env lean Hadwiger/Audit.lean`.

**Recorded results** (2026-10-31, `certs/`): `CERTCHECK_FAST_VERDICT.txt`, `CERTCHECK_STRICT_VERDICT.txt` and
`CERTCHECK_EXT_VERDICT.txt` read COMPLETE AND VALID; `CERTCHECK_STRICT_EXT.log` (2026-11-01) reads STRICT VALID
for the 87,828 extension leaves; `replay_strict_SUMMARY.json` records 7,561,583 leaves,
0 failures, 39,950 of 39,950 start boxes tiled, 0 orbits without a certified box.

## Trust base

* **The certificate data** (on Zenodo; hashes here). Every leaf is checked by two exact checkers with separate
  code: `newmethod/certcheck_strict.py` and `replay/hcheck.py`. Neither uses the certifier's code.
* **The checkers.** `certcheck_fast.py`, `certcheck_strict.py` and `certcheck_ext.py` share `certcheck.py`'s
  translation from a stored structure to constraints, and `hw.py`. `hcheck.py` was written from the specification
  only, ignores the stored structures, and shares no code with them.
* **The link from data to Lean.** The Lean axiom `all_leaves_pass` states that every certified start box is tiled
  by leaves satisfying `LeafOK`. That this is what `certcheck_strict.py` checks is established by reading both
  (index conventions are listed in `lean/README.md`).
* **Lean.** Besides `all_leaves_pass`, the main theorem uses one `native_decide` computation (finite orbit
  coverage of the 1/10 grid, also checked in Python twice) and Lean's three standard axioms. Lemma M,
  Papadoperakis's Lemmas 3(a) and 4, existence of the minimal box, Prymak's Lemma 2.4(i), leaf soundness and the
  symmetry reduction are all proved.

## Limits

* No human has refereed this work.
* The JSON-to-Lean correspondence is checked by reading code, not formally.
* Papadoperakis 1999 has not been read in the original. Nothing in the proof depends on it (both lemmas are
  re-proved), but the attribution of those lemmas rests on Prymak's account.
* Lemma M is elementary and may be known; we did not find it.
* The novelty search (`NOVELTY.md`) cannot rule out unindexed work.

## Files

* `PROOF.md`: the argument.
* `NOVELTY.md`: the searches.
* `CHECKLIST.md`: release checklist.
* `code/`: the checkers (`newmethod/certcheck*.py`, `hw.py`), the independent replay (`replay/`), the
  counterexample scripts (`corner_example.py`, `symmetry_gap.py`), `newmethod/dstar_domain.py`, the certifier
  for reproducibility (`newmethod/certify2.py`, `newmethod/fast/certify2f.py` and their modules), and
  `check_sample.sh` / `check_all.sh`. The checker files are byte-identical to the ones that produced the recorded
  verdicts (their SHA-256 is part of each ledger entry), so they still contain some local paths from the machine
  that ran them; these only take effect if those paths exist.
* `lean/`: the Lean 4 project.
* `certs/`: SHA-256 sums of the archives and of every chunk file, the start-box lists, the verdict files, the
  production policies, the replay summary, and `sample/`.
* `zenodo_note/`: LaTeX source and PDF of the technical note, doi:10.5281/zenodo.XXXXXXXX.

## References

* I. Papadoperakis, An estimate for the problem of illumination of the boundary of a convex body in E³,
  Geom. Dedicata 75 (1999) 275–285, doi:10.1023/A:1005056207406.
* A. Prymak, A new bound for Hadwiger's covering problem in E³, SIAM J. Discrete Math. 37 (2023) 17–24,
  doi:10.1137/22M1490314, arXiv:2112.10698.
* A. Prymak and V. Shepelska, On the Hadwiger covering problem in low dimensions, J. Geom. 111 (2020) 42,
  arXiv:1811.08962.
* A. Arman, A. Bondarenko and A. Prymak, On Hadwiger's covering problem in small dimensions, Canad. Math. Bull.
  (2025), arXiv:2404.00547.
* M. Lassak, Covering a three-dimensional convex body by smaller homothetic copies, Beiträge Algebra Geom. 39
  (1998) 259–262.
* K. Bezdek and M. A. Khan, The geometry of homothetic covering and illumination, in Discrete Geometry and
  Symmetry, Springer (2018) 1–30.

## AI disclosure

This work was generated with AI (Claude, from Anthropic), directed by Lee Marshall. It was checked by two
separately written exact checkers and a Lean formalisation, and attacked by a hostile AI referee. It has not been
through human peer review. Please check the computations yourself with the checkers provided.

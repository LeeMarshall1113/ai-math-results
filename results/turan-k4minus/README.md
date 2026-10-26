# Turan density of K4-: pi(K4-) < 0.28674302

**Status: exact certificate, refereed (hostile AI referee, independent reproduction). No Lean check.**

## Claim
pi(K4-) <= 14337150771693633 / (5 * 10^16) = 0.28674302..., where K4- is the 3-graph on 4 vertices with 3 edges.
- Previous best upper bound: 0.286889, Falgas-Ravry and Vaughan (2012, flag algebras).
- Conjectured value: 2/7 = 0.285714 (Frankl-Furedi construction).
- Improvement 1.46e-4, about 12% of the gap to 2/7.

## Prior work
Falgas-Ravry-Vaughan 2012 (flag algebras). The novelty sweep of 2026-10-24 (arXiv, the authors' listings of arXiv 2609.27495, web, GitHub) found 0.2868 as the best bound everywhere; see `notes/NOVELTY.md`. Searches can miss things.

## Method
1. A plain 7-vertex flag-algebra certificate reproduces the old bound exactly: pi <= 0.2868885176 (`certs/cert_plain.json.gz`).
2. Degree regularity: an extremal limit object can be taken with all vertex degrees equal to the edge density (a min-degree deletion lemma). The resulting constraints are linearised on the density interval [0.2867, 0.28689] and the SDP is solved with a floor (`certs/cert_C.json.gz`). That certificate shows no density in the interval exceeds 0.28674302; the plain one excludes densities above 0.28689.
3. Integer Gram factors and exact rational arithmetic.

## How to verify
From this folder (pure Python 3, standard library only; about 25 minutes on one thread, peak about 0.85 GB):
```
python check_bound.py 0.2867433 certs/cert_plain.json.gz certs/cert_C.json.gz
```
It checks completeness of the catalogue of the 8157 K4--free 7-vertex 3-graphs (cross-checked against brute force at 5 and 6 vertices), positive semidefiniteness of the integer Gram factors, every graph's inequality, and that the intervals cover the range. Expected last line: `BOUND 2867433/10000000 ESTABLISHED`. The recorded run is `certs/check_bound.out`.
`check_bound.py` has no dry-run mode.

## Files
- `REPORT.md`: the full report. `notes/REFEREE.md`: hostile referee report. `notes/NOVELTY.md`: novelty search.
- `notes/indep_check.py`, `notes/indep_catalogue.py` (need numpy): the referee's independent re-implementation, output in `notes/referee_indep.out`; `notes/referee_rerun_check_bound.out` is the referee's rerun of the checker.
- `k4enum.py`, `flaglib.py`, `classes.pkl`: catalogue generator and flag-algebra library from the search. They are not needed by `check_bound.py`, which rebuilds and checks the catalogue itself. Both scripts need numpy.
- The solver (`sdp_ipm.py`) and its scratch output are not included.

## Caveats
- No formal (Lean) verification. Correctness rests on the exact checker, the degree-regularity lemma (standard, argued in the referee report), and the referee's independent reproduction.
- The bound could gain roughly 3e-5 by bisecting the interval further; not done.
- Still far from the conjectured 2/7.

# Turán density of K₄⁻: upper bound < 0.28674302 (candidate)

Written 2026-10-24 by the coordinator from the lane agent's report; commits 609e86f, d1846f5 and 99879aa. **Status: REFEREED (SURVIVES, notes/REFEREE.md; independent reproduction of both exact maxima); novelty NEW (notes/NOVELTY.md).** There is no Lean check.

## Claim
π(K₄⁻) ≤ 14337150771693633 / 5·10¹⁶ = 0.28674301543…

K₄⁻ is the 3-graph on 4 vertices with 3 edges.
- **Previous best upper bound:** 0.286889, Falgas-Ravry–Vaughan (2012, flag algebras).
- **Conjecture:** 2/7 ≈ 0.285714, Frankl–Füredi construction.
- **Improvement:** 1.46·10⁻⁴, which closes about 12% of the gap to 2/7.

## Method
1. **Plain 7-vertex flag-algebra certificate.** This reproduces the previous bound in exact form: π ≤ 0.2868885176 (certs/cert_plain.json.gz).
2. **Degree regularity.** By a min-degree deletion lemma, the limit object of an extremal sequence can be taken with every vertex degree equal to the edge density. We add the resulting linear constraints, using 7 regularity rows. The constraint is linearised on the density interval [0.2867, 0.28689], and the SDP is solved with a λ floor (`sdp_ipm.py --floor`).
   - The certificate on that interval (certs/cert_C.json.gz) shows that no density in the interval is attainable above 0.28674302.
   - The plain certificate excludes densities above 0.28689.
   - Together these give the claim.
3. **Exact rounding.** Integer Gram factors and exact rational arithmetic.

## How to check it
```
python check_bound.py 0.2867433 certs/cert_plain.json.gz certs/cert_C.json.gz
```
This is pure Python with Fractions: about 25 min on 1 thread, peak 0.85 GB. The checker verifies:
- completeness of the catalogue of the 8157 K₄⁻-free 7-vertex 3-graphs (cross-checked by brute force at 5 and 6 vertices);
- that the Gram factors are positive semidefinite (integer factors);
- every graph's inequality;
- that the intervals cover the range.

The recorded output is certs/check_bound.out, which ends "BOUND 2867433/10000000 ESTABLISHED".

## Novelty
Clear as of 2026-10-24:
- arXiv;
- the listings of the authors of the 2026 tetrahedron paper (arXiv 2609.27495);
- the web and GitHub.

All of them still show 0.2868 as the best bound.

## Not done
- A Lean check. The tetrahedron paper's Lean build needed 12.5 h, and one module needed 36.7 GiB, which is over this machine's approval line.
- Further bisection of the interval, which might gain about 3·10⁻⁵.

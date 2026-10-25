# Novelty check: pi(K4-) <= 0.2867430154 (checked 2026-10-24, web search only)

Verdict: **NEW, with caveats** (no published bound below 0.286889 found; see "Not reached").

K4- = {123,124,134}. Ordinary Turan density. Not codegree (pi_2) / min-degree (pi_1) variants.

## Sources and what they say about the K4- upper bound
| Source | Says |
|---|---|
| Razborov, "On 3-hypergraphs with forbidden 4-vertex configurations" (SIAM J. Discrete Math 2010) | pi(K4-) <= 0.2978 (as cited by later papers) |
| Baber-Talbot, "Hypergraphs do jump" https://arxiv.org/pdf/1004.3733 | pi(K4-) <= 0.2871 |
| Baber-Talbot, "New Turan densities for 3-graphs" https://arxiv.org/abs/1110.4287 | Same group; 0.2871 is the figure later papers credit to Baber-Talbot |
| Falgas-Ravry-Vaughan, "On applications of Razborov's flag algebra calculus to extremal 3-graph theory" https://arxiv.org/abs/1110.1623 (Flagmatic) | Source credited for 0.286889 (not machine-read here; PDF was compressed) |
| Falgas-Ravry-Vaughan, "Turan H-densities for 3-graphs" https://arxiv.org/abs/1201.4326 (EJC 19(3) P40) | Other K4- results (pi(K4-,F3,2)=5/18 etc.); no better ordinary bound found |
| "Four-vertex traces of finite sets" https://arxiv.org/abs/2301.05830 | 2/7 <= pi(K4-) <= 0.286889 (FRV, flag algebras) |
| "Turan densities of hypergraph augmentations" https://arxiv.org/abs/2610.11757 (2026) | Best known pi(H_3^3 = K4-) <= 0.2868..., FRV; mentions no later improvement. Lower bound 2/7 (Frankl-Furedi, Giraud), conjectured sharp |
| Jeong-Park-Im-Lee-Yang, tetrahedron paper https://arxiv.org/abs/2609.27495 | Tetrahedron only (<0.557789). 7-vertex certificate with "degree-stationarity" from Razborov's differential method, Lean 4 proof, GitHub code. Never mentions K4-. Notes Baber's 0.5615 used degree-regularity constraints (tetrahedron). |
| Pikhurko et al. codegree K4 paper https://wrap.warwick.ac.uk/id/eprint/174832/ | Codegree variant; different quantity |
| Lean flag algebras https://arxiv.org/abs/2607.23500 | Search hit only; not read for K4- |

## Findings
- Every post-2012 source reachable (incl. 2026 ones) still quotes 0.286889 / "0.2868..." as the record. No source gives a K4- bound below it.
- No hit for min-degree / degree-regularity applied to ordinary pi(K4-). The only degree-regularity use found is Baber's for the tetrahedron, and the 2026 tetrahedron paper's degree-stationarity (same idea, other graph).
- Technique is not new in general (Razborov differential method; Baber thesis for K4). The *result* for K4- looks new.
- Searches for 0.2867/0.28674 returned nothing.

## Not reached (residual risk)
- Baber's UCL thesis full text (2011): not fetched. It could contain K4- computations with degree constraints; check before claiming. Baber-Talbot's published 0.2871 suggests it did not beat 0.286889 (thesis predates FRV).
- Flagmatic docs/examples (maths.qmul.ac.uk/~ev/flagmatic) and GitHub code search: no repo hit; not browsed.
- Google Scholar citation list of FRV 2012, Keevash survey full text, Glebov-Kral-Volec papers: not reached directly. PDFs of 1110.1623 / 2301.05830 were unreadable via fetch (compressed); attribution of 0.286889 rests on secondary text.
- Recommend: pull 1110.1623 and Baber's thesis locally and grep for "K_4^-" and "regular" before publishing.

## Coordinator addendum (2026-10-25)
- **Baber's thesis** ("Some results in extremal combinatorics", UCL 2011; discovery.ucl.ac.uk/1306175) predates Falgas-Ravry–Vaughan (FRV). arXiv 2301.05830 (2023) and arXiv 2610.11757 (2026) both cite FRV's 0.286889 as the record. A lower thesis bound would have been cited, so the thesis is very unlikely to beat it.
- **arXiv 1110.1623 (FRV / Flagmatic)** gives K₄⁻ results only together with extra forbidden graphs (π(K₄⁻, C₅, F₃,₂) = 12/49; π(K₄⁻, F₃,₂) = 5/18), not a sharper plain bound.
- **Verdict: NEW**, on the strength of the 2026 citation of 0.2868 as current. Residual risk: an unindexed note or talk.

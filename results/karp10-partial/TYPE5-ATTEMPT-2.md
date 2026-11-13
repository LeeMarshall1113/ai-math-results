# Karp n = 10, type 5: second attempt (2026-11-11, Eastern)

Author: Claude (Fable 5.2) for Lee.  One bounded fresh attempt.  Local commits only; nothing pushed; nobody contacted; no sub-agents; no GPU; no cloud.
Code and logs: `stage5/`.  Status words used below: PROVED (paper argument), CERTIFIED (replayable computation), SAMPLE (checked on a sample), SOLVER (solver/heuristic fact, not replayed), CONJ (belief).

## 0. Reading and due diligence (before any computation)

* Read `PARTIAL-THEOREM.md`, `NOTES.md` (stages 1-4), `note/NOTES-FOR-LEE.md`, the Karp entries of `math-comms/COMMS.md`.
* Re-ran `python replay/replay_all.py` on the present tree: manifest OK, 23 rows, **ALL VALID** (`stage5/logs/replay_all_check.log`).  CERTIFIED (unchanged).
* Re-derived by hand: the type weights (+4, -2, -4, +4, -1 and negatives) and that each of the nine listed ideals has weight sum 1; the type-4 argument for the order-192 group (containments K_{1,1,8} <= G2bar, K_{4,4} u K_2 <= K_{5,5}, 2C_5 <= K_2 v K_{4,4}, 2K_5 <= K_2 v 2K_4, equation 2 + x - y = 1 with y <= x).  No error found in the eight exclusions.  (Not a full independent re-implementation.)

## 1. Triage of the five directions (thinking, before computing)

The log below records, per direction, the size estimate, the stopping point fixed in advance, and the outcome.

### Outcome in one paragraph

**Type 5 is NOT killed, and no counterexample candidate was found.**  New in this attempt: (i) a small CERTIFIED lemma showing that the certified facts about P give an adversary nothing (section 4); (ii) a first test of explicit closed-form candidates for a near-counterexample, with a negative SOLVER result (section 5); (iii) a quantified statement of the two walls (sections 2-4) and of the one rung that is still affordable but was not run (section 6).  Compute used: about 0.4 CPU-hours, all in jobs under 1 GB and under 10 minutes (no lease needed; the broker budget was fully leased, 15.0 of 15.0 GB, and commit charge was at 87 % when I checked at about 11:30 ET, so nothing larger could have started for about 2.5 hours anyway).

## 2. Directions 1 and 2 (more groups; exact linear algebra instead of clauses)

**What the method can ever see (standard theory, cited not re-proved).**  For a Z-acyclic S_10-complex the fixed-point constraints in use are: F_p-acyclicity of the fixed complex of every p-subgroup (Smith), chi = 1 for G^1, chi = 1 mod q for G^q (Oliver), and the same for the minors lk(e), dl(e) that the decision tree is known to visit.  By Oliver's realisation theorems (as used in `NOTES.md` 12.8) nothing else follows for abstract complexes.  The fixed-point programme is therefore a finite ladder indexed by the subgroup classes, and the bottom rung (trivial subgroup) is the statement "Delta_P is Z-acyclic" itself.  So the programme can end in only two ways: a contradiction on some rung, or a Z-acyclic S_10-invariant complex between Delta_D and Delta_D*, after which no fixed-point argument exists and evasiveness of type 5 needs input that is not topological.

**Where the ladder stands (from stages 1-4, re-read, not re-run).**  Every rung with at most 16 edge orbits passes jointly in all four Z5 branches (Euler, link, deletion, relative Smith, depth-2 adjacent, Euler of eight further 2-groups).  The Euler defect of D is nonzero for only 24 classes (uniform cycle types), all already imposed except the trivial group.  Every remaining rung has zero defect.

**Groups not yet tried, with sizes (the estimate asked for before any push).**  The untested Smith rungs are the p-groups with many orbits.  Cyclic ones, with k = number of edge orbits and 2^k = number of invariant graphs: 2^5 (sigma) k = 25; 2^4 1^2 k = 25; 2^3 1^4 k = 27; 3 1^7 k = 29; 2^2 1^6 k = 31; transposition k = 37; trivial k = 45.  Non-cyclic: class 8 (3^2 1^4, k = 19), the eight 2-groups with k = 17-19, 14 more p-groups with k <= 19.  Non-p-groups mixing 3 and 5, or groups acting at depth 2-3, add only Euler equations or congruences; those all have zero defect, and stage 2 found the X-part defect of all 7,632 depth-2 and 4,928 depth-3 equations to be 0, so I did not rebuild them.

* sigma rung, cost.  Measured in stage 3: one full F_2-homology evaluation of the sigma complex takes 128 s on one thread; delta debugging needs a few hundred evaluations per clause; the library needed about 5,500 clauses for 1,049 complexes of at most 2^16 faces.  With at least 10^4 clauses for sigma: about 10^4 x 300 x 30 s = 25,000 CPU-hours by delta debugging; about 10^4 x 2 x 128 s = 700 CPU-hours if every clause comes from one cycle/cocycle pair (design in section 6).  Central estimate 2,000 CPU-hours, range 700-25,000.
* The rungs 2^2 1^6, 3 1^7, transposition and trivial have 2^29 to 2^45 faces.  One homology evaluation of a complex with 10^11 cells does not fit in this machine (31 GB), and the S_10-isotypic reduction of the trivial rung still has boundary matrices with about 10^8-10^9 rows (roughly 10^6 classes per edge count times irreducible dimensions up to 768).  Beyond 10^5 CPU-hours.  **WALL 1.**
* Expected outcome of the affordable rungs: satisfiable.  Reason (heuristic, CONJ): the lower the rung, the larger the share of faces whose stabiliser is exactly the group in question; these enter only that rung, so each rung has its own free generators to absorb its conditions.  This is what happened on every rung so far (the solver converges in 10-25 iterations once the rank-profile inequalities are present).

**Decision.**  Not run: the sigma rung and the Smith conditions of the eight 2-groups.  The command for the latter (branch 0, about 20 CPU-hours, 5 GB) is in section 6.  I did not queue it: my estimate of an empty core there is 2-3 %, the machine was saturated, and "satisfiable" there changes nothing in the conclusion.

## 3. Direction 3 (direct cube-and-conquer encoding)

Size estimate, made before building anything.  "There is a monotone P in the type-5 ideal with a decision tree of depth <= 44" quantifies over (a) the down-sets X of the poset of about 11.7 million CC classes and (b) a strategy on the positions (yes-graph, no-graph), of which there are 3^45 / 10! = 8.1 x 10^14 up to isomorphism.  Even for ONE fixed P the exact game solve needs about 10^15 position evaluations (more than 2 x 10^5 CPU-hours at 1 microsecond each, and 10^14 bytes of table).  "Collapsible" needs an acyclic matching on up to 2^45 faces; "Q-acyclic" is a rank condition on the same faces, not a CNF.  The only parts of the problem that fit a SAT solver are the fixed-complex systems of stages 2-4, and those are SATISFIABLE, so there is nothing for a cube cover to refute.  The pg39 pipeline is the right tool only after some rung turns UNSAT (then: CaDiCaL + LRAT + cake_lpr on the CNF of monotonicity and sandwich clauses, with each sandwich clause certified separately by a rank computation).  **WALL 2** (shared with direction 4).  Nothing built.

## 4. Direction 4 (adversary), with one new certified lemma

An adversary proves evasiveness by ending the game at a critical pair (Y, Y+e), Y in P, Y+e not in P.  What is certified about P is IN = D u (subgraphs of the nine graphs) and OUT = (graphs with disconnected complement) u (supergraphs of the nine complements).

**Lemma (CERTIFIED, `stage5/critical_pairs.py`, standard library only, 3 s, own embedding code, log `stage5/logs/critical_pairs.log`).**  IN and OUT are disjoint, and if Y is in IN and Y+e is in OUT then e = hv with v isolated in Y and h adjacent in Y to the other eight vertices.  (Checks: every generator O of OUT is connected and embeds in none of the nine graphs; K_{1,9} is the only generator with a bridge; none of the 3,195 triples (O, edge e of O, nine-graph H) has O - e inside H.)

Consequence (PROVED, short): an adversary strategy that does not depend on P must end at a star pair, and the algorithm avoids every star pair after two queries: ask {0,2} and {1,2}; if the answers agree, keep {0,1} for the end (its endpoints both have a present edge, or both have an absent edge to the third vertex 2); if they differ, keep {2,3} for the end (vertex 2 has a present and an absent edge, vertex 3 is not on the absent one, so 2 is neither the isolated vertex nor the hub; and if 3 were the isolated vertex, 2 would have to be the hub).  So the nine graphs add nothing for an adversary, and any adversary proof must depend on X, that is, must be carried out for each P of a family that is known only through necessary conditions.  With the game size of section 3 this is WALL 2.  (Stage 1 had the n = 4, 5, 6 version of the same failure for the pair (D, D*).)

## 5. Direction 5 (is type 5 a near-counterexample?), about 15 % of the effort

Stopping point fixed in advance: one closed-form family, three hill-climb seeds of at most 10 minutes, one exact solver question; stop unless a candidate reaches 0 violated library equations.

* Observation: the certified facts fit a cut picture.  With phi(G) = min and Phi(G) = max over bipartitions of e_G(A,B)/(|A||B|), the dual pair P_< = {phi + Phi < 1}, P_<= = {phi + Phi <= 1} is monotone, lies between D and D*, and puts the seven named Z5 graphs in and their complements out (`cand1b.py`).  No property of this kind had been tested before (stage 1 tried edge-count and degree thresholds: 1,247-2,901 violated equations).
* `cand1.py`: P_< violates 401 of the 4,164 library equations and its Z5 trace (16 classes) is not one of the 8 admissible patterns.  Not a candidate.
* `cand2.py`: hill-climb over the 200-parameter family "min_a alpha_a(emin_a(G)) < min_a beta_a(a(10-a) - emax_a(G))".  Best: 45 violated equations (11 of them equations of sigma) and a Z5 trace at Hamming distance 1 from an admissible pattern (seed 1, 9 minutes, best reached after 5.5; seeds 2 and 3 ran 90 s only and reached 56 and 92; `logs/cand2_seed*.log`).  Heuristic.
* `prof_z3.py ext`: exact question "is there ANY monotone X, containing the certified units, whose membership depends only on the extreme-cut profile (emin_a, emax_a), a = 1..5, and which satisfies the library Euler system?"  z3: **UNSAT in 1 s** (24,494 profile variables, 698 distinct equations).  SOLVER fact, not replayed.  The full cut spectrum separates 40,620 of the 41,134 classes, so it is no restriction (`prof_z3.py full`, stopped after 110 s without an answer, not informative).
* Conclusion: no closed-form near-counterexample; the stage 2-4 witnesses remain unstructured solver output defined on 41,134 classes only.  A non-evasive P: still below 1 % in my estimate.  A Z-acyclic P (which would end the topological method without contradicting Karp): about 60 % (range 30-85), CONJ; nothing here tests it.

## 6. What would move it, in order

1. **The sigma rung with one-shot clauses** (about 700-2,000 CPU-hours, one agent-week).  From the persistence reduction keep a representative cycle z of an excess class of H_s(Rel^sigma(X)) and a cocycle phi with phi(z) != 0.  Sound clause: NOT x_c for the maximal classes c of supp(z), OR x_d for the minimal classes d that own an (s+1)-cell t with phi(dt) != 0 (phi extended by 0 to new s-cells).  One reduction and one co-reduction per clause instead of several hundred.  My probability that this rung is contradictory: 5 % (range 2-10).  If it is, certification = LRAT proof of the propositional part (pg39 chain) plus one rank certificate per clause used.
2. **Smith of the eight 2-groups, branch 0** (about 20 CPU-hours, 5 GB), from Git Bash in `karp10/` (assembled from `stage4/run_branches.sh` and the `lns4g` log; not test-run in this session):
   `export PY=/c/Python314/python PYTHONUNBUFFERED=1 NUMBA_NUM_THREADS=1 EXTRAS=0 THREADS=4 NCORES=4 MAXDD=40 MAXVIOL=6 ZTIMEOUT=600000 PATTERN=0 EULER8=22,24,70,80,18,20,76,17 OMIT=22,24,70,80,18,20,76,17 D2=a LOAD=lnsb_0,lns4f,i0,i2,i3 SEEDTAG=lnsb_0 VERBOSE=1`
   then `lease-run.sh --project karp10 --min-gb 5 --minutes 480 --eta-minutes 300 --match "run_lns.py" ... -- $PY stage4/run_lns.py lns5_b0 > stage5/logs/lns5_b0.log 2>&1`.  Probability of an empty core: 2-3 %.
3. **A non-topological necessary condition for non-evasiveness** that can be evaluated on orbit data.  I looked for one and found none: degree, sensitivity and the (1+t)-divisibility identities all reduce to Euler characteristics of the minors on the tree.
4. **Turn the question round**: prove that the Euler/Oliver method cannot decide n = 10 by exhibiting a full monotone property that passes every Euler equation of every subgroup class.  Route: P = D u (down-closure of a stage-4 witness) is a genuine property that keeps every equation of the imposed classes (their invariant graphs are all universe classes); the other classes have to be repaired on graphs outside the universe.  Days of engineering (Hasse diagram on an interval of the 12 million classes); a negative result, not the theorem.

## 7. Estimates

* Full ten-vertex theorem within three months with what is now known: **about 2 % (range 0.5-6)**.  It would take either a contradiction on the sigma rung (5 %) followed by a certified refutation, or a new idea that uses the decision tree itself rather than acyclicity.
* The eight exclusions and the partial theorem are unaffected; I found no error in them.

## 8. Files

`stage5/critical_pairs.py` (certified lemma, stdlib), `cand_common.py`, `cand1.py`, `cand1b.py`, `cand2.py`, `cand2_an.py`, `prof_z3.py` (candidate tests; need numpy, scipy, z3 and the regenerable `stage2/data2`, `stage3/data3`), `stage5/logs/`.

# ENCODER-SEMANTICS: what the CNFs of the pg39 proof say, and why

Written 2026-11-11 in answer to an AI referee's request for a proof of the encoder semantics.  Scope: every clause family that occurs in the base CNFs `F(case)`, the plane CNFs `G(L)` (level-1 and level-2 completeness CNFs), the
planar frame CNFs (`planar_complete_*`), and the cube/child unit clauses.  The encoder is the function `replay.sum_in_set` (file `replay/replay.py`, lines 153-198), used through `replay/rebuild_all.py`
(`Case.build_base`, lines 138-152; `Case.plane_cnf`, lines 155-180; `planar_cnf`, lines 296-338).  `code/cnfenc.py` (the module that produced the files originally) contains the same function; the
byte-for-byte equality of the 1,424 rebuilt CNFs with the manifest hashes shows that the two produce the same clauses wherever they are used.

Every claim is tagged **[proved]** (argued in this text), **[exhaustive]** (decided by `replay/encoder_exhaustive.py` for every input of a stated finite set) or **[sampled]** (random inputs; size stated).

## 0. What actually needs to be true for the theorem

The theorem is a non-existence statement proved by UNSAT.  Let `F` be one of the CNFs.  If a genuine object (a (2,5)-set K with the normalisations) yields a model of `F`, and `F` is UNSAT, there is no such object.
So the logical chain of PROOF-CHAIN uses only **completeness of the encoding**: *every genuine indicator vector extends, through auxiliary variables, to a model of the CNF* (Lemma 1(a), Lemma 2, the unit-clause table in section 3).
**Soundness** (every model restricts to a genuine indicator) is *not* needed for any refutation; it is needed (i) to be able to say that the SAT control and the planar solution lists consist of genuine sets
(this is checked directly, geometrically, by `replay/planar_validity.py`, not through the CNF) and (ii) for the referee's request of an exact "iff".  Both directions are proved below.

## 1. The DP gadget `sum_in_set(cnf, lits, weights, allowed, amo)`

### 1.1 The code (quoted from `replay/replay.py`, lines 153-198)

```python
153 def sum_in_set(cnf, lits, weights, allowed, amo=True):
154     n = len(lits)
155     allowed = set(allowed)
156     top = max(allowed)
157     suffix = [None] * (n + 1)
158     suffix[n] = {0}
159     for i in range(n - 1, -1, -1):
160         w = weights[i]
161         suffix[i] = suffix[i + 1] | {s + w for s in suffix[i + 1] if s + w <= top}
162     layer = [None] * (n + 1)
163     layer[0] = {0} if any((a - 0) in suffix[0] for a in allowed) else set()
164     for i in range(n):
165         w = weights[i]
166         nxt = set()
167         for s in layer[i]:
168             for t in (s, s + w):
169                 if t <= top and any((a - t) in suffix[i + 1] for a in allowed):
170                     nxt.add(t)
171         layer[i + 1] = nxt
172     layer[n] = {s for s in layer[n] if s in allowed}
173     var = [dict() for _ in range(n + 1)]
174     for i in range(n + 1):
175         for s in sorted(layer[i]):
176             var[i][s] = cnf.new_var()
177     if not layer[0]:
178         cnf.add()
179         return {}
180     cnf.add(var[0][0])
181     for i in range(n):
182         z, w = lits[i], weights[i]
183         for s, v in var[i].items():
184             t0 = var[i + 1].get(s)
185             t1 = var[i + 1].get(s + w)
186             cnf.add(*([-v, z] + ([t0] if t0 else [])))
187             cnf.add(*([-v, -z] + ([t1] if t1 else [])))
188         for t, v in var[i + 1].items():
189             p0 = var[i].get(t)
190             p1 = var[i].get(t - w)
191             cnf.add(*([-v, z] + ([p0] if p0 else [])))
192             cnf.add(*([-v, -z] + ([p1] if p1 else [])))
193         if amo and i + 1 < n + 1:
194             vs = list(var[i + 1].values())
195             for a in range(len(vs)):
196                 for b in range(a + 1, len(vs)):
197                     cnf.add(-vs[a], -vs[b])
198     return dict(var[n])
```

In words.  Let the input literals be z_0..z_{n-1} (the value of z_i is 0 or 1) and all weights w_i = 1 in every use in this proof (the lemma is proved for positive integer weights).
`layer[i]` is the set of *states* (i, s) kept: s is a possible prefix sum after i items that is at most `top` and that can still be completed to an element of `allowed` (the pruning: `suffix[i+1]` holds
all subset sums of the remaining items that are `<= top`).  `var[i][s]` is a fresh variable for each kept state.  The clauses are:

* (U) the unit `var[0][0]` (line 180);
* (F0) for each state (i,s): `(-var[i][s] \/ z_i \/ var[i+1][s])`, and if `(i+1, s)` was pruned, simply `(-var[i][s] \/ z_i)` (line 186);
* (F1) `(-var[i][s] \/ -z_i \/ var[i+1][s+w_i])`, resp. `(-var[i][s] \/ -z_i)` if pruned (line 187);
* (B0) for each state (i+1,t): `(-var[i+1][t] \/ z_i \/ var[i][t])`, resp. without the last literal if (i,t) is not a kept state (line 191);
* (B1) `(-var[i+1][t] \/ -z_i \/ var[i][t-w_i])`, resp. without the last literal (line 192);
* (A) when `amo=True`, pairwise "at most one" among the states of each layer 1..n (lines 193-197);
* the returned dictionary maps each allowed sum `s` to the variable of the final state `(n, s)` (line 198).  If no allowed sum is reachable at all the empty clause is added (line 178), so the CNF is UNSAT, as it should be.

Notation: for an input vector z in {0,1}^n let `s_i(z) = sum_{j<i} w_j z_j` (prefix sums) and `sigma = s_n(z)`.

### 1.2 Lemma 1 (exactness)  [proved]; [exhaustive] tests in section 4

> **Lemma 1.**  Let all `w_i` be positive integers and `allowed` a finite set of non-negative integers.
> (a) *Completeness.* If `sigma(z) in allowed`, there is an assignment of the auxiliary variables that, together with z, satisfies every clause (U),(F),(B),(A).  It sets `var[i][s]` true iff `s = s_i(z)`.
> (b) *Soundness.* If z together with some assignment of the auxiliary variables satisfies (U),(F),(B), then `sigma(z) in allowed`.
> (c) *Functionality.* In every model of (U),(F),(B), `var[i][s]` is true iff `s = s_i(z)`; consequently the model is unique given z, `var[n][s]` is true iff `sigma = s` (so the returned `fin[5]` is exactly the indicator [sum = 5], which is how the `y` variables of lines and the local-count variables are used), and (A) is implied (redundant but harmless).

*Proof.*  Pruning invariant (P): if `sigma(z) in allowed` then `s_i(z) in layer[i]` for every i.  Indeed `s_0 = 0`; every subset sum of the items `i..n-1` that is `<= top` lies in `suffix[i]` by the recursion on
line 161 (induction from i = n, using that weights are positive so every sub-sum of a sum `<= top` is `<= top`); `s_i <= sigma <= top` and `sigma - s_i` is such a subset sum of items `i..n-1`, so `sigma = s_i + (sigma-s_i)` with
`(sigma - s_i) in suffix[i]`; applying this at `i+1` gives the test on line 169 for `t = s_{i+1}`, and `s_{i+1}` is one of `s_i, s_i+w_i`, so `s_{i+1} in layer[i+1]` by induction; finally `sigma in allowed` survives line 172.

(a) Put `var[i][s] = [s = s_i(z)]`; by (P) `var[i][s_i]` exists for every i.  (U): `s_0 = 0`.  (F0),(F1): the clause is satisfied unless `var[i][s]` is true, i.e. `s = s_i`; then if `z_i = 0` we have `s_{i+1} = s`
and `var[i+1][s]` exists and is true by (P), and if `z_i = 1` we have `s_{i+1} = s + w_i` and `var[i+1][s+w_i]` is true.  (B0),(B1): the clause can only fail if `var[i+1][t]` is true, i.e. `t = s_{i+1}`; if `z_i = 0` this
means `t = s_i` so `var[i][t]` is true; if `z_i = 1` it means `t - w_i = s_i`, so `var[i][t-w_i]` is true.  (A): exactly one variable per layer is true.

(c) and (b).  Take any model.  We show by induction on i that the set of true variables of layer i is exactly `{ var[i][s_i] }` and that this variable exists.  i = 0: layer 0 is `{0}` (line 163) and (U) makes it true.
Step: (forward) `var[i][s_i]` is true; if `z_i = 0` then (F0) forces `var[i+1][s_i]` to exist and to be true (if it does not exist the clause `(-var[i][s_i] \/ z_i)` is violated), and `s_i = s_{i+1}`;
if `z_i = 1`, (F1) forces `var[i+1][s_i + w_i] = var[i+1][s_{i+1}]` likewise.  (backward) if some `var[i+1][t]` is true, then if `z_i = 0` (B0) forces `var[i][t]` true (it must exist), so `t = s_i = s_{i+1}`;
if `z_i = 1`, (B1) forces `var[i][t-w_i]` true, so `t - w_i = s_i`, `t = s_{i+1}`.  Hence layer i+1 has exactly the true variable `var[i+1][s_{i+1}]`.
At i = n the existing variables are the states in `layer[n]`, a subset of `allowed` (line 172), and `var[n][sigma]` exists and is true, so `sigma in allowed`: this is (b).  The unique-extension statement and the meaning of `fin[s]` follow.
(A) holds in the model because exactly one variable per layer is true.  If `layer[0]` is empty no allowed sum is reachable, no `z` is valid, and the empty clause is the correct output.  QED

Remarks.  (1) The pairwise (A) clauses are redundant (by (c)); `amo=True` is used for the 10-point lines and `amo=False` for the sums and local counts; this does not change the semantics.  (2) The `any(...)` test on
lines 163 and 169 uses `suffix[i+1]` capped at `top = max(allowed)`; it is exactly the co-reachability pruning of (P).  (3) Positive weights are essential for (P); all weights used are 1.

## 2. The remaining clause families: what each says and which proof step implies it

Notation: x_p = [p in K] is a variable `x[p]` (variables 1..820 in the base CNF, 1..91 in plane/planar CNFs: the point variables are allocated first); `y_l = fin[5]` of line l's gadget, i.e. `y_l = [|l n K| = 5]` by Lemma 1(c).

| # | clause family | code | meaning | implied by |
|---|---|---|---|---|
| 1 | line gadgets, allowed {2,5}, 10 literals `x[p]`, p in the line | `rebuild_all.py` 143-145 (base), 164-166 (plane), 301-303 (planar) | every line meets K (resp. S) in 2 or 5 points, and `y_l = [meets in 5]` | the hypothesis (two-intersection set of type (2,5)); for a section S = K n pi: the lines of pi are lines of PG(3,9) (PROOF-CHAIN 2.1) |
| 2 | size gadget, allowed {k} resp. {m}, 820 resp. 91 literals | 146, 168, 304 | sum of all x = k (resp. m) | k in {260,287}: PROOF-CHAIN 1.1; m in {26,35} and m = 26 iff `|Sigma n ell| = 5`: 1.2 (planes), 1.3 (planar sets) |
| 3 | local counts | 169-174 (plane), 305-310 (planar) | for each point p of the plane: `sum of y over the 10 lines through p` in {d, d+3}, and `x_p <-> [sum = d+3]` (`d = (m-20)/3`) | PROOF-CHAIN 2.1: p in S lies on `a` 2-secants and `b` 5-secants with `a+b = 10`, `a+4b = m-1`, so `b = (m-11)/3 = d+3`; p outside: `a+b = 10`, `2a+5b = m`, `b = (m-20)/3 = d`.  Both counts are in {d, d+3}, and the equivalence `x_p <-> (b = d+3)` is the dichotomy just proved (checked on the 1,345 real planar sets by `replay/lemma_checks.py` and again by `replay/planar_validity.py`) |
| 4 | point units P0 and the three height-0 points `(v_i:0)` (base and plane CNFs; in a plane CNF only those lying in it) | 147-148 (base), 178-180 (plane; loop `for p, v in {**self.pt_units, **fixed}.items()`) | these points are in K | PROOF-CHAIN 1.4 (P0 in K; elations normalise one point on each line of direction v_i), verified on the actual geometry by `replay/lemma_checks.py` |
| 5 | line units on the lines through P0 | 149-150 (base), 175-177 (plane) | `y_l = 0` (two-secant) iff the direction of l lies in Sigma, `y_l = 1` otherwise | definition of Sigma (1.2) and the choice of the representative (1.4) |
| 6 | cube units (92 clauses: the 4 base units + the 88 free points of pi1) | `rebuild_all.cube_units` | the section of K on pi1 equals the representative | the cube is one case of the level-1 cover; by 2.1/2.2 some homothety image of K has this section |
| 7 | child units (80 further units) | `children` in `run_star` | the section of K on pi2 equals the child | level-2 cover (PROOF-CHAIN 3) |
| 8 | planar frame units: `x_p` for p in R, `-x_p` for the other points of l0; `x[origin]` | 312-314 | K n l0 = R, origin in K | PROOF-CHAIN 1.3 step 1: a 5-secant moved to l0 by PGammaL(3,9), its trace in one of the two orbit representatives; translations of the affine plane are regular on the 81 affine points |
| 9 | lex-leader clauses (flag variables `e`) | 316-335 | `x <=_lex x o h` for every h in H_R | Lemma 3 below; H_R is constructed independently and shown to be the full stabiliser in `replay/planar_validity.py` |
| 10 | blocking clauses | 336-338 (planar), 428-431 (level 1), 450-453 (level 2) | the assignment to the free (resp. non-l0) points differs from a listed one in at least one point | Lemma 4 below |

### 2.1 Lemma 3 (the lex-leader clause set)  [proved] + [exhaustive on small instances]

The loop is lines 316-335 of `rebuild_all.py`:

```python
    for h in H:
        if h != ident:
            prev = None
            for p in order:
                a, b = x[p], x[h[p]]
                if a == b:
                    continue
                if prev is None:
                    cnf.add(-a, b)
                else:
                    cnf.add(-prev, -a, b)
                e = cnf.new_var()
                if prev is None:
                    cnf.add(a, b, e)
                    cnf.add(-a, -b, e)
                else:
                    cnf.add(-prev, a, b, e)
                    cnf.add(-prev, -a, -b, e)
                prev = e
```

`order` = the points off l0 in increasing index.  Let `x_h(p) = x[h[p]]`.  Positions with `h[p] = p` are skipped (they never differ).  For the remaining positions p_1 < p_2 < ... let `e_j` be the flag created at p_j
(and `prev_j = e_{j-1}`, `prev_1 = true`).  Clauses at position j: `(prev_j -> (x(p_j) <= x_h(p_j)))`, and `prev_j /\ (x(p_j) = x_h(p_j)) -> e_j` (the two clauses with `a,b,e`).
*Soundness for lex-least:* if `x <=_lex x_h` (0 < 1), set `e_j := [x and x_h agree on p_1..p_j]`; then `prev_j` true means agreement before j, so `x(p_j) <= x_h(p_j)` (else x would be lex-larger at the first difference), and the implications defining `e_j` hold.
So the lex-least element of every H_R-orbit, which satisfies `x <=_lex x_h` for every h, extends to a model.  *Soundness of the clauses (converse):* in a model, by induction `e_j` is true whenever x and x_h agree on p_1..p_j, so at the first
differing position j the clause `(prev_j -> x(p_j) <= x_h(p_j))` forces `x(p_j) = 0 < 1 = x_h(p_j)`; thus `x <=_lex x_h`.  The flags `e_j` after the first difference are left unconstrained in one direction (not functions of x), which is harmless.
[exhaustive] `E6`: 120 random instances (n in {8,10,12}, 1-4 random permutations, partial compare sets), all 2^n inputs: accepted inputs = `{x : x <=_lex x o h for all h}`, where a complete DPLL search decides the extension.
The text of the loop in `rebuild_all.py` is compared with the copy tested.

Because H_R is a group, `x <=_lex x o h for all h in H_R` says exactly that x is the lex-least element of its orbit `{x o h}`; every orbit has one such element, so no orbit is lost (completeness of the symmetry break).

### 2.2 Lemma 4 (blocking clauses)  [proved]

A blocking clause over a set `Fr` of variables with signs of a listed assignment `a` is `\/_{p in Fr} (x_p != a_p)`.  Suppose `F` is a CNF containing units that fix all the points outside `Fr`, and `F /\ blocking(a) for all a in List` is UNSAT.
Then every model of `F` restricts to an assignment of `Fr` that equals some `a in List`; i.e. the unit-fixed points together with `a` give the full assignment.  Cover logic: for level 1, `F = G(L1)` over the 91 points of pi1 with `Fr` the 88 points not fixed by units
(`rebuild_all.py` 428-431), so every section (a model of `G(L1)` by section 2 rows 1-5) is one of the blocked ones.  Level 2: `Fr` = the 80 points of pi2 not on the common line and not fixed (450-453).  Planar: `Fr = order` (81 points) and the 10 points of l0 are unit-fixed (336-338).
`replay/cover_asserts.py` computes `Fr` from the geometry and shows that the blocked set is exactly the homothety expansion of the cube list (level 1) resp. the children (level 2), cube by cube.

## 3. Assembly: the complete argument that the CNFs contain every genuine object

For a hypothetical K with the normalisations of PROOF-CHAIN 1.4 and the cube section chosen as in 2.1-2.2: Lemma 1(a) applied to each gadget in the base CNF (rows 1, 2 of the table) gives auxiliary values satisfying all gadget clauses
(the gadgets share only the input variables x_p and read the `y` outputs of other gadgets as inputs, e.g. in rows 3 and 5; Lemma 1(a) applies to those too because `y_l` is itself the true indicator [sum = 5] by Lemma 1(c),
and the local-count clause `x_p <-> fin[d+3]` is true by row 3); rows 4, 5, 6 and 7 are true by the cited steps; hence K is a model of the cube CNF and of the child CNF.
For the plane CNFs, the section S = K n pi satisfies rows 1, 2, 3, 4, 5 (and for level 2, row 4's `fixed` units, which are the cube section on the common line).  By Lemma 4, the cover CNFs show S is one of the listed sections.
The planar CNF is analogous with rows 1, 2, 3, 8, 9 (by 1.3 step 1 and Lemma 3: the lex-least member of the H_R orbit of K is a model).

Things that are **not** covered by this document: that the GF(9)/PG(3,9) objects fed to the encoder are the intended geometry (checked by independent re-implementation: `replay/rebuild_all.py` against `code/`, byte-for-byte; and by `replay/lemma_checks.py`),
and cake_lpr's correctness.

## 4. The tests: `replay/encoder_exhaustive.py` (full log `logs/encoder_exhaustive.log`)

Method.  For an input vector z, the CNF (with its units) has an extension iff unit propagation (UP) from the inputs ends without conflict in a total assignment satisfying every clause.  A UP conflict refutes every extension (UP is sound),
so each input is decided with nothing assumed; if UP neither conflicts nor completes the test fails (for the lex gadget, whose flags are not functions of x, a complete DPLL search is used instead).
The engine is independent of the encoder (separate code, no shared state), and a second, independent DPLL model counter (no input fixed) is used on the small shapes.

Headline (run of 2026-11-11, 153 s, one thread, `logs/encoder_exhaustive.log`): `ENCODER TESTS: ALL PASS`.  Results, one line per test:

```
E1/E2: production shapes, all 2^10 inputs, plus independent model counting
  line shape {2,5} amo=True             :  297 valid inputs of 1024 accepted, all others rejected; models of the CNF counted = 297 (= valid inputs); 55 variables, 262 clauses
  line shape {2,5} amo=False            :  297 valid inputs of 1024 accepted, all others rejected; models of the CNF counted = 297 (= valid inputs); 55 variables, 175 clauses
  local count m=26 {d, d+3} = {2,5}     :  297 valid inputs of 1024 accepted, all others rejected; models of the CNF counted = 297 (= valid inputs); 55 variables, 175 clauses
  local count m=35 {d, d+3} = {5,8}     :  297 valid inputs of 1024 accepted, all others rejected; models of the CNF counted = 297 (= valid inputs); 55 variables, 175 clauses
E3: random shapes and unit-weight shapes up to n = 16, ALL inputs
  300 random weighted shapes (n <= 12, weights 1..3, 1-3 allowed sums, amo on/off): all 261034 inputs correct; model count = valid inputs for the n <= 7 ones (13s)
  unit weights n=14 allowed [5], amo=False: all 16384 inputs correct (2002 valid; 74 vars, 237 clauses; 1s)
  unit weights n=14 allowed [7], amo=False: all 16384 inputs correct (3432 valid; 78 vars, 253 clauses; 1s)
  unit weights n=14 allowed [2, 5], amo=False: all 16384 inputs correct (2093 valid; 83 vars, 271 clauses; 1s)
  unit weights n=14 allowed [5, 8], amo=False: all 16384 inputs correct (5005 valid; 95 vars, 319 clauses; 1s)
  unit weights n=16 allowed [2, 5], amo=False: all 65536 inputs correct (4488 valid; 97 vars, 319 clauses; 4s)
E4: real planar CNF (rebuild_all.planar_cnf), q = 9
  m=26 (10063 vars, 47075 clauses): 74 real solutions all accepted; 296 perturbed sets (swap/delete/add) all decided as the direct geometric check decides (0 of them were genuine (2,5)-sets and accepted) (3s)
  m=35 (10333 vars, 48155 clauses): 1271 real solutions all accepted; 5084 perturbed sets (swap/delete/add) all decided as the direct geometric check decides (0 of them were genuine (2,5)-sets and accepted) (50s)
E5: the plane pipeline (line shapes, size, local counts) on PG(2,3) and PG(2,4): brute force vs CNF
  replica clause stream == rebuild_all.planar_cnf (PG(2,9), m=26, 47075 clauses): identical
  PG(2,3) type (1,4), m=4 (the lines)                 : 8192 subsets tested (all subsets of size 0..13, plus random neighbours for q = 4): CNF == brute force for every one; 13 genuine sets (1s)
  PG(2,4) type (0,2), m=6 (hyperovals)                : 60264 subsets tested (all subsets of size [6], plus random neighbours for q = 4): CNF == brute force for every one; 168 genuine sets (30s)
  PG(2,4) type (3,5), m=15 (complements of hyperovals): 60264 subsets tested (all subsets of size [15], plus random neighbours for q = 4): CNF == brute force for every one; 168 genuine sets (25s)
E6: lex-leader clause gadget
  120 random instances (n in {8,10,12}, 1-4 permutations, partial compare order), ALL 2^n inputs: accepted inputs == {x : x <=_lex x o h for all h}
E7: plain sum constraints at the real sizes (SAMPLED)
  [1]*91 = 26 (1873 vars, 7125 clauses): 100 random inputs of weight 26 extend (unique total model satisfying every clause), 402 inputs of weight 26+-1,+-2, 0, 91 rejected (1s)
  [1]*91 = 35 (2143 vars, 8205 clauses): 100 random inputs of weight 35 extend (unique total model satisfying every clause), 402 inputs of weight 35+-1,+-2, 0, 91 rejected (1s)
  [1]*820 = 260 (147241 vars, 585681 clauses): 3 random inputs of weight 260 extend (unique total model satisfying every clause), 14 inputs of weight 260+-1,+-2, 0, 820 rejected (3s)
  [1]*820 = 287 (154612 vars, 615165 clauses): 3 random inputs of weight 287 extend (unique total model satisfying every clause), 14 inputs of weight 287+-1,+-2, 0, 820 rejected (3s)
ENCODER TESTS: ALL PASS (153s)
```

Reading guide.  E1/E2: for each of the four production shapes all 1,024 inputs were decided and the independent model counter found exactly 297 = 45 + 252 models (= C(10,2) + C(10,5)) in the CNFs of the line shape (with and without the pairwise at-most-one clauses),
so every valid input has exactly one extension and the CNF has no other model: this is Lemma 1 (a),(b),(c) for the real line constraint, exhaustively.  E5 (replica) is checked clause for clause identical to `rebuild_all.planar_cnf` on PG(2,9), m = 26,
and then run on planes where all subsets can be enumerated; the counts of genuine sets (13 lines of PG(2,3); 168 hyperovals and 168 complements in PG(2,4)) are the known ones.  E4: no random perturbation of a real solution happened to be a genuine set
(the 'genuine and accepted' counters are 0), so the E4 negatives only show that near-misses are rejected; the sharp single-constraint negatives are the exhaustive E1-E3 runs.

What is exhaustive and what is not, plainly.

* **Exhaustive** over every input: the production line shape (10 points, {2,5}), the production local-count shapes (10 inputs, {2,5} and {5,8}), 300 random weighted shapes up to n = 12, unit-weight shapes n = 14, 16, the plane pipeline on PG(2,3) (all 2^13 subsets) and,
  for PG(2,4), every subset of size m (hyperovals, and complements), and the lex gadget on 120 random small instances.  With model counting for the small ones: the number of models equals the number of valid inputs (no spurious auxiliary freedom).
* **Sampled**: the real plane/planar CNF at q = 9 (all 1,345 real solutions as positives, plus about 4 random perturbations each as negatives, decided against the direct geometric check) and the plain sum constraints at n = 91, 820 (random inputs of weight k, k+-1, k+-2, 0, n).
  The number of inputs at these sizes is astronomical (C(91,26) is about 1e23), so an exhaustive or model-counting check is not possible; the argument at these sizes is the proof of section 1 (which is independent of n), supported by the exact identity of the replica pipeline with `planar_cnf`.
* **Not tested by this script**: the base CNF as a whole at q = 9 (no positive example exists: that is the theorem; the Hermitian SAT control of PROOF-CHAIN section 4 exercises the same encoder with allowed set {1,4,10} and passes), the cube/child unit clauses (they are plain unit clauses; their assignment is checked by `replay/cover_asserts.py`).

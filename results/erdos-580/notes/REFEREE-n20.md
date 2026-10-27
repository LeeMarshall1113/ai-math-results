# Hostile referee pass: Erdos #580 vertex form, n = 20 (m = 10), hence n = 21
Date 2026-10-26. Verdict: SURVIVES (minor documentation gaps, none affecting soundness).

## 1. Completeness
- Independent tree enumeration (own AHU canon, not lks.py): 106 trees on 10 vertices. OK.
- Independent re-classification (brute force over all paths, two-star test, no networkx): D 58, P 11, C 6, E 31, identical to classification_m10.json tree by tree (0 mismatches).
- results/m10.jsonl: 152 rows; the set of trees with a verified:true row is exactly the 31 E-indices. Each final row has audit OK, cadical_rc 20, gate_rc 0, cake_rc 0 "s VERIFIED UNSAT" (tree 45 see 5).
- 58+11+6+31 = 106. No tree missing.

## 2. Cited theorems (read Piguet-Stein arXiv 0712.3382 directly)
- Thm 2: any k, n, >= n/2 vertices of degree >= k, diameter <= 5, <= k edges. Applies (k = 9, host has 10 vertices of degree >= 10).
- Thm 3: l >= c, and (c even or l + c >= n/2-ish). All 6 C-trees have c in {2,4}, i.e. c even, so the garbled/misquoted odd-c bound is never used at m = 10. l >= c holds in all 6 (l = 4,c = 4; l = 3,c = 2; l = 2,c = 2). T(k,l,c) = path of length k-l plus two stars at distance c; k = 9 matches.
- Path-plus-star: Bazgan-Li-Wozniak as described in the PS intro (paths, and a path with a star centred on a path vertex). 11 P-trees fit.
- Hypotheses "degree >= k" with k = 9 are weaker than the vertex form's >= 10, so the theorems apply.

## 3. Reduction
- Edge-minimal host rechecked from scratch: |L| = m (else delete an edge that drops a degree-m vertex); an L-vertex of degree > m would have all neighbours in L (no S neighbour, else the edge is deletable), giving degree <= m-1, contradiction. So all L have degree exactly m, S independent. OK.
- Trees with < m vertices embed in an m-vertex tree: OK.
- Leaf-extension: a marked image in L has >= r+1 unused neighbours for r deleted leaves; greedy works. OK.
- n = 21 from n = 20: delete any vertex; at least 10 vertices keep degree >= 10, trees have <= 10 vertices. OK.

## 4. Encoding
- Exhaustive check at m = 4 (n = 8): all 16,618 hosts satisfying the base constraints (L exact degree 4, S independent, deg <= 3) and both trees: whenever a rooted-core embedding (marked parents in L) exists, the tree embeds. Zero violations, so the avoidance clauses are implied by tree-avoidance.
- Base, audit and clause validator logic in lks.py read; consistent with the stated semantics.

## 5. Tree 45 (S(2,3,4))
- work/certs/m10_t045.cnf sha256 e323aa7e... equals the cnf_sha256 in both logged generations; header "p cnf 1685 14823490" = 5890 + 14,817,600; actual clause-line count 14,823,490. Header comment: core vertices [0,1,2,3,4,5,7], marked [3,5,7] = spider legs 2,3,4. OK.
- m10_t045.lrat is 649,828,266 bytes, sha 6a199fd0... equal to the lrat_sha256 recorded by the CaDiCaL run; the d73 script runs cake_lpr on exactly those two files (heap 8500). Log: s VERIFIED UNSAT, CAKE_RC=0.

## Gaps (non-fatal)
- The tree-45 final row records no hashes; I matched them by hand above. Suggest recording them.
- Tree 47 and others use keep/keepS variants; soundness relies on the audit (AUDIT OK recorded), which I did not re-run.
- Trusted base unchanged: cake_lpr, lks.py, cited theorems, the two lemmas (re-derived above).

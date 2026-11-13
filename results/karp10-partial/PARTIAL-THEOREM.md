# Karp / Aanderaa-Rosenberg for n = 10: the partial theorem

Status: **no exclusion of type 5**.  This note states exactly what is proved (modulo the cited theorems) about a non-evasive monotone graph property on 10 vertices, what the replayable certificates are, and what is only computational evidence. 

## 0. Setting and conventions

* V = {0..9}, E = the 45 edges of K_10.  P is a family of graphs on V that is closed under isomorphism and under taking subgraphs (a monotone decreasing property; the increasing version is the complement), nonempty, and not all graphs.  Delta_P is the simplicial complex on E whose faces are the edge sets of graphs in P.  P is **non-evasive** if some decision tree for "G in P" has depth <= 44.
* D = disconnected graphs; D* = graphs with connected complement (= D u CC); CC = connected graphs with connected complement; N = connected graphs with disconnected complement.  P* = {G : complement(G) not in P} is again non-evasive.
* Cited theorems (not re-proved): Kahn-Saks-Sturtevant (non-evasive => collapsible => Z-acyclic => chi(Delta_P) = 1); Smith theory (p-group fixed sets of an F_p-acyclic complex are F_p-acyclic); Oliver / Angel-Borja Thm 4.1 / Engstrom Thm 4, 5 (fixed-point Euler characteristics for groups in Oliver's classes G^1, G^q); S_10 is transitive on the 45 edges, so lk(e) and dl(e) of any edge are non-evasive.  An AI referee read the argument and the source and judged the method sound; it did not run the replay (see README).

## 1. The theorem (what is certified)

Let P be non-evasive as above.  Then:

1. **Type.**  P is of Angel-Borja type 5: P meets the ten special classes {empty, K_10, G1..G5 and complements} in {empty, G1, G2, G3, G4, G5}.  Types 1, 2, 3, 4, 6, 7, 8, 9 are impossible.  Certificates `certs/type{1,2,3,4,6,7,8,9}*.json` (14 files; types 1, 9: groups of order 8/16; 3, 7: order 50; 4, 6: order 16 and order 192; 2, 8: the pair of groups of order 112 and 192).
2. **Sandwich.**  D <= P <= D*, i.e. P = D u X with X a down-set of CC graphs; in particular K_a u K_{10-a} in P for a = 1..4 and no graph with disconnected complement (e.g. the star K_{1,9}) lies in P.  Certificates `certs/type5_force_K{1,2,3,4}_K{9,8,7,6}.json`.
3. **Nine graphs.**  P contains these nine graphs, and none of their complements (so, by monotonicity, every subgraph of one of them is in P and every supergraph of a complement is not): C_10 (10 edges); the pentagon with five pendant edges (10 edges, degrees 3^5 1^5); the Petersen graph and the pentagonal prism C_5 x K_2 (15 edges, 3-regular); K_5 with five pendant edges (15 edges, 5^5 1^5); a 15-edge graph with degrees 4^5 2^5; a 20-edge 4-regular graph; (K_5 u C_5) + perfect matching (20 edges, 5^5 3^5); K_5 x K_2 (25 edges, 5-regular).  Certificates `certs/lemma_z5.json` (first seven) and `certs/lemma_z5_euler.json` (all nine, with the Euler equations of two further groups of order 80 and 160). 
4. **Z5 trace.**  The trace of P on the 32 CC isomorphism classes that occur as unions of orbits of Q = <(01234)(56789)> is one of **8** explicit patterns (sizes 9, 11, 11, 13, 19, 21, 21, 23 of the 32 classes in X).  Of the 2,539 down-sets of the class poset, 16 pass the relative-homology test of `replay/replay_smith.py` (a necessary condition for F_5-acyclicity of the fixed complex, not an equivalent one); `replay/check_z5_absolute.py` computes the F_5-homology of the whole fixed complex for every down-set and selects exactly the same 16, and then the same 8 after the Euler equations of the two groups of order 80 and 160.  These counts are pinned in the certificates (key `expect`) and replayed by `replay_all.py`.  The 8 patterns form four pairs under P <-> P*.  The nine graphs are identified with named graphs by explicit isomorphism tests in `replay/check_nine_graphs.py`.
5. **P is never self-dual.**  P != P*.  (`replay/replay_nonselfdual.py certs/lemma_z5_euler.json`, standard library only: no graph on 10 vertices is self-complementary (45 is odd), so the complement map is a fixed-point-free involution on the 32 classes; a self-dual trace would contain exactly one class of each of the 16 pairs, i.e. 16 classes, and the 8 patterns have 9, 11, 11, 13, 19, 21, 21, 23.  The replay also checks that the 8 patterns are closed under duality and none is invariant, which is consistent with but not needed for the argument.)

Replay everything: `python replay/replay_all.py` (about a minute and a half, standard library only; it first checks the SHA-256 manifest of the certificates and scripts) prints ALL VALID for the 20 certificates, the non-self-duality check in 5, the absolute-acyclicity check of the Z5 patterns and the identification of the nine graphs.

## 2. What is NOT claimed

* Nothing says that a type-5 property exists or does not exist.  Section 3 lists necessary conditions that a type-5 property must satisfy and that are *jointly satisfiable* so far; these are solutions of a system of necessary conditions on the indicator vector of X on a few thousand library classes, not graph properties.
* The F_p-Smith conditions are only necessary (they use acyclicity of fixed complexes, not collapsibility).  By Oliver's theory any consistent family of fixed-point data comes from an abstract acyclic S_10-complex; a contradiction must use that the faces are graphs (stabilisers are automorphism groups, the order is the subgraph order) or non-topological input.

## 3. Not included

The authors also ran solver searches for a type-5 property (necessary conditions on a library of about 1,500 subgroup classes). They were inconclusive: no exclusion and no counterexample. Those logs and models are not part of this release; the README gives a summary.

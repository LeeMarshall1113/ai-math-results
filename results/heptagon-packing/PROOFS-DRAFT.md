# Heptagon packing: working proofs of the structural lemmas (DRAFT, 2026-11-04)

Internal working proofs, written for the next hostile check. This is not paper prose.
Sources: NOTES.md (lemma skeleton L1-L6, the stated reduction, facts 1-9), AUDIT-2026-11-04.md
(sections 1, 2h and 5, and "What remains before Lean"), `src/bnb2.py` (**current certified code: hash 9f4fb628**, run a1 of 2026-11-04; earlier runs used 7ed516ab),
`src/local_tube_arb.py`, `src/chart_excl.py`, `runs/D5_*`, `runs/charts_w030*.json`.
Nothing in `src/`, `runs/` or NOTES.md was changed.

Markers used below:
- **[FLAG]**: a step I am unsure of, or a definition in NOTES that is ambiguous.
- **[GAP]**: something the proof needs that is not yet certified or proved.
- **[CITE]**: a standard result that is cited rather than reproved, with how far I checked it.

**Main structural change against NOTES.** L2 and L5 are proved here with a separate **finite**
Delaunay triangulation in each window. No triangulation of the infinite point set is used. This
removes the need for the existence and local finiteness theorem for infinite Delaunay
subdivisions, and for any argument that degenerate cells are triangulated "consistently" across
the plane. Each window gets its own finite triangulation, and the density bound in each window
uses only that one. The counting is done by angle sums, not by Euler's formula, which gives the
inequality we need directly (the Euler form is noted in a remark).

---

## 0. Definitions and constants

- u = cos(pi/7) and A = (7/2) sin(2pi/7) = 2.7364101886..., the area of a regular heptagon of
  circumradius 1.
- delta_DL = (2/97)(-111 + 492u - 356u^2) = 0.892690686126508...
- acrit := A/(2 delta_DL) = 1.5326754446782212... (mpmath, 40 digits).
- tau := the double 1.5326754446782214, the B&B target. tau - acrit = +1.53e-16 > 0 (mpmath, on
  the binary value).
- RTH := 2u sqrt2 = 2.5483247845270692...
- rb0 := the double 2.5483247845270696. rb0 - RTH = +2.75e-16 (mpmath).
- Heptagon H(c, t) := conv{ c + (cos(t + 2pi k/7), sin(t + 2pi k/7)) : k = 0..6 }. Then
  H(c, t) = H(c, t + 2pi/7), so orientations are taken mod 2pi/7. H(c, t) lies inside the
  closed disc B(c, 1) and contains the closed disc B(c, u).
- Packing: a family of such heptagons with pairwise disjoint interiors.

### 0.1 The pair transfer f

The edge configuration of an ordered pair (H(c_i, t_i), H(c_j, t_j)) is (r, a, b), where
r = |c_j - c_i|, phi = arg(c_j - c_i), a = t_i - phi and b = t_j - phi (mod 2pi/7).

f is the function `FTrans.f` of `bnb2.py` with the coefficients `runs/f_w030.npy` (DR = 3,
DF = 2; 100 coefficients, |c|_1 = 0.060875). Coefficient 9 is replaced by c_9 + mu*, where mu* is
the exact exactification value. mu* lies in the arb ball [-1.168485374e-12 +/- 1.9e-22]
(audit 2h), which lies inside the B&B's interval HEPT_MU.

    g(r, a, b) = sum_{k, ja, jb} c_{k, ja, jb} T_k(x(r)) F_ja(7a) F_jb(7b),
        x(r) = (2r - 2u - 4)/(4 - 2u),   F = (1, cos, sin, cos 2., sin 2.)
    f = g - g o R + g o m - g o S,
    R(r, a, b) = (r, b + pi/7, a + pi/7),   m(r, a, b) = (r, pi/7 - b, pi/7 - a),
    S(r, a, b) = (r, -a, -b).

**Lemma 0.1 (symmetries).** Mod 2pi/7 in a and b:
- R^2 = m^2 = S^2 = id, m o R = S and S o R = m;
- f o R = -f, f o m = f, f o S = -f;
- f = 0 on R-fixed configurations (b = a + pi/7, the point-reflected pairs) and on S-fixed ones
  (a, b in {0, pi/7}).
- |f| <= 4(|c|_1 + |mu*|) < 0.2436 for r in [2u, 4].

*Proof.* The compositions are direct substitutions using -pi/7 = pi/7 (mod 2pi/7). For example,
m(R(a, b)) = m(b + pi/7, a + pi/7) = (-a, -b).

Then
- f o R = g o R - g o R^2 + g o mR - g o SR = g o R - g + g o S - g o m = -f;
- f o m = g o m - g o Rm + g - g o Sm, with Rm = S and Sm = R, so f o m = f;
- f o S = g o S - g o RS + g o mS - g, with RS = m and mS = R, so f o S = -f.

At a fixed point p of R, f(p) = f(Rp) = -f(p), so f(p) = 0; the same holds for S.

For the bound: |T_k(x)| <= 1 for |x| <= 1, x([2u, 4]) = [-1, 1], and each F has modulus <= 1.
So |g| <= |c|_1 + |mu*|, and f is a sum of four copies of g. QED.

(Float check with a random 2pi/7-periodic g, 1000 points: largest violation 3.5e-14.
Scratch: `scratchpad/hp_proofs/chk_sym.py`.)

**Lemma 0.2 (reversal and mirror).**
- The reversed edge j -> i has configuration R(r, a, b).
- Reflect the configuration in a line and relabel it counterclockwise. Each ccw edge of the
  image carries m(the configuration of the corresponding original ccw edge).

*Proof.*
- Reversal: phi_ji = phi_ij + pi, so the reversed edge has a' = t_j - phi - pi and
  b' = t_i - phi - pi. Since -pi = pi/7 (mod 2pi/7), (a', b') = (b + pi/7, a + pi/7).
- Mirror: reflect in the y-axis (any reflection is this one up to a rigid motion). A vertex angle
  t goes to pi - t, and an edge direction phi goes to pi - phi. So the image of edge i -> j has
  (a, b) -> (-a, -b), and that image edge is clockwise. The ccw edge is its reverse, with
  configuration R(r, -a, -b) = (r, pi/7 - b, pi/7 - a) = m(r, a, b). QED.

### 0.2 Triangle functionals

A *configuration* T is a ccw triple (H_0, H_1, H_2) with centres c_q and orientations t_q. Its
edges are 01, 12 and 20, and the angle at c_k is opposite edge e_k.

- **b(T)** := area(c_0 c_1 c_2) + f(e_01) + f(e_12) + f(e_20), each a ccw edge configuration.
- **D5(r, c)**, defined on [rb0, 4] x [-1, 0], is the real bilinear interpolant of the float
  nodes `runs/D5_nodes.npy` (30 x 21, sha 87d669a7...) on the float grid `runs/D5_bins.json`
  (sha afa2c4ac...). Cells are closed. On a shared cell edge both cells give the same edge-linear
  interpolant, so D5 is single-valued and continuous. Set D5(r, c) := 0 for r < rb0.
  - Checked: node rows 0 and 1 are exactly 0, the column c = 0 is exactly 0, all nodes lie in
    [0, 0.23347], rb[-1] = 4 and cb runs over [-1, 0]. So 0 <= D5 <= Dmax := 0.2335, and D5 is
    continuous on [0, 4] x [-1, 0].
- **Outgoing transfer.** Dout(T) := D5(l_T, cos gamma_T) if T is obtuse, where gamma_T > pi/2 is
  its obtuse angle and l_T is the (unique, Lemma 3.1) longest edge. Otherwise Dout(T) := 0.
- **Charge function.** For r > 0 and c in [-1, 1]:

      Dhat(r, c) := max{ D5(r, c') : c' in [-c, 0] }   if c > 0,     Dhat(r, c) := 0 if c <= 0.

- **Charge bound.** Chat(T) := sum of Dhat(r_e, cos alpha_e) over the edges e of T other than the
  long edge of an obtuse T. Here alpha_e is T's own angle opposite e.

**[CLOSED 2026-11-05] NOTES fact 9 "bins are half-open"** (NOTES now corrected). That was true of the step tables D2/D3. For D4
and D5 (continuous bilinear) the cells are closed and the value does not depend on the cell
chosen. The proofs below use only the continuous D5.

---

## 1. The interface to L4: exactly what the sweeps certify

Call a configuration T **admissible** when
- (i) its three heptagons have pairwise disjoint interiors, and
- (ii) its circumradius satisfies R(T) <= 2.

(i) implies every edge is >= 2u (disjoint inscribed discs), and (ii) implies every edge is <= 4.

**Statement (M), the main sweep plus the chart lemma.** For every admissible, nonobtuse T:

    b(T) >= tau (> acrit).

**Statement (O), the obtuse-mode sweep.** For every admissible T whose longest edge is >= 2.548:

    b(T) + Dout(T) - Chat(T) >= tau.

Each human lemma uses (M) only for triangles with longest edge < RTH, and (O) only for triangles
with longest edge >= RTH.

**What the code literally checks, and why it implies (M) and (O).**

*Canonical coordinates.* c_0 = 0, c_1 = (d, 0), c_2 = (x, y), plus (t_0, t_1, t_2). The sweep
regions are:
- main: d in [1.8019377358048383, 4], x in [0, 2], y in [1.2, 3.5], t_0 in [0, 0.897598] (union
  of 45 slices), t_1, t_2 in [0, 0.8975979010256552];
- obtuse: d in [2.548, 4], x in [0, 2], y in [0.8, 3.5], t_0 in [0, 0.8975979010256552],
  t_1, t_2 the same.

A box is discarded only if one of the following holds:
- it is certainly infeasible: an edge < 2u, c_0c_1 not longest, x > d/2, R > 2, overlap on the
  whole box by SAT, or certainly obtuse in main mode;
- its interval lower bound is >= the box target;
- (main slice 39 only) it lies inside a certified chart region (Section 6).

The box targets are:
- main: tau;
- obtuse: `up(tau + up(extra))`, where `extra = obtuse_extra(B)` (`bnb2.py:330-345`).

For every point T of a box B (audit 1b/1d/2d, rederived):

    extra(B) >= Chat(T) - Dout(T).

The reasons:
- `dhat_ub(r, ca)` is >= the max of D5 over r' in the box and c' in [-ca.hi, 0], which is
  >= Dhat(r_e, cos alpha_e) at each point;
- `d_lb` is <= D5(r_01, cos gamma) at each obtuse point;
- the straddling case takes the max of the two branches.

The main-sweep pruning of obtuse boxes is harmless: (M) is only claimed for nonobtuse T.

**Lemma 1.1 (canonical reduction).** b, Dout and Chat are invariant under:
- rigid motions;
- t_q -> t_q + 2pi/7;
- cyclic relabelling;
- reflection followed by ccw relabelling.

So every admissible T is equivalent to one with:
- c_0c_1 a longest edge, c_0 = 0, c_1 = (d, 0) and y > 0;
- 0 <= x <= d/2 and t_q in [0, 2pi/7).

Such a T then satisfies:
- 2u <= d <= 4;
- y <= d sqrt3/2 <= 3.4641 < 3.5;
- y >= u^2 = 0.8117 > 0.8;
- if T is nonobtuse, y >= sqrt2 u = 1.2742 > 1.2.

*Proof of the invariances.*
- b depends only on lengths, area and the differences t - phi.
- Dout and Chat depend only on lengths and angles.
- Cyclic relabelling permutes the same ccw edges.
- Reflection: area and lengths are preserved, and every f term is replaced by f o m = f
  (Lemma 0.2 with Lemma 0.1).

*Proof of the reduction.* Choose a longest edge and rotate it onto the x-axis, ccw. If x > d/2,
reflect in the perpendicular bisector of c_0c_1 and relabel.
- x >= 0, because the angle at c_0 is <= the angle at c_2 (it is opposite a shorter or equal
  edge), so it is <= pi/2.
- x <= d/2 <= 2.
- y <= d sqrt3/2 because |c_0c_2|, |c_1c_2| <= d.
- y = |c_0c_2||c_1c_2|/(2R) >= (2u)^2/4 = u^2.
- Nonobtuse at c_2 means y^2 >= x(d - x) >= x^2, using x <= d/2. Also y^2 >= 4u^2 - x^2. So
  y^2 >= 2u^2.

QED (this is NOTES facts 3-5, checked independently in audit 1a/2e.)

**Lemma 1.2 (the d = 2u sliver is empty).** No admissible T has its longest edge in
[2u, 2u + eps) with eps = 4e-17. The same holds for any eps < pi/(7 * 27.2), so this closes the
d-part of FLAW-2.

*Proof.* Take two heptagons H_i, H_j with disjoint interiors. A separating line with unit normal
n (pointing from H_i to H_j) gives

    r_ij >= n.(c_j - c_i) >= h_i(n) + h_j(-n),

where h is the support function about the centre. For H(c, t), h(theta) = cos(dist(theta - t,
(2pi/7)Z)), which is >= u. By concavity of cos on [0, pi/7], h - u >= (1 - u)(pi/7 - delta)/(pi/7)
with delta = dist(...).

So if r_ij < 2u + eps, then theta_n - t_i and theta_n + pi - t_j are both within
kappa eps of pi/7 (mod 2pi/7), with kappa = (pi/7)/(1 - u) = 4.532. Hence
t_j - t_i = pi/7 +/- 2 kappa eps (mod 2pi/7). The relation is symmetric under i <-> j, since
-pi/7 = pi/7 (mod 2pi/7).

If all three edges are < 2u + eps, summing over the three pairs gives 0 = 3pi/7 = pi/7
(mod 2pi/7), up to 6 kappa eps = 27.2 eps. That is impossible, since pi/7 is at distance pi/7
from (2pi/7)Z. QED

**[CLOSED 2026-11-05: re-run a1 uses outward bounds t in [-1e-12, fl(2pi/7) + 1e-12]; AUDIT2 FLAW-2 SOUND] FLAW-2 orientation sliver.** The upper bound 0.8975979010256552 used for t_1, t_2,
and for t_0 in obtuse slice 44, lies below 2pi/7 by 3.5e-17. Orientations in that gap are not
covered.
- I found no symmetry that maps the gap into the covered region while keeping x <= d/2.
  Reflection sends t to pi/7 - t, which is covered, but it also sends x to d - x.
- So (M) and (O) are certified only on the covered set until there is a re-run with outward
  orientation bounds, as the audit recommends.

**[CLOSED 2026-11-05: `bnb2.half_w`, provenance block in every log; AUDIT2 SOUND] FLAW-1 (half-widths not outward-rounded) and UNCLEAR-1 (main-sweep provenance).** These
affect the truth of the certification claims (M) and (O), not their statements. They stay open
until the re-run.

**[FLAG] Only pairwise non-overlap is used.** (M) and (O) are claimed for all admissible
configurations: only the three heptagons of T, pairwise. Neither emptiness of the circumdisc nor
the other heptagons is used. That makes them stronger than needed, which is fine. The lemmas
below supply the Delaunay information only through L3 (the angle inequality).

---

## 2. L1: saturation and the density definition

**Definition.** W_L(p) := p + [-L, L]^2 and N_P(W) := #{centres of P in W}. The upper density of
a packing P is

    delta+(P) := limsup_{L -> oo} sup_p (2L)^{-2} sum_{H in P} area(H ∩ W_L(p)),

and delta_hept := sup_P delta+(P).

The usual variants have no sup over p (windows centred at the origin), or use discs as windows.
Each is <= this one, or is handled by the same proof with discs. The bound proved in L5 holds
for every window centre with uniform constants, so it bounds every variant.

**[CLOSED 2026-11-05]** This is THE density definition of the project (NOTES now uses the same
one). The result covers the other standard variants because L5 is uniform in p.

**Lemma L1.** Every packing P is contained in a saturated packing P', meaning no heptagon (in any
position and orientation) can be added to P' without overlapping an interior. Then
delta+(P) <= delta+(P').

*Proof.*
- Order the packings containing P by inclusion.
- The union of a chain is a packing: any two of its members already lie in one element of the
  chain, so their interiors are disjoint.
- By Zorn's lemma there is a maximal element P', and maximal means saturated.
- Every window sum in the definition of delta+ is monotone under adding heptagons, so
  delta+(P) <= delta+(P'). QED

**Lemma L1'.** Let C be the set of centres of a saturated packing. Then:
- (a) C is 2u-separated: |c - c'| >= 2u for distinct c, c' in C;
- (b) every open disc of radius 2 meets C.

*Proof.*
- (a) The inscribed discs B(c, u) lie inside the heptagons, so their interiors are disjoint.
- (b) Suppose B(z, 2) ∩ C is empty, so |z - c| >= 2 for all c in C. Then the closed discs
  B(z, 1) and B(c, 1) have disjoint interiors. Any heptagon centred at z lies in B(z, 1), so it
  could be added, contradicting saturation. QED

So it suffices to prove delta+(P) <= delta_DL for saturated P. From here on P is saturated and C
is its set of centres.

**[CITE]** HK (arXiv:1602.07220v2, section 2.1) make the same two remarks without proof; they are
standard. The Zorn argument above is complete as written.

---

## 3. L2: the Delaunay facts used

Fix L >= 1 and a window centre p; take p = 0 without loss of generality.

- C_L := C ∩ W_{L+20}. It is finite by (a), and it is not collinear by (b).
- **D_L** := any Delaunay triangulation of the finite set C_L, in the sense of de Berg et al.: any
  triangulation of conv(C_L) with vertex set C_L that refines the Delaunay subdivision.
  Degenerate (cocircular) Delaunay faces are triangulated arbitrarily, using their own vertices.

**[CITE] Finite Delaunay theory.** M. de Berg, O. Cheong, M. van Kreveld, M. Overmars,
*Computational Geometry: Algorithms and Applications*, 3rd ed., Springer 2008, chapter 9:
- the Delaunay graph is plane;
- "a triangulation T of P is a Delaunay triangulation iff the circumcircle of any triangle of T
  contains no point of P in its interior";
- "T is legal iff T is a Delaunay triangulation".

These statements are confirmed verbatim in the Utrecht course slides for that chapter
(ics-websites.science.uu.nl/docs/vakken/ga/2024/slides/slides9.pdf). We rely only on the three statements as quoted above, not on any theorem numbering.

The degenerate case is covered by these statements, because de Berg defines the Delaunay
triangulation as an arbitrary refinement of the Delaunay subdivision.

We use only one consequence: **(E0)** every triangle of D_L has an open circumdisc containing no
point of C_L.

**Lemma 2.1 (triangles near the window).** Let T in D_L meet W_{L+4}. Then:
- (a) R(T) < 2, so every edge of T is in [2u, 4), diam T < 4, and the three heptagons at its
  vertices form an admissible configuration;
- (b) the open circumdisc of T contains no point of the whole of C, so T is a Delaunay triangle
  of C;
- (c) every edge of T is interior to conv(C_L), so it is shared by exactly two triangles of D_L.

*Proof.*
- (a) Take x in T ∩ W_{L+4}, and let D = B(z, R) be the open circumdisc, so |x - z| <= R.
  Suppose R >= 2.
  - If |x - z| <= R - 2, put z' = x.
  - Otherwise put z' = z + (R - 2)(x - z)/|x - z|.
  - In both cases B(z', 2) ⊂ D and |z' - x| <= 2, so B(z', 2) ⊂ W_{L+8} ⊂ W_{L+20}.
  - By L1'(b) the disc B(z', 2) contains a point of C, which lies in C_L ∩ D. That contradicts
    (E0).

  So R < 2. Edges are chords, hence < 4, and they are >= 2u by L1'(a). The heptagons of a packing
  are pairwise interior-disjoint.
- (b) D ⊂ B(x, 4) ⊂ W_{L+8}, so D ∩ C = D ∩ C_L, which is empty.
- (c) Suppose a point q of T lay on the boundary of conv(C_L). A supporting line through q gives
  an open half-plane with no point of C_L in it. That half-plane contains the disc
  B(q + 2.5n, 2), where n is the outward normal, and this disc lies in W_{L+4+4.5} ⊂ W_{L+20}.
  By L1'(b) the disc contains a point of C, and so of C_L: a contradiction.

  So T lies in the interior of conv(C_L), and each of its edges is shared by two triangles. QED

**Lemma 2.2 (local Delaunay angle condition).** Let pq be an edge of D_L that is interior to
conv(C_L), with adjacent triangles pqa and pqb. Then angle(paq) + angle(pbq) <= pi, with equality
iff a, b, p, q are cocircular.

*Proof.* a and b lie on opposite sides of the line pq. By (E0), b is not in the open circumdisc of
pqa.

By the inscribed angle theorem, a point x on b's side of the line sees pq at angle:
- exactly pi - angle(paq) if x is on the circle;
- more than that if x is inside the disc;
- less than that if x is outside.

So angle(pbq) <= pi - angle(paq). Equality forces b onto the circle. QED

(This is de Berg's legal-edge criterion; the slides cited above state the inscribed angle
theorem in this form.)

**Corollary 2.3.**
- (a) Two triangles sharing an edge are not both obtuse at the angles opposite it.
- (b) If T' is obtuse with obtuse angle gamma' opposite the shared edge e, and T is the other
  triangle on e with angle alpha opposite e, then cos gamma' >= -cos alpha. In particular
  cos alpha > 0.

*Proof.*
- (a) Two angles above pi/2 would sum to more than pi.
- (b) gamma' <= pi - alpha, and cos is decreasing on [0, pi]. QED

**Degenerate cells.** Inside a cocircular face, each diagonal has opposite angles summing to
exactly pi (Lemma 2.2, equality case). All the facts above hold with no further hypothesis.
"Consistency" across cells is automatic, because a single triangulation D_L is fixed per window.
The per-triangle statements (M) and (O) do not care which triangulation of a face was chosen:
they are statements about single admissible triangles, together with Corollary 2.3 for each
neighbour.

**[CITE] Remark (not needed).** For the infinite Delone set C, the Delaunay subdivision exists
and is a locally finite, face-to-face tiling by convex polygons with at most 6 vertices. This goes
back to B. Delaunay, "Sur la sphère vide", *Izv. Akad. Nauk SSSR, Otd. Mat. Estestv. Nauk* 7
(1934), 793-800 (bibliographic data as cited in Edelsbrunner, Garber and Saghafian,
arXiv:2410.21204). **[FLAG]** I did not verify Delaunay's exact hypotheses from the original text.
The window construction above avoids the theorem, so nothing rests on it.

---

## 4. L3: the bookkeeping

Work in D_L, and let T_L* be the set of triangles of D_L that meet W_{L+4}. By Lemma 2.1 these
are admissible, and each of their edges has exactly two sides.

**Lemma 3.1 (obtuse triangles).** Let T in T_L* be obtuse at gamma. Then:
- its longest edge l_T is the edge opposite gamma, and it is the unique longest edge;
- |l_T| > RTH;
- the triangle N(T) on the other side of l_T exists and is unique.

*Proof.* With a and b the other two edges:
- c^2 = a^2 + b^2 - 2ab cos gamma > a^2 + b^2 >= max(a, b)^2 + (2u)^2, so c is the unique
  longest edge;
- c^2 > 2(2u)^2, so c > RTH.

N(T) exists and is unique by Lemma 2.1(c). QED

**Definition (the transfer, exact).** For every obtuse T in T_L*, move the amount

    tr(T) := D5(|l_T|, cos gamma_T)

from N(T) to T. This uses the real bilinear interpolant of Section 0.2, with D5 := 0 for
|l_T| < rb0. Then for T in T_L*:

    b*(T) := b(T) + Dout(T) - Din(T),   where Dout(T) = tr(T) [0 if T is not obtuse],
                                        Din(T) = sum of tr(T') over obtuse T' with N(T') = T.

Din(T) involves only triangles T' that share an edge with T. Those meet W_{L+4} whenever T
meets W_L, so Din(T) is defined for every T meeting W_L.

**Lemma 3.2 (each receiver is charged at most Chat).** For every T in T_L* meeting W_L,
Din(T) <= Chat(T).

*Proof.*
- If T' donates to T, then e = l_{T'} is an edge of T.
- By Corollary 2.3(a), e is not the long edge of T if T is itself obtuse. So e is one of the edges
  that Chat sums over.
- Each edge of T has at most one donor: the single triangle on its other side.
- For that donor, |e| = |l_{T'}| and cos gamma_{T'} lies in [-cos alpha_e, 0), by Corollary
  2.3(b). So tr(T') = D5(|e|, cos gamma_{T'}) <= Dhat(|e|, cos alpha_e).
- Sum over the edges. QED

**Lemma 3.3 (per-triangle bound, consuming L4).** For every T in D_L meeting W_L,
b*(T) >= tau >= acrit.

*Proof.* T is admissible (Lemma 2.1).
- **Case: longest edge < RTH.** By Lemma 3.1, T is not obtuse, so Dout(T) = 0. A donor T' would
  have l_{T'} > RTH as an edge of T, longer than T's longest edge; so Din(T) = 0. Then
  b*(T) = b(T) >= tau by (M).
- **Case: longest edge >= RTH > 2.548.** By (O) and Lemma 3.2:

      b*(T) = b(T) + Dout(T) - Din(T) >= b(T) + Dout(T) - Chat(T) >= tau.

QED

**Lemma 3.4 (transfers cancel).** Let F be a finite set of triangles of D_L meeting W_L, and
call an edge of a triangle in F a *boundary edge of F* when the triangle on its other side is not
in F. Then

    | sum_{T in F} (b*(T) - area(T)) | <= (Fmax + Dmax) * #(boundary edges of F),

with Fmax = 0.2436 and Dmax = 0.2335.

*Proof.* Expand sum (b* - area) as a sum over pairs (T, e) with e an edge of T:
- the f term f(e as a ccw edge of T);
- plus tr(T) if e = l_T and T is obtuse;
- minus tr(T') if T' is across e, T' is obtuse and l_{T'} = e.

Take an edge e with both sides T_1, T_2 in F:
- The ccw orientations from T_1 and from T_2 are opposite, so their f terms are f(p) and
  f(Rp) = -f(p) (Lemmas 0.1 and 0.2), and they sum to 0.
- At most one of T_1, T_2 is obtuse opposite e (Corollary 2.3(a)). If one is, its tr appears
  once with + (Dout) and once with - (Din of the other), and they cancel.

Every edge with only one side in F contributes one f term (|f| <= Fmax), plus at most one tr
term (0 <= tr <= Dmax). QED

**Dropped from NOTES L3: the "chains absorbed by recursion over r" remark.** It is not used and is
not needed.
- The proof never orders triangles and never follows chains.
- An obtuse triangle can donate across its long edge while receiving across one or both short
  edges (short edges can be >= rb0). (O) charges this directly, because `inc_short` is added in
  the obtuse branch too (`bnb2.py:343`).
- A nonobtuse triangle can receive across all three edges, and (O) charges all three
  (`bnb2.py:339, 344`).
- How D2-D5 were built (sampling, recursion, padding) does not matter for soundness. Only the
  final table D5 enters, through Lemma 3.2 and statement (O).

---

## 5. L5: averaging

**Theorem.** For a saturated packing P and any window centre p,

    (2s)^{-2} sum_H area(H ∩ W_s(p)) <= A/(2 acrit) + O(1/s),

with an O-constant independent of p and of P. So delta+(P) <= A/(2 acrit) = delta_DL, and by L1
the same bound holds for every packing.

*Proof.* Put p = 0 and L := s + 5, and let N(t) := #(C ∩ W_t).

1. **Area side.** A heptagon meeting W_s has its centre in W_{s+1}, because H ⊂ B(c, 1). So

       sum_H area(H ∩ W_s) <= A N(s + 1) = A N(L - 4).

2. **The window family.** Let F_L := {T in D_L : T ⊂ W_L}. These triangles have disjoint
   interiors inside W_L, so sum_{F_L} area(T) <= 4L^2.
3. **Angle count (the Euler step).** Take v in C ∩ W_{L-4}.
   - By Lemma 2.1(c) and its proof, v is interior to conv(C_L), so the triangles of D_L at v
     surround it, and their angles at v sum to 2pi.
   - Each of these triangles contains v, so it meets W_L, has R < 2, and lies in the open disc of
     radius 4 about v, hence in W_L. So each lies in F_L.
   - Summing all angles of all triangles in F_L therefore gives

         pi #F_L >= 2pi N(L - 4),   i.e.   #F_L >= 2 N(L - 4).

4. **Boundary edges.** Let e be a boundary edge of F_L, with T in F_L on one side and T' not in
   F_L on the other.
   - T' meets W_L (it contains e), so diam T' < 4.
   - T' is not inside W_L, so it has a point q with |q|_inf > L.
   - Both endpoints of e are within distance < 4 of q, so they lie in the strip
     S := W_L \ int W_{L-4}.
   - Points of C are 2u-separated, so #(C ∩ S) <= area(S + B(0, u))/(pi u^2)
     <= (8 + 4u)(4L)/(pi u^2).
   - The degree in D_L of a vertex, counting edges of triangles meeting W_L, is at most
     #(C ∩ B(v, 4)) - 1 <= ((4 + u)/u)^2 - 1 < 29, because the discs B(c, u) are disjoint.
   - So #(boundary edges) <= 29 (8 + 4u)(4L)/(pi u^2) =: kappa L.
5. **Combine.** Lemmas 3.3 and 3.4 with F = F_L give

       acrit #F_L <= sum_{F_L} b*(T) <= 4L^2 + (Fmax + Dmax) kappa L.

   Then by step 1:

       sum_H area(H ∩ W_s) <= A N(L - 4) <= A #F_L / 2 <= A (4L^2 + c' L)/(2 acrit).

   Dividing by 4s^2 with L = s + 5 gives A/(2 acrit) + O(1/s). The constants depend only on u, A,
   Fmax and Dmax. QED

**Remark (Euler form).** In a finite triangulated disc, #triangles = 2 #(interior vertices)
+ #(boundary vertices) - 2. The angle-sum inequality in step 3 is the half of this that we need,
and it needs no assumption about the topology of the union of F_L, which can have holes or
pinches. So the angle-sum form is preferred to NOTES' "Euler, + O(L)".

**[FLAG] Where acrit enters.** acrit is defined as A/(2 delta_DL), and the target tau is >= acrit
(Section 0). L5 never needs area(DL triangle) = acrit. That identity is needed only for the chart
lemma (Section 6, hypothesis H-DL), where beta(0) must equal acrit.

---

## 6. The chart/tube lemma in convexity form (UNCLEAR-3)

### 6.1 Setting

- **Chart frame** (`local_lemma_arb.py`): H0 = H(0, 0). Ht1 is centred at P1 + (z_1, z_2) with
  orientation z_3, where P1 = (1 + u, 0). H2p2 is centred at P2 + (z_4, z_5) with orientation
  pi/7 + z_6, where P2 is the KK point (runs/dl_ref.json). z is in R^6.
- **Contacts:** g_k(z), k = 1..5, are the vertex-to-edge-line functions of the five DL contacts
  (CONTACTS = (0,2,4,1), (1,0,0,3), (1,2,5,2), (2,0,1,4), (2,1,2,5)).
- **Coordinates:** Phi(z) := (g_1, ..., g_5, v.z), where v is the exact null vector of DG(0).
  b(z) is b of this configuration, with f at mu = mu*. All of these are real-analytic near z = 0.

**Hypothesis H-DL (symbolic; NOT proved here; audit "remains" item 2).**
- g(0) = 0 exactly;
- rank DG(0) = 5, with ker DG(0) = R v;
- b(0) = acrit;
- grad b(0) . v = 0 at mu = mu*.

Notes on H-DL:
- b(0) = acrit follows from area(DL triangle) = A/(2 delta_DL) plus f = 0 on the three DL edge
  configurations.
  - The area identity holds because all Delaunay triangles of the DL are congruent and there are
    2 per heptagon.
  - f vanishes on the two point-reflected edges (R-fixed) and on the translation edge (S-fixed).
- **[FLAG]** The translation edge has configuration (1 + u, 0, 0) only because the DL translation
  vector is exactly parallel to a vertex direction of H0. That is part of the symbolic check. For
  a general a = b, f(r, a, a) is NOT zero (float check: f(0.3, 0.3) = 2.04 for a random g), so
  this cannot be dropped.
- grad b(0) . v = 0 is the exactification. mu* exists and lies in the ball, given that v is exact
  and the arb solve enclosed it.
- **[FLAG]** The code's V is an arb ball from a linear solve. The proof needs the exact v to lie
  in that ball, and needs every chart quantity (J_c, the windows) to be read with the exact v.
  That holds under ball semantics if DG(0) is enclosed correctly. It should be stated in the
  symbolic note (`src/dl_exact.py` appears to be in progress elsewhere).

**Certified per chart c** (`local_tube_arb.chart`, arb at 128 bits; read with ball semantics, so
every statement holds for the exact matrix M_c = J_c^{-1}, J_c = DPhi(z_c)):
- K_c := { z : |J_c (z - z_c)|_inf <= rho_c };
- **(A)** On K_c, for each contact, the touching vertex lies strictly inside the other six
  half-planes of the owner. So g_k(z) < 0 puts that vertex in the interior of the owner, which is
  an overlap. Hence a non-overlapping z in K_c has g(z) >= 0.
- **(B)** ||DPhi(z) J_c^{-1} - I||_inf <= q_c < 1 on K_c.
- **(C)** w(z) := DPhi(z)^{-T} grad b(z) has w_1, ..., w_5 > 0 on K_c.
- **(D)** On K_c, the second derivative along the contact curve,

      beta2(z) := xi^T Hess b(z) xi + grad b(z) . xi2,
      xi = DPhi(z)^{-1} e_6,   xi2 = -DPhi(z)^{-1} (D^2 g(z)[xi, xi], 0),

  is >= B2_c > 0 (the `bpp` of the code). Recorded minimum: 0.0135 (tube b) and 0.0723 (tube a).

**Chart region used by the B&B** (`chart_excl.inside`):

    E_c := { z in K_c : g_k(z) <= eta_c (k = 1..5), |v.(z - z_c)| <= eta_c },
    with the s-window I_c := [v.z_c - eta_c, v.z_c + eta_c].

**Extra conditions the proof needs** (all CERTIFIED as of 2026-11-05; see below):
- **(E) offset, strict:** eta_c + off_c < (1 - q_c) rho_c, where off_c := |g(z_c)|_inf.
  **[CLOSED 2026-11-05]** The re-run tube (`runs/charts_w030c.json`) uses eta_c = dn((1 - q_c) rho_c - 2 off_c); exact slack >= 2.98e-21 (AUDIT2). FLAW-4 (exported data rounded, not outward) is closed by `src/tube_recheck_arb.py` (log `runs/tube_recheck_w030c.log`): every exported number widened by 2 ulp, (E') checked on the inflated cube K'_c = {|J_c(z - z_c)| <= rho_c (1 + 5e-10)} (Lemma 6.0 below). Historical note: the old code had eta_c = (1 - q_c) rho_c exactly, so (E) failed formally. Recorded
  margin (float): min over charts of (1 - q) rho - eta - off = -6.7e-18 (tube b), -1.8e-18
  (tube a).
  - Cheapest fix: re-certify each chart with rho'_c = rho_c (1 + 1e-9), keeping the exported
    eta_c and the chart_excl test (which uses rho_c, not rho'_c) unchanged.
  - Then (E) holds with rho'_c in place of rho_c, as long as the re-certified q' stays below
    1 - (eta_c + off_c)/rho'_c, which leaves about 1e-9 rho of slack. All the arguments below run
    on the larger cube K'_c ⊃ K_c.
  - The B&B does not need re-running for this.
- **(G) gluing:** for consecutive charts c -> c' in each walk of each tube file, the curve point
  gamma_{c'}(v.z_{c'}) lies in K_c. It is enough that

      |J_c (z_{c'} - z_c)|_inf + ||J_c J_{c'}^{-1}||_inf tiny_{c'} <= rho_c.

  **[CLOSED 2026-11-05: certified in arb by `src/tube_glue_arb.py` (log `runs/tube_glue_w030c.log`) and again with 2-ulp widened data by `src/tube_recheck_arb.py`; 1529 pairs, 0 failures, worst ratio 0.8744; uses 2 off_{c'} >= off_{c'}/(1 - q) in place of tiny_{c'}]** Original float evidence (`scratchpad/hp_proofs/chk_glue.py`, mid-point matrices):
  the left side is <= 0.874 rho_c in both tubes, and the first five components are <= 4e-14 rho_c.
  So an arb check will almost certainly pass with a wide margin. It needs only the JSON data plus
  one interval matrix product per pair.
- **(Ov) windows overlap:** I_c ∩ I_{c'} is nonempty for consecutive charts.
  **[CLOSED 2026-11-05: certified in arb, worst gap upper bound -6.39e-5; the window check s* = v.z_{c'} in I_c (worst |v.(z_{c'} - z_c)|/eta_c = 0.900) also certified]** Original: the worst gap is -6.4e-5 (an overlap) in both tubes, which
  agrees with the audit. Redo it with exact endpoints. The margin is 1e11 times the float error,
  so this is a formality.

### 6.0 Lemma (inflated cube; closes AUDIT2 FLAW-4)

Let n = 2 and h = 2 rho / n = rho. `local_tube_arb.chart` certifies (A)-(D) on the 2^6 arb boxes
with centres -rho + h(i + 1/2) (i = 0, 1) and radius r = (h/2)(1 + 1e-9) in each coordinate of zeta.
Their union is the cube {|zeta|_inf <= rho/2 + r} = {|zeta|_inf <= rho (1 + 5e-10)} =: K'_c,
because in each coordinate the two intervals [-rho/2 - r, -rho/2 + r] and [rho/2 - r, rho/2 + r]
overlap (r > rho/2) and together cover [-rho/2 - r, rho/2 + r]. So (A)-(D) hold on K'_c with the
EXACT chart matrix J_c (the arb computation encloses it), and Steps 1-5 below run on K'_c.
The B&B (`chart_excl`) tests membership with the exported J'_c, |J'_c - J_c| <= 2 ulp entrywise.
For every excluded z: |J_c(z - z_c)|_inf <= rho + d and |v.(z - z_c)| <= eta + d, where
d = 6 max|J_c - J'_c| ||J'^{-1}_c|| rho <= ~1e-18. `src/tube_recheck_arb.py` certifies, for all
1530 charts, rho + d <= rho (1 + 5e-10) and eta + max(off+, d) <= (1 - q+) rho (1 + 5e-10) with
q+, off+ the exported values raised by 2 ulp (worst ratio 0.9999999995). So every excluded z is
in K'_c and every Banach target used in Steps 1 and 5 is within reach. QED

### 6.2 The argument

**Step 1 (Banach chart).** Put Psi_c(zeta) := Phi(z_c + M_c zeta) on the cube
Q := {|zeta|_inf <= rho_c}. By (B), ||DPsi_c - I|| <= q_c on Q. Then:
- Psi_c is injective on Q. Indeed, |Psi(zeta) - Psi(zeta') - (zeta - zeta')| <= q|zeta - zeta'|
  by the mean-value inequality on the convex set Q.
- For every p with |p - Psi_c(0)|_inf <= (1 - q_c) rho_c, there is a zeta in Q with Psi_c(zeta) = p.
  Proof: zeta -> zeta - (Psi_c(zeta) - p) maps Q into Q and is a q-contraction.
- That zeta satisfies |zeta| <= |p - Psi_c(0)|/(1 - q_c). So if the inequality on p is strict,
  zeta lies in int Q, Psi_c^{-1} is C^2 near p (inverse function theorem), and DPhi is invertible
  at the preimage.

**Step 2 (curve pieces).** For s in I_c, the target p = (0, s) has
|p - Phi(z_c)| <= max(off_c, eta_c) < (1 - q_c) rho_c, by (E). So
gamma_c(s) := z_c + M_c Psi_c^{-1}(0, s) is the unique point of K_c with Phi = (0, s).
gamma_c is C^2 on I_c, and differentiating Phi(gamma_c(s)) = (0, s) gives

    gamma_c'(s) = DPhi(gamma_c(s))^{-1} e_6 =: F(gamma_c(s)).

**Step 3 (one curve; uses (G) and (Ov)).**
- DPhi is invertible on the open set U ⊃ ∪_c K_c, so F is C^1, hence locally Lipschitz, on U.
- Let Gamma be the maximal solution of z' = F(z) in U with Gamma(0) = 0.
- For the central chart c_0: Phi(0) = (0, 0) by H-DL and 0 is in K_{c_0}, so gamma_{c_0}(0) = 0
  by uniqueness. By Picard-Lindelöf uniqueness, gamma_{c_0} = Gamma on I_{c_0}.
- Induction along each walk: (G) gives a point s* in I_c ∩ I_{c'} (namely s* = v.z_{c'}, which is
  in I_c by the tube spacing 0.9 eta) with gamma_{c'}(s*) in K_c. Uniqueness in K_c then gives
  gamma_{c'}(s*) = gamma_c(s*) = Gamma(s*). ODE uniqueness gives gamma_{c'} = Gamma on I_{c'}.
- So Gamma is defined on J := ∪_c I_c, which is an interval by (Ov), and gamma_c = Gamma on I_c for
  every chart c.
- This step is what makes "beta" one function. Without it, two charts could in principle follow
  different branches of {Phi = (0, s)}. **[FLAG]** NOTES and the audit both assume the curve is a
  single one without saying why.

**Step 4 (convexity).** Put beta := b o Gamma on J.
- beta is C^2, and beta''(s) = beta2(Gamma(s)) >= B2_c > 0 whenever s is in I_c, because
  Gamma(s) lies in K_c.
- So beta' is increasing on J. Also beta'(0) = grad b(0) . v/|v|^2 = 0 (H-DL), using
  Gamma'(0) = DPhi(0)^{-1} e_6 = v/|v|^2.
- Hence beta' <= 0 on J ∩ (-oo, 0] and beta' >= 0 on J ∩ [0, oo), which gives
  **beta(s) >= beta(0) = acrit for all s in J**.
- **The per-chart bounds m_c (`excess_lower`, negative for 2 charts in each file) are not used.**

**Step 5 (monotone path off the curve; uses (A), (C), (E)).** Let z be in E_c and
non-overlapping, and put s = v.z (which is in I_c).
- By (A), g(z) is in [0, eta_c]^5.
- For t in [0, 1] the targets p(t) = (t g(z), s) satisfy
  |p(t) - Phi(z_c)| <= max(eta_c + off_c, eta_c) < (1 - q_c) rho_c, by (E).
- So z(t) := z_c + M_c Psi_c^{-1}(p(t)) is a C^1 path in K_c, with z(1) = z (injectivity) and
  z(0) = Gamma(s) (Step 3).
- Differentiating Phi(z(t)) = p(t) gives DPhi z' = (g(z), 0), so

      d/dt b(z(t)) = grad b . DPhi^{-1} (g(z), 0) = sum_{k <= 5} w_k(z(t)) g_k(z) >= 0,

  by (C).
- Hence **b(z) >= b(Gamma(s)) = beta(s) >= acrit**.

**Theorem 6.1.** Assume H-DL and, for each tube file, (A)-(D) for every chart used, together with
(E), (G) and (Ov). Then every non-overlapping configuration in the union of the regions E_c has
b >= acrit.

Combined with the rigorous canonical-box -> z map of `chart_excl.to_z` (audit 2e: the same f,
the cyclic relabelling and the KROT rotation all match), every canonical box the B&B discards by
charts has b >= acrit everywhere on its feasible part. That is the "chart lemma" part of (M).

**Note on targets.** The B&B discards chart boxes without comparing against tau. The chart part
of (M) is therefore b >= acrit, not b >= tau. That is harmless, since L5 needs only acrit: state
(M) as "b >= tau off the chart regions, b >= acrit on them". Lemma 3.3 holds verbatim with acrit.

**[RESOLVED by run a1] Which tube was used where.** Run a1 (the certificate of record) uses only `charts_w030c.json` (1530 charts), logged with its hash in every main-sweep log. Historical: the old DL-slice chain used `charts_w030.json` (459 charts) in
its early steps and `charts_w030b.json` (1530 charts) after the extension. Each file is its own
tube from the same central chart, and Theorem 6.1 applies to each separately. The proof does not
need the two to agree. The replay should record which file each step of the k39 chain loaded.

**[FLAG] FALSE-AS-STATED in isolation (code docstrings, not the mathematics).**
- `chart_excl.py`'s docstring says each certified chart region has b >= acrit "by the chart
  lemma". For a single non-central chart that does not follow: its own bound m_c can be negative
  (2 charts). It is true only through Steps 3-4, which need (G) and (Ov).
- NOTES' "beta is convex on the union with beta'(0) = 0" likewise silently assumes the
  single-curve fact (G).

---

## 7. Status

| item | status | what it rests on |
|---|---|---|
| Lemma 0.1/0.2 (symmetries, reversal, mirror) | PROVED | algebra above |
| Lemma 1.1 (canonical reduction, domain bounds) | PROVED | |
| Lemma 1.2 (d = 2u sliver empty) | PROVED | support-function argument; closes the d-part of FLAW-2 |
| L1 (saturation, density) | PROVED | Zorn; density definition fixed in section 2 (shared with NOTES) |
| L2 (Delaunay facts, degenerate cells) | PROVED | finite Delaunay theory (de Berg et al. ch. 9, statements quoted in section 3; no theorem numbers relied on); windowed construction, no infinite triangulation |
| L3 (bookkeeping, D5 transfer, unique neighbour, cancellation) | PROVED, given (O) | |
| L5 (averaging) | PROVED, given L3 and (M)/(O) | angle-sum count |
| Chart/tube (convexity form) | PROVED | (E)/(E') certified (charts_w030c + tube_recheck_arb, Lemma 6.0), (G) and (Ov) certified in arb, H-DL symbolic (`src/dl_exact.py`; null vector via arb enclosure, AUDIT2 section 3) |
| Interface (M), (O) as certified | stated precisely; certified by run a1 (bnb2 9f4fb628), pending UNCLEAR-4 replays | FLAW-1, FLAW-2, UNCLEAR-1 closed; UNCLEAR-4 (six checkpoint-resumed slices) being replayed from scratch |

Nothing in L1-L3 or L5 is false as stated in NOTES. Corrections to NOTES (ALL APPLIED 2026-11-04/05):
- the "recursion over r" remark should be dropped (L3);
- "bins are half-open" is stale for D5;
- the claim "b >= acrit on the whole tube" needs the extra hypotheses (E), (G) and (Ov);
- the density definition was never fixed.

Scratch checks (outside the repo): `<scratch directory outside the repo>`
and `...\hp_proofs\chk_sym.py`. Both are float only, under 1 s and < 100 MB.

# Heptagon packing lane: running log

Target: every packing of congruent regular heptagons has density <= delta_DL =
(2/97)(-111+492u-356u^2), u = cos(pi/7), = 0.892690686126508 (Kuperberg double lattice).
Phase goal (3 days from 2026-10-31): kill test, not the proof.

All numerics below are floating-point and multistart local optimisation: **unverified**
(no interval arithmetic) unless marked otherwise.

## 2026-10-31 (Sat), Step 0: collision check

Searches run 2026-10-31 13:00-13:40 EDT.

- arXiv API: all:heptagons (50 newest), heptagon AND packing, "double lattice" AND packing,
  abs:pessimal, au:Kusner, au:Kallus AND packing. Only the known papers turned up
  (1305.0289 Kallus, 1509.02241 Kallus-Kusner, 1602.07220 Hales-Kusner, 2207.08959 Torda
  plane-group packings, 2212.13849 RSA). No claim.
- Zenodo API: heptagon packing, "regular heptagons", heptagon densest, title:heptagon*,
  title:heptagons, creators Tomar / Tejasvi, "heptagon" AND "packing density", Kallus AND
  Kusner, "double lattice" AND polygon, heptagon AND Kuperberg, pessimal packing,
  double-lattice packing. No claim. Near-misses, all unrelated: 10.5281/zenodo.23066520
  (anonymous, 2026-09-30, local minimality of heptagon/octagon for area x Dirichlet eigenvalue),
  23824544 / 23824216 (regular-polygon nesting/containment, 2026-10-30), 7577671 (circles in
  a heptagon). No Tomar record on packings found.
- GitHub: repos "heptagon packing", "heptagon density", "double lattice packing",
  "polygon packing kuperberg" -> nothing. PRs and issues for the same terms -> nothing
  relevant (one PR "Publish certified local minima for the regular heptagon and octagon",
  ipitchford/evidence-press-site, 2026-09-30 = the eigenvalue record above).
- **Crowding signal:** elder-plinius/V3R1T4S (AI "known unknowns" bank, pushed 2026-10-30)
  lists `known-unknowns/packings-physics/regular-heptagon-densest-packing.md` as a
  "moonshot", tractability 1/5, framed as a hunt for a packing DENSER than the double lattice
  (the disproof direction), "open as of 2026-10-30". Not a claim, but AI agents are being
  pointed at this problem from this week.
- AI catalogues: MColbrook/AIM (tree + code search: no heptagon/double-lattice item; nearest
  are 238 tetrahedron packing, 349 Ulam convex-solid packing); openai/math (tree grep:
  nothing on heptagons; only a planar circle-packing Fourier challenge); Epoch FrontierMath
  open problems (10 of 49 titles visible on the page, none relevant; full list not visible:
  partly unverified); erdosproblems: site search returned nothing (not an Erdos problem).
- Web search (2 queries): only the 2014-2016 sources; all say open.

Verdict Step 0: **no claim, nobody visibly close on the proof direction.** Proceed.

## Step 1: reduction design and numerical kill tests

### Basic constants (float; closed form reproduced to 1e-15)

- Heptagon circumradius 1: area A = 2.7364101886, inradius u = 0.9009688679.
- KK half-length parallelogram: a = 0.12651019814, b = 0.31174490093, area 1.5326754447.
- delta_DL = A/(2*1.5326754447) = 0.892690686126508 (matches closed form).
- **All Delaunay triangles of the double lattice are congruent**, area
  acrit = A/(2 delta) = 1.5326754447, edges 1.83132 / 1.90097 / 1.91498,
  angles 57.36 / 60.94 / 61.71 deg, circumradius 1.0874. Every heptagon touches 6 neighbours.
  Edge types: 1.83132 and 1.91498 join point-reflected heptagons, 1.90097 joins translates.
- Saturation gives Delaunay circumradius <= 2 (same argument as Hales-Kusner).
  Minimum centre distance 2u = 1.80194 (edge-to-edge contact, needs opposite parity).

### Single triangles (src/min_tri.py)

| quantity | heptagon | pentagon (HK) |
|---|---|---|
| min area, nonobtuse Delaunay triangle | 1.484980 (0.9689 acrit) | 1.23719 (0.9588 acrit) |
| min area, any triangle with R <= 2 | 1.376913 (0.8984 acrit), obtuse, R = 2 | ~0.98 (0.76 acrit) |

Minimal nonobtuse triangle: edges 1.8031 / 1.9265 / 1.8329; two point-reflected edges and
one translation edge (same parity pattern as the DL triangle). The single-triangle gaps are
no worse than the pentagon's.

### Dimers (src/min_dimer.py)

- Two adjacent nonobtuse Delaunay triangles: min total 2.969960 = 2 x 1.48498
  (0.9689 x 2acrit). Two minimal triangles glue along an edge. **Dimers alone give nothing.**
  (Pentagon: the ice-ray dimer is a local minimum and pseudo-dimers fall below only by ~1e-5.)
- Consistent with Kallus-Kusner sec. 3.2: for heptagons the linearised cell (dimer) problem
  is unbounded; they need auxiliary transfer terms mu_r g_r on the 4 outer contact edges,
  cancelling between neighbouring cells, even for LOCAL optimality.

### Vertex stars (src/min_star.py)

- Star of 6 triangles around one heptagon, started at the DL star (area 6 acrit = 9.19605):
  local descent reaches 9.13805 (0.99369). **The DL star is not a local minimum of area.**
  So plain-area clusters (dimer, star) cannot give a sharp bound; a transfer/redistribution
  rule is mandatory, as in KK (local) and HK (b(T) with N, M corrections).

### Obstruction test for "triangle + antisymmetric pair transfer" (src/balanced_dimer.py)

Certificate class: b(T) = area(T) + sum over ccw edges (i->j) of f(pair config i->j), with
f(i->j) = -f(j->i) a function of the two heptagons only. Pair config (r, a, b), r = centre
distance, a, b = orientations relative to the edge direction mod 2pi/7.
Reverse map: (r, a, b) -> (r, b + pi/7, a + pi/7). Point-reflected pairs are fixed by it, so
**f = 0 on every point-reflected pair**. A dimer sharing a translation edge whose four outer
edges are point-reflected is therefore a balanced 2-triangle measure: any such f sums to 0
over it, so min b <= average area.

- Allowing obtuse triangles: min total 2.819382 (0.91976 x 2acrit): two obtuse triangles
  (edges 1.901 / 3.285 / 1.806, R = 2). **Kills pure pair transfer** once obtuse triangles
  are included.
- Nonobtuse only: 300/300 starts converge to exactly the DL dimer, 3.0653508893 = 2acrit.
  **This obstruction does not occur among nonobtuse triangles.** Same split as Hales-Kusner:
  obtuse triangles need separate handling (their N correction: the long-edge neighbour of an
  obtuse triangle has surplus).

### LP test of pair-transfer certificates, nonobtuse triangles (src/lp_transfer.py) - IN PROGRESS

- Run 1 (DR=3 Chebyshev in r, DF=2 Fourier in 7a, 7b, 100 antisymmetrised basis functions,
  box 1.0, 40k sampled triangles + DL triangles, 19 cutting-plane rounds, 40 min, lease
  heptagon-4383): the LP value stays at exactly acrit (it is pinned by the DL triangle), but
  local search with f fixed keeps finding triangles with b far below acrit (from -0.18 up to
  +1.24). **Not converged and not informative yet**: the LP picks an extreme point of a
  degenerate optimal face, so f is huge wherever sampling is thin.
- Fix written: phase 2 keeps t >= t* - ETA and minimises ||c||_1 (smallest transfer); box 0.3.
  Run 2 (`lp_transfer.py 3 2 60 2500 nonobtuse 0.3 0.002`, lease heptagon-4388) was
  **killed at startup for Lee's reboot (2026-10-31 ~15:00 EDT)**, before any round finished.
  Lease released. Nothing auto-relaunches.

### Pilot interval B&B (src/iv_pilot.py) - written, NOT yet run

Lemma: nonobtuse triangle area >= TARGET, 6-dim box space, vectorised numpy intervals with
nextafter outward rounding (cost pilot, not a rigorous certificate). Plan: run TARGET = 1.40,
1.45, 1.47, 1.48 on a t0 sub-range to get box count vs margin, then extrapolate.

## Re-arm after reboot

1. Check GAMING flag and `broker.py status` (1 thread if GAMING; 2 after 18:00 EDT otherwise).
2. Re-run LP run 2 under lease (about 45 min):
   `PY=/c/Python314/python lease-run.sh --project heptagon --min-gb 1 --minutes 75 --eta-minutes 42 --match lp_transfer.py -- /c/Python314/python src/lp_transfer.py 3 2 60 2500 nonobtuse 0.3 0.002`
3. Iv pilot: short runs (< 15 min each) or leased.

## 2026-10-31 after the reboot (resumed about 15:00 EDT; no GAMING flag, up to 2 threads)

### LP run 2: pair-transfer certificate for nonobtuse triangles CONVERGES (unverified, float)

`lp_transfer.py 3 2 60 2500 nonobtuse 0.3 0.002`, lease heptagon-4390, 42 min, 11 rounds.
Phase 1 LP value stays at exactly acrit through every round, so the cuts never found a
triangle that forces t below acrit. Local search with f fixed (min b / acrit):
0.9783, 0.9945, 0.9968, 0.9985, 0.99868, 0.99850, 0.99870, 0.99870, 0.99871, 0.99868, 0.99869.
That is the floor acrit - ETA = 0.99870 exactly: the phase-2 L1 step deliberately spends the
allowed ETA = 0.002 slack.

Independent check (`src/validate_f.py`, seed 777, LP used 11): fresh samples give min
b = 1.00358 acrit (min area in that sample 0.97494 acrit). Local search from the 60 worst
plus 40 random samples gives min b = 0.998695 acrit, attained at d = 1.90372,
angles (0.50258, 0.50257, 0.05378). That has the double-lattice parity pattern and sits next
to the DL triangle. So **apart from the DL neighbourhood, b >= acrit with this 100-function
basis.** The only place the slack is used is the tangency at the double lattice, which is
what Kallus-Kusner's local theory covers (null direction = KK parallelogram flex).

### Obtuse triangles (`src/obtuse_pair.py`)

Obtuse triangle plus its neighbour across the longest edge (Delaunay, R <= 2): min total
3.41679 = **1.1146 x 2acrit**, an 11% surplus. The HK-style "N" rule (the obtuse triangle
draws from its long-edge neighbour) has plenty of room. Not yet handled: a neighbour shared
by several obtuse triangles, and chains of them (HK Lemmas 76-85 analogues).

### Interval pilot (`src/iv_pilot.py`, numpy intervals with nextafter rounding: a COST PILOT, not rigorous)

Lemma: every nonobtuse triangle of three non-overlapping heptagons has area >= T.
True minimum (float) 1.484980.

| T | margin | boxes | seconds (1 thread) | result |
|---|---|---|---|---|
| 1.45 | 0.035 | 238,241 | 8.3 | certified |
| 1.47 | 0.015 | 990,919 | 37 | certified |
| 1.48 | 0.005 | 3,643,589 | 138 | certified |
| 1.484 | 0.001 | > 20,793,453 | 780 (stopped, 33,636 boxes open, 0 stuck) | nearly done |
| 1.486 | false | 7,993,453 | 300 (stopped) | **fails as it should**: 275 boxes stuck at width 1e-7 |

Box count grows about like margin^-1.1 to margin^-1.4: polynomial, no blow-up near
degenerate shapes in this lemma. Throughput is about 27k boxes/s per numpy thread.

### Cost extrapolation for the main nonobtuse inequality (estimate)

b(T) = area + sum of f >= acrit over the 6-dim nonobtuse domain, with a small ball around
the DL triangle cut out and handled by a local second-order argument.
- At the ball boundary the margin is about c rho^2. Take rho = 0.01 and margin about 1e-4.
  Scaling from 21M boxes at 1e-3 with exponent 1.1-1.4 gives 0.25-0.5 x 10^9 boxes.
- Evaluating f costs about 10x a plain box: 3 edges x 100 basis terms.
- That is about 30-100 CPU-hours in numpy, or about 1-5 CPU-hours in C or arb.
  Obtuse clusters (MITM as in HK) are probably the larger part: HK took about 60 h on a
  2016 quad-core.
- **Total: CPU-days, under 2 GB.** Not the bottleneck.

## VERDICT (2026-10-31): GO, conditional

The kill condition was "needs clusters of more than 6 triangles". It is not met.
- A single-triangle main inequality with an antisymmetric pair transfer works numerically for
  nonobtuse triangles.
- Obtuse triangles have an 11% surplus through their long-edge neighbour.
- Interval cost scales polynomially.

Hard parts left, all mathematical:
1. **Tangency at the DL.** f must satisfy the Kallus-Kusner local LP conditions exactly
   (antisymmetric transfers whose gradients reproduce their mu), and the null flex direction
   needs a second-order bound. That means the fitted f has to be projected onto those exact
   conditions, which is research.
2. **Obtuse bookkeeping.** Analogues of the HK cluster lemmas: shared neighbours, chains,
   at most k triangles per cluster.
3. **Rigour.** Port to arb or C with directed rounding, then Lean.

## 2026-11-01: exact local lemma at the double lattice (Lee D83 = A: continue)

### Infinitesimal local lemma: CERTIFIED (ball arithmetic, 300 bits), `src/local_lemma_arb.py`

Setting: one DL triangle s2 = (H0, Ht1, H2p2), H0 fixed, 6 variables. b = area + sum of f over
the 3 oriented edges. f is any C^2 pair transfer that is antisymmetric (f o R = -f) and
mirror-symmetric (f o m = f, with m(r,a,b) = (r, pi/7-b, pi/7-a)). Mirror symmetry is WLOG:
the problem is mirror-invariant and the certificate conditions are linear in f.

Symmetry forces, at the s2 edge configurations:
- translation edge (r*, 0, 0): f = 0, grad f = (0, ga, gb). This uses S = R o m:
  (r,a,b) -> (r,-a,-b), with f o S = -f.
- each reflected edge: f = 0, grad f = alpha_k (0, 1, -1).
- admissible Hessian terms are dt*(.) on reflected edges and dr*(da, db) on the translation
  edge.

Certified facts:
1. **Five active contacts**, each vertex strictly inside the edge it touches:
   - H0-H2p2 segment contact: two constraints;
   - Ht1-H2p2 segment contact: two constraints;
   - H0-Ht1 point contact: vertex k0 of H0 at the midpoint of edge 3 of Ht1.

   The next-nearest vertex/edge pair has gap 0.188. The contacts are exactly zero because
   H2p2 = Ref_p2 H0, Ht1 = Ref_p2 Ref_p3 H0, p2 lies on edge k1k2 and p3 on edge k2k3, and
   p4 is the midpoint of edge k3k4 (the KK construction).
2. **First order.** With alpha2 = 0 and lambda1 = lambda3 = 1/20, the system grad b = G^T lambda
   has the unique solution ga = -0.3178997, gb = -0.2618162, alpha1 = 0.1714933, with
   lambda = (0.05, 0.6375924, 0.05, 0.7079860, 0.2205221), all > 0.
   Plain area (f = 0) is NOT first-order stable: area decreases along the flex.
3. **Second order.** G has rank 5. Its null vector (the flex, pure translations)
   v = (0, 1, 0, -0.9749279, 0.2225209, 0) gives v^T Hess(L) v = 0.9749279 > 0.
   Along v, dr(e01) = dt(e12) = dt(e20) = 0, so **no admissible Hessian of f can change this**.

Conclusion: for every C^2 antisymmetric, mirror-symmetric f with those three gradient values,
the DL triangle is a strict local minimiser of b (standard second-order sufficient
condition: strict complementarity plus a positive-definite Lagrangian Hessian on the critical
cone). s6 follows by mirror symmetry. Unlike plain area, this needs **no cluster**: one
triangle suffices.

### Still open for the local part
- Quantitative radius rho (explicit neighbourhood) via interval Taylor bounds, to glue
  with the global B&B.
- Constrain the global LP basis f to these jets: f = 0 at the edge configs, the gradients
  above, or the whole feasible family of (ga, gb, alpha1, alpha2) as a linear constraint
  set. Then re-run the cutting plane with ETA -> 0 outside the rho-ball.

### Explicit neighbourhood (ball arithmetic, 128 bits): `src/local_radius_arb.py`, `src/local_tube_arb.py`

Jets-only f (linear near the s2 edges). Gauge: ga = -1/2, lambda4 = lambda5 = 0.13. Then
gb = -0.0797159, alpha1 = 0.0383694, alpha2 = -0.3090046,
lambda = (0.627986, 0.637592, 0.140522, 0.13, 0.13). The first-order solution set is 3-dim.
The smallest multiplier can reach at most 0.1353 (LP), whatever the gauge.

- Contact coordinates: Phi(z) = (g_1..g_5, v.z), with J0 = DPhi(0). J0 is ill-conditioned
  (singular values down to 0.105), so a z-box would be wasteful; work in zeta = J0 z.
- **Central chart certified:** the zeta-box of radius 0.001 (2^6 sub-boxes) gives Banach
  contraction q = 0.475, half-plane margin -0.16, w_k >= 0.029, beta'' >= 0.1885.
  This yields D0 = {0 <= g_k <= 5.25e-4, |s| <= 5.25e-4}, with z-radius below 0.0128.
  Radius 0.0015 fails (w4 goes negative; that is enclosure wrapping, and true failure is
  around 0.002).
- **Tube along the flex:** 123 certified charts cover s in [-0.0401, +0.0142]. In every one,
  w_k > 0 and beta'' >= 0.172. beta is convex on the union with beta'(0) = 0, so
  **b >= acrit on the whole tube, with equality only at the DL triangle**.
  The +s side stops at s = 0.0147: one multiplier reaches 0 there (with this gauge), meaning
  a contact wants to open. beta - acrit there is 2.4e-5.
- Margins the global B&B must certify on the tube boundary: about 1e-6 to 4e-5 on the
  contact-opening faces (w_min x eta + beta - acrit); 1.9e-4 at s = -0.040; 2.4e-5 at
  s = +0.014. Small, so near the tube the B&B needs second-order Taylor bounds; plain
  intervals will not do.
- Caveat: the tube used jets-only f. It must be re-run with the final global f, whose
  Hessian and higher terms differ. The code takes any f evaluable on arb series.

### Collision check rerun, 2026-11-01

Zenodo (title:heptagon*, creators Tomar, heptagon AND packing, double lattice AND heptagon,
pessimal AND heptagon): the only new item is 10.5281/zenodo.23686812, "Improved packings of
equal circles in regular polygons", 2026-10-27, which is unrelated.
arXiv: the API was rate-limited (HTTP 429), so I used the arxiv.org search pages instead:
"heptagon packing" returns only 2212.13849 and 1305.0289; "double lattice" has nothing new
that is relevant.
GitHub repos, code, issues and PRs: only the V3R1T4S bank entry (no new commits since
2026-10-30) and a Wikipedia mirror. MColbrook/AIM and openai/math: nothing.
**Clear.**

### LP with local-lemma jets (2026-11-01)

- `lp_transfer2.py` (mirror-symmetric basis, jets FIXED at the hand-picked gauge
  ga = -1/2, ...): LP value drops to **0.99034 acrit**. That gauge is bad, not the method.
- `lp_transfer_sym.py` (symmetric basis, no jets, zero slack): converges to acrit
  (local-search min 1.00000 acrit). So mirror symmetry costs nothing.
- `lp_transfer3.py` (jets allowed to range over the admissible family, i.e.
  Fm D c - G^T lambda = -grad area with lambda >= LAM_MIN, zero slack): **feasible at acrit**
  - LAM_MIN = 0.02: jets p = (-0.3504, -0.2293, 0.1227, -0.1404),
    lambda = (0.541, 0.638, 0.02, 0.217, 0.251).
    Independent check (seed 777): fresh-sample min 1.0058 acrit. Local search min
    acrit - 4.35e-7 at a triangle with edges 1.835 / 1.876 / 1.863, which is NOT the DL
    triangle. So the f on the LP's zero-slack optimal face is razor-thin in places.
  - LAM_MIN = 0.1: also feasible at acrit (local-search min 0.99999 acrit).
- Exactification (`src/f_arb.py`): f + mu * basis_9 with mu in a ball, so that
  grad b . v = 0 exactly. Before the correction it was -2.4e-11.
- Central chart with the actual global f (not jets only): certified at zeta radius
  1e-4 (LAM_MIN 0.02 f) and 2e-4 (LAM_MIN 0.1 f: eta = 1.8e-4, w >= 0.079, beta'' >= 0.16).
  The radius is limited by f's large higher derivatives (frequencies up to 14 in a, b).
- Robust LP (`lp_transfer4.py`): maximise s with b >= acrit + s*w(T), where
  w = min(1, d_DL^2 / 0.05^2). The sample margin settles at s = **1.58e-3**, but without
  regularisation the cutting plane does not converge in 20 rounds (local search finds
  0.983 acrit). Rerunning with an L1 penalty, 200 rounds, 95 min (`runs/lp4b.log`).
- **Robust LP converged** (`runs/lp4b.log`, L1 penalty 1e-3, 25 rounds, 93 min, leased).
  The LP value stays at acrit, and local search has found nothing below acrit since round 7.
  Sample margin s = 5.8e-4: b >= acrit + 5.8e-4 * min(1, d_DL^2 / 0.05^2) on all cuts.
  |c|_1 = 0.032. Jets p = (-0.3854, -0.1943, 0.1039, -0.1420),
  lambda = (0.389, 0.638, 0.05, 0.369, 0.221). Coefficients are in `runs/lp4b_c.npy`.
- **Independent validation** (seed 777): fresh samples give min b = 1.01059 acrit. Local search
  gives min b = acrit + 5.5e-7, at the DL triangle (translation edge 1.900969). **No violation
  found** (float, unverified).
- Central chart with this f (exactified, mu = -1.0e-12 on basis 9): certified at zeta radius
  2e-4 (eta = 1.8e-4, w >= 0.014, beta'' >= 0.114). It fails at 5e-4: f's curvature is large
  (frequencies up to 14).
- Next: (a) tube along the flex with this f; (b) rigorous B&B of the main nonobtuse inequality
  using second-order Taylor bounds near the DL, with a cost pilot first; (c) obtuse clusters.

## 2026-11-01 (cont.): steps (a), (b), (c)

### (a) Tube along the flex with the real f (lp4b, exactified): DONE (ball arithmetic)
422 of 423 charts are certified (zeta radius 2e-4, 2^6 sub-boxes each, 177 s in total). They
cover s in [-0.0501, +0.0179], with w > 0 (min 4.7e-6 at the far +s end) and beta'' >= 0.098.
The +s end stops where a multiplier reaches 0.
beta - acrit at the ends: 3.1e-4 (s = -0.050) and 4.0e-5 (s = +0.018).

### (b) Interval pilot for the main nonobtuse inequality (`src/bnb_main.py`, cost pilot, not rigorous)
- Bounds: interval forward-mode AD, taking the larger of the natural extension and the
  mean-value form. Transfer f evaluated in factorised form.
  Cross-check against the float b: agreement to 1e-14.
- **Measured rate: 2,400-2,900 boxes/s per thread.**
- Slice theta0 in [0.40, 0.42] (positions d <= 2.4, y in [1.2, 2.1]): **done**, 1.98M boxes,
  818 s, nothing stuck. So b >= acrit there (pilot grade).
- Slice [0, 0.02]: not finished after 5.65M boxes / 40 min (12.8k boxes still open).
  Slice [0.78, 0.80], which contains the canonical DL triangle (d = 1.915, x = 0.890,
  theta = (0.794, 0.345, 0.345)), with target acrit - 1e-5: not finished after 5.7M.
- The open boxes are dominated by the **feasibility boundary**: boxes straddling a contact
  surface, where the infeasible side has b < acrit, at widths around 0.007-0.01. The test
  "a vertex is strictly inside for the whole box" only prunes once boxes are small.
  The margin of b is not what limits the count.
- Running now: 4 slices to completion (2.5 h cap each) to measure the box-count distribution.

### (c) Obtuse clusters with f (`src/obtuse_pair_f.py`)
Obtuse triangle plus its neighbour across the longest edge, with b including f: min average
**1.0797 acrit** (33 starts, all converging to the same value). The binding case is a right
angle, i.e. the edge of the nonobtuse domain. Still to test: a neighbour shared by 2-3
obtuse triangles, and chains of them.

### (b) cont.: interval SAT pruning and the full sweep (pilot grade)
- The new overlap test (interval separating-axis: all 14 axis gaps < 0 on the whole box) is
  sound by SAT. Slice [0, 0.02] went from > 6M boxes unfinished to **114k boxes, 46 s**.
- **Full sweep** over all canonical nonobtuse triangles:
  - domain: d in [2u, 4], x in [0, 2], y in [1.2, 3.5], all orientations; y >= 1.2 and
    the other bounds are derived in the notes above;
  - 45 theta0 slices, 5 threads, lease heptagon-4432;
  - 133.9M boxes, 47,271 thread-s, **2,832 boxes/s/thread**, nothing stuck.
  - 43 of 45 slices done. That includes the DL slice [0.78, 0.80], run with target
    acrit - 1e-6.
  - Slices 12-13 (theta0 in [0.24, 0.28]) hit the 1.5 h cap with 25k / 12.6k boxes open.
    Almost all open boxes have infeasible centres (contact-surface straddlers at width
    ~1e-3). Feasible open centres have b - acrit >= 5.7e-3, so it is just slow pruning, not a
    near-violation. Re-running both slices with a 3.9 h cap.
- **Measured cost: about 1.5-2e8 boxes, about 15-20 thread-hours in numpy** (the estimate
  was ~1e8). A C/arb port would be 30-100x faster, which matters for the rigorous version.

### (c) cont.: clusters of obtuse triangles (float, few starts)
Central triangle receiving K obtuse triangles across its edges, minimum average b / acrit:
K = 1: 1.0797; K = 2: 1.2175; K = 3: 1.2926. Chains (the neighbour itself obtuse) are not
yet tested.
- Slice 13 finished in the long run: 27.3M boxes, 9,646 s.
- Slice 12 was not done after 42.3M boxes / 14,000 s (16.2k open). The open boxes sit around
  a double/triple-contact region where b is within ~1e-3 of acrit: median centre
  b - acrit = 2.4e-4, and feasible centres are at >= 1.05e-3. Pruning that region is slow.
- Added Lagrangian **dual bounding**. When exactly one SAT axis can still separate a pair on
  the box, feasibility forces the vertex-on-edge inequality h >= 0. Then
  LB = min over the box of (b - sum mu h), mean-value form, with mu >= 0 from the projected
  centre gradient. On 4000 open boxes it discards 15% of those that survive the plain bound.
- Fixed a note in the SAT code: the "reversed side" gap was identical to the forward gap.
  Testing outward normals of both polygons is complete for SAT (a separating line can be
  taken through an edge of one polygon), so the pruning was sound throughout.
- Obtuse chain T -> T0 (obtuse) -> T2: min average b = 1.3897 acrit (5 feasible starts).
- Resuming slice 12 from its dump with dual bounds: lease, 2 h cap.
- Zenodo has blocked this network for "unusual traffic" (coordinator, 2026-11-01): no
  Zenodo API calls for the rest of 2026-11-01.
- Slice 12 resume with dual bounds (2 h, lease heptagon-4434): 12.6M more boxes, open count
  **grew** to 30.4k. Float local minima started from the open boxes (`src/localmin_open.py`)
  give b - acrit = 6.5e-4, 7.0e-4, 8.1e-4, ... at theta0 ~ 0.17-0.25, double/triple contact.
  (A 4.1e-4 one is the mirror DL with x > d/2, outside the canonical domain.)
  So slice 12 holds a **broad near-tight zone with margins of 4e-4 to 1e-3**: the flat part of
  the LP optimum, s = 5.8e-4. First-order interval bounds then need boxes of about 1e-4 there.
  This is slow, but no violation has been found.
- **Pilot verdict for (b):** 44 of 45 slices certified at pilot grade (161M boxes,
  2,800 boxes/s/thread, about 17 thread-hours). The remaining slice is limited by thin margins.
  For the rigorous run:
  (i) re-fit f for a larger margin: bigger basis, margin weight in the LP;
  (ii) second-order (Taylor-model) bounds plus the dual bounding;
  (iii) port to C/arb (30-100x).
  With (i)-(iii), a full rigorous sweep should take CPU-hours, not days.

## 2026-11-01 evening: fixes 1-3

### Fix 1: refit for margin
- DR3 DF3 (196 coefficients): converged, but margin s = 5.9e-4. No better than DR3 DF2.
- DR4 DF3 (245 coefficients): s = **2.46e-3** (4x better), but not converged in 2 h (local
  search 0.9976 acrit at round 15).
- Running DR4 DF3 and DR5 DF2 to convergence (6 h cap, lease).

### Fixes 2 and 3: rigorous bounds (`src/rig.py`, `src/bnb2.py`)
- No C compiler on Windows. Ubuntu WSL has gcc but no flint/arb/mpfr dev headers, and
  installing them needs sudo. So "C/arb port" became **rigorous vectorised numpy**:
  - outward rounding after every IEEE op: dn(a) = a - (|a| 2^-51 + 2^-1000), which is
    >= 2 ulp, so round-to-nearest cannot land above the true value;
  - **own rigorous sin/cos/atan**: argument reduction with a 1-ulp enclosure of pi/2,
    Taylor polynomials with explicit remainders; libm only picks branches;
  - self-test against mpmath (40 digits): 20,001 sin/cos points, 10,004 atan points,
    3,000 random cos intervals, all enclosed;
  - forward AD of order 2 (value, gradient, Hessian) over intervals; Hessians checked
    against finite differences;
  - bounding ladder: geometric + SAT prune -> natural/mean-value (order 1) -> second-order
    Taylor (order 2) -> second-order Lagrangian dual.
- Benchmark on slice [0.40, 0.42], 120 s: 138k boxes (1,130 boxes/s). The old pilot handled
  2,800 boxes/s but needs about 10x more boxes for the same open front, so the net is about
  4x faster, and now rigorous.
  Replacing nextafter by the multiplicative widening gave 4.4x by itself.
- Hard-slice pilot (theta0 0.24-0.26) with bnb2 and the lp4b f: running (lease, 1.5 h cap).
- **Hard-slice pilot with bnb2** (theta0 0.24-0.26, lp4b f): 8.0M boxes in 90 min,
  1,480 boxes/s, nothing stuck. Not finished, but the open front stays at 10-15k and shrinks,
  where the old method's front grew to 30k after 55M boxes.
  Breakdown: 3.09M infeasible, 0.51M first-order, 64k second-order, 318k dual.
- **Long refits converge to the same margin:** DR4 DF3 (43 rounds) s = 5.83e-4; DR5 DF2
  (68 rounds) s = 5.81e-4. The 2.46e-3 seen earlier was before convergence.
  **s = 5.8e-4 is a cap set by the weight function, not by the basis.** With
  w = min(1, d^2 / 0.05^2), near-DL triangles at d ~ 0.05 (along the flex, where b - acrit is
  fixed by the geometry, ~0.1 s^2) cap s.
  The slice-12 near-tight triangles have DL parity (one translation pair, two reflected pairs)
  and are rotated by ~0.017: they lie in the valley of DL-like triangles.
- Testing wider weight radii (WDELTA = 0.15, 0.30) to see whether the valley margin grows.

## 2026-11-02

- **Wider margin weights fix the margin.** Weight radius 0.15: s = 5.2e-3. Radius 0.30:
  **s = 1.25e-2**. Both converged (local search reaches acrit).
  f_w030 = `runs/f_w030.npy` (DR3 DF2, |c|_1 = 0.061). Independent validation: fresh samples
  1.0096 acrit; local-search min acrit + 2.4e-4.
- **Hard slice (theta0 0.24-0.26) with f_w030 and the rigorous bnb2: DONE**, 1,276,873 boxes,
  1,208 s (20 min). Breakdown: 254k infeasible, 45k first-order, 107k second-order,
  232k dual. With the old f it had not finished after 55M boxes.
- Chart radius with f_w030 is smaller (f bends more): central chart certified at zeta radius
  6e-5 (eta = 5.8e-5). The 6e-5 tube has **459 charts, all certified**, covering
  s in [-0.012, 0.012] (w >= 0.0134, beta'' >= 0.072).
  The exactification mu = -1.16848e-12 on basis 9 enters the B&B as an interval coefficient.
- `src/chart_excl.py`: a rigorous map from canonical boxes to the chart frame. The canonical
  DL labelling is (H2p2, H0, Ht1), with rotation ph = -t1 - 2(2pi/7). A box is discarded if
  its image lies in a certified chart region. Tested: a 1e-6 box at the DL is inside, a
  1e-5 box is not.
- Running: the DL slice [0.78, 0.80] with chart exclusion at target exactly acrit, and slice
  [0.40, 0.42], both with f_w030 (lease, 1.5 h).
- **Rigorous sweep with f_w030** (`runs/sweepR/`, lease heptagon-4441, 5 threads):
  **43 of the 44 non-DL slices DONE**, 56.8M boxes, 65,803 thread-s (~18.3 thread-hours),
  864 boxes/s/thread, nothing stuck. Slice 38 [0.76, 0.78], next to the DL slice, was still
  running at 5,200 s (cap 10,800 s).
- Ordinary slice [0.40, 0.42]: done in 1.72M boxes, 33 min.
- DL slice [0.78, 0.80] with chart exclusion: 6.6M + 11.4M boxes, 116k boxes excluded by
  charts, not finished (17k open). The open boxes near the DL sit at distance ~1e-3 with
  widths 1e-5 to 1.5e-5.
  Cause: dual v1 needs a pair to have a *unique* separating axis, but at a segment contact
  two near-antiparallel axes both have gap ~0, so v1 never fires there.
- **Dual v2** (`lb_dual2`): the vertex-on-edge inequality h_{k,m} >= 0 is valid whenever
  vertex m of the other heptagon is strictly inside the owner's other six half-planes on the
  whole box. Then h < 0 puts the vertex strictly inside the owner, which is an overlap.
  So both inequalities of a segment contact are usable. On 1,121 near-DL boxes that survive
  the second-order bound: v1 decides 283, v2 590, either 709.
  Resuming the DL slice with v2 (lease, 4 h cap).
- DL slice resume with dual v2 (4 h): 9.1M more boxes, 461k excluded by charts, still 21k
  open. The open near-DL boxes lie along the flex at |s| up to 0.021, i.e. **outside the
  ±0.012 tube**, where beta - acrit is only ~1e-5.
- **Extended tube for f_w030:** 1,530 certified charts covering s in [-0.0301, +0.0500].
  The -s end is where a multiplier reaches 0; beta - acrit there is 1.2e-4.
  Saved as `runs/charts_w030b.json`.
  Chart pre-filter widened to 0.12, with a float pre-selection of candidate charts by the
  flex coordinate (the selection is heuristic only; membership is still certified by
  intervals).
- Slice 38 resume (dual v2, old ±0.012 charts) is still running. DL slice relaunched with
  the extended tube.

## 2026-11-02: MILESTONE: rigorous sweep of the main nonobtuse inequality COMPLETE

**Statement certified** (computer-assisted; code not yet independently audited):
for every nonobtuse triangle of three pairwise non-overlapping regular heptagons
(circumradius 1) with circumradius <= 2,
    b(T) = area(T) + f(e01) + f(e12) + f(e20) >= acrit = 1.5326754446782214...,
where f = f_w030 (`runs/f_w030.npy`, DR3 DF2, mirror-symmetric, antisymmetric; coefficient 9
shifted by the exactification interval mu in [-1.1684855e-12, -1.1684853e-12]).
Equality holds only at the DL triangle. The local part (inside the certified charts) uses the
ball-arithmetic chart lemma: 1,530 + 459 charts, `runs/charts_w030b.json`,
`runs/charts_w030.json`.

- **All 45 theta0 slices done.**
  - Slices 0-44 except 38/39: single runs (`runs/sweepR/`).
  - Slice 38: run plus resume.
  - Slice 39 (DL): chain open -> open2 -> open3 -> open4 -> 4 shards, all finished.
    Two runs I killed are not part of the chain; each relaunch resumed from the last
    *dumped* open set.
- Totals: **124.6M boxes, 39.1 thread-hours, 928k boxes excluded by charts, nothing stuck**.
- Engine: `src/rig.py`, `src/bnb2.py`, `src/chart_excl.py`, `src/local_tube_arb.py`,
  `src/local_lemma_arb.py`, `src/f_arb.py`.

**What this does NOT yet give:**
- the obtuse-cluster inequalities (rigorous);
- the structural lemmas: saturation gives R <= 2; Delaunay triangulation; the cluster
  bookkeeping and finiteness; the averaging argument (HK Lemma 11 analogue);
- independent replay or audit of the interval code; Lean.

**Domain-reduction facts used** (re-check in the audit):
- b is invariant under cyclic relabelling and under mirror images, because f o m = f.
  Hence c0c1 can be taken as a longest edge and x <= d/2.
- Nonobtuse with c0c1 longest gives y >= 1.2 and y <= d sqrt(3)/2 <= 3.47.
- R <= 2 with nonobtuse gives d <= 4.
- The SAT overlap test is complete using outward normals only.
- The vertex-in-slab inequality is valid whenever the touching vertex lies strictly inside
  the other six half-planes.

## 2026-11-02: audit packaging (coordinator: keep the interval layer small and isolated)

### Interval-layer contract (`src/rig.py`, the only homemade arithmetic)
- `R(lo, hi)`: arrays of closed intervals. Every operation (+, -, *, sq, inv, sqrt) returns
  an interval containing the exact real result for all real inputs in the argument intervals.
  Implementation: IEEE round-to-nearest, then widening by dn/up, where
  dn(a) = a - (|a| 2^-51 + 2^-1000) and up is symmetric.
  Claim: for a computed value v with |v - exact| <= ulp(v)/2, dn(v) <= exact <= up(v).
- `rcos(a), rsin(a)`: interval enclosures for any interval a. Argument reduction uses
  PIO2 = [1.5707963267948966, next double up], which contains pi/2. The degree-24 Taylor
  polynomials carry explicit remainders |r|^25/25! and |r|^24/24! for |r| <= 0.8.
  Interior extrema are detected conservatively (any doubt includes +-1).
- `ratan(a)` for a.lo >= 0: atan(t) = pi/2 - atan(1/t) when t > 1, two half-angle reductions,
  a 14-term alternating series and its remainder. Monotone, so endpoints suffice.
- **Swap-in test for the replay:** replace `R` and these three functions by mpmath.iv or
  python-flint arb balls with the same interface (lo/hi arrays), re-run `bnb2.run` on sampled
  boxes, and check that every certified box is still certified.
- **Rounding helpers used outside rig.py** (also to audit): `bnb2.py` lines ~181, 188
  (pruning thresholds), 205-216 (mean-value and Taylor bound assembly), 353 (per-box targets).
- **Non-rigorous floats used ONLY for heuristic choices** (soundness never depends on them):
  `bnb2._min_vertex` and `_cand_float` (which vertex/edge inequality to try), the chart
  pre-selection in `chart_excl.inside`, and the box-splitting rule.
- Ball arithmetic (python-flint arb) is used separately for the chart lemmas
  (`local_lemma_arb.py`, `local_radius_arb.py`, `local_tube_arb.py`). Chart data is exported
  to the B&B as float lo/hi bounds (`runs/charts_*.json`).

### Domain-reduction facts the sweeps ASSUME (each needs a written proof in the paper)
1. Saturation means every Delaunay triangle has circumradius R <= 2. If a circumcircle had
   radius > 2, its centre is > 2 from every heptagon centre, so a heptagon fits there.
2. Non-overlap gives centre distance >= 2u = 2 cos(pi/7) (inscribed discs).
3. b is invariant under cyclic relabelling (it is a sum over ccw edges). b is invariant under
   reflection because f o m = f and the mirror reverses orientation. Hence c0c1 is a longest
   edge, c1 = (d, 0) and x <= d/2.
4. Bounds: d <= 4 (R <= 2); y <= d sqrt(3)/2 <= 3.47 (c0c1 longest); y >= 1.2 for
   nonobtuse triangles; y >= (2u)^2/4 = 0.81 for all triangles with R <= 2 (area >= abc/8).
5. Orientations can be reduced mod 2 pi / 7 (heptagon symmetry); f is 2 pi/7-periodic in a, b.
6. SAT is complete with outward edge normals of both polygons: a separating line can be
   taken through an edge of one polygon.
7. Vertex-in-slab: if vertex q of H_j lies strictly inside the six half-planes l != k of H_i,
   then h_{k,q} = n_k.(q - c_i) - u < 0 implies q is in int H_i (overlap). So non-overlap
   implies h >= 0.
8. Chart lemma: inside each certified chart, overlap holds whenever some g_k < 0, and
   b >= acrit on Phi^{-1}(D_c) (ball arithmetic).
9. Obtuse rule:
   - an obtuse triangle's longest edge is > 2u sqrt(2) = RTH;
   - two triangles sharing an edge cannot both be obtuse at the angles opposite it
     (Delaunay: the opposite angles sum to <= pi);
   - a receiver's angle alpha opposite the shared edge satisfies cos gamma >= -cos alpha;
   - bins are half-open, [rb_i, rb_{i+1}) x [cb_j, cb_{j+1}).

### Structural part: proof SKELETON (Lee writes the prose)

Constants: A = (7/2) sin(2pi/7) (heptagon area, circumradius 1); acrit = A/(2 delta_DL).

- **L1 (saturation).** We may assume the packing is saturated; adding heptagons does not lower
  the density.
- **L2 (Delaunay facts).** Take the Delaunay triangulation of the centres; in a degenerate
  (cocircular) cell, any triangulation of that cell.
  - Every triangle has circumradius <= 2 (fact 1).
  - Every edge is >= 2u (fact 2) and <= 4.
  - For any interior edge, the two opposite angles sum to <= pi (empty circles).
- **L3 (no cluster finiteness needed).**
  - Transfers are one step: f across every edge (antisymmetric, so it cancels between the two
    sides), and D from each obtuse triangle to its long-edge neighbour.
  - Chains of obtuse triangles are absorbed inside D by the recursion over increasing r: an
    obtuse triangle's incoming edges are shorter than its long edge.
  - So the per-triangle inequality b*(T) >= acrit is all we need.
  - [Check in the audit: the table D2 was built by that recursion from samples; the rigorous
    sweep verifies the per-triangle inequality directly, so the way D2 was built does not
    matter for soundness.]
- **L4 (per-triangle inequality).** b*(T) = area(T) + sum_edges f + D_out(T) - sum D_in(T) >= acrit
  for every Delaunay triangle T. This holds by:
  - main sweep (nonobtuse, longest edge < RTH: D terms are 0);
  - obtuse-mode sweep (longest edge >= RTH, all triangles);
  - chart lemma at the DL.
- **L5 (averaging; HK Lemma 11 analogue).** Take a window of side L.
  - #triangles = 2 #centres + O(L) (Euler; the boundary has O(L) vertices because centres
    are 2u-separated and edges are <= 4).
  - The sum of b* over triangles in the window equals the sum of their areas + O(L), because
    every interior transfer cancels and |f|, D are bounded.
  - Hence N A / L^2 <= A/(2 acrit) + O(1/L) -> delta_DL.
- **L6 (optional: uniqueness).** Equality forces b*(T) = acrit for almost all T, i.e. T is
  (up to congruence) the DL triangle by the chart lemma's strictness. Then the packing is
  the DL. Needs a rigidity argument; not required for the density bound.

## 2026-11-03: obtuse sweep with D3

- First D3 run (detached, lease 4450): **invalid environment**. A WMI-detached job does not
  reliably see the user site-packages. 18 workers failed with "no numpy". The rest
  misbehaved: pruning counters froze after ~200 s and boxes split until minw, giving
  millions of "stuck" boxes. Discarded (`runs/sweepO_broken/`).
- Re-run of slices 44 and 18 in a clean environment, same code (git-identical except stuck
  instrumentation) and same D3: **done**. 621k boxes / 740 s and 703k boxes / 871 s,
  0 stuck. So the stuck boxes were an environment artefact, not a gap in D3.
- Launcher fixed (`runs/sweepO4.sh`): PYTHONPATH pinned to the user site, numpy version and
  path logged by each worker, no --detach. Relaunched the 43 remaining slices on 7 threads,
  lease heptagon-4463, ETA 120 min.

## 2026-11-03/04: table gap found; D4

- **Stuck boxes in obtuse slices 9-13 = a real gap in table D3.** Of 321k captured stuck boxes,
  the feasible ones are obtuse isosceles triangles with two edge contacts: short edges about
  1.802 (= 2u), long edge 2.796-2.799, cos gamma about -0.20, R about 1.43.
  b - acrit goes down to -0.0067, but their cell (r in [2.749, 2.799), c in [-0.25, -0.20))
  had D = 0. The sampled table missed this family. This is not a counterexample to the
  theorem: the sweep catches it as non-certifiable, which is exactly what it should do.
- Per-cell maximum deficits by multistart local maximisation (`src/cell_deficit.py`).
  The first run lacked y > 0 and gave bogus values; rerun with y > 0 and 0 <= x <= d/2.
  Row maxima peak at 0.183 (r ~ 3.30). The gap cell reads +0.0077.
- **D4** (`runs/D4_nodes.npy`, `src/build_D4.py`): **continuous** bilinear table on nodes.
  - Cell target = max(deficit, 0) + 0.03, plus 0.02 within one cell of any near-tight cell
    (deficit > -0.01). Node = max over adjacent cells.
  - Max 0.233, 92 nonzero nodes.
  - **sha256 65eb79fab861b9007488e65775546020ef5a0ced46761ae624ec94a682ae0e13**.
    Bins: `runs/D4_bins.json`, sha256 afa2c4ac0d7330f73feb17f0cb16cb6426a7b6889ceef306db626ace6fbd9b85.
  - B&B (`ObtTable4`): the obtuse lower bound is exact bilinear corner evaluation inside a
    cell (bilinear extremes over a rectangle are at its corners; minus 1e-12), or the min of
    nodes across cells. The receiver upper bound is the max of nodes over compatible cells.
- **Which certified runs are affected by D3 -> D4:**
  - The obtuse-mode sweep (all 45 slices) is redone with D4. It is the rigorous receiver
    check: every triangle with longest edge >= 2u sqrt2, charged by the table's upper bound.
  - Unaffected: the non-obtuse main sweep (it never uses D; it covers K = 0, where no D
    terms exist); the chart lemma and tube (f only); L1-L3 and L5 (they use only that D is
    bounded). L4 is exactly the two sweeps.
- D4 pilot k9/k11 (45 min each): **0 stuck** (D3: ~600k), open front ~25k and falling.
- Full obtuse sweep with D4 on **Modal** (one container per slice; python 3.13.3, numpy 2.5.4;
  file hashes printed in every log and matching local: bnb2.py 2ac34278..., rig.py
  1ad0261b..., f_w030 ced5a7aa..., D4 65eb79fa...). Logs go to `runs/sweepM/`.
  Local replay pair: k9 and k11, resumed from the pilot dumps.

### D4 bounds: soundness argument (for the hostile audit, check explicitly)
`ObtTable4._rect(rlo, rhi, clo, chi, want_max)` bounds the bilinear table over a rectangle
Q = [rlo, rhi] x [clo, chi].

1. **Covering.** Cells C_ij = [rb_i, rb_{i+1}] x [cb_j, cb_{j+1}] are CLOSED and cover the
   table domain [RTH, 4] x [-1, 0]. A cell is visited iff it meets Q in the closed sense
   (rhi >= rb_i and rlo <= rb_{i+1}, the same in c). So every point of Q (clipped to the
   domain) lies in some visited cell. There is no gap at a cell boundary: a point on a
   shared edge belongs to both cells, and the bilinear interpolant is continuous across
   edges (both cells give the edge-linear interpolant of the two shared nodes).
2. **Clipping is harmless.**
   - Obtuse lower bound: an obtuse triangle has r >= RTH (fact 9) and c in [-1, 0).
   - Receiver upper bound: compatible c lies in [-ca.hi, 0] and r in [RTH, 4]; edges are
     <= 4 since R <= 2.
3. **Exact extremes.** On Q ∩ C_ij the interpolant is bilinear in local coordinates (t, s),
   linear in each separately. So its min and max over the sub-rectangle are attained at the
   sub-rectangle's 4 corners. The code evaluates exactly those 4 corners.
4. **Float error.**
   - Corners are computed as t = (clip(r) - rb_i)/(rb_{i+1} - rb_i), with error <= a few ulps
     (<= 1e-15).
   - Node values are in [0, 0.24], and the bilinear form has 4 products of numbers in [0, 1].
     The total value error is <= 0.24 x (|dt| + |ds| + 8 ulp) < 1e-14.
   - The code subtracts 1e-12 from lower bounds and adds 1e-12 to upper bounds, which covers
     this with margin.
   - [Auditor: re-derive this bound, or replace the evaluation with rig.R arithmetic.]
5. **Proof-side definition.** The transfer from an obtuse T is D4(r_T, cos gamma_T), the
   bilinear interpolant (continuous, single-valued). The receiver is charged exactly that
   amount, and it is <= sup over compatible (r, c) of D4, which is what `dhat_ub` bounds.
6. **Straddling boxes** (cos gamma ranges over 0): extra = max(obtuse-part bound with
   cg in [cg.lo, 0], non-obtuse-part bound with cg in [0, cg.hi]).

### Modal graceful stop
`python -m modal run src/modal_obt.py::request_stop` sets the Dict flag heptagon-ctl[stop].
Containers poll it every 30 s and SIGTERM bnb2, which dumps its open boxes (STOP_REQ) and
returns TIMEOUT. Tested: slice 20 stopped after 122 s and returned its dump.

## 2026-11-04: table D5, padding-rule history, and the stated reduction (for the auditor)

### Table history (the padding rule was changed after seeing stuck runs; legitimate, since the
sweep certifies whatever table it is given, but recorded here)
- D2: sampled cell maxima + 0.02. Stuck at the r ~ 3.27 cell, deficit not covered.
- D3: D2 widened to neighbouring cells + 0.04. Stuck in slices 9-13: the two-contact isosceles
  family at r in [2.749, 2.799) had deficit 0.0077 in a D = 0 cell.
- D4: continuous bilinear, padding +0.03 (+0.02 near tight cells) on every cell with local-search
  deficit > -0.05. Stuck in slice 13: near-right ACUTE receivers (cos gamma ~ +0.025, long edge
  2.549, b - acrit ~ 0.030) were charged 0.03 by padding put on cells whose deficit is ~ -0.02.
- **D5** = D4 rule with padding only where deficit > -0.01 (ZTOL 0.01). 17 cells zeroed;
  D5 <= D4 at every node. `runs/D5_nodes.npy`
  **sha256 87d669a7419b037b801f753c173da484833a49ab2e00726375a8448d6ff1c020**.
  Bins identical to D4 (`runs/D5_bins.json`).

### Stated reduction: one table, two sweeps
- Transfers in the proof: f across every edge, and D5(r_T, cos gamma_T) from each obtuse T to
  its long-edge neighbour.
- **Main (non-obtuse) sweep** (`runs/sweepR`, `runs/sweep`): certifies b >= acrit, with NO table
  (HEPT_OBT/HEPT_OBT4 unset). It is used only for triangles with longest edge < 2u sqrt2.
  Such a triangle cannot be obtuse (obtuse means longest^2 > 2 (2u)^2) and cannot receive
  (the shared edge would have to be >= 2u sqrt2). So D_out = D_in = 0 and b* = b.
- **Obtuse-mode sweep** (D5): domain = all canonical triangles with longest edge d >= 2.548
  (< 2u sqrt2), obtuse or not. It certifies b - acrit - extra >= 0, where extra is an upper
  bound of (incoming D) - (outgoing D).
  **All three edges are charged** (`src/bnb2.py` `obtuse_extra`, lines ~339-345):
  - `inc_short` = dhat_ub over edges 12 and 20, for every triangle;
  - edge 01 (the longest) is charged as incoming (`extra_non`, + dhat_ub(r01, ...)) unless the
    triangle is certainly obtuse;
  - for an obtuse triangle edge 01 is outgoing (`extra_obt`, - d_lb(r01, cos gamma));
  - boxes straddling the right angle take the max of both.
  So a triangle with two or three edges >= 2u sqrt2 is charged on each of them.
  `dhat_ub` returns 0 only when the edge is certainly < RTH or its opposite angle is certainly
  >= 90 deg (no obtuse neighbour can then be Delaunay-adjacent).
- **Boundary-leak check (exact):** the threshold node rb[0] = 2.5483247845270696 exceeds
  2u sqrt2 = 2.548324784527069289... by 2.7e-16 (mpmath, 40 digits). The first two node rows
  of D5 are exactly 0.0, so D5 vanishes identically on r in [rb[0], rb[1]]. Defining D5 = 0 for
  r < rb[0] is continuous, and no donor with r in (2u sqrt2, rb[1]) carries a nonzero amount.
  Receivers are charged whenever r.hi >= rb[0], which covers all r where D5 can be nonzero.

### Modal driver (`src/modal_obt.py`)
- Checkpoints every 10 min to a Volume; a preempted container resumes from its checkpoint
  (logged).
- On STUCK the driver sets the graceful-stop flag and keeps collecting every other slice's dump.
- The table is selected by HEPT_TABLE / --table (D5).

## 2026-11-04: MILESTONE: obtuse-mode sweep with D5 COMPLETE (Modal)

- **45/45 slices DONE, 0 stuck, 0 timeout, 0 fail.** 35,586,755 boxes, 20.96 container-hours
  (Modal, python 3.13.3, numpy 2.5.4 in all 45 logs).
  Logs in `runs/sweepM_d5/`, driver log `runs/sweepM_d5_driver.log`.
- Every log has identical hashes: bnb2.py 7ed516ab8d6535cc..., D5_nodes 87d669a7419b037b...,
  D5_bins afa2c4ac0d7330f7....
- **Replay pair (local Windows numpy 2.5.4 vs Modal Linux): identical box counts and verdicts.**
  k13: 1,083,199 vs 1,083,199. k30: 459,705 vs 459,705. Both DONE, 0 stuck.
- Slices 29 and 40 were preempted and resumed from their 10-min Volume checkpoints (logged).
  The other preemptions restarted before their first checkpoint, i.e. from scratch.
- **Status of the computer-assisted chain:**
  - main nonobtuse sweep (no table; needed for longest edge < 2u sqrt2) 45/45;
  - obtuse-mode sweep (D5; longest edge >= 2u sqrt2) 45/45;
  - chart lemma and tube for the DL neighbourhood;
  - the rest is the structural skeleton L1-L6 (written proofs needed).
- Next, per the coordinator: fresh-agent replay and hostile audit (assigned separately), then
  Lean.

## 2026-11-04: fixes for the hostile audit (AUDIT-2026-11-04.md) and PROOFS-DRAFT corrections

### Corrections to earlier NOTES (per PROOFS-DRAFT.md)
- **Drop the "chains absorbed by the recursion over r" remark (L3).** It is irrelevant to
  soundness: the obtuse-mode sweep charges all three edges directly against the fixed table D5.
  How D2-D5 were built does not matter.
- **D5 cells are CLOSED** (bilinear, continuous). The old "bins are half-open" (fact 9) applied
  to the piecewise-constant D2/D3 only.
- **Density definition (shared with PROOFS-DRAFT section 2):** delta+(P) = limsup_{L -> inf}
  sup_p (2L)^{-2} sum_{H in P} area(H ∩ (p + [-L, L]^2)). Other standard variants are <= this one,
  or follow by the same proof, since L5 is uniform in p.
- "beta convex on the union, beta'(0) = 0" needs (E), (G), (Ov) (below). The tube theorem is
  PROOFS-DRAFT 6.2 / Thm 6.1. The `chart_excl.py` docstring is corrected accordingly: a single
  non-central chart's own bound m_c can be negative.

### Fixes
- **FLAW-1:** `bnb2.half_w` gives h = up(max(up(hi - c), up(c - lo))) about the float centre.
  It is used at all four former `(hi - lo)/2` sites (lb_first, lb_second, lb_dual, lb_dual2).
  `chart_excl` and the arb codes have no such pattern.
- **FLAW-2:**
  - orientations t1, t2 in [-1e-12, fl(2pi/7) + 1e-12];
  - t0: slice 0 starts at -1e-12 and the last slice ends at fl(2pi/7) + 1e-12;
  - main sweep d starts at 1.8019377358048 (< 2u). PROOFS-DRAFT Lemma 1.2 closes the d = 2u
    sliver anyway.
  - Regions are generated in `src/modal_sweep.py` (`jobs`, `region`).
- **FLAW-3 / (E):** I chose **eta_c = dn((1 - q_c) rho_c - 2 off_c)**, with the chart window
  centred at Phi(z_c): g in [0, eta], |v.(z - z_c)| <= eta. Every target is then within
  eta + off <= (1 - q) rho of Phi(z_c).
  Sub-box radii are inflated by (1 + 1e-9), so there are no seams. The tube was re-run:
  `runs/charts_w030c.json` (sha256 3237cf54...). 1,530 certified charts, s in [-0.0301, +0.0500],
  max off 7.2e-18, eta >= 5.81e-5, B2 >= 0.0135, w >= 1.73e-4.
- **(G) gluing and (Ov) overlap: CERTIFIED in arb (200 bits)** by `src/tube_glue_arb.py` on
  `charts_w030c.json`. 1,529 consecutive pairs (both walks from the central chart), 0 failures.
  - (G): |J_c(z_{c'} - z_c)|_inf + ||J_c J_{c'}^{-1}||_inf (2 off_{c'}) <= rho_c. Worst ratio
    0.8744. Here 2 off >= off/(1 - q) bounds the distance to the exact curve point (q < 0.033).
  - (Ov): the windows [v.z - eta, v.z + eta] intersect. Worst gap upper bound -6.39e-5, i.e.
    they overlap.
- **UNCLEAR-1 (provenance):** every bnb2 run now prints, at the start: python version and
  platform, numpy version and file, the sha256 of rig.py, bnb2.py, chart_excl.py, the f file,
  dl_ref.json, the table and chart files and any resume file, all HEPT_* env vars, and argv.
- **UNCLEAR-2:** Modal checkpoints are kept: the resumed-from checkpoint is copied to
  used_*.npy, and a finished one is renamed *.final.npy on the Volume heptagon-ckpt.
- **H-DL (exactness), `src/dl_exact.py` (sympy, Q(c, s), with negative controls):**
  - the KK parallelogram identities;
  - the 5 contacts are exactly 0;
  - area(s2) * 2 delta_DL = A, i.e. area(s2) = A/(2 delta_DL) = acrit exactly;
  - the touching vertices lie strictly inside their edges (margins >= 0.2196);
  - the translation edge is exactly (1 + c, 0, 0);
  - [running] exact null vector: G v = 0.
  - The exact mu is the real root of grad b . v = 0, which is affine in mu with nonzero slope.
    Its arb enclosure [-1.168485374e-12 +- 1.9e-22] lies inside HEPT_MU (audit 2h).
- **Per-box certificate logging (item 6): FUTURE WORK.**
  - Logging (rule id + box) is ~100 bytes per box. At ~160M boxes per full run that is ~16 GB,
    over the few-GB budget.
  - A compact encoding (split-tree path bits + rule id, ~6-20 bytes/box) needs path tracking in
    the vectorised stack. Not done now, to keep these reruns unchanged in logic. Note this for
    the Lean/checker phase.
- **Re-runs (Modal, launched 2026-11-04):**
  - main sweep: 45 slices, DL slice split into 10;
  - obtuse D5 sweep: 45 slices;
  - all with fixed code: bnb2 9f4fb628..., rig 1ad0261b..., chart_excl 304d8743...,
    charts_w030c 3237cf54..., D5 87d669a7....
  - Local spot replays: main k5 and k39c, obtuse k13 and k30.
- H-DL, null vector: the sympy check passed rows 0 and 1 of G v = 0. The remaining rows had
  not finished when the 50-minute timeout killed the run: LUsolve expression swell
  (`runs/dl_exact.log`). **Not needed for rigour:**
  - the proof uses the exact null vector n of the exact matrix G (exact by H-DL);
  - the arb solve in `local_lemma_arb.py` gives a ball V containing n;
  - mu* is DEFINED by grad b(0) . n = 0, which is affine in mu with nonzero slope;
  - exactify's arb computation with the ball V gives a ball containing mu*, and that ball is
    inside HEPT_MU (audit 2h). The B&B carries coefficient 9 as that interval.

  So beta'(0) = grad b . Gamma'(0) = 0 holds exactly for the f with mu*. A cleaner symbolic n
  (elimination with reduction mod the relations at each step) is optional future work.

## 2026-11-05: MILESTONE: audit-fix reruns of BOTH sweeps complete (Modal)

- **Obtuse-mode sweep (D5), run a1: 45/45 DONE, 0 stuck.**
- **Main (nonobtuse) sweep, run a1: 54/54 jobs DONE, 0 stuck.** That is 44 slices, plus the DL
  slice split into 10 sub-slices (39a-39j).
  - Slice 38 reached TIMEOUT twice and returned a dump each time (a1: 16,658 open; a1r resumed
    from dump sha 0e58e965..., 19,980 open).
  - The a1r dump was split into 6 shards; their multiset equals the dump exactly. All shards
    DONE, 0 stuck (run a1s).
- Totals for **both sweeps combined** (post-resume work only for the resumed slices): main 91,520,252
  boxes / 48.4 h, obtuse 34,047,698 / 20.6 h; together 125,567,950 boxes / 68.9 container-hours
  (AUDIT2 2a).
- **Every log** starts with the provenance block. All 106 bnb2 logs share
  bnb2.py 9f4fb628..., rig.py 1ad0261b..., chart_excl.py 304d8743..., f_w030 ced5a7aa...,
  plus charts_w030c 3237cf54... (main) or D5 87d669a7... (obtuse), and numpy 2.5.4.
- Container restarts resumed from Volume checkpoints: main k2, k5, k39f, k39g; obtuse k7, k26,
  k44. Checkpoints are kept on the Volume (used_*.npy, *.final.npy).
- **Local spot replays** (Windows, numpy 2.5.4, lease, 2 threads):
  - obtuse k13 1,066,185 and k30 452,191 are identical to Modal;
  - main k5 and 39c are DONE with 0 stuck. Modal main k5 was resumed from a checkpoint, so only
    the verdict is comparable.
- This supersedes the earlier sweeps (runs/sweepR, runs/sweep, runs/sweepM_d5) for the
  certificate. Those ran before the FLAW-1/FLAW-2 fixes.

## 2026-11-05: second audit (AUDIT2-2026-11-05.md) follow-up

- **FLAW-4 (exported chart data rounded to nearest, not outward): CLOSED without a re-run.**
  `src/tube_recheck_arb.py` (log `runs/tube_recheck_w030c.log`, exit 0) widens every exported
  number by 2 ulp in the safe direction and rechecks, in arb (200 bits), for all 1530 charts:
  (A)-(D); K' (the excluded set lies in the inflated cube); E' (Banach reach on K');
  (G); (Ov); and the window check W (s* = v.z_{c'} lies in I_c). 1529 pairs, 0 failures.
  Worst ratios: E' 0.9999999995, G 0.8744, W 0.900. Worst Ov gap -6.39e-5.
  The stated argument is PROOFS-DRAFT Lemma 6.0, the inflated cube: n = 2 sub-boxes of radius
  (h/2)(1 + 1e-9) cover |zeta| <= rho (1 + 5e-10).
  Since (a) passed, main 39a-39j are not re-run.
- **UNCLEAR-5:** the (G)/(Ov) output is now committed (`runs/tube_glue_w030c.log`), and the
  window check is part of `tube_recheck_arb.py`.
- **UNCLEAR-4:** the six checkpoint-resumed slices (main k2, 39f, 39g; obtuse k7, k26, k44)
  are being replayed FROM SCRATCH on Modal (run rp1, --nockpt 1, cap 7.5 h, function timeout
  8 h). Every attempt's log is kept on the Volume (/ckpt/rp1/log_*).
- PROOFS-DRAFT updated: closed gaps marked, current bnb2 hash 9f4fb628 cited, the shared
  density definition, the de Berg statements quoted (no theorem numbers relied on), and
  Lemma 6.0 added.
- **UNCLEAR-4 CLOSED (2026-11-05):** all six checkpoint-resumed slices were replayed FROM SCRATCH
  on Modal (run rp1, no checkpointing, same bnb2 9f4fb628). **All DONE, 0 stuck:**
  - main: k2 499,793 boxes; 39f 2,367,945; 39g 8,430,765 (one preemption, restarted from
    scratch; both attempts' logs are on the Volume, /ckpt/rp1/log_main_39g_*);
  - obtuse: k7 1,267,429; k26 505,137; k44 855,327.
  Logs are in `runs/sweep_main_rp1/`, `runs/sweep_obt_rp1/`. The original logs record
  post-resume work only, so box counts are not comparable; the verdicts agree.
  **Every slice of both sweeps is now certified by a from-scratch run or a verified resume chain**
  (main 38: two dump resumes plus 6 shards, multiset-checked).

## 2026-11-05 referee pass (math-comms/ASTRA-HEPTAGON-REFEREE-2026-11-05.md)
- FLAW (statement gap) fixed: (M) split into (M-off) b >= tau off the chart set and (M-on) b >= acrit on it;
  Lemma 3.3 short-edge branch now b* >= acrit; averaging uses acrit only. Verified from logs/code: charts are
  enabled only in the main sweep (54/54 main logs HEPT_CHARTS, 0 HEPT_OBT; 45/45 obtuse logs HEPT_OBT4, 0 charts),
  and every chart region is within ~0.026 of the DL triangle (edges 1.83-1.92 < RTH, acute), so the short-edge
  branch with Din = Dout = 0 applies. PROOFS-DRAFT sections 1, 4, 6, 6a, 7.
- UNCLEAR (dl_exact) addressed: `src/dl_exact_cert.py` + `runs/dl_exact_cert.log` (arb, ALL CERTIFIED: c,s isolating
  intervals, all denominators nonzero, det(Asub) != 0, five contact interiorities strict, arb null vector).
  `dl_exact.py` only gained recording lists; its final symbolic LUsolve section still does not finish (not needed).
- Hypothesis-to-script linkage and UNCLEAR-4 status: PROOFS-DRAFT section 6a (slice 38 remains a resume chain).

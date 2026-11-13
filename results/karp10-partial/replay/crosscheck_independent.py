"""crosscheck_independent.py -- compare the replay's graph and homology routines with the SECOND implementations (indep_graph.py, indep_homology.py).

NOT part of replay_all.py (it needs numpy, sympy and z3 as well; the replay itself is standard library only).  Prints one line per section and
ends with CROSSCHECK: ALL AGREE (exit 0) or CROSSCHECK: DISAGREEMENT (exit 1).   usage:  python replay/crosscheck_independent.py  (about 5 minutes)

A. every (G, H) embedding query that the replay of all 20 certificates and the three further checks actually makes (harvested by wrapping
   replay_lib.find_embedding): result and witness compared with indep_graph.embeds (packing search); for equal edge counts the replay's isomorphism
   verdict is compared with equality of canonical forms; a sample is also given to z3.
B. the whole Z5 computation redone independently: orbits of Q, the 32 classes by canonical form, the poset by packing embedding, the 2,539 down-sets by
   antichain enumeration, absolute and relative homology of every down-set by Smith normal form over Z (and by dense F_5 rank), the 16 acyclic patterns,
   the Euler equations of the two extra groups, the 8 survivors, the nine named graphs, non-self-duality; compared with the replay's own sets and with
   note/z5_patterns.json and with the appendix table note/z5_table.tex (membership columns and e, degrees, t, q, complement of every class).
C. group witnesses of every certificate group re-checked with sympy.combinatorics (order, normality of N, |N| a power of p, <g,N> = Gamma) and edge
   orbits by breadth-first search.
D. randomised cross-tests: embedding (replay vs packing vs z3), isomorphism (replay vs canonical form), relative and absolute homology of random
   complexes (replay_smith.homology_by_size vs SNF / UCT vs dense numpy rank) for p = 2, 3, 5.
"""
import glob
import io
import itertools
import json
import os
import random
import re
import sys
import time
import contextlib

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, ".."))
sys.path.insert(0, here)
import replay_lib as RL                    # noqa: E402  (the implementation under test)
import replay_smith as RS                  # noqa: E402
import replay_cert as RC                   # noqa: E402
import indep_graph as IG                   # noqa: E402  (the second implementation)
import indep_homology as IH                # noqa: E402

dis = []                                   # disagreements, reported at the end
t0 = time.time()


def say(*a):
    print(*a, flush=True)


def note_dis(msg):
    dis.append(msg)
    say("   DISAGREEMENT:", msg)


# independent decoding of the 45-bit masks: bit k <-> the k-th pair (i<j) in lexicographic order
PAIRS = [(i, j) for i in range(10) for j in range(i + 1, 10)]


def rows_of_mask(mask):
    return IG.rows_from_edges([PAIRS[k] for k in range(45) if (mask >> k) & 1])


def mask_of_rows(rows):
    return sum(1 << PAIRS.index(e) for e in IG.edges_from_rows(rows))


# ----------------------------------------------------------------------------------------------------------------- A
say("== A. every embedding query of the replay, against the packing search and canonical forms")
calls = {}
orig = RL.find_embedding


def spy(g, h):
    r = orig(g, h)
    calls[(g, h)] = r
    return r


RL.find_embedding = spy
quiet = io.StringIO()
cert_files = sorted(glob.glob(os.path.join(root, "certs", "*.json")))
for f in cert_files:
    with contextlib.redirect_stdout(quiet):
        if "prime" in json.load(open(f)):
            rc = RS.main(f)
        else:
            nm = os.path.basename(f)
            decl = None
            if nm.startswith("type5_force"):
                a_, b_ = nm[len("type5_force_K"):-len(".json")].split("_K")
                decl = (int(a_), int(b_))
            rc = RC.main(f, False, decl)
    assert rc == 0, f
import replay_nonselfdual as RN            # noqa: E402
import check_z5_absolute as RA             # noqa: E402
import check_nine_graphs as RG             # noqa: E402
z5e = os.path.join(root, "certs", "lemma_z5_euler.json")
for fn in (lambda: RN.main(z5e), lambda: RA.main(z5e), lambda: RG.main()):
    with contextlib.redirect_stdout(quiet):
        assert fn() == 0
RL.find_embedding = orig
say(f"   harvested {len(calls)} distinct (G, H) embedding queries made while replaying all 23 rows "
    f"({sum(1 for v in calls.values() if v is not None)} answered 'embeds', {sum(1 for v in calls.values() if v is None)} 'does not embed')")
bad_a = 0
canon_cache = {}


def canon_of(mask):
    if mask not in canon_cache:
        canon_cache[mask] = IG.canonical(rows_of_mask(mask))
    return canon_cache[mask]


n_iso_pairs = 0
for (g, h), phi in calls.items():
    rg, rh = rows_of_mask(g), rows_of_mask(h)
    mine = IG.embeds(rg, rh)
    if mine and not IG.check_witness(IG.embedding(rg, rh), rg, rh):
        note_dis('own witness invalid')
    if mine != (phi is not None):
        bad_a += 1
        note_dis(f"embedding verdict differs for masks {g}, {h}: replay {phi is not None}, packing {mine}")
        continue
    if phi is not None:
        if not IG.check_witness(list(phi), rg, rh):
            bad_a += 1
            note_dis(f"the replay's witness is not an embedding by the independent check: {g}, {h}")
    if bin(g).count("1") == bin(h).count("1"):
        n_iso_pairs += 1
        iso_canon = canon_of(g) == canon_of(h)
        if iso_canon != (phi is not None):
            bad_a += 1
            note_dis(f"isomorphism verdict differs (canonical form) for {g}, {h}")
say(f"   all {len(calls)} verdicts agree with the packing search (witnesses re-checked); {n_iso_pairs} equal-edge-count pairs agree with canonical-form equality"
    if not bad_a else f"   {bad_a} problems")
# z3 on a sample of the harvested queries
try:
    import z3
    rnd = random.Random(11)
    sample = rnd.sample(sorted(calls), min(300, len(calls)))
    bad_z = 0
    for g, h in sample:
        rg, rh = rows_of_mask(g), rows_of_mask(h)
        s = z3.Solver()
        P = [z3.Int(f"p{i}") for i in range(10)]
        s.add(z3.Distinct(P))
        for i in range(10):
            s.add(P[i] >= 0, P[i] < 10)
        hadj = [[(rh[a] >> b) & 1 for b in range(10)] for a in range(10)]
        for a, b in IG.edges_from_rows(rg):
            # P[a], P[b] must be an edge of H: forbid every non-edge pair
            for x in range(10):
                for y in range(10):
                    if x != y and not hadj[x][y]:
                        s.add(z3.Not(z3.And(P[a] == x, P[b] == y)))
        res = s.check() == z3.sat
        if res != (calls[(g, h)] is not None):
            bad_z += 1
            note_dis(f"z3 disagrees on masks {g}, {h}")
    say(f"   z3 (third method) on a random sample of {len(sample)} harvested queries: {'all agree' if not bad_z else bad_z}")
except ImportError:
    say("   z3 not installed: sample skipped")

# ----------------------------------------------------------------------------------------------------------------- B
say("== B. the Z5 computation redone independently")
cert = json.load(open(z5e))
qgens = [tuple(g) for g in cert["gens"]]


def own_orbits(gens):
    """orbits on the 45 edges by breadth-first search under the generators"""
    seen = {}
    orbs = []
    for e in PAIRS:
        if e in seen:
            continue
        comp = [e]
        seen[e] = len(orbs)
        i = 0
        while i < len(comp):
            a, b = comp[i]
            i += 1
            for g in gens:
                x, y = g[a], g[b]
                f = (x, y) if x < y else (y, x)
                if f not in seen:
                    seen[f] = len(orbs)
                    comp.append(f)
        orbs.append(comp)
    return orbs


orbs = own_orbits(qgens)
k = len(orbs)
assert k == 9 and all(len(o) == 5 for o in orbs)
rep_orbs = [sorted(RL.edges_of(o)) for o in RL.edge_orbits(qgens)]
if sorted(sorted(o) for o in orbs) != sorted(rep_orbs):
    note_dis("edge orbits of Q differ between the two implementations")
union_rows = {}
for T in range(1, 1 << k):
    es = [e for i in range(k) if (T >> i) & 1 for e in orbs[i]]
    union_rows[T] = IG.rows_from_edges(es)
kind = {}
for T, r in union_rows.items():
    kind[T] = "D" if not IG.is_connected(r) else ("N" if not IG.is_connected(IG.complement_rows(r)) else "CC")
cnt = {s: sum(1 for v in kind.values() if v == s) for s in ("D", "N", "CC")}
say(f"   nonempty unions: {cnt} (note: 20 disconnected, 21 connected with disconnected complement, 470 in CC)")
if cnt != {"D": 20, "N": 21, "CC": 470}:
    note_dis(f"union statuses differ from the note: {cnt}")
canon_class = {}
cls_rows = []
cls_of = {}
for T in range(1, 1 << k):
    if kind[T] != "CC":
        continue
    c = IG.canonical(union_rows[T])
    if c not in canon_class:
        canon_class[c] = len(cls_rows)
        cls_rows.append(union_rows[T])
    cls_of[T] = canon_class[c]
nC = len(cls_rows)
say(f"   isomorphism classes of CC unions (canonical forms): {nC} (note: 32)")
if nC != 32:
    note_dis("number of classes is not 32")
# replay's classes, for comparison
r_orbs = RL.edge_orbits(qgens)
r_unions = {}
for T in range(1 << k):
    m = 0
    for i in range(k):
        if (T >> i) & 1:
            m |= r_orbs[i]
    r_unions[T] = m
# the replay's orbit order may differ from mine: compare class multisets by canonical form of each union mask
r_classes = {}
for T in range(1, 1 << k):
    if RS.connected(r_unions[T]) and RS.connected(RL.complement(r_unions[T])):
        r_classes.setdefault(canon_of(r_unions[T]), 0)
        r_classes[canon_of(r_unions[T])] += 1
if set(r_classes) != set(canon_class):
    note_dis("the replay's 32 classes are not the same graphs as the canonical-form classes")
else:
    say("   the replay's CC union classes are exactly the same 32 graphs (compared as canonical forms)")
sizes_mine = {}
for T, c in cls_of.items():
    sizes_mine[c] = sizes_mine.get(c, 0) + 1
# complement map
cmap = {}
for c in range(nC):
    cc = IG.canonical(IG.complement_rows(cls_rows[c]))
    cmap[c] = canon_class[cc]
assert all(cmap[cmap[c]] == c and cmap[c] != c for c in range(nC))
say("   complement map is a fixed-point-free involution on the 32 classes (canonical forms)")
# the poset by packing embedding
below = [[j for j in range(nC) if j != i and IG.n_edges(cls_rows[j]) < IG.n_edges(cls_rows[i]) and IG.embeds(cls_rows[j], cls_rows[i])] for i in range(nC)]
nrel = sum(len(b) for b in below)
say(f"   poset: {nrel} strict order relations among the 32 classes (packing search)")
# down-sets by antichain enumeration
above = [[i for i in range(nC) if j in below[i]] for j in range(nC)]
comparable = [[(j in below[i]) or (i in below[j]) or i == j for j in range(nC)] for i in range(nC)]
downs = set()


def rec_ac(start, chosen):
    ideal = set()
    for c in chosen:
        ideal.add(c)
        ideal.update(below[c])
    # closure under 'below' of below (transitive): below is already the full strict relation (embedding is transitive), but close anyway
    ch = True
    while ch:
        ch = False
        for c in list(ideal):
            for d in below[c]:
                if d not in ideal:
                    ideal.add(d)
                    ch = True
    downs.add(frozenset(ideal))
    for c in range(start, nC):
        if all(not comparable[c][d] for d in chosen):
            rec_ac(c + 1, chosen + [c])


rec_ac(0, [])
say(f"   down-sets by antichain enumeration: {len(downs)} (replay: 2,539)")
if len(downs) != 2539:
    note_dis(f"number of down-sets {len(downs)} != 2539")
# homology of every down-set
D_faces = {0} | {T for T in range(1, 1 << k) if kind[T] == "D"}
fbc = {}
for T, c in cls_of.items():
    fbc.setdefault(c, []).append(T)
hD = IH.homology(D_faces, primes=(5,))
say(f"   reduced homology of K_D: F_5 dims by face size {{s: dim}} = { {s: v for s, v in hD['fp'][5].items() if v} }; integral Betti {{s: b}} = "
    f"{ {s: v for s, v in hD['betti'].items() if v} }, torsion = { {s: t for s, t in hD['torsion'].items() if t} }")
expected_rel = {s + 1: v for s, v in hD["fp"][5].items() if v}
acyclic, relok = [], []
nonclosed = 0
zinfo = {}
t1 = time.time()
for X in sorted(downs, key=lambda s: (len(s), sorted(s))):
    faces_x = {T for c in X for T in fbc[c]}
    K = D_faces | faces_x
    habs = IH.homology(K, primes=(5,))
    if not habs["closed"]:
        nonclosed += 1
    hrel = IH.homology(K, L=D_faces, primes=(5,))
    absz = all(v == 0 for v in habs["fp"][5].values())
    relz = {s: v for s, v in hrel["fp"][5].items() if v} == expected_rel
    if absz:
        acyclic.append(X)
        zinfo[X] = (habs["betti"], habs["torsion"], habs["closed"])
    if relz:
        relok.append(X)
say(f"   all {len(downs)} down-sets: {len(acyclic)} have vanishing ABSOLUTE reduced F_5-homology (SNF/UCT = dense-rank, asserted for every complex), "
    f"{len(relok)} pass the relative test ({time.time()-t1:.0f} s)")
if set(acyclic) != set(relok):
    note_dis("absolute and relative criteria select different patterns (independent homology)")
if len(acyclic) != 16:
    note_dis(f"{len(acyclic)} acyclic patterns, replay says 16")
say(f"   of the {len(downs)} down-sets, {nonclosed} give a K that is not closed under taking faces (a union of orbits inside a CC union is of type N, not in P by (b));"
    f" among the 16 acyclic ones: {sum(1 for X in acyclic if not zinfo[X][2])} not closed")
say("   the 16 acyclic patterns over Z: " + "; ".join(f"|X|={len(X)}: betti {sum(zinfo[X][0].values())} torsion {[t for t in zinfo[X][1].values() if t]}" for X in acyclic[:16]))
# Euler equations of the two extra groups
survivors = list(acyclic)
for ent in cert["euler"]:
    gam = [tuple(x) for x in ent["gens"]]
    gorbs = own_orbits(gam)
    gk = len(gorbs)
    tab = {}
    for T in range(1, 1 << gk):
        es = [e for i in range(gk) if (T >> i) & 1 for e in gorbs[i]]
        r = IG.rows_from_edges(es)
        if not IG.is_connected(r):
            tab[T] = ("D", None)
        elif not IG.is_connected(IG.complement_rows(r)):
            tab[T] = ("N", None)
        else:
            tab[T] = ("CC", canon_class[IG.canonical(r)])
    keep = []
    for X in survivors:
        val = sum((-1) ** (bin(T).count("1") - 1) * (1 if s == "D" else (0 if s == "N" else (1 if c in X else 0))) for T, (s, c) in tab.items())
        if val == 1:
            keep.append(X)
    say(f"   Euler equation of the group of order {ent['order']} ({gk} orbits): {len(keep)} of {len(survivors)} survive")
    survivors = keep
if len(survivors) != 8:
    note_dis(f"{len(survivors)} survivors, replay says 8")
sizes = sorted(len(X) for X in survivors)
say(f"   survivors have sizes {sizes} (note: 9, 11, 11, 13, 19, 21, 21, 23); a self-dual trace would need 16: {'present' if 16 in sizes else 'absent'}")
if sizes != [9, 11, 11, 13, 19, 21, 21, 23]:
    note_dis("pattern sizes differ from the note")
dual = lambda X: frozenset(c for c in range(nC) if cmap[c] not in X)
if not all(dual(X) in set(survivors) and dual(X) != X for X in survivors):
    note_dis("survivors not closed under duality or a self-dual one exists")
else:
    say("   survivors closed under duality, none self-dual")
# the nine graphs
ins = [canon_class[IG.canonical(IG.rows_from_edges([tuple(e) for e in g]))] for g in cert["claims_in"]]
outs = [canon_class[IG.canonical(IG.rows_from_edges([tuple(e) for e in g]))] for g in cert["claims_out"]]
ok9 = all(c in X for X in survivors for c in ins) and all(c not in X for X in survivors for c in outs) and sorted(cmap[c] for c in ins) == sorted(outs)
say(f"   nine claimed-in graphs are in all 8 patterns, nine claimed-out in none, and they are the complements of each other: {ok9}")
if not ok9:
    note_dis("nine graphs claim not reproduced")
# the replay's own sets (re-run its recursion, import its functions) compared with mine, as sets of canonical forms
r_cc = {}
for T in range(1, 1 << k):
    if RS.connected(r_unions[T]) and RS.connected(RL.complement(r_unions[T])):
        r_cc[T] = canon_of(r_unions[T])
mine_as_canon = {frozenset(IG.canonical(cls_rows[c]) for c in X) for X in survivors}
zp = json.load(open(os.path.join(root, "note", "z5_patterns.json")))
label_canon = {c["label"]: IG.canonical(IG.rows_from_edges([tuple(e) for e in c["rep_edges"]])) for c in zp["classes"]}
json_sets = {frozenset(label_canon[l] for l in labs) for labs in zp["patterns"].values()}
if mine_as_canon != json_sets:
    note_dis("the 8 independent patterns differ from note/z5_patterns.json")
else:
    say("   the 8 independent patterns equal the 8 patterns of note/z5_patterns.json (as sets of canonical forms)")
if len(set(label_canon.values())) != 32 or set(label_canon.values()) != set(canon_class):
    note_dis("z5_patterns.json classes are not the 32 classes")
# replay's 16 patterns recomputed with replay's own functions, compared with my 16
r_below = None
r_classes_list = []
r_cls_of = {}
for T in range(1, 1 << k):
    if RS.connected(r_unions[T]) and RS.connected(RL.complement(r_unions[T])):
        for ci, rep in enumerate(r_classes_list):
            if RL.isomorphic(rep, r_unions[T]):
                r_cls_of[T] = ci
                break
        else:
            r_classes_list.append(r_unions[T])
            r_cls_of[T] = len(r_classes_list) - 1
r_nC = len(r_classes_list)
r_below = [[j for j in range(r_nC) if j != i and RL.popcount(r_classes_list[j]) < RL.popcount(r_classes_list[i]) and RL.embeds(r_classes_list[j], r_classes_list[i])] for i in range(r_nC)]
r_dfaces = {0: [0]}
for T in range(1, 1 << k):
    if not RS.connected(r_unions[T]):
        r_dfaces.setdefault(RL.popcount(T), []).append(T)
r_hD = {s: v for s, v in RS.homology_by_size(r_dfaces, 5).items() if v}
r_exp = {s + 1: v for s, v in r_hD.items()}
r_fbc = {}
for T, c in r_cls_of.items():
    r_fbc.setdefault(c, []).append(T)
r_good = []


def r_rec(pos, cur, order):
    if pos == r_nC:
        faces = {}
        for c in cur:
            for T in r_fbc[c]:
                faces.setdefault(RL.popcount(T), []).append(T)
        H = {s: v for s, v in RS.homology_by_size(faces, 5).items() if v} if faces else {}
        if H == r_exp:
            r_good.append(frozenset(cur))
        return
    i = order[pos]
    r_rec(pos + 1, cur, order)
    if all(j in cur for j in r_below[i]):
        r_rec(pos + 1, cur | {i}, order)


r_rec(0, frozenset(), sorted(range(r_nC), key=lambda i: RL.popcount(r_classes_list[i])))
r_as_canon = {frozenset(canon_of(r_classes_list[c]) for c in X) for X in r_good}
mine16 = {frozenset(IG.canonical(cls_rows[c]) for c in X) for X in acyclic}
if r_as_canon != mine16:
    note_dis("the replay's 16 patterns differ from the independent 16")
else:
    say(f"   the replay's {len(r_good)} filter survivors equal the independent {len(acyclic)} absolutely acyclic patterns (sets of canonical forms)")

# the appendix table (tex) against recomputed invariants and the json
tex = open(os.path.join(root, "note", "z5_table.tex"), encoding="utf-8").read()
rows_tex = [l for l in tex.splitlines() if "&" in l and not l.startswith("class") and "\\toprule" not in l and r"\midrule" not in l]
rows_tex = [l for l in rows_tex if re.match(r"\s*\d+", l)]
cells_all = []
for l in rows_tex:
    cs = [c.strip() for c in l.replace("\\\\", "").split("&")]
    cells_all.append(cs)
tab_rows = {}
for cs in cells_all:
    for half in (cs[:14], cs[14:28]):
        if len(half) < 14 or not half[0]:
            continue
        lab = int(re.match(r"\d+", half[0]).group())
        tab_rows[lab] = dict(star="*" in half[0], e=int(half[1]), deg=half[2], t=int(half[3]), q=int(half[4]), cpl=int(half[5]), pats=[j + 1 for j in range(8) if "bullet" in half[6 + j]])
say(f"   the appendix table has {len(tab_rows)} class rows")
bad_t = 0
for c in zp["classes"]:
    row = tab_rows.get(c["label"])
    r = IG.rows_from_edges([tuple(e) for e in c["rep_edges"]])
    degs = sorted((bin(x).count("1") for x in r), reverse=True)
    tri = sum(1 for a, b, d in itertools.combinations(range(10), 3) if (r[a] >> b) & 1 and (r[a] >> d) & 1 and (r[b] >> d) & 1)
    quad = 0
    for a, b in itertools.combinations(range(10), 2):
        common = bin(r[a] & r[b] & ~((1 << a) | (1 << b))).count("1")
        quad += common * (common - 1) // 2
    quad //= 2
    cplab = next(cc["label"] for cc in zp["classes"] if label_canon[cc["label"]] == IG.canonical(IG.complement_rows(r)))
    pats = [j for j in range(1, 9) if c["label"] in zp["patterns"][str(j)]]
    if row is None:
        bad_t += 1
        note_dis(f"class {c['label']} missing from the table")
        continue
    from collections import Counter
    ds = Counter(degs)
    want_deg = "".join(f"{d}^{{{ds[d]}}}" for d in sorted(ds, reverse=True))
    got_deg = row["deg"].replace("$", "").replace("\\,", "")
    if not (row["e"] == IG.n_edges(r) and got_deg == want_deg and row["t"] == tri and row["q"] == quad and row["cpl"] == cplab and row["pats"] == pats):
        bad_t += 1
        note_dis(f"table row of class {c['label']} differs: tex {row} vs recomputed e={IG.n_edges(r)} deg={want_deg} t={tri} q={quad} cpl={cplab} pats={pats}")
say("   every table row (e, degrees, triangles, 4-cycles, complement, membership in the 8 patterns, star marks pending) matches the recomputation" if not bad_t else f"   {bad_t} table rows differ")
stars = sorted(l for l, r in tab_rows.items() if r["star"])
nine_labels = sorted(label_canon_l for label_canon_l in [next(l for l, cc in label_canon.items() if cc == IG.canonical(IG.rows_from_edges([tuple(e) for e in g]))) for g in cert["claims_in"]])
if stars != nine_labels:
    note_dis(f"starred classes {stars} are not the nine graphs {nine_labels}")
else:
    say(f"   the starred classes of the table {stars} are exactly the nine graphs of the certificate")
# statements of the appendix text
pj = {int(j): set(v) for j, v in zp["patterns"].items()}
stmt = (pj[1] == set(nine_labels)) and (pj[2] - pj[1] == {4, 9}) and (pj[3] - pj[1] == {4, 15}) and (pj[4] - pj[1] == {4, 9, 15, 20}) \
    and all(pj[9 - j] == {l for l in range(1, 33) if next(c for c in zp['classes'] if c['label'] == l)['complement'] not in pj[j]} for j in range(1, 9))
say(f"   appendix text: pattern 1 = the nine graphs; 2, 3 add {{4,9}}, {{4,15}}; 4 adds {{4,9,15,20}}; pattern j dual to 9-j: {stmt}")
if not stmt:
    note_dis("an appendix sentence about the patterns is not reproduced")

# ----------------------------------------------------------------------------------------------------------------- C
say("== C. group witnesses of all certificate groups with sympy.combinatorics")
from sympy.combinatorics import Permutation, PermutationGroup


def PG(gl):
    return PermutationGroup([Permutation(list(g)) for g in gl])


ngroups = 0
bad_c = 0
for f in cert_files:
    c = json.load(open(f))
    ents = []
    if "groups" in c:
        ents = [(f"{os.path.basename(f)}:{n}", g) for n, g in c["groups"].items()]
    for i, e in enumerate(c.get("euler", [])):
        ents.append((f"{os.path.basename(f)}:euler{i}", e))
    for name, g in ents:
        ngroups += 1
        G = PG(g["gens"])
        Nn = PG(g["N_gens"]) if g["N_gens"] else None
        p = g["p"]
        ok = G.order() == g["order"] and Nn is not None and all(G.contains(Permutation(list(x))) for x in g["N_gens"]) and Nn.is_normal(G)
        no = Nn.order()
        t = no
        while t % p == 0:
            t //= p
        ok = ok and t == 1 and PG(g["N_gens"] + [g["g"]]).order() == G.order()
        # quotient order and cyclicity: |G/N| = order/|N|; G = <g>N means G/N is generated by the image of g
        if not ok:
            bad_c += 1
            note_dis(f"group witness fails with sympy: {name}")
        # fixed edges for lk/dl equations: orbits by BFS agree with the replay
        ob = own_orbits([tuple(x) for x in g["gens"]])
        rb = [sorted(RL.edges_of(o)) for o in RL.edge_orbits([tuple(x) for x in g["gens"]])]
        if sorted(sorted(o) for o in ob) != sorted(rb):
            bad_c += 1
            note_dis(f"edge orbits differ for {name}")
say(f"   {ngroups} group witnesses (Gamma, N, g): order of Gamma, N normal in Gamma, |N| a power of p, Gamma = <g, N>, edge orbits by BFS: "
    f"{'all verified' if not bad_c else bad_c}")
# fixed edges: for every (lk)/(dl) equation the named edge is a fixed edge (orbit of size one)
for f in cert_files:
    c = json.load(open(f))
    for eq in c.get("equations", []):
        if eq["kind"] != "fix":
            g = c["groups"][eq["group"]]
            a, b = eq["fixed_edge"]
            e = (a, b) if a < b else (b, a)
            if not all(((x[a], x[b]) if x[a] < x[b] else (x[b], x[a])) == e for x in g["gens"]):
                note_dis(f"{os.path.basename(f)}: {eq['group']} does not fix {e}")
say("   every (lk)/(dl) equation names an edge fixed by all generators of its group")

# ----------------------------------------------------------------------------------------------------------------- D
say("== D. randomised cross-tests")
rnd = random.Random(20261111)


def rand_graph(dens):
    return [e for e in PAIRS if rnd.random() < dens]


def perm_edges(E, p):
    return [(p[a], p[b]) for a, b in E]


n_emb = 0
bad_d = 0
t1 = time.time()
for trial in range(6000):
    dens = rnd.choice([0.1, 0.2, 0.3, 0.45, 0.6, 0.8])
    H = rand_graph(dens)
    p = list(range(10))
    rnd.shuffle(p)
    mode = trial % 4
    if mode == 0:                     # a true subgraph of a relabelled H
        G = [e for e in perm_edges(H, p) if rnd.random() < 0.8]
    elif mode == 1:                   # a true subgraph plus one extra edge (often false)
        G = [e for e in perm_edges(H, p) if rnd.random() < 0.9]
        extra = rnd.choice(PAIRS)
        if extra not in G and (extra[1], extra[0]) not in G:
            G.append(extra)
    elif mode == 2:                   # unrelated graph with similar density
        G = rand_graph(dens * rnd.choice([0.5, 0.8, 1.0]))
    else:                             # structured: a double edge switch of a relabelled H (same degree sequence) as G, or H itself
        G = perm_edges(H, p)
        if len(G) >= 2:
            (a, b), (c, d) = rnd.sample(G, 2)
            if len({a, b, c, d}) == 4 and (min(a, c), max(a, c)) not in [(min(x), max(x)) for x in G] and (min(b, d), max(b, d)) not in [(min(x), max(x)) for x in G]:
                G = [e for e in G if e not in ((a, b), (c, d))] + [(a, c), (b, d)]
    gm = RL.mask_of([tuple(e) for e in G])
    hm = RL.mask_of([tuple(e) for e in H])
    mine = IG.embeds(rows_of_mask(gm), rows_of_mask(hm))
    theirs = RL.embeds(gm, hm)
    n_emb += 1
    if mine != theirs:
        bad_d += 1
        note_dis(f"random embedding query disagrees: {gm} {hm}")
    if RL.popcount(gm) == RL.popcount(hm):
        if (IG.canonical(rows_of_mask(gm)) == IG.canonical(rows_of_mask(hm))) != RL.isomorphic(gm, hm):
            bad_d += 1
            note_dis(f"random isomorphism query disagrees: {gm} {hm}")
say(f"   {n_emb} random embedding queries (about {n_emb//4} each: true subgraph, subgraph + 1 edge, unrelated, degree-preserving switch), embeddings replay vs packing and isomorphism replay vs canonical form (equal edge counts): "
    f"{'all agree' if not bad_d else bad_d}  ({time.time()-t1:.0f} s)")
# random pairs of regular graphs / vertex-transitive-like: circulants with random connection sets
bad_circ = 0
n_circ = 0
for _ in range(400):
    S1 = sorted(rnd.sample([1, 2, 3, 4, 5], rnd.randint(1, 4)))
    S2 = sorted(rnd.sample([1, 2, 3, 4, 5], rnd.randint(1, 4)))
    def circ(S):
        return [(i, (i + s) % 10) for i in range(10) for s in S if i != (i + s) % 10]
    g1 = RL.mask_of(circ(S1))
    g2 = RL.mask_of(circ(S2))
    n_circ += 1
    if (IG.canonical(rows_of_mask(g1)) == IG.canonical(rows_of_mask(g2))) != RL.isomorphic(g1, g2) or IG.embeds(rows_of_mask(g1), rows_of_mask(g2)) != RL.embeds(g1, g2):
        bad_circ += 1
        note_dis(f"circulant pair disagrees {S1} {S2}")
say(f"   {n_circ} pairs of random circulant graphs (many isomorphic-looking, large automorphism groups): {'all agree' if not bad_circ else bad_circ}")
# z3 on random embedding queries
try:
    import z3
    bad_z = 0
    for trial in range(250):
        H = rand_graph(rnd.choice([0.3, 0.5, 0.7]))
        G = [e for e in perm_edges(H, rnd.sample(range(10), 10)) if rnd.random() < 0.85]
        if trial % 2:
            extra = rnd.choice(PAIRS)
            if extra not in G:
                G.append(extra)
        gm, hm = RL.mask_of([tuple(e) for e in G]), RL.mask_of([tuple(e) for e in H])
        rg, rh = rows_of_mask(gm), rows_of_mask(hm)
        s = z3.Solver()
        P = [z3.Int(f"q{i}") for i in range(10)]
        s.add(z3.Distinct(P))
        for i in range(10):
            s.add(P[i] >= 0, P[i] < 10)
        hadj = [[(rh[a] >> b) & 1 for b in range(10)] for a in range(10)]
        for a, b in IG.edges_from_rows(rg):
            for x in range(10):
                for y in range(10):
                    if x != y and not hadj[x][y]:
                        s.add(z3.Not(z3.And(P[a] == x, P[b] == y)))
        if (s.check() == z3.sat) != RL.embeds(gm, hm) or (s.check() == z3.sat) != IG.embeds(rg, rh):
            bad_z += 1
            note_dis("z3 disagrees on a random embedding query")
    say(f"   250 random embedding queries, z3 vs replay vs packing: {'all agree' if not bad_z else bad_z}")
except ImportError:
    pass
# homology of random complexes
bad_h = 0
nh = 0
t1 = time.time()
for trial in range(700):
    nv = rnd.randint(4, 8)
    nfac = rnd.randint(2, 9)
    facets = [rnd.randrange(1, 1 << nv) & rnd.randrange(1, 1 << nv) | (1 << rnd.randrange(nv)) for _ in range(nfac)]
    K = set()
    for f in facets:
        sub = f
        while True:
            K.add(sub)
            if sub == 0:
                break
            sub = (sub - 1) & f
    K.add(0)
    # subcomplex L: the down-closure of a random subset of faces of K (possibly empty-face only)
    L = set()
    for f in rnd.sample(sorted(K), rnd.randint(0, min(4, len(K)))):
        sub = f
        while True:
            L.add(sub)
            if sub == 0:
                break
            sub = (sub - 1) & f
    for primes in ((2, 3, 5),):
        Hrel = IH.homology(K, L=L, primes=primes)           # asserts SNF/UCT == dense rank internally
        for p in primes:
            faces = {}
            for t in K - L:
                faces.setdefault(bin(t).count("1"), []).append(t)
            theirs = RS.homology_by_size(faces, p) if faces else {}
            mine = {s: v for s, v in Hrel["fp"][p].items() if v}
            theirs = {s: v for s, v in theirs.items() if v}
            nh += 1
            if mine != theirs:
                bad_h += 1
                note_dis(f"random pair homology disagrees p={p}: replay {theirs} independent {mine}")
say(f"   {nh} (complex, subcomplex, prime) homology comparisons, replay F_p elimination vs SNF/UCT vs dense numpy rank: {'all agree' if not bad_h else bad_h}  ({time.time()-t1:.0f} s)")
# known torsion: RP^2 (6-vertex triangulation) has Z/2 torsion in H_1: F_2 and F_3 dimensions must differ
rp2 = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 5, 1), (1, 2, 4), (2, 3, 5), (3, 4, 1), (4, 5, 2), (5, 1, 3)]
K = set()
for f in rp2:
    m = sum(1 << v for v in f)
    sub = m
    while True:
        K.add(sub)
        if sub == 0:
            break
        sub = (sub - 1) & m
hr = IH.homology(K, primes=(2, 3, 5))
say(f"   control, 6-vertex RP^2: integral reduced Betti {hr['betti']} torsion {hr['torsion']}; F_2 dims {hr['fp'][2]}, F_3 dims {hr['fp'][3]}")
if not (hr["torsion"][2] == [2] and hr["fp"][2] != hr["fp"][3]):
    note_dis("RP^2 torsion control failed")
faces = {}
for t in K:
    faces.setdefault(bin(t).count("1"), []).append(t)
if {s: v for s, v in RS.homology_by_size(faces, 2).items() if v} != {s: v for s, v in hr["fp"][2].items() if v} or \
        {s: v for s, v in RS.homology_by_size(faces, 3).items() if v} != {s: v for s, v in hr["fp"][3].items() if v}:
    note_dis("RP^2: replay homology differs")

say(f"== done in {time.time()-t0:.0f} s")
if dis:
    say("CROSSCHECK: DISAGREEMENT (%d)" % len(dis))
    sys.exit(1)
say("CROSSCHECK: ALL AGREE")

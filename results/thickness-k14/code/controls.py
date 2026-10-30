"""Positive/negative controls for renc.py.
  table   : small-n table (K4/1, K5/1, K8/2, K9/2, K12/3) with the full encoding; "table13" adds
            K13/3 (SAT; about 5 min with pysat CaDiCaL 1.9.5).
  real    : real drawings (Seamons K13 witness, random integer drawings coloured by SAT) must
            satisfy the full encoding (incl. normal form + Euler) after the referee's own
            normalisation, checked by unit-propagating the induced X and C into a solver.
"""
import itertools, json, random, sys, time, os
from fractions import Fraction
from pysat.solvers import Cadical195
from renc import Enc

HERE = os.path.dirname(os.path.abspath(__file__))


def orient(P, a, b, c):
    (ax, ay), (bx, by), (cx, cy) = P[a], P[b], P[c]
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def crosses(P, a, b, c, d):
    return orient(P, a, b, c) * orient(P, a, b, d) < 0 and orient(P, c, d, a) * orient(P, c, d, b) < 0


def check_drawing(P, layer):
    """exact check: distinct points, no collinear triple, every layer crossing-free"""
    m = len(P)
    assert len(set(map(tuple, P))) == m
    for a, b, c in itertools.combinations(range(m), 3):
        assert orient(P, a, b, c) != 0, "collinear"
    E = list(itertools.combinations(range(m), 2))
    assert set(layer) == set(E)
    bad = 0
    for e, f in itertools.combinations(E, 2):
        if len(set(e) | set(f)) == 4 and layer[e] == layer[f] and crosses(P, *e, *f):
            bad += 1
    return bad


def normalise(P):
    """Referee's normal form: p = lowest (then leftmost) point, others sorted ccw around p.
    Returns label order: new label i+1 -> original index order[i]."""
    m = len(P)
    p = min(range(m), key=lambda i: (P[i][1], P[i][0]))
    rest = [i for i in range(m) if i != p]
    import functools
    def cmp(i, j):
        o = orient(P, p, i, j)
        return -1 if o > 0 else 1
    rest.sort(key=functools.cmp_to_key(cmp))
    return [p] + rest


def induced_units(enc, P, order, layer, sign=1):
    n = enc.n
    units = []
    for (a, b, c), v in enc.X.items():
        o = orient(P, order[a - 1], order[b - 1], order[c - 1]) * sign
        units.append(v if o > 0 else -v)
    # layers: permute layer names so that value precedence holds (first-use order along lex edges)
    ren = {}
    for (a, b) in enc.edges:
        e = tuple(sorted((order[a - 1], order[b - 1])))
        L = layer[e]
        if L not in ren:
            ren[L] = len(ren)
    for (a, b) in enc.edges:
        e = tuple(sorted((order[a - 1], order[b - 1])))
        for k in range(enc.t):
            v = enc.C[((a, b), k)]
            units.append(v if ren[layer[e]] == k else -v)
    return units


def accept(P, layer, t, **kw):
    n = len(P)
    order = normalise(P)
    enc = Enc(n, t, **kw)
    units = induced_units(enc, P, order, layer)
    with Cadical195(bootstrap_with=enc.cl) as s:
        return s.solve(assumptions=units)


def colour_points(P, t, seed=0):
    """find a t-layer colouring of the straight-line K_m on P (graph colouring of crossing graph)"""
    m = len(P)
    E = list(itertools.combinations(range(m), 2))
    idx = {e: i for i, e in enumerate(E)}
    var = lambda e, k: idx[e] * t + k + 1
    cl = [[var(e, k) for k in range(t)] for e in E]
    for e, f in itertools.combinations(E, 2):
        if len(set(e) | set(f)) == 4 and crosses(P, *e, *f):
            cl += [[-var(e, k), -var(f, k)] for k in range(t)]
    with Cadical195(bootstrap_with=cl) as s:
        if not s.solve():
            return None
        M = set(x for x in s.get_model() if x > 0)
    return {e: next(k for k in range(t) if var(e, k) in M) for e in E}


def table(k13=False):
    rows = []
    cases = [(4, 1, True), (5, 1, False), (8, 2, True), (9, 2, False), (12, 3, True)]
    for n, t, exp in cases + ([(13, 3, True)] if k13 else []):
        for euler in (False, True):
            enc = Enc(n, t, euler=euler)
            t0 = time.time()
            with Cadical195(bootstrap_with=enc.cl) as s:
                r = s.solve()
            rows.append((n, t, euler, r, exp, round(time.time() - t0, 1)))
            print("K%d t=%d euler=%s -> %s (expected %s) %.1fs" % rows[-1], flush=True)
            assert r == exp
    return rows


def real():
    W = json.load(open(os.path.join(HERE, "..", "certs", "witnesses", "k13_seamons.json")))
    P = [tuple(p) for p in W["points"]]
    layer = {(min(u, v), max(u, v)): L for u, v, L in W["edges"]}
    bad = check_drawing(P, layer)
    print("Seamons K13: exact crossings within layers =", bad)
    assert bad == 0
    for euler in (False, True):
        r = accept(P, layer, 3, euler=euler)
        print("Seamons K13 accepted by referee encoding (nf, prec, euler=%s): %s" % (euler, r))
        assert r
    # random real drawings: n points, t = smallest colourable from 3 up
    rnd = random.Random(2026)
    cnt = 0
    for trial in range(40):
        n = rnd.choice([9, 10, 11, 12])
        P = [(rnd.randint(-10**6, 10**6), rnd.randint(-10**6, 10**6)) for _ in range(n)]
        if trial % 2:  # small hull: big triangle + interior points
            P = [(-10**7, -10**7), (10**7, -10**7 + 3), (17, 10**7)] + P[3:]
        if any(orient(P, *tr) == 0 for tr in itertools.combinations(range(n), 3)):
            continue
        for t in (3, 4, 5, 6, 7, 8):
            lay = colour_points(P, t)
            if lay:
                break
        assert check_drawing(P, lay) == 0
        r = accept(P, lay, t)
        hull = sum(1 for v in range(n) if not any(
            len({v, a, b, c}) == 4 and orient(P, a, b, c) * orient(P, a, b, v) > 0 and
            orient(P, a, b, c) * orient(P, b, c, v) > 0 and orient(P, a, b, c) * orient(P, c, a, v) > 0
            for a, b, c in itertools.combinations(range(n), 3)))
        print("random K%d h=%d t=%d accepted=%s" % (n, hull, t, r), flush=True)
        assert r
        cnt += 1
    print("random real drawings accepted:", cnt)


if __name__ == "__main__":
    {"table": table, "table13": lambda: table(True), "real": real}[sys.argv[1]]()

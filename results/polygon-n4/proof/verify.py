"""Independent checker for the n = 4 polygon-conjecture proof.  See ../PROOF.md.

Usage:  python verify.py

Exact arithmetic only: sympy for the identity and for the face polynomials, python-flint
for polynomial arithmetic in the proof-tree replay, python integers for Bernstein
coefficients.  It does not import any of the search code (faces.py, phiface.py,
bern.py, prover2.py); it only reads the proof trees they produced.
"""
import itertools, json, glob, sys
from math import factorial, comb
import sympy as sp
import flint

# ---------- 1. the algebraic criterion ---------------------------------------------
def check_identity():
    a, b, c, r, x, y, t = sp.symbols('a b c r x y t', real=True)
    p = t**3 - a*t**2 + b*t - c
    disc = sp.discriminant(p, t)
    Phi = -disc - 4*(a**2 - 3*b + 4)**2
    sub = {a: r + 2*x, b: 2*r*x + x**2 + y**2, c: r*(x**2 + y**2)}
    lhs = sp.expand(Phi.subs(sub, simultaneous=True))
    d = r - x
    rhs = sp.expand(4*(y**2 - 1)*(d**2 + (y + 2)**2)*(d**2 + (y - 2)**2))
    assert sp.expand(lhs - rhs) == 0
    return Phi, (a, b, c)

# ---------- 2. faces -----------------------------------------------------------------
OFF = [(i, j) for i in range(4) for j in range(4) if i != j]
IDX = {e: k for k, e in enumerate(OFF)}

def faces_all():
    out = []
    for f in itertools.product('01*', repeat=12):
        if f.count('*') != 6:
            continue
        if any(sum(f[IDX[(i, j)]] == '*' for j in range(4) if j != i) > 2 for i in range(4)):
            continue
        out.append(''.join(f))
    return out

def canon(f):
    best = None
    for p in itertools.permutations(range(4)):
        for comp in (0, 1):
            g = [None]*12
            for k, (i, j) in enumerate(OFF):
                v = f[k]
                if comp and v != '*':
                    v = '1' if v == '0' else '0'
                g[IDX[(p[i], p[j])]] = v
            g = ''.join(g)
            if best is None or g < best:
                best = g
    return best

# ---------- 3. face polynomial  -Phi(a,b,c)  (sympy) ---------------------------------
def face_poly(face, Phi, abc):
    ts = sp.symbols('t0:6')
    M = sp.zeros(4, 4)
    c = 0
    for k, (i, j) in enumerate(OFF):
        v = face[k]
        if v == '*':
            M[i, j] = -ts[c]; c += 1
        elif v == '1':
            M[i, j] = -1
    for i in range(4):
        M[i, i] = -sum(M[i, j] for j in range(4) if j != i)
    A = M.trace()
    B = sum(M.extract([i, j], [i, j]).det() for i, j in itertools.combinations(range(4), 2))
    C = sum(M.extract(list(S), list(S)).det() for S in itertools.combinations(range(4), 3))
    # det(tI-M) = t^4 - A t^3 + B t^2 - C t + det M, and det M = 0
    assert sp.expand(M.det()) == 0
    N = sp.expand(-Phi.subs(dict(zip(abc, (A, B, C))), simultaneous=True))
    return sp.Poly(N, *ts)

# ---------- 4. Bernstein check (own implementation) -----------------------------------
def bern_coeffs(terms, k):
    """Bernstein coefficients on [0,1]^k, each axis scaled by d_l! (a positive constant)."""
    d = [max([m[l] for m in terms] + [0]) for l in range(k)]
    arr = dict(terms)
    for l in range(k):
        dl = d[l]; new = {}
        for m, cf in arr.items():
            j = m[l]
            for i in range(j, dl + 1):
                key = m[:l] + (i,) + m[l+1:]
                new[key] = new.get(key, 0) + cf*comb(i, j)*factorial(j)*factorial(dl - j)
        arr = new
    full = 1
    for x in d:
        full *= x + 1
    return arr, full, d

def bern_nonneg(P, ctx, depth=0, maxdepth=14, budget=None):
    """True only if P >= 0 on [0,1]^k is certified by Bernstein coefficients,
    possibly after bisections (each half re-expanded from scratch)."""
    if budget is None:
        budget = [0]
    budget[0] += 1
    if budget[0] > 6000:
        return False
    k = ctx.nvars()
    raw = P.to_dict()
    if not raw:
        return True
    L = 1
    for v in raw.values():
        L = sp.ilcm(L, int(v.q))
    terms = {tuple(int(e) for e in m): int(v*L) for m, v in raw.items()}
    arr, full, deg = bern_coeffs(terms, k)
    vals = list(arr.values())
    if len(arr) < full:
        vals.append(0)
    if min(vals) >= 0:
        return True
    if depth >= maxdepth:
        return False
    act = [l for l in range(k) if deg[l] > 0]
    if not act:
        return False  # negative constant: nothing to bisect
    def spread(l):
        s = 0
        for m, cf in arr.items():
            if m[l] < deg[l]:
                m2 = m[:l] + (m[l] + 1,) + m[l+1:]
                s = max(s, abs(arr.get(m2, 0) - cf))
        return s
    l = max(act, key=spread)
    g = ctx.gens()
    for side in (0, 1):
        args = list(g)
        args[l] = g[l]/2 if side == 0 else (g[l] + 1)/2
        if not bern_nonneg(P.compose(*args), ctx, depth + 1, maxdepth, budget):
            return False
    return True

# ---------- 5. proof-tree replay ---------------------------------------------------
def replay(P, node, ctx):
    g = ctx.gens(); names = ctx.names(); t = node['t']; k = ctx.nvars()
    if t == 'zero':
        return P.is_zero()
    if t == 'const':
        d = P.to_dict()
        return all(m == (0,)*k for m in d) and d.get((0,)*k, 0) >= 0
    if t == 'bern':
        return bern_nonneg(P, ctx)
    if t == 'factor':
        c, facs = P.factor()
        prod = ctx.from_dict({(0,)*k: c})
        for f, e in facs:
            prod = prod * f**e
        if prod != P:
            return False
        sgn = 1 if c > 0 else -1
        kids = {kd['f']: kd for kd in node['odd']}
        for f, e in facs:
            if e % 2 == 0:
                continue
            kd = kids.get(str(f))
            if kd is None or not replay(kd['s']*f, kd['p'], ctx):
                return False
            sgn *= kd['s']
        return sgn > 0
    if t == 'mono':
        l = node['var']; at = node['at']
        D = P.derivative(names[l])
        if not bern_nonneg(D if at == 0 else -D, ctx):
            return False
        return replay(P.subs({names[l]: at}), node['p'], ctx)
    if t == 'split':
        l = node['var']
        for side, kd in zip((0, 1), node['kids']):
            args = list(g)
            args[l] = g[l]/2 if side == 0 else (g[l] + 1)/2
            if not replay(P.compose(*args), kd, ctx):
                return False
        return len(node['kids']) == 2
    if t == 'simplex':
        # cube in local coords x = |t - corner| over the 3 variables = simplex {sum x <= 1}
        # union 4 tetrahedra coned from (1,1,1); each piece = Duffy image of the box
        vs = node['vars']; c = {int(j): v for j, v in node['corner'].items()}
        E3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        tets = [((0, 0, 0), E3[0], E3[1], E3[2]), ((1, 1, 1), E3[0], E3[1], E3[2]),
                ((1, 1, 1), (0, 1, 0), (0, 0, 1), (0, 1, 1)),
                ((1, 1, 1), (1, 0, 0), (0, 0, 1), (1, 0, 1)),
                ((1, 1, 1), (1, 0, 0), (0, 1, 0), (1, 1, 0))]
        if len(vs) != 3 or len(node['kids']) != 5 or sorted(c) != sorted(vs):
            return False
        u = [g[v] for v in vs]
        for tet, kd in zip(tets, node['kids']):
            w0 = 1 - u[0]; w1 = u[0]*(1 - u[1]); w2 = u[0]*u[1]*(1 - u[2]); w3 = u[0]*u[1]*u[2]
            args = list(g)
            for i, v in enumerate(vs):
                xv = w0*tet[0][i] + w1*tet[1][i] + w2*tet[2][i] + w3*tet[3][i]
                args[v] = xv if c[v] == 0 else 1 - xv
            if not replay(P.compose(*args), kd, ctx):
                return False
        return True
    if t == 'splitat':
        l = node['var']; r = flint.fmpq(*node['r'])
        if not (0 < r < 1) or len(node['kids']) != 2:
            return False
        for side, kd in zip((0, 1), node['kids']):
            args = list(g)
            args[l] = r*g[l] if side == 0 else r + (1 - r)*g[l]
            if not replay(P.compose(*args), kd, ctx):
                return False
        return True
    if t == 'blowup':
        S = node['S']; c = {int(j): v for j, v in node['corner'].items()}
        if sorted(kd['m'] for kd in node['kids']) != sorted(S):
            return False
        for kd in node['kids']:
            m = kd['m']; args = list(g)
            for j in S:
                x = g[m] if j == m else g[m]*g[j]
                args[j] = x if c[j] == 0 else 1 - x
            if not replay(P.compose(*args), kd['p'], ctx):
                return False
        return True
    return False

# ---------- 6. certificate input validation and exact coverage ----------------------
TUBE_FILES = ['tube_t158.jsonl', 'tube_tubeA.jsonl', 'tube_tubeB.jsonl', 'tube_tubeC.jsonl', 'tube_tubeD.jsonl']
EXPECTED_CANON = 836

def _isint(x):
    return isinstance(x, int) and not isinstance(x, bool)

def _var(x):
    assert _isint(x) and 0 <= x < 6, ('bad variable index', x)

def _bit(x):
    assert _isint(x) and x in (0, 1), ('value not in {0,1}', x)

def validate_face(f):
    assert isinstance(f, str) and len(f) == 12 and set(f) <= set('01*') and f.count('*') == 6, ('bad face key', f)

def validate_tree(n):
    """Type/range checks on a replay() proof tree."""
    assert isinstance(n, dict), ('node not a dict', n)
    t = n.get('t')
    if t in ('zero', 'const', 'bern'):
        return
    if t == 'factor':
        assert isinstance(n['odd'], list)
        for kd in n['odd']:
            assert isinstance(kd['f'], str) and kd['s'] in (1, -1)
            validate_tree(kd['p'])
    elif t == 'mono':
        _var(n['var']); _bit(n['at']); validate_tree(n['p'])
    elif t == 'split':
        _var(n['var']); assert isinstance(n['kids'], list) and len(n['kids']) == 2
        for kd in n['kids']: validate_tree(kd)
    elif t == 'simplex':
        assert isinstance(n['vars'], list) and len(n['vars']) == 3 and len(set(n['vars'])) == 3
        for v in n['vars']: _var(v)
        assert isinstance(n['corner'], dict)
        for j, v in n['corner'].items(): _var(int(j)); _bit(v)
        assert isinstance(n['kids'], list) and len(n['kids']) == 5
        for kd in n['kids']: validate_tree(kd)
    elif t == 'splitat':
        _var(n['var']); r = n['r']
        assert isinstance(r, list) and len(r) == 2 and all(_isint(x) for x in r) and 0 < r[0] < r[1], ('bad r', r)
        assert isinstance(n['kids'], list) and len(n['kids']) == 2
        for kd in n['kids']: validate_tree(kd)
    elif t == 'blowup':
        S = n['S']; assert isinstance(S, list) and S
        for j in S: _var(j)
        assert len(set(S)) == len(S)
        for j, v in n['corner'].items(): _var(int(j)); _bit(v)
        assert sorted(int(j) for j in n['corner']) == sorted(S), 'corner keys must equal S'
        for kd in n['kids']:
            _var(kd['m']); assert kd['m'] in S; validate_tree(kd['p'])
    else:
        raise AssertionError(('unknown node type', t))

def validate_kd(n):
    """Checks on the k-d tree of a tube certificate."""
    assert isinstance(n, dict)
    t = n.get('t')
    if t == 'split':
        _var(n['var']); assert isinstance(n['kids'], list) and len(n['kids']) == 2
        for kd in n['kids']: validate_kd(kd)
    elif t == 'bern':
        pass
    elif t == 'tree':
        validate_tree(n['p'])
    elif t == 'reduce':
        S = n['S']; assert isinstance(S, list) and S
        assert len(set(l for l, v in S)) == len(S)
        for l, v in S: _var(l); _bit(v)
    else:
        raise AssertionError(('unknown k-d node type', t))

def validate_tube(r):
    validate_face(r['face']); validate_kd(r['tree'])
    for key, sub in r['TG'].items():
        S = json.loads(key); assert isinstance(S, list) and S
        for l, v in S: _var(l); _bit(v)
        if isinstance(sub, dict) and sub.get('t') == 'tube':
            validate_tube({'face': r['face'], 'tree': sub['tree'], 'TG': sub['TG']})
        else:
            validate_tree(sub)

def load_certificates(canonical=None, quiet=False):
    """Return (easy, trees, tubes): easy = set of Bernstein faces, trees = {face: tree},
    tubes = {face: record}.  Validates every input, and ASSERTS that the three families are
    pairwise disjoint, that every face key is canonical, and that their union is exactly the
    836 canonical faces."""
    if canonical is None:
        canonical = sorted(set(canon(f) for f in faces_all()))
    can = set(canonical)
    assert len(can) == EXPECTED_CANON, len(can)
    easy = set()
    for line in open('result_0_836.jsonl'):
        r = json.loads(line)
        validate_face(r['face'])
        if r['ok'] is True:
            easy.add(r['face'])
    trees = {}
    for fn in sorted(glob.glob('hard[0-9]_*.jsonl') + glob.glob('kd_*.jsonl') + glob.glob('iso*_*.jsonl')):
        for line in open(fn):
            r = json.loads(line)
            if r.get('ok') and r.get('tree'):
                validate_face(r['face']); validate_tree(r['tree'])
                trees[r['face']] = r['tree']
    tubes = {}
    for fn in TUBE_FILES:
        for line in open(fn):
            r = json.loads(line)
            if r.get('ok'):
                validate_tube(r); tubes[r['face']] = r
    for name, S in (('easy', easy), ('trees', set(trees)), ('tubes', set(tubes))):
        assert S <= can, (name, 'non-canonical faces', sorted(S - can)[:5])
    assert not (easy & set(trees)) and not (easy & set(tubes)) and not (set(trees) & set(tubes)), 'families overlap'
    assert easy | set(trees) | set(tubes) == can, ('coverage gap', sorted(can - (easy | set(trees) | set(tubes)))[:5])
    assert len(easy) + len(trees) + len(tubes) == EXPECTED_CANON
    if not quiet:
        print('certificate inputs valid; coverage exact and disjoint: %d Bernstein + %d trees + %d tubes = %d canonical faces'
              % (len(easy), len(trees), len(tubes), len(can)), flush=True)
    return easy, trees, tubes

def main(only_with_trees=False):
    Phi, abc = check_identity()
    from fractions import Fraction as Fr
    tv = [((0,0,0),(1,0,0),(0,1,0),(0,0,1)),((1,1,1),(1,0,0),(0,1,0),(0,0,1)),((1,1,1),(0,1,0),(0,0,1),(0,1,1)),((1,1,1),(1,0,0),(0,0,1),(1,0,1)),((1,1,1),(1,0,0),(0,1,0),(1,1,0))]
    vol = sum(abs(sp.Matrix([[b[i]-a[i] for i in range(3)] for b in (p,q,r)]).det()) for a,p,q,r in tv)
    assert vol == 6  # five tetrahedra of total volume 1 inside the unit cube
    print('simplex split: 5 tetrahedra, total volume 1: OK', flush=True)
    print('identity  -Disc-4(Q+4)^2 = 4(y^2-1)(d^2+(y+2)^2)(d^2+(y-2)^2): OK', flush=True)
    allf = faces_all()
    print('faces with 6 free coordinates, <=2 per row:', len(allf), flush=True)
    can = sorted(set(canon(f) for f in allf))
    print('canonical faces (S4 x complement):', len(can), flush=True)
    easy, trees0, tubes = load_certificates(can)
    trees = dict(trees0)
    for f in easy:
        trees[f] = {'t': 'bern'}
    missing = [f for f in can if f not in trees and f not in tubes]
    print('canonical faces without a certificate:', len(missing), flush=True)
    ctx = flint.fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(6)), 'lex')
    bad = []
    for n, f in enumerate(can):
        if f not in trees:
            continue  # tube faces are replayed by verify_tube.py
        Psym = face_poly(f, Phi, abc)
        P = ctx.from_dict({tuple(m): flint.fmpq(int(sp.Rational(v).p), int(sp.Rational(v).q))
                           for m, v in Psym.terms()})
        ok = replay(P, trees[f], ctx)
        if not ok:
            bad.append(f)
        print(n, f, 'ok' if ok else 'FAIL', flush=True)
    print('checked', len(can) - len(missing) - len(tubes), 'faces (tube faces: verify_tube.py); failures:', len(bad), bad[:10])
    print('ALL OK' if not bad and not missing else 'NOT COMPLETE')
    if bad or missing:
        sys.exit(1)

if __name__ == '__main__':
    main()

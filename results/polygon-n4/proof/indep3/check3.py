"""Third independent checker for the n = 4 polygon-conjecture proof (../../PROOF.md).

Python standard library only: no python-flint, no sympy, no numpy, and no code from the
project.  Inputs: the JSON certificate files in ../ (read-only).

  python check3.py all            # everything, prints a summary and ALL OK / FAILED
  python check3.py algebra        # (1) Phi and identity (*)
  python check3.py faces          # (2) 38016 faces, 836 orbits, coverage, symmetry, N_F
  python check3.py easy|trees|tubes [--limit n] [--only FACE]

Exit status 0 iff every requested check passed.
"""
import itertools, json, os, random, sys, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import poly as PL
from poly import add, sub, mul, pw, scale, neg, var, const, ONE
from bernstein import bern_nonneg, STATS

DATA = os.path.normpath(os.path.join(HERE, '..'))
EASY_FILE = 'result_0_836.jsonl'
TREE_FILES = ['hard3_0_20.jsonl', 'hard3_20_146.jsonl', 'hard5_0_200.jsonl', 'hard5_0_3_long.jsonl',
              'hard6_0_200.jsonl', 'iso2_one158.jsonl', 'iso2_rem2.jsonl', 'iso2_rem3.jsonl',
              'iso2_resA.jsonl', 'iso2_resA2.jsonl', 'iso2_resB.jsonl', 'iso_one32.jsonl',
              'iso_one593.jsonl', 'iso_remA.jsonl', 'iso_remB.jsonl', 'iso_remC.jsonl',
              'kd_one32.txt.jsonl', 'kd_remaining2.txt.jsonl']
TUBE_FILES = ['tube_t158.jsonl', 'tube_tubeA.jsonl', 'tube_tubeB.jsonl', 'tube_tubeC.jsonl',
              'tube_tubeD.jsonl']


def log(*a):
    print(*a, flush=True)


# ===================== (1) the algebraic criterion ====================================

def det(Mx):
    """Determinant of a square matrix of polynomials, Laplace expansion along the first
    row with memoisation on the set of remaining columns (no division)."""
    n = len(Mx)
    memo = {}

    def rec(r, cols):
        if r == n:
            return ONE
        key = (r, cols)
        if key in memo:
            return memo[key]
        tot = {}
        sgn = 1
        for idx, cidx in enumerate(cols):
            e = Mx[r][cidx]
            if e:
                minor = rec(r + 1, cols[:idx] + cols[idx + 1:])
                tot = add(tot, scale(mul(e, minor), sgn))
            sgn = -sgn
        memo[key] = tot
        return tot
    return rec(0, tuple(range(n)))


def sylvester(f, g):
    """f, g: coefficient lists, highest degree first (entries are polynomials)."""
    m, n = len(f) - 1, len(g) - 1
    rows = []
    for i in range(n):
        rows.append([{}] * i + list(f) + [{}] * (n - 1 - i))
    for i in range(m):
        rows.append([{}] * i + list(g) + [{}] * (m - 1 - i))
    return rows


def reduce_i(P, li):
    """Reduce modulo i^2 = -1, where i is variable li."""
    R = {}
    s = PL.W * li
    for k, c in P.items():
        e = (k >> s) & PL.FM
        k2 = k - (e << s) + ((e % 2) << s)
        c2 = c if e % 4 in (0, 1) else -c
        R[k2] = R.get(k2, 0) + c2
    return {k: v for k, v in R.items() if v}


def phi_from_text():
    """The expansion printed in PROOF.md, typed in by hand (variables a,b,c = 0,1,2)."""
    a, b, c = var(0), var(1), var(2)
    terms = [(-4, [4, 0, 0]), (4, [3, 0, 1]), (-1, [2, 2, 0]), (24, [2, 1, 0]), (-32, [2, 0, 0]),
             (-18, [1, 1, 1]), (4, [0, 3, 0]), (-36, [0, 2, 0]), (96, [0, 1, 0]), (27, [0, 0, 2]),
             (-64, [0, 0, 0])]
    return {PL.mono(e + [0, 0, 0]): cf for cf, e in terms}


def star_identity(Phi):
    """Identity (*): with roots r, x +- iy, i.e. a = r+2x, b = 2rx+x^2+y^2, c = r(x^2+y^2),
    Phi = 4 (y^2-1) (d^2+(y+2)^2) (d^2+(y-2)^2), d = r - x."""
    r, x, y = var(0), var(1), var(2)
    x2y2 = add(mul(x, x), mul(y, y))
    A = add(r, scale(x, 2))
    B = add(scale(mul(r, x), 2), x2y2)
    C = mul(r, x2y2)
    lhs = PL.subst(Phi, {0: A, 1: B, 2: C})
    d = sub(r, x)
    d2 = mul(d, d)
    rhs = scale(mul(mul(sub(mul(y, y), ONE), add(d2, pw(add(y, const(2)), 2))),
                    add(d2, pw(sub(y, const(2)), 2))), 4)
    return lhs == rhs


def check_algebra():
    a, b, c, t, I = var(0), var(1), var(2), var(3), var(4)
    # p(t) = t^3 - a t^2 + b t - c  as coefficient list in t
    pc = [ONE, neg(a), b, neg(c)]
    dpc = [const(3), scale(a, -2), b]
    res = det(sylvester(pc, dpc))
    disc = neg(res)                      # Disc = (-1)^{n(n-1)/2} Res(p,p') / lc, n = 3
    hand = {}
    for cf, e in [(18, [1, 1, 1]), (-4, [3, 0, 1]), (1, [2, 2, 0]), (-4, [0, 3, 0]), (-27, [0, 0, 2])]:
        hand[PL.mono(e + [0, 0, 0])] = cf
    assert disc == hand, 'discriminant of the cubic'
    Q = add(sub(mul(a, a), scale(b, 3)), const(4))
    Phi = sub(neg(disc), scale(mul(Q, Q), 4))
    assert Phi == phi_from_text(), 'Phi expansion in PROOF.md'
    # Phi = Res_t(p(t+i), p(t-i)) / (8i)
    def shifted(sg):
        u = add(t, scale(I, sg))          # t + sg*i
        P = sub(add(sub(pw(u, 3), mul(a, pw(u, 2))), mul(b, u)), c)
        P = reduce_i(P, 4)
        co = []
        for e in range(3, -1, -1):
            co.append({k - (e << (PL.W * 3)): v for k, v in P.items() if PL.ex(k, 3) == e})
        return co
    r2 = reduce_i(det(sylvester(shifted(1), shifted(-1))), 4)
    assert r2 == reduce_i(scale(mul(I, Phi), 8), 4), 'Res(p(t+i),p(t-i)) = 8i Phi'
    assert star_identity(Phi), 'identity (*)'
    # Phi = -Disc - 4(...)^2 <= 0 when all roots real: Disc >= 0 then; nothing to compute.
    return Phi


# ===================== (2) faces and orbits ===========================================

OFF = [(i, j) for i in range(4) for j in range(4) if i != j]   # face string position k
POS = {e: k for k, e in enumerate(OFF)}


def all_faces():
    """6 free coordinates ('*'), at most 2 per row, the other 6 fixed in {0,1}.
    Built row by row (each row: 3 entries)."""
    rows = []
    for i in range(4):
        rows.append([k for k, (a, b) in enumerate(OFF) if a == i])
    out = []
    for frees in itertools.product(*[[c for n in range(3) for c in itertools.combinations(r, n)]
                                     for r in rows]):
        fr = [k for part in frees for k in part]
        if len(fr) != 6:
            continue
        fixed = [k for k in range(12) if k not in fr]
        for bits in itertools.product('01', repeat=6):
            s = ['*'] * 12
            for k, bt in zip(fixed, bits):
                s[k] = bt
            out.append(''.join(s))
    return out


GROUP = [(p, cp) for p in itertools.permutations(range(4)) for cp in (0, 1)]


def act(g, f):
    """Image of face f under relabelling p (entry (i,j) -> (p i, p j)) and complement."""
    p, cp = g
    s = [None] * 12
    for k, (i, j) in enumerate(OFF):
        v = f[k]
        if cp and v != '*':
            v = '1' if v == '0' else '0'
        s[POS[(p[i], p[j])]] = v
    return ''.join(s)


def orbits(faces):
    fs = set(faces)
    orb = {}
    n = 0
    for f in faces:
        if f in orb:
            continue
        for g in GROUP:
            h = act(g, f)
            assert h in fs, 'group does not preserve the face set'
            orb[h] = n
        n += 1
    return orb, n


# ===================== (3) face polynomial N_F = -Phi(a,b,c) ==========================

def face_matrix(face):
    """M = D - W as a 4x4 matrix of polynomials in t0..t5 (free entries in string order)."""
    Wm = [[{} for _ in range(4)] for _ in range(4)]
    c = 0
    for k, (i, j) in enumerate(OFF):
        v = face[k]
        if v == '*':
            Wm[i][j] = var(c)
            c += 1
        elif v == '1':
            Wm[i][j] = ONE
    M = [[neg(Wm[i][j]) for j in range(4)] for i in range(4)]
    for i in range(4):
        rs = {}
        for j in range(4):
            if j != i:
                rs = add(rs, Wm[i][j])
        M[i][i] = rs
    return M


def leibniz(M, idx):
    tot = {}
    for perm in itertools.permutations(range(len(idx))):
        inv = sum(1 for x in range(len(perm)) for y in range(x + 1, len(perm)) if perm[x] > perm[y])
        term = ONE
        for r, cc in enumerate(perm):
            term = mul(term, M[idx[r]][idx[cc]])
            if not term:
                break
        if term:
            tot = add(tot, term if inv % 2 == 0 else neg(term))
    return tot


def face_poly(face, Phi):
    M = face_matrix(face)
    A = {}
    for i in range(4):
        A = add(A, M[i][i])
    B = {}
    for S in itertools.combinations(range(4), 2):
        B = add(B, leibniz(M, S))
    C = {}
    for S in itertools.combinations(range(4), 3):
        C = add(C, leibniz(M, S))
    assert leibniz(M, (0, 1, 2, 3)) == {}, 'det M != 0'
    return neg(PL.subst(Phi, {0: A, 1: B, 2: C}))


def numeric_charpoly_check(face, N, Phi, rng, npts=2):
    """Faddeev-LeVerrier on a rational point (a different route to a, b, c)."""
    for _ in range(npts):
        pt = [Fr(rng.randint(0, 40), 40) for _ in range(6)]
        Wn = [[Fr(0)] * 4 for _ in range(4)]
        c = 0
        for k, (i, j) in enumerate(OFF):
            v = face[k]
            if v == '*':
                Wn[i][j] = pt[c]; c += 1
            elif v == '1':
                Wn[i][j] = Fr(1)
        Mn = [[(-Wn[i][j] if i != j else sum(Wn[i])) for j in range(4)] for i in range(4)]
        n = 4
        Mk = [[Fr(0)] * n for _ in range(n)]
        ck = [Fr(1)]
        for k in range(1, n + 1):
            prev = [[sum(Mn[i][l] * Mk[l][j] for l in range(n)) + (ck[-1] if i == j else 0)
                     for j in range(n)] for i in range(n)]
            Mk = prev
            AM = [[sum(Mn[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
            ck.append(-sum(AM[i][i] for i in range(n)) / k)
        assert ck[4] == 0
        a, b, cc = -ck[1], ck[2], -ck[3]
        ph = PL.evalq(Phi, [a, b, cc, 0, 0, 0])
        assert PL.evalq(N, pt) == -ph, ('N_F mismatch', face)


def symmetry_check(face, g, Phi):
    """N_{g F} composed with the variable correspondence equals N_F (exactly)."""
    p, cp = g
    G = act(g, face)
    vf = {}; c = 0
    for k in range(12):
        if face[k] == '*':
            vf[k] = c; c += 1
    vg = {}; c = 0
    for k in range(12):
        if G[k] == '*':
            vg[k] = c; c += 1
    mp = {}
    for k, (i, j) in enumerate(OFF):
        if face[k] == '*':
            tf = var(vf[k])
            mp[vg[POS[(p[i], p[j])]]] = sub(ONE, tf) if cp else tf
    NG = face_poly(G, Phi)
    NF = face_poly(face, Phi)
    assert PL.subst(NG, mp) == NF, ('symmetry', face, g)


# ===================== (4) certificate replay =========================================

class Bad(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Bad(msg)


def isint(x):
    return isinstance(x, int) and not isinstance(x, bool)


def chkvar(x):
    need(isint(x) and 0 <= x < 6, ('bad variable', x))


def chkbit(x):
    need(isint(x) and x in (0, 1), ('bad 0/1 value', x))


TETS = [((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)), ((1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
        ((1, 1, 1), (0, 1, 0), (0, 0, 1), (0, 1, 1)), ((1, 1, 1), (1, 0, 0), (0, 0, 1), (1, 0, 1)),
        ((1, 1, 1), (1, 0, 0), (0, 1, 0), (1, 1, 0))]


def check_tets():
    """The 5 tetrahedra tile the unit cube: volumes sum to 1, and a lattice of sample points
    (avoiding the shared faces) lies in exactly one of them."""
    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    vol = sum(abs(det3([[q[i] - tt[0][i] for i in range(3)] for q in tt[1:]])) for tt in TETS)
    need(vol == 6, 'tetrahedra volume')
    for pt in itertools.product([Fr(2 * i + 1, 14) for i in range(7)], repeat=3):
        if sum(pt) == 1 or sum(pt) == 2:
            continue
        inside = 0
        for tt in TETS:
            # barycentric coordinates
            m = [[tt[k][i] - tt[0][i] for k in (1, 2, 3)] for i in range(3)]
            D = det3(m)
            rhs = [pt[i] - tt[0][i] for i in range(3)]
            lam = []
            for col in range(3):
                mm = [row[:] for row in m]
                for i in range(3):
                    mm[i][col] = rhs[i]
                lam.append(Fr(det3(mm), D))
            if all(x >= 0 for x in lam) and sum(lam) <= 1:
                inside += 1
        need(inside == 1, ('tetrahedra overlap/gap', pt, inside))


def duffy(tet, i, u):
    """i-th coordinate of the Duffy map of the box (u0,u1,u2) onto tetrahedron tet."""
    one = ONE
    w0 = sub(one, u[0])
    w1 = sub(u[0], mul(u[0], u[1]))
    w2 = sub(mul(u[0], u[1]), mul(mul(u[0], u[1]), u[2]))
    w3 = mul(mul(u[0], u[1]), u[2])
    x = {}
    for wk, vk in zip((w0, w1, w2, w3), tet):
        if vk[i]:
            x = add(x, wk)
    return x


NODES = {}


def replay(P, node):
    """True iff the proof tree `node` certifies P >= 0 on [0,1]^6.  Raises Bad on a
    malformed node.  P may be any positive multiple of the polynomial in question."""
    need(isinstance(node, dict), 'node not a dict')
    t = node.get('t')
    NODES[t] = NODES.get(t, 0) + 1
    if t == 'zero':
        return not P
    if t == 'const':
        return all(k == 0 for k in P) and P.get(0, 0) >= 0
    if t == 'bern':
        return bern_nonneg(P)
    if t == 'mono':
        l, at = node['var'], node['at']
        chkvar(l); chkbit(at)
        D = PL.deriv(P, l)
        if not bern_nonneg(D if at == 0 else neg(D)):
            return False
        return replay(PL.prim(PL.subs01(P, l, at)), node['p'])
    if t == 'split':
        l = node['var']; chkvar(l)
        need(isinstance(node['kids'], list) and len(node['kids']) == 2, 'split kids')
        return (replay(PL.affine(P, l, 0, 1, 2), node['kids'][0]) and
                replay(PL.affine(P, l, 1, 1, 2), node['kids'][1]))
    if t == 'splitat':
        l = node['var']; chkvar(l)
        r = node['r']
        need(isinstance(r, list) and len(r) == 2 and all(isint(x) for x in r), 'splitat r')
        p, q = r
        need(0 < p < q, ('splitat r not in (0,1)', r))
        need(isinstance(node['kids'], list) and len(node['kids']) == 2, 'splitat kids')
        return (replay(PL.affine(P, l, 0, p, q), node['kids'][0]) and
                replay(PL.affine(P, l, p, q - p, q), node['kids'][1]))
    if t == 'blowup':
        S = node['S']
        need(isinstance(S, list) and S and len(set(S)) == len(S), 'blowup S')
        for j in S:
            chkvar(j)
        cor = {}
        for j, v in node['corner'].items():
            cor[int(j)] = v; chkbit(v)
        need(sorted(cor) == sorted(S), 'blowup corner keys != S')
        kids = node['kids']
        need(sorted(kd['m'] for kd in kids) == sorted(S), 'blowup cones != S')
        P1 = P
        for j in S:
            if cor[j] == 1:
                P1 = PL.reflect(P1, j)
        for kd in kids:
            if not replay(PL.prim(PL.blowup_map(P1, kd['m'], S)), kd['p']):
                return False
        return True
    if t == 'simplex':
        vs = node['vars']
        need(isinstance(vs, list) and len(vs) == 3 and len(set(vs)) == 3, 'simplex vars')
        for v in vs:
            chkvar(v)
        cor = {}
        for j, v in node['corner'].items():
            cor[int(j)] = v; chkbit(v)
        need(sorted(cor) == sorted(vs), 'simplex corner keys')
        need(isinstance(node['kids'], list) and len(node['kids']) == 5, 'simplex kids')
        P1 = P
        for v in vs:
            if cor[v] == 1:
                P1 = PL.reflect(P1, v)
        u = [var(v) for v in vs]
        for tet, kd in zip(TETS, node['kids']):
            mp = {vs[i]: duffy(tet, i, u) for i in range(3)}
            if not replay(PL.prim(PL.subst(P1, mp)), kd):
                return False
        return True
    if t == 'factor':
        odd = node['odd']
        need(isinstance(odd, list), 'factor odd')
        R = P
        sign = 1
        for kd in odd:
            need(isinstance(kd['f'], str) and kd['s'] in (1, -1), 'factor entry')
            f = PL.parse(kd['f'])
            need(bool(f), 'zero factor')
            f = PL.prim(f)                       # positive multiple
            e = 0
            while True:
                Q = PL.exact_div(R, f)
                if Q is None:
                    break
                R = Q; e += 1
                need(e <= 64, 'factor multiplicity')
            if e == 0:
                return False                     # listed factor does not divide
            if e % 2:
                if not replay(scale(f, kd['s']), kd['p']):
                    return False
                sign *= kd['s']
        # what is left must be (rational constant) * square
        R = PL.prim(R)
        if R[max(R)] < 0:
            R = neg(R); sign = -sign
        if PL.sqrt_poly(R) is None:
            return False
        return sign > 0
    raise Bad(('unknown node type', t))


def check_tube(N, tree, TG):
    """Tube certificate (see PROOF.md Step 4b).  N: face polynomial (or a restriction)."""
    need(isinstance(TG, dict), 'TG')
    tg = {}
    for key, subc in TG.items():
        S = json.loads(key)
        need(isinstance(S, list) and S, 'TG key')
        for pin in S:
            need(isinstance(pin, list) and len(pin) == 2, 'pin'); chkvar(pin[0]); chkbit(pin[1])
        need(len(set(l for l, v in S)) == len(S), 'TG key repeats a variable')
        fz = frozenset((l, v) for l, v in S)
        need(fz not in tg, 'duplicate TG key')
        PG = N
        for l, v in S:
            PG = PL.subs01(PG, l, v)
        PG = PL.prim(PG)
        if isinstance(subc, dict) and subc.get('t') == 'tube':
            ok = check_tube(PG, subc['tree'], subc['TG'])
        else:
            ok = replay(PG, subc)
        if not ok:
            log('   TG failed', key)
            return False
        tg[fz] = True
    vol = [Fr(0)]

    def walk(node, P, lo, hi):
        need(isinstance(node, dict), 'kd node')
        t = node.get('t')
        if t == 'split':
            l = node['var']; chkvar(l)
            need(isinstance(node['kids'], list) and len(node['kids']) == 2, 'kd split kids')
            mid = (lo[l] + hi[l]) / 2
            return (walk(node['kids'][0], PL.affine(P, l, 0, 1, 2), lo, hi[:l] + [mid] + hi[l + 1:]) and
                    walk(node['kids'][1], PL.affine(P, l, 1, 1, 2), lo[:l] + [mid] + lo[l + 1:], hi))
        v = Fr(1)
        for a, b in zip(lo, hi):
            v *= b - a
        vol[0] += v
        if t == 'bern':
            return bern_nonneg(P)
        if t == 'tree':
            return replay(P, node['p'])
        if t == 'reduce':
            S = node['S']
            need(isinstance(S, list) and S, 'reduce S')
            for pin in S:
                need(isinstance(pin, list) and len(pin) == 2, 'pin'); chkvar(pin[0]); chkbit(pin[1])
            need(len(set(l for l, _ in S)) == len(S), 'reduce repeats a variable')
            if frozenset((l, v) for l, v in S) not in tg:
                return False
            Q = P
            for l, v in S:
                if (v == 0 and lo[l] != 0) or (v == 1 and hi[l] != 1):
                    return False                 # box does not touch the pinned face
                D = PL.deriv(Q, l)
                if not bern_nonneg(D if v == 0 else neg(D)):
                    return False
                Q = PL.prim(PL.subs01(Q, l, v))
            return True
        raise Bad(('unknown k-d node', t))
    ok = walk(tree, N, [Fr(0)] * 6, [Fr(1)] * 6)
    need(vol[0] == 1, ('k-d leaves do not have total volume 1', vol[0]))
    return ok


# ===================== loading and coverage ===========================================

def jl(fn):
    with open(os.path.join(DATA, fn)) as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def chkface(f):
    need(isinstance(f, str) and len(f) == 12 and set(f) <= set('01*') and f.count('*') == 6, ('bad face', f))


def load():
    easy = []
    for r in jl(EASY_FILE):
        chkface(r['face'])
        if r['ok'] is True:
            easy.append(r['face'])
    trees = {}
    for fn in TREE_FILES:
        for r in jl(fn):
            if r.get('ok') and r.get('tree'):
                chkface(r['face'])
                trees.setdefault(r['face'], []).append((fn, r['tree']))
    tubes = {}
    for fn in TUBE_FILES:
        for r in jl(fn):
            if r.get('ok'):
                chkface(r['face'])
                tubes.setdefault(r['face'], []).append((fn, r))
    return easy, trees, tubes


def coverage(orb, norb, easy, trees, tubes):
    need(len(set(easy)) == len(easy), 'duplicate easy face')
    E, T, U = set(easy), set(trees), set(tubes)
    need(not (E & T) and not (E & U) and not (T & U), 'certificate families overlap')
    hit = {}
    for f in E | T | U:
        need(f in orb, ('certified face is not a face', f))
        o = orb[f]
        need(o not in hit, ('two certified faces in one orbit', f, hit.get(o)))
        hit[o] = f
    need(len(hit) == norb, ('orbits without a certificate', norb - len(hit)))
    return len(E), len(T), len(U)


# ===================== driver =========================================================

def main(argv):
    what = argv[0] if argv else 'all'
    limit = None; only = None
    if '--limit' in argv:
        limit = int(argv[argv.index('--limit') + 1])
    if '--only' in argv:
        only = argv[argv.index('--only') + 1]
    t0 = time.time()
    fails = []
    Phi = check_algebra()
    log('(1) Disc via own Sylvester determinant; Phi matches PROOF.md; Phi = Res(p(t+i),p(t-i))/(8i); identity (*): OK')
    check_tets()
    log('    simplex move: 5 tetrahedra tile the cube (volume 1, lattice test): OK')
    faces = all_faces()
    orb, norb = orbits(faces)
    log('(2) faces:', len(faces), ' distinct:', len(set(faces)), ' orbits under S4 x complement:', norb)
    need(len(faces) == 38016 and len(set(faces)) == 38016 and norb == 836, 'face/orbit counts')
    easy, trees, tubes = load()
    nE, nT, nU = coverage(orb, norb, easy, trees, tubes)
    log('(5) coverage: %d Bernstein + %d trees + %d tubes = %d, one certified face in each of the %d orbits, families disjoint: OK'
        % (nE, nT, nU, nE + nT + nU, norb))
    canon_txt = sorted(l.strip() for l in open(os.path.join(DATA, 'faces_canon.txt')) if l.strip())
    log('    (info) certified faces == faces_canon.txt:', sorted(set(easy) | set(trees) | set(tubes)) == canon_txt)
    rng = random.Random(20261028)
    if what in ('all', 'faces'):
        cert = sorted(set(easy) | set(trees) | set(tubes))
        ts = time.time()
        for f in cert:
            numeric_charpoly_check(f, face_poly(f, Phi), Phi, rng)
        log('(3) N_F vs Faddeev-LeVerrier at 2 rational points, all 836 certified faces: OK (%.0fs)' % (time.time() - ts))
        ts = time.time()
        sample = rng.sample(faces, 150)
        for f in sample:
            symmetry_check(f, GROUP[rng.randrange(len(GROUP))], Phi)
        for f in cert[:20]:
            for g in GROUP:
                symmetry_check(f, g, Phi)
        log('    symmetry: N_{gF} o (relabel/complement) == N_F exactly, 150 random (F,g) + 20 faces x all 48 g: OK (%.0fs)'
            % (time.time() - ts))
    def items(d):
        out = sorted(d.items())
        if only:
            out = [x for x in out if x[0] == only]
        return out[:limit] if limit else out
    if what in ('all', 'easy'):
        ts = time.time(); n = 0
        for f in (sorted(easy) if not only else [x for x in easy if x == only])[:limit]:
            N = face_poly(f, Phi)
            ok = bern_nonneg(N)
            n += 1
            if not ok:
                fails.append(('easy', f)); log('FAIL easy', f)
        log('(4) Bernstein faces: %d/%d OK (%.0fs)' % (n - sum(1 for x in fails if x[0] == 'easy'), n, time.time() - ts))
    if what in ('all', 'trees'):
        ts = time.time(); n = 0; nf = 0
        for f, lst in items(trees):
            N = face_poly(f, Phi)
            for fn, tr in lst:
                t1 = time.time()
                try:
                    ok = replay(N, tr)
                except (Bad, KeyError, TypeError, ValueError, AssertionError) as e:
                    ok = False; log('   malformed:', repr(e)[:200])
                n += 1
                if not ok:
                    nf += 1; fails.append(('tree', f, fn))
                log('tree', f, fn, 'OK' if ok else 'FAIL', '%.1fs' % (time.time() - t1))
        log('(4) proof trees: %d/%d OK (%.0fs)' % (n - nf, n, time.time() - ts))
    if what in ('all', 'tubes'):
        ts = time.time(); n = 0; nf = 0
        for f, lst in items(tubes):
            N = face_poly(f, Phi)
            for fn, r in lst:
                t1 = time.time()
                try:
                    ok = check_tube(N, r['tree'], r['TG'])
                except (Bad, KeyError, TypeError, ValueError, AssertionError) as e:
                    ok = False; log('   malformed:', repr(e)[:200])
                n += 1
                if not ok:
                    nf += 1; fails.append(('tube', f, fn))
                log('tube', f, fn, 'OK' if ok else 'FAIL', '%.1fs' % (time.time() - t1))
        log('(4) tube certificates: %d/%d OK (%.0fs)' % (n - nf, n, time.time() - ts))
    log('node types replayed:', dict(sorted(NODES.items())), ' bern calls:', STATS['calls'], ' Bernstein boxes:', STATS['boxes'])
    log('total %.0fs' % (time.time() - t0))
    log('ALL OK' if not fails else 'FAILED: %s' % fails[:20])
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

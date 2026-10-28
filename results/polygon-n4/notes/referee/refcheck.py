"""Referee's independent checker (does not import verify.py or any project code).
Usage: refcheck.py PROOFDIR [easy|trees|tubes|all]
"""
import sys, json, glob, random, time, itertools, math
from math import comb
import flint

PD = sys.argv[1]; WHAT = sys.argv[2] if len(sys.argv) > 2 else 'all'
NAMES = tuple(f't{i}' for i in range(6))
CTX = flint.fmpq_mpoly_ctx.get(NAMES, 'lex')
G = CTX.gens()
ZERO = CTX.from_dict({})
ONE = CTX.from_dict({(0,)*6: 1})
OFF = [(i, j) for i in range(4) for j in range(4) if i != j]

# ---- face polynomial: own a,b,c via explicit minors, hand-expanded Phi from PROOF.md ----
def face_N(face):
    W = [[ZERO for _ in range(4)] for _ in range(4)]
    k = 0
    for (i, j), ch in zip(OFF, face):
        if ch == '*': W[i][j] = G[k]; k += 1
        elif ch == '1': W[i][j] = ONE
    assert k == 6
    M = [[None]*4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            M[i][j] = sum((W[i][q] for q in range(4) if q != i), ZERO) if i == j else -W[i][j]
    a = M[0][0] + M[1][1] + M[2][2] + M[3][3]
    b = sum((M[i][i]*M[j][j] - M[i][j]*M[j][i] for i, j in itertools.combinations(range(4), 2)), ZERO)
    def det3(S):
        p, q, r = S
        m = lambda x, y: M[x][y]
        return (m(p,p)*(m(q,q)*m(r,r)-m(q,r)*m(r,q)) - m(p,q)*(m(q,p)*m(r,r)-m(q,r)*m(r,p))
                + m(p,r)*(m(q,p)*m(r,q)-m(q,q)*m(r,p)))
    c = sum((det3(S) for S in itertools.combinations(range(4), 3)), ZERO)
    Phi = (-4*a**4 + 4*a**3*c - a**2*b**2 + 24*a**2*b - 32*a**2 - 18*a*b*c + 4*b**3
           - 36*b**2 + 96*b + 27*c**2 - 64)
    return -Phi

# ---- Bernstein sign test via (1+x)^d p(x/(1+x)) = sum_i C(d,i) b_i x^i ----
def to_int_terms(P):
    d = P.to_dict()
    L = 1
    for v in d.values():
        q = int(v.q); L = L*q//math.gcd(L, q)
    return {tuple(int(e) for e in m): int(v*L) for m, v in d.items()}

def bern_info(P):
    T = to_int_terms(P)
    if not T: return True, None
    deg = [max(m[l] for m in T) for l in range(6)]
    for l in range(6):
        d = deg[l]
        if d == 0: continue
        new = {}
        for m, cf in T.items():
            j = m[l]
            for kk in range(d - j + 1):
                key = m[:l] + (j+kk,) + m[l+1:]
                new[key] = new.get(key, 0) + cf*comb(d-j, kk)
        T = new
    full = 1
    for d in deg: full *= d+1
    vals = list(T.values())
    if len(T) < full: vals.append(0)
    return min(vals) >= 0, (T, deg)

def bern(P, depth=0, budget=None, maxdepth=16):
    if budget is None: budget = [0]
    budget[0] += 1
    if budget[0] > 20000: return False
    ok, info = bern_info(P)
    if ok: return True
    if depth >= maxdepth: return False
    T, deg = info
    act = [l for l in range(6) if deg[l] > 0]
    if not act: return False
    def spread(l):  # compare normalised Bernstein coefficients b_i = T_i / C(d,i)
        s = 0
        for m, cf in T.items():
            if m[l] < deg[l]:
                m2 = m[:l] + (m[l]+1,) + m[l+1:]
                s = max(s, abs(T.get(m2, 0)/comb(deg[l], m[l]+1) - cf/comb(deg[l], m[l])))
        return s
    l = max(act, key=spread)
    for lo in (0, 1):
        args = list(G); args[l] = (G[l] + lo)/2
        if not bern(P.compose(*args), depth+1, budget, maxdepth): return False
    return True

# ---- proof-tree replay, own implementation of the node semantics in PROOF.md ----
TETS = [((0,0,0),(1,0,0),(0,1,0),(0,0,1)), ((1,1,1),(1,0,0),(0,1,0),(0,0,1)),
        ((1,1,1),(0,1,0),(0,0,1),(0,1,1)), ((1,1,1),(1,0,0),(0,0,1),(1,0,1)),
        ((1,1,1),(1,0,0),(0,1,0),(1,1,0))]
STATS = {}
def rp(P, n):
    t = n['t']; STATS[t] = STATS.get(t, 0) + 1
    if t == 'zero': return P.is_zero()
    if t == 'const':
        d = P.to_dict(); return set(d) <= {(0,)*6} and d.get((0,)*6, 0) >= 0
    if t == 'bern': return bern(P)
    if t == 'split':
        l = n['var']; assert 0 <= l < 6 and len(n['kids']) == 2
        return all(rp(P.compose(*[(G[i] + s)/2 if i == l else G[i] for i in range(6)]), kd)
                   for s, kd in zip((0, 1), n['kids']))
    if t == 'splitat':
        l = n['var']; r = flint.fmpq(*n['r']); assert 0 < r < 1 and len(n['kids']) == 2 and 0 <= l < 6
        m0 = [r*G[i] if i == l else G[i] for i in range(6)]
        m1 = [r + (1-r)*G[i] if i == l else G[i] for i in range(6)]
        return rp(P.compose(*m0), n['kids'][0]) and rp(P.compose(*m1), n['kids'][1])
    if t == 'mono':
        l = n['var']; at = n['at']; assert at in (0, 1) and 0 <= l < 6
        D = P.derivative(NAMES[l])
        if not bern(D if at == 0 else -D): return False
        return rp(P.subs({NAMES[l]: at}), n['p'])
    if t == 'factor':
        c, facs = P.factor()
        prod = ONE*c
        for f, e in facs: prod *= f**e
        assert prod == P
        odd = {kd['f']: kd for kd in n['odd']}
        sgn = 1 if c > 0 else (-1 if c < 0 else 0)
        for f, e in facs:
            if e % 2:
                kd = odd[str(f)]; s = kd['s']; assert s in (1, -1)
                if not rp(s*f, kd['p']): return False
                sgn *= s
        return sgn > 0
    if t == 'blowup':
        S = list(n['S']); cor = {int(k): v for k, v in n['corner'].items()}
        assert len(set(S)) == len(S) and sorted(kd['m'] for kd in n['kids']) == sorted(S)
        for kd in n['kids']:
            m = kd['m']; args = list(G)
            for j in S:
                x = G[m] if j == m else G[m]*G[j]
                args[j] = x if cor[j] == 0 else 1 - x
            if not rp(P.compose(*args), kd['p']): return False
        return True
    if t == 'simplex':
        vs = n['vars']; cor = {int(k): v for k, v in n['corner'].items()}
        assert len(vs) == 3 and len(set(vs)) == 3 and len(n['kids']) == 5
        u = [G[v] for v in vs]
        bw = [1-u[0], u[0]*(1-u[1]), u[0]*u[1]*(1-u[2]), u[0]*u[1]*u[2]]
        for tet, kd in zip(TETS, n['kids']):
            args = list(G)
            for i, v in enumerate(vs):
                x = sum((bw[q]*tet[q][i] for q in range(4)), ZERO)
                args[v] = x if cor[v] == 0 else 1 - x
            if not rp(P.compose(*args), kd): return False
        return True
    raise ValueError('unknown node ' + str(t))

# ---- tube certificate ----
def tube(N, tree, TG, lvl=0):
    for key, sub in TG.items():
        S = json.loads(key)
        assert all(v in (0, 1) and 0 <= l < 6 for l, v in S)
        PG = N.subs({NAMES[l]: v for l, v in S})
        ok = tube(PG, sub['tree'], sub['TG'], lvl+1) if sub.get('t') == 'tube' else rp(PG, sub)
        if not ok: print('  TG fail', key); return False
    vol = [flint.fmpq(0)]; cnt = {}
    def walk(n, lo, hi):
        if n['t'] == 'split':
            l = n['var']; assert 0 <= l < 6 and len(n['kids']) == 2
            mid = (lo[l] + hi[l])/2
            return (walk(n['kids'][0], lo, hi[:l]+[mid]+hi[l+1:]) and
                    walk(n['kids'][1], lo[:l]+[mid]+lo[l+1:], hi))
        v = flint.fmpq(1)
        for i in range(6): v *= hi[i]-lo[i]
        vol[0] += v; cnt[n['t']] = cnt.get(n['t'], 0) + 1
        P = N.compose(*[lo[i] + (hi[i]-lo[i])*G[i] for i in range(6)])
        if n['t'] == 'bern': return bern(P)
        if n['t'] == 'tree': return rp(P, n['p'])
        if n['t'] == 'reduce':
            S = [tuple(x) for x in n['S']]
            key = json.dumps(sorted([list(x) for x in S]))
            if key not in TG: return False
            Q = P
            for l, val in S:
                assert val in (0, 1)
                if val == 0 and lo[l] != 0: return False
                if val == 1 and hi[l] != 1: return False
                D = Q.derivative(NAMES[l])
                if not bern(D if val == 0 else -D): return False
                Q = Q.subs({NAMES[l]: val})
            return True
        raise ValueError(n['t'])
    ok = walk(tree, [flint.fmpq(0)]*6, [flint.fmpq(1)]*6)
    assert vol[0] == 1, vol[0]
    print('  ' * (lvl+1) + 'leaves', cnt, 'volume', vol[0], flush=True)
    return ok

def main():
    random.seed(20261028)
    easy = [json.loads(l)['face'] for l in open(PD+'/result_0_836.jsonl') if json.loads(l)['ok'] is True]
    trees = {}
    for fn in sorted(glob.glob(PD+'/hard[0-9]_*.jsonl') + glob.glob(PD+'/kd_*.jsonl') + glob.glob(PD+'/iso*_*.jsonl')):
        for l in open(fn):
            r = json.loads(l)
            if r.get('ok') and r.get('tree'): trees[r['face']] = r['tree']
    if WHAT in ('easy', 'all'):
        for f in random.sample(sorted(easy), 20):
            t0 = time.time(); ok = bern(face_N(f)); print('easy', f, ok, round(time.time()-t0, 1), flush=True)
    if WHAT in ('trees', 'all'):
        slow = ['**0**0**1110', '**0**0*011*0', '**0*0*0*0*01', '**0**0*0*110']
        rest = random.sample(sorted(set(trees) - set(slow)), 6)
        for f in slow + rest:
            t0 = time.time(); STATS.clear(); ok = rp(face_N(f), trees[f])
            print('tree', f, ok, round(time.time()-t0, 1), dict(STATS), flush=True)
    if WHAT in ('tubes', 'all'):
        for fn in ['tube_t158.jsonl', 'tube_tubeA.jsonl', 'tube_tubeB.jsonl', 'tube_tubeC.jsonl', 'tube_tubeD.jsonl']:
            for l in open(PD + '/' + fn):
                r = json.loads(l)
                if not r.get('ok'): continue
                t0 = time.time(); ok = tube(face_N(r['face']), r['tree'], r['TG'])
                print('tube', fn, r['face'], ok, round(time.time()-t0, 1), flush=True)

main()

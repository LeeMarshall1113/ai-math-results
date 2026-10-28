"""Referee sanity checks: identity, faces/orbits, polynomial vs numpy, numerical search, mutation tests."""
import sys, itertools, json, random
import numpy as np, sympy as sp
sys.argv = [sys.argv[0], sys.argv[1], 'none']
PD = sys.argv[1]
import importlib.util
spec = importlib.util.spec_from_file_location('rc', __file__.replace('sanity.py', 'refcheck.py'))
rc = importlib.util.module_from_spec(spec); spec.loader.exec_module(rc)

# 1. identity: Phi formula, discriminant of cubic by the textbook formula, factorisation
a, b, c, r, x, y, t = sp.symbols('a b c r x y t')
disc = 18*a*b*c - 4*a**3*c + a**2*b**2 - 4*b**3 - 27*c**2      # disc of t^3 - a t^2 + b t - c
assert sp.expand(disc - sp.discriminant(t**3 - a*t**2 + b*t - c, t)) == 0
Phi = -4*a**4 + 4*a**3*c - a**2*b**2 + 24*a**2*b - 32*a**2 - 18*a*b*c + 4*b**3 - 36*b**2 + 96*b + 27*c**2 - 64
assert sp.expand(Phi - (-disc - 4*(a**2 - 3*b + 4)**2)) == 0
sub = {a: r + 2*x, b: 2*r*x + x**2 + y**2, c: r*(x**2 + y**2)}
d = r - x
assert sp.expand(Phi.subs(sub, simultaneous=True) - 4*(y**2-1)*(d**2+(y+2)**2)*(d**2+(y-2)**2)) == 0
p = t**3 - a*t**2 + b*t - c
res = sp.resultant(sp.expand(p.subs(t, t + sp.I)), sp.expand(p.subs(t, t - sp.I)), t)
print('resultant/(8i) == Phi:', sp.expand(res/(8*sp.I) - Phi) == 0)
print('identity OK')

# 2. faces and orbits, independent enumeration via union-find
OFF = rc.OFF; IDX = {e: k for k, e in enumerate(OFF)}
faces = []
for free in itertools.combinations(range(12), 6):
    if any(sum(1 for k in free if OFF[k][0] == i) > 2 for i in range(4)): continue
    rest = [k for k in range(12) if k not in free]
    for bits in itertools.product('01', repeat=6):
        s = ['*']*12
        for k, bt in zip(rest, bits): s[k] = bt
        faces.append(''.join(s))
print('faces', len(faces), len(set(faces)))
fs = set(faces); parent = {f: f for f in faces}
def find(u):
    while parent[u] != u: parent[u] = parent[parent[u]]; u = parent[u]
    return u
for f in faces:
    for perm in itertools.permutations(range(4)):
        for comp in (0, 1):
            g = [None]*12
            for k, (i, j) in enumerate(OFF):
                v = f[k]
                if comp and v != '*': v = '10'[int(v)]
                g[IDX[(perm[i], perm[j])]] = v
            g = ''.join(g); assert g in fs
            parent[find(g)] = find(f)
orbits = {}
for f in faces: orbits.setdefault(find(f), []).append(f)
canon = set(l.strip() for l in open(PD + '/faces_canon.txt') if l.strip())
print('orbits', len(orbits), 'canon file', len(canon),
      'one canon per orbit:', all(sum(1 for f in o if f in canon) == 1 for o in orbits.values()))

# 3. face polynomial vs numpy characteristic polynomial and eigenvalues
rng = np.random.default_rng(1)
worst = 0
for f in random.Random(2).sample(sorted(canon), 60):
    N = rc.face_N(f)
    for _ in range(5):
        tv = rng.random(6)
        W = np.zeros((4, 4)); k = 0
        for (i, j), ch in zip(OFF, f):
            if ch == '*': W[i, j] = tv[k]; k += 1
            elif ch == '1': W[i, j] = 1
        M = np.diag(W.sum(1)) - W
        cp = np.poly(M)  # t^4 + c1 t^3 + ...
        A, B, C = -cp[1], cp[2], -cp[3]
        PhiNum = (-4*A**4 + 4*A**3*C - A**2*B**2 + 24*A**2*B - 32*A**2 - 18*A*B*C + 4*B**3
                  - 36*B**2 + 96*B + 27*C**2 - 64)
        Nv = float(N(*[rc.flint.fmpq(int(round(z*2**20)), 2**20) for z in tv]))
        # recompute at the rounded point
        tv2 = np.array([round(z*2**20)/2**20 for z in tv]); W = np.zeros((4, 4)); k = 0
        for (i, j), ch in zip(OFF, f):
            if ch == '*': W[i, j] = tv2[k]; k += 1
            elif ch == '1': W[i, j] = 1
        M = np.diag(W.sum(1)) - W; cp = np.poly(M); A, B, C = -cp[1], cp[2], -cp[3]
        PhiNum = (-4*A**4 + 4*A**3*C - A**2*B**2 + 24*A**2*B - 32*A**2 - 18*A*B*C + 4*B**3
                  - 36*B**2 + 96*B + 27*C**2 - 64)
        worst = max(worst, abs(Nv + PhiNum) / (1 + abs(Nv)))
        ev = np.linalg.eigvals(M); ym = max(abs(ev.imag))
        if ym > 1 + 1e-9: print('VIOLATION', f, tv2, ev)
        if (ym > 1 + 1e-9) != (PhiNum > 1e-7) and ym > 1e-6 and abs(ym - 1) > 1e-6:
            print('criterion mismatch', f, ym, PhiNum)
print('max rel diff N_F vs numpy charpoly:', worst)

# 4. numerical attack: maximise max|Im eig| over [0,1]^12 (projected random ascent)
def maxim(w):
    W = np.zeros((4, 4))
    for k, (i, j) in enumerate(OFF): W[i, j] = w[k]
    return max(abs(np.linalg.eigvals(np.diag(W.sum(1)) - W).imag))
best = 0; rs = np.random.default_rng(7)
for s in range(400):
    w = rs.random(12) if s % 2 else (rs.random(12) < 0.5).astype(float)
    v = maxim(w); step = 0.2
    for it in range(1500):
        w2 = np.clip(w + step*rs.normal(size=12)*(rs.random(12) < 0.4), 0, 1)
        v2 = maxim(w2)
        if v2 >= v: w, v = w2, v2
        else: step *= 0.997
        if step < 1e-6: break
    best = max(best, v)
print('numerical max |Im mu| over 400 runs:', best)

# 5. mutation tests: the checker must reject wrong claims
f = '**0**0**1110'  # tight face
N = rc.face_N(f)
print('bern(N) on tight face (expect False):', rc.bern(N, maxdepth=6))
print('bern(N - 1/10^6) (expect False):', rc.bern(N - rc.flint.fmpq(1, 10**6), maxdepth=8))
trees = {}
import glob
for fn in sorted(glob.glob(PD+'/hard[0-9]_*.jsonl') + glob.glob(PD+'/kd_*.jsonl') + glob.glob(PD+'/iso*_*.jsonl')):
    for l in open(fn):
        r_ = json.loads(l)
        if r_.get('ok') and r_.get('tree'): trees[r_['face']] = r_['tree']
small = '**0**0*0111*'
tr = trees[small]
print('tree replay (expect True):', rc.rp(rc.face_N(small), tr))
print('tree replay on N - 1/10^6 (expect False):', rc.rp(rc.face_N(small) - rc.flint.fmpq(1, 10**6), tr))
def flip_mono(n):
    if isinstance(n, dict):
        if n.get('t') == 'mono': n['at'] = 1 - n['at']; return True
        return any(flip_mono(v) for v in n.values())
    if isinstance(n, list): return any(flip_mono(v) for v in n)
    return False
import copy
for f2, tr2 in trees.items():
    tt = copy.deepcopy(tr2)
    if flip_mono(tt):
        print('flipped mono in', f2, '(expect False):', rc.rp(rc.face_N(f2), tt)); break

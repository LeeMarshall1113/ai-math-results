"""Self-tests of the replay checker: agreement fast vs exact, 8- vs 64-vertex Q_P,
tamper detection, tiling detection, and a hand-made sanity case."""
import gzip, json, random, copy, sys, os
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hcheck as H

FC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'newmethod', 'fullcert')
random.seed(11)


def load(ci, nrec=3):
    out = []
    with gzip.open(os.path.join(FC, 'chunk_%05d.jsonl.gz' % ci), 'rt') as f:
        for _ in range(nrec):
            out.append(json.loads(next(f)))
    return out


def fast1(L):
    lo, hi = H.leaf_box(L)
    LO = np.array([[int(x * H.G) for x in lo]]); HI = np.array([[int(x * H.G) for x in hi]])
    return bool(H.fast_plain_batch(LO, HI, np.array(L['c']['q'][:39]).reshape(1, 13, 3))[0])


leaves = []
for ci in (0, 37, 77, 120, 152, 175, 190):
    for r in load(ci, 2):
        leaves += [L for L in r['leaves'] if L['k'] == 'p']
sample = random.sample(leaves, 60)
agree = 0; v64 = 0
for L in sample:
    e, _ = H.exact_leaf(L); f = fast1(L); e64, _ = H.exact_leaf(L, use64=True)
    assert e == e64, 'Q_P 8 vs 64 disagree'
    v64 += 1
    if f:
        assert e, 'FAST ACCEPTED BUT EXACT REJECTED'
    agree += (e == f)
print('fast/exact agree on %d/%d genuine leaves; 8-vs-64 vertex Q_P identical verdict on %d' % (agree, len(sample), v64))

# tamper: shift one translate by 1e-3 .. 5e-2 in one coordinate
rej_e = rej_f = 0; N = 0
for L in sample[:40]:
    T = copy.deepcopy(L)
    m = random.randrange(39)
    T['c']['q'][m] += random.choice((-1, 1)) * random.choice((10**6, 10**7, 5 * 10**7))
    e, _ = H.exact_leaf(T); f = fast1(T)
    assert not (f and not e)
    N += 1; rej_e += (not e); rej_f += (not f)
print('tampered translates: exact rejects %d/%d, fast rejects %d/%d (a shift need not break a cover)' % (rej_e, N, rej_f, N))

# tamper: enlarge leaf box (claims cover for a bigger box) -> should usually fail
rej = 0
for L in sample[:20]:
    T = copy.deepcopy(L)
    i = random.randrange(6)
    T['d'][i] = str(Fr(T['d'][i]) * 4)
    e, _ = H.exact_leaf(T); rej += (not e)
print('leaf box enlarged 4x in one direction: exact rejects %d/20' % rej)

# drop a translate to the far away: replaces one of the 13 bodies
rej = 0
for L in sample[:20]:
    T = copy.deepcopy(L)
    m = random.randrange(13)
    T['c']['q'][3 * m:3 * m + 3] = [10**9 * 1, 10**9 * 1, -10**9]
    e, _ = H.exact_leaf(T); rej += (not e)
print('one translate removed (moved far away): exact rejects %d/20' % rej)

# tiling: remove a leaf / duplicate a leaf / shift a leaf
r = load(0, 1)[0]
LO = np.array([[int(Fr(x) * 640) for x in L['a']] for L in r['leaves']])
DD = np.array([[int(Fr(x) * 640) for x in L['d']] for L in r['leaves']])
sb = np.array(r['sb']) * 64
assert H.tiling_ok(sb, sb + 64, LO, LO + DD)[0]
print('tiling genuine ok; drop leaf ->', H.tiling_ok(sb, sb + 64, LO[1:], LO[1:] + DD[1:]))
print('tiling dup leaf ->', H.tiling_ok(sb, sb + 64, np.vstack([LO, LO[:1]]), np.vstack([LO + DD, LO[:1] + DD[:1]])))
LO2 = LO.copy(); LO2[3, 0] += 0  # widen leaf 3 to overlap neighbour
HI2 = LO + DD; HI2 = HI2.copy(); HI2[3, np.argmax(DD[3])] += DD[3].max()
print('tiling overlapping leaf ->', H.tiling_ok(sb, sb + 64, LO2, HI2))

# chord tamper: flip the axis of a chord leaf; corrupt one sign
with gzip.open(os.path.join(FC, 'chunk_00196.jsonl.gz'), 'rt') as f:
    ch = []
    for line in f:
        ch += [L for L in json.loads(line)['leaves'] if L['k'] == 'c']
        if len(ch) > 30:
            break
ch = random.sample(ch, 6)
print('genuine chord leaves:', [H.exact_leaf(L)[0] for L in ch])
res = []
for L in ch:
    T = copy.deepcopy(L); T['ax'] = (T['ax'] + 1) % 3
    res.append(H.exact_leaf(T)[0])
print('chord axis changed -> accepted:', res)
res = []
for L in ch:
    T = copy.deepcopy(L); T['minus'] = copy.deepcopy(T['plus'])  # use the plus solution for minus
    res.append(H.exact_leaf(T)[0])
print('minus certificate replaced by plus certificate -> accepted:', res)
res = []
for L in ch:
    T = {'k': 'p', 'c': L['plus'], 'a': L['a'], 'd': L['d']}  # chord leaf treated as plain (no chord)
    res.append(H.exact_leaf(T)[0])
print('chord leaf plus-solution checked WITHOUT chord body -> accepted:', res)

# STRICT mode unit tests (boundary cases)
from fractions import Fraction as Fr
H.STRICT = False
assert H.covers_unit_exact([(Fr(-1), Fr(1, 2)), (Fr(1, 2), Fr(2))])
assert H.covers_unit_exact([(Fr(0), Fr(1))])
H.STRICT = True
assert not H.covers_unit_exact([(Fr(-1), Fr(1, 2)), (Fr(1, 2), Fr(2))])   # 1/2 only on boundaries
assert not H.covers_unit_exact([(Fr(0), Fr(2))])                          # 0 on boundary
assert not H.covers_unit_exact([(Fr(-1), Fr(1))])                         # 1 on boundary
assert H.covers_unit_exact([(Fr(-1), Fr(1, 2)), (Fr(1, 3), Fr(2))])
cube = [((1, 0, 0), Fr(1)), ((-1, 0, 0), Fr(0)), ((0, 1, 0), Fr(1)), ((0, -1, 0), Fr(0)), ((0, 0, 1), Fr(1)), ((0, 0, -1), Fr(0))]
z = (Fr(0),) * 3
assert not H.point_in_exact(cube, Fr(1), z, (Fr(1), Fr(1, 2), Fr(1, 2)))
assert H.point_in_exact(cube, Fr(1), z, (Fr(1, 2),) * 3)
assert H.interval_exact(cube, Fr(1), z, (0, 0, 0), 0) is None              # edge on the boundary of the body
H.STRICT = False
assert H.interval_exact(cube, Fr(1), z, (0, 0, 0), 0) == (Fr(0), Fr(1))
print('strict-mode unit tests passed')

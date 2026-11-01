"""Full-region certification for H3 <= 13 via R1' (fast version).  Resumable, chunked, multi-process.

Region: Prymak's U^D_{2M} with M = 10 (boxes of side 1/20, p1..p3 in [0,1/2], p4..p6 in [0,1],
k2+k4 <= k1+k6+1 <= k3+k5+2), which contains D (Prymak sec. 3).  Start boxes: the 1/10-grid boxes that contain at
least one U^D box.  Each start box is bisected (longest side) until every leaf is certified (plain or chord, exact
integer check against the full Q_P H-rep), or a leaf of side < MINW fails -> FAIL record, run stops.

Usage: python certify2.py OUTDIR [--workers 3] [--chunk 200]
Output: OUTDIR/chunk_XXXXX.jsonl.gz (one per chunk of start boxes, written atomically), OUTDIR/FAIL (if any)."""
import sys, os, json, gzip, time, itertools, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
_v = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_vendor')
if os.path.isdir(_v) and _v not in sys.path: sys.path.insert(0, _v)
_u = r'C:\Users\hackathon\AppData\Roaming\Python\Python314\site-packages'
if os.path.isdir(_u) and _u not in sys.path: sys.path.append(_u)
from fractions import Fraction as F
import numpy as np
import hw
import lpfast as LF
import cover13 as C13
from qpfast import QP_hrep_fast

MINW = F(1, 640)          # finest allowed side
MAXEVALS = 4000           # per start box; exceeding it is reported as FAIL (never hidden)
CHORD_R = 0.1             # chord bodies only for boxes whose centre is within this L_inf distance of c
GOOD = 0.002
NSHORT = 40
H20 = F(1, 20)

def in_UD(k):
    k1, k2, k3, k4, k5, k6 = k
    return k2 + k4 <= k1 + k6 + 1 <= k3 + k5 + 2

def start_boxes():
    """1/10-grid boxes containing at least one U^D_{20} box, in a fixed order."""
    out = []
    for c in itertools.product(range(5), range(5), range(5), range(10), range(10), range(10)):
        ok = any(in_UD([2 * c[i] + e[i] for i in range(6)]) for e in itertools.product((0, 1), repeat=6))
        if ok: out.append(c)
    return out

ALL = []
for f0 in range(6):
    for v0 in LF.face_vertices(f0):
        others = [f for f in range(6) if f != f0]
        for combo in itertools.product(*[LF.face_edges(f) + [None] for f in others]):
            es = [e for e in combo if e is not None]
            if len(set(es)) == len(es): ALL.append((f0, v0, combo))
import random as _r; _r.Random(1).shuffle(ALL)

SHORT = {}

def box_data(a, d):
    Hint = QP_hrep_fast(a, d)
    A = np.array([[float(h[1]), float(h[2]), float(h[3])] for h in Hint]); b = np.array([float(h[0]) for h in Hint])
    nr = np.linalg.norm(A, axis=1); A = A / nr[:, None]; b = b / nr
    from scipy.optimize import linprog
    from scipy.spatial import HalfspaceIntersection
    r = linprog([0, 0, 0, -1], A_ub=np.hstack([-A, np.ones((len(A), 1))]), b_ub=b,
                bounds=[(None, None)] * 3 + [(0, None)], method='highs')
    if r.status != 0 or r.x[3] < 1e-9: return None
    V = HalfspaceIntersection(np.hstack([-A, -b[:, None]]), r.x[:3]).intersections
    tight = np.array([np.sum(np.abs(V @ row + bb) < 1e-9) >= 3 for row, bb in zip(A, b)])
    # dedupe near-identical rows
    Ar, br = A[tight], b[tight]
    key = np.round(np.hstack([Ar, br[:, None]]), 9)
    _, ix = np.unique(key, axis=0, return_index=True)
    return Hint, Ar[ix], br[ix], V

def try_struct(R, st, axis, sigma, Hint, Nf, cf):
    f0, v0, combo = st
    merged = axis is not None
    items, nv, nmu = LF.items_for(R, f0, v0, combo, merged)
    s_vec = LF.chord_dir(axis, sigma) if merged else None
    m, x = LF.float_solve(items, nv, nmu, Nf, cf, s_vec)
    if m is None or m <= 1e-9: return None
    q = LF.round_solution(x, nv, nmu)
    if not LF.exact_ok(items, nv, nmu, q, Hint, s_vec): return None
    return {'st': [f0, v0, list(combo)], 'q': q, 'm': round(float(m), 7)}

def rank(p, V, R, axis, sigma, budget=6000, top=4):
    H = C13.hrep(V); Hm = None
    if axis is not None:
        s = np.array(LF.chord_dir(axis, sigma), float)
        Hm = C13.hrep(np.vstack([(1 - float(LF.LAM)) * V, (1 - float(LF.LAM)) * V + float(LF.LAM * LF.T) * s]))
    fp = [[tuple(float(x) for x in pt) for pt in r] for r in R]
    out = []
    for st in ALL[:budget]:
        f0, v0, combo = st
        others = [f for f in range(6) if f != f0]
        m = C13.solve(p, f0, v0, dict(zip(others, combo)), H, Hm, fp)
        if m is not None and m > 1e-7:
            out.append((m, st))
            if sum(1 for x in out if x[0] >= GOOD) >= top: break
    out.sort(key=lambda x: x[0], reverse=True)  # ties: never compare st (may hold None)
    return [st for m, st in out[:top]]

def attempt(R, p, V, Hint, Nf, cf, axis, sigma, allow_rank):
    key = (axis, sigma)
    sl = SHORT.get(key, [])
    for st in sl[:NSHORT]:
        c = try_struct(R, st, axis, sigma, Hint, Nf, cf)
        if c: SHORT[key] = [st] + [x for x in sl if x != st][:60]; return c
    if not allow_rank: return None
    for st in rank(p, V, R, axis, sigma):
        if st in sl[:NSHORT]: continue
        c = try_struct(R, st, axis, sigma, Hint, Nf, cf)
        if c: SHORT[key] = [st] + [x for x in SHORT.get(key, []) if x != st][:60]; return c
    return None

def certify_box(a, d):
    bd = box_data(a, d)
    if bd is None: return None
    Hint, Nf, cf, V = bd
    R = hw.rects(a, d)
    p = [float(a[i] + d[i] / 2) for i in range(6)]
    # rank (expensive) only on boxes of side <= 1/40, or when the short-list is empty
    small = max(d) <= F(1, 40)
    c = attempt(R, p, V, Hint, Nf, cf, None, 1, small or not SHORT.get((None, 1)))
    if c: return {'k': 'p', 'c': c}
    if max(abs(x - 0.5) for x in p) <= CHORD_R + float(max(d)):
        for axis in range(3):
            g = {}
            for sigma in (1, -1):
                r = attempt(R, p, V, Hint, Nf, cf, axis, sigma, small or not SHORT.get((axis, sigma)))
                if not r: break
                g[sigma] = r
            if len(g) == 2: return {'k': 'c', 'ax': axis, 'plus': g[1], 'minus': g[-1]}
    return None

def do_start(sb):
    """certify one start box; returns (sid, leaves, fails, seconds)."""
    t0 = time.time()
    a0 = [F(sb[i], 10) for i in range(6)]; d0 = [F(1, 10)] * 6
    leaves = []; fails = []
    stack = [(a0, d0)]; evals = 0
    while stack:
        a, d = stack.pop()
        r = certify_box(a, d)
        if r is not None:
            r['a'] = [str(x) for x in a]; r['d'] = [str(x) for x in d]; leaves.append(r); continue
        evals += 1
        if max(d) <= MINW or evals > MAXEVALS:
            fails.append({'a': [str(x) for x in a], 'd': [str(x) for x in d]}); continue
        i = max(range(6), key=lambda k: d[k])
        d2 = list(d); d2[i] = d[i] / 2; a2 = list(a); a2[i] = a[i] + d2[i]
        stack.append((a2, d2)); stack.append((a, d2))
    return sb, leaves, fails, time.time() - t0

def main():
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    args = sys.argv[2:]
    def opt(name, default, typ=int):
        return typ(args[args.index(name) + 1]) if name in args else default
    workers = opt('--workers', 3); chunk = opt('--chunk', 200); maxh = opt('--max-hours', 1e9, float)
    wk = opt('--worker', None); nw = opt('--nworkers', 1)
    SB = start_boxes()
    nchunks = (len(SB) + chunk - 1) // chunk
    if wk is None:
        # driver: manifest + independent worker processes (each command line contains certify2.py)
        with open(os.path.join(out, 'manifest.json'), 'w') as fh:
            json.dump({'n_start': len(SB), 'chunk': chunk, 'n_chunks': nchunks, 'grid': '1/10',
                       'region': 'U^D_20 (M=10)', 'lam': str(LF.LAM), 't': str(LF.T), 'den': LF.DEN,
                       'minw': str(MINW), 'maxevals': MAXEVALS}, fh)
        import subprocess
        ps = [subprocess.Popen([sys.executable, os.path.abspath(__file__), out, '--chunk', str(chunk),
                                '--max-hours', str(maxh), '--worker', str(k), '--nworkers', str(workers)])
              for k in range(workers)]
        for p in ps: p.wait()
        done = sum(1 for ci in range(nchunks) if os.path.exists(os.path.join(out, 'chunk_%05d.jsonl.gz' % ci)))
        fail = os.path.exists(os.path.join(out, 'FAIL'))
        print('driver: chunks done', done, 'of', nchunks, 'FAIL' if fail else '', flush=True)
        with open(os.path.join(out, 'STATUS'), 'w') as fh:
            fh.write('FAIL\n' if fail else ('DONE\n' if done == nchunks else 'PARTIAL %d/%d\n' % (done, nchunks)))
        return
    t_start = time.time()
    todo = [ci for ci in range(nchunks) if ci % nw == wk and not os.path.exists(os.path.join(out, 'chunk_%05d.jsonl.gz' % ci))]
    print('worker', wk, 'todo chunks', len(todo), flush=True)
    for ci in todo:
        if os.path.exists(os.path.join(out, 'FAIL')): print('FAIL present, stopping', flush=True); return
        if time.time() - t_start > maxh * 3600: print('worker', wk, 'time budget reached', flush=True); return
        # Lee's rule (2026-10-25): while the GAMING flag exists, at most 1 CPU thread per lane -> only worker 0 runs
        gflag = os.environ.get('GAMING_FLAG', r'C:\Users\hackathon\.compute-broker\GAMING')
        while wk > 0 and os.path.exists(gflag):
            if time.time() - t_start > maxh * 3600: return
            time.sleep(60)
        t0 = time.time(); recs = []; nl = 0
        for sb in SB[ci * chunk:(ci + 1) * chunk]:
            while wk > 0 and os.path.exists(gflag): time.sleep(60)  # broker 10-28: also pause mid-chunk while GAMING
            sb_, leaves, fails, sec = do_start(sb)
            recs.append({'sb': list(sb), 'leaves': leaves, 'sec': round(sec, 2)}); nl += len(leaves)
            if fails:
                with open(os.path.join(out, 'FAIL'), 'a') as fh:
                    fh.write(json.dumps({'sb': list(sb), 'fails': fails}) + "\n")
                print('FAIL in start box', sb, flush=True); return
        tmp = os.path.join(out, 'chunk_%05d.tmp' % ci)
        with gzip.open(tmp, 'wt') as fh:
            for r in recs: fh.write(json.dumps(r) + "\n")
        os.replace(tmp, os.path.join(out, 'chunk_%05d.jsonl.gz' % ci))
        print('w%d chunk %d leaves %d sec %.0f %s' % (wk, ci, nl, time.time() - t0, time.strftime('%m-%d %H:%M:%S')), flush=True)

if __name__ == "__main__":
    main()

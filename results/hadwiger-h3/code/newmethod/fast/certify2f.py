"""FAST variant of certify2.py (see EFFICIENCY notes in fast/DEPLOY.md).  Same region, same start boxes, same
bisection rule, same exact integer check and same output format; LPs go straight to HiGHS and the rank() LPs are
built vectorised.

Full-region certification for H3 <= 13 via R1' (fast version).  Resumable, chunked, multi-process.

Region: Prymak's U^D_{2M} with M = 10 (boxes of side 1/20, p1..p3 in [0,1/2], p4..p6 in [0,1],
k2+k4 <= k1+k6+1 <= k3+k5+2), which contains D (Prymak sec. 3).  Start boxes: the 1/10-grid boxes that contain at
least one U^D box.  Each start box is bisected (longest side) until every leaf is certified (plain or chord, exact
integer check against the full Q_P H-rep), or a leaf of side < MINW fails -> FAIL record, run stops.

Usage: python certify2.py OUTDIR [--workers 3] [--chunk 200]
Output: OUTDIR/chunk_XXXXX.jsonl.gz (one per chunk of start boxes, written atomically), OUTDIR/FAIL (if any)."""
import sys, os, json, gzip, time, itertools, math
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', '..'))   # hadwiger3d (hw)
sys.path.insert(0, os.path.join(_here, '..'))         # newmethod (lpfast, cover13, qpfast)
sys.path.insert(0, _here)
_v = os.path.join(_here, '..', '_vendor')
if os.path.isdir(_v) and _v not in sys.path: sys.path.insert(0, _v)
_u = r'C:\Users\hackathon\AppData\Roaming\Python\Python314\site-packages'
if os.path.isdir(_u) and _u not in sys.path: sys.path.append(_u)
from fractions import Fraction as F
import numpy as np
import hw
import lpfastf as LF
import cover13 as C13
import cover13f as C13F
import hlp
from qpfastf import QP_hrep_fast

MINW = F(1, 640)          # finest allowed side
MAXEVALS = 4000           # per start box; exceeding it is reported as FAIL (never hidden)
CHORD_R = 0.1             # chord bodies only for boxes whose centre is within this L_inf distance of c
GOOD = 0.002
NSHORT = int(os.environ.get('CERT_NSHORT', 40))          # short-list tries when rank() may follow
NSHORT_BIG = int(os.environ.get('CERT_NSHORT_BIG', NSHORT))  # short-list tries on boxes where rank() is not allowed
# L3 knob (2026-10-29): CERT_CHORD_RANK_BIG=0 -> chord rank() only on boxes of side <= 1/40.  Default (1) is the old
# rule "rank() also when that chord key's short-list is empty", which in the chord zone re-ran three full 6000-structure
# rank() scans on EVERY failing node (no chord leaf had ever been found, so the chord short-lists never filled).
CHORD_RANK_BIG = os.environ.get('CERT_CHORD_RANK_BIG', '1') != '0'
POLICY_TAG = os.environ.get('CERT_POLICY_TAG')        # if set, stored as 'pol' in every start-box record
FALLBACK = os.environ.get('CERT_FALLBACK') == '1'     # a start box that FAILs under the policy is redone with caps 40
                                                      # and the old chord rule before FAIL is declared
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
if os.environ.get('CERT_PRESOLVE') == 'off': hlp._OPT.presolve = 'off'
_C4 = np.array([0., 0., 0., -1.]); _LB4 = np.array([-hlp.INF] * 3 + [0.]); _UB4 = np.array([hlp.INF] * 4)
# optional search order for rank() (does not change what is certified, only which structures are tried first)
ORDER = os.environ.get('CERT_ORDER')
if ORDER:
    _freq = json.load(open(ORDER))
    _key = lambda st: _freq.get(json.dumps([st[0], st[1], list(st[2])]), 0)
    _pos = {st: i for i, st in enumerate(ALL)}
    _top = sorted([st for st in ALL if _key(st) > 0], key=lambda st: (-_key(st), _pos[st]))
    _rest = [st for st in ALL if _key(st) == 0]
    ALL[:] = _top + _rest
# optional: before a full rank(), try the historically used structures directly (cheap: ~50 LPs instead of up to 6000)
FREQ_FIRST = [st for st in ALL[:len(_top)]] if (ORDER and os.environ.get('CERT_FREQ_FIRST') == '1') else []
RANK_BUDGET = int(os.environ.get('CERT_RANK_BUDGET', 6000))

def box_data(a, d):
    Hint = QP_hrep_fast(a, d)
    Hf = np.array(Hint, dtype=float)                   # same values as float(h[i]) (all |entries| < 2**53)
    A = Hf[:, 1:]; b = Hf[:, 0]
    nr = np.linalg.norm(A, axis=1); A = A / nr[:, None]; b = b / nr
    from scipy.spatial import HalfspaceIntersection
    ok, x = hlp.solve_ub(_C4, np.hstack([-A, np.ones((len(A), 1))]), b, _LB4, _UB4)
    if not ok or x[3] < 1e-9: return None
    V = HalfspaceIntersection(np.hstack([-A, -b[:, None]]), x[:3]).intersections
    tight = (np.abs(A @ V.T + b[:, None]) < 1e-9).sum(1) >= 3   # = the per-row test (values are ~1e-15 or >> 1e-9)
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
    C13F.clear_cache()
    out = []
    for st in ALL[:budget]:
        f0, v0, combo = st
        others = [f for f in range(6) if f != f0]
        m = C13F.solve(p, f0, v0, dict(zip(others, combo)), H, Hm, fp)
        if m is not None and m > 1e-7:
            out.append((m, st))
            if sum(1 for x in out if x[0] >= GOOD) >= top: break
    out.sort(key=lambda x: x[0], reverse=True)  # ties: never compare st (may hold None)
    return [st for m, st in out[:top]]

def attempt(R, p, V, Hint, Nf, cf, axis, sigma, allow_rank):
    key = (axis, sigma)
    sl = SHORT.get(key, [])
    for st in sl[:(NSHORT if allow_rank else NSHORT_BIG)]:
        c = try_struct(R, st, axis, sigma, Hint, Nf, cf)
        if c: SHORT[key] = [st] + [x for x in sl if x != st][:60]; return c
    if not allow_rank: return None
    lim = NSHORT if allow_rank else NSHORT_BIG
    if FREQ_FIRST:
        tried = set(sl[:lim])
        for st in FREQ_FIRST:
            if st in tried: continue
            c = try_struct(R, st, axis, sigma, Hint, Nf, cf)
            if c: SHORT[key] = [st] + [x for x in SHORT.get(key, []) if x != st][:60]; return c
    for st in rank(p, V, R, axis, sigma, budget=RANK_BUDGET):
        if FREQ_FIRST and st in FREQ_FIRST: continue
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
                r = attempt(R, p, V, Hint, Nf, cf, axis, sigma,
                            small or (CHORD_RANK_BIG and not SHORT.get((axis, sigma))))
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
    # F1 repair (2026-10-31): CERT_SB_FILE=<json with key 'boxes'> replaces the U^D_20 start-box list (same 1/10 grid,
    # same bisection, same leaf checks, same chunk format).  Used for newmethod/fullcert_ext/ only.
    if os.environ.get('CERT_SB_FILE'):
        SB = [tuple(c) for c in json.load(open(os.environ['CERT_SB_FILE']))['boxes']]
        print('start boxes from', os.environ['CERT_SB_FILE'], len(SB), flush=True)
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
    # --chunks a,b,c: explicit chunk list (overrides the ci % nw == wk stripe); --sbrange lo:hi: only start boxes
    # SB[ci*chunk+lo : ci*chunk+hi] of each listed chunk (a split: several workers can share one chunk).
    # --parts: checkpoint finished start boxes to OUT/parts/ every CERT_CKPT_SEC (default 1200 s) and skip start boxes
    # already there; when all start boxes of a chunk are in parts, the chunk file is written in SB order (the same
    # format certcheck reads) and the parts move to OUT/parts/merged/.  certcheck reads only OUT/*.jsonl.gz.
    chunks_opt = opt('--chunks', None, str); sbr = opt('--sbrange', None, str); use_parts = '--parts' in args or '--queue' in args or sbr is not None
    lo, hi = (int(x) for x in sbr.split(':')) if sbr else (0, chunk)
    pdir = os.path.join(out, 'parts'); ckpt_sec = float(os.environ.get('CERT_CKPT_SEC', 1200))
    if use_parts: os.makedirs(os.path.join(pdir, 'merged'), exist_ok=True)
    final = lambda ci: os.path.join(out, 'chunk_%05d.jsonl.gz' % ci)
    if chunks_opt: todo = [int(x) for x in chunks_opt.split(',') if x and not os.path.exists(final(int(x)))]
    else: todo = [ci for ci in range(nchunks) if ci % nw == wk and not os.path.exists(final(ci))]
    print('worker', wk, 'todo chunks', len(todo), todo[:8], 'sbrange %d:%d' % (lo, hi), 'policy', POLICY_TAG, flush=True)

    def part_files(ci):
        pre = 'chunk_%05d.' % ci
        return sorted(os.path.join(pdir, f) for f in os.listdir(pdir) if f.startswith(pre) and f.endswith('.part.gz'))
    def load_parts(ci):
        got = {}
        if not use_parts: return got
        for fn in part_files(ci):
            try:
                for line in gzip.open(fn, 'rt'):
                    r = json.loads(line); got[tuple(r['sb'])] = r
            except (OSError, EOFError, ValueError) as e: print('unreadable part', fn, e, flush=True)
        return got
    CLAIM = [None]
    def write_part(ci, recs, tag):
        if CLAIM[0] and os.path.exists(CLAIM[0]): os.utime(CLAIM[0])
        if not recs: return
        p = os.path.join(pdir, 'chunk_%05d.%s.%d.%d.part.gz' % (ci, tag, os.getpid(), int(time.time() * 1000)))
        with gzip.open(p + '.tmp', 'wt') as fh:
            for r in recs: fh.write(json.dumps(r) + "\n")
        os.replace(p + '.tmp', p)
    def try_merge(ci):
        """write the chunk file if parts hold every start box of chunk ci (in SB order, as certify2 writes it)."""
        if os.path.exists(final(ci)): return True
        got = load_parts(ci); want = [tuple(s) for s in SB[ci * chunk:(ci + 1) * chunk]]
        if any(s not in got for s in want): return False
        tmp = os.path.join(out, 'chunk_%05d.tmp%d' % (ci, os.getpid()))
        with gzip.open(tmp, 'wt') as fh:
            for s in want: fh.write(json.dumps(got[s]) + "\n")
        if not os.path.exists(final(ci)): os.replace(tmp, final(ci))
        else: os.remove(tmp)
        for fn in part_files(ci):
            os.replace(fn, os.path.join(pdir, 'merged', os.path.basename(fn)))
        return True

    def run_unit(ci, lo, hi):
        if os.path.exists(os.path.join(out, 'FAIL')): print('FAIL present, stopping', flush=True); return 'stop'
        if time.time() - t_start > maxh * 3600: print('worker', wk, 'time budget reached', flush=True); return 'stop'
        # Lee's rule (2026-10-25): while the GAMING flag exists, at most 1 CPU thread per lane -> only worker 0 runs
        gflag = os.environ.get('GAMING_FLAG', r'C:\Users\hackathon\.compute-broker\GAMING')
        while wk > 0 and os.path.exists(gflag):
            if time.time() - t_start > maxh * 3600: return 'stop'
            time.sleep(60)
        if use_parts and try_merge(ci): print('w%d chunk %d complete from parts' % (wk, ci), flush=True); return None
        done = load_parts(ci)
        t0 = time.time(); recs = []; pend = []; nl = 0; t_ck = time.time()
        for sb in SB[ci * chunk + lo:ci * chunk + min(hi, chunk)]:
            if tuple(sb) in done: continue
            while wk > 0 and os.path.exists(gflag): time.sleep(60)  # broker 10-28: also pause mid-chunk while GAMING
            if time.time() - t_start > maxh * 3600 and use_parts:
                write_part(ci, pend, 'w%d' % wk); print('worker', wk, 'time budget reached mid-chunk (parts saved)', flush=True); return 'stop'
            sb_, leaves, fails, sec = do_start(sb); pol = POLICY_TAG
            if fails and FALLBACK:
                g = globals(); keep = (g['NSHORT'], g['NSHORT_BIG'], g['CHORD_RANK_BIG'])
                g['NSHORT'], g['NSHORT_BIG'], g['CHORD_RANK_BIG'] = 40, 40, True
                print('w%d sb %s failed under %s; fallback' % (wk, sb, POLICY_TAG), flush=True)
                try: sb_, leaves, fails, sec2 = do_start(sb)
                finally: g['NSHORT'], g['NSHORT_BIG'], g['CHORD_RANK_BIG'] = keep
                sec += sec2; pol = (POLICY_TAG or '') + '+fallback40'
            rec = {'sb': list(sb), 'leaves': leaves, 'sec': round(sec, 2)}
            if pol: rec['pol'] = pol
            recs.append(rec); pend.append(rec); nl += len(leaves)
            if fails:
                write_part(ci, pend, 'w%d' % wk) if use_parts else None
                with open(os.path.join(out, 'FAIL'), 'a') as fh:
                    fh.write(json.dumps({'sb': list(sb), 'fails': fails}) + "\n")
                print('FAIL in start box', sb, flush=True); return 'stop'
            if use_parts and time.time() - t_ck > ckpt_sec:
                write_part(ci, pend, 'w%d' % wk); pend = []; t_ck = time.time()
        if use_parts:
            write_part(ci, pend, 'w%d' % wk)
            ok = try_merge(ci)
            print('w%d chunk %d sb %d:%d new %d leaves %d sec %.0f %s%s' % (wk, ci, lo, hi, len(recs), nl, time.time() - t0,
                  time.strftime('%m-%d %H:%M:%S'), ' MERGED' if ok else ' (other parts pending)'), flush=True)
            return None
        tmp = os.path.join(out, 'chunk_%05d.tmp' % ci)
        with gzip.open(tmp, 'wt') as fh:
            for r in recs: fh.write(json.dumps(r) + "\n")
        os.replace(tmp, final(ci))
        print('w%d chunk %d leaves %d sec %.0f %s' % (wk, ci, nl, time.time() - t0, time.strftime('%m-%d %H:%M:%S')), flush=True)
        return None

    if '--queue' in args:
        # dynamic queue (2026-10-29): units = (chunk, start-box range); chunks with chord-zone start boxes are cut into
        # pieces of CERT_PIECE (default 25) start boxes, heaviest first.  A unit is claimed by creating
        # parts/claim_<ci>_<lo>_<hi> exclusively; claims untouched for CERT_STALE_H (default 8) hours are retaken.
        piece = int(os.environ.get('CERT_PIECE', 25)); stale = float(os.environ.get('CERT_STALE_H', 8)) * 3600
        cz = lambda c: max(abs(c[i] / 10 + 0.05 - 0.5) for i in range(6)) <= 0.1500001
        units = []
        for ci in range(nchunks):
            if os.path.exists(final(ci)): continue
            s_ = SB[ci * chunk:(ci + 1) * chunk]
            if any(cz(c) for c in s_):
                for l0 in range(0, len(s_), piece):
                    units.append((-sum(cz(c) for c in s_[l0:l0 + piece]), ci, l0, l0 + piece))
            else: units.append((0, ci, 0, chunk))
        units.sort()
        print('worker', wk, 'queue units', len(units), flush=True)
        for _, ci, l0, h0 in units:
            if os.path.exists(final(ci)): continue
            want = [tuple(c) for c in SB[ci * chunk + l0:ci * chunk + min(h0, chunk)]]
            if all(c in load_parts(ci) for c in want): try_merge(ci); continue
            cl = os.path.join(pdir, 'claim_%05d_%03d_%03d' % (ci, l0, h0))
            try:
                fd = os.open(cl, os.O_CREAT | os.O_EXCL | os.O_WRONLY); os.write(fd, b'w%d pid %d' % (wk, os.getpid())); os.close(fd)
            except FileExistsError:
                if time.time() - os.path.getmtime(cl) < stale: continue
                os.utime(cl)
            CLAIM[0] = cl
            if run_unit(ci, l0, h0) == 'stop':
                if not os.path.exists(os.path.join(out, 'FAIL')) and os.path.exists(cl): os.remove(cl)  # release claim
                return
            os.replace(cl, os.path.join(pdir, 'merged', os.path.basename(cl) + '.done'))
            CLAIM[0] = None
        return
    for ci in todo:
        if run_unit(ci, lo, hi) == 'stop': return

if __name__ == "__main__":
    main()

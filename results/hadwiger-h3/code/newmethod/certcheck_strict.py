"""certcheck_strict: certcheck_fast with STRICT inequalities (added 2026-10-31 for the Lean formalisation).

Identical to certcheck_fast.py except: every containment test `c0*w + n.z >= 0` is replaced by `> 0`, i.e. each
covered point must lie in the INTERIOR of its translate of Q_P (resp. (1-lam)Q_P).  The Lean predicate `LeafOK`
(hadwiger3d/lean/Hadwiger/Leaf.lean) requires exactly this; the non-strict check is not enough because Prymak's
Lemma 2.4(ii) (a translate of Q_P inside int O_p) fails at some box corners.  Ledgers: certcheck_strict_done*.jsonl
(separate from certcheck_fast's), keyed by sha256(certcheck_strict.py + certcheck.py + hw.py).
"""
"""certcheck_fast: the same checks as certcheck.py, with the exact rational arithmetic replaced by exact INTEGER arithmetic.

certcheck.py is left untouched and stays the reference.  This module reuses from it only the parts that are not
arithmetic-heavy and carry the meaning of the certificate: constraints() (structure -> constraint list), chord_dir(),
tiles() (bisection-tree coverage), in_UD() and the constants LAM, T, DEN.  It does NOT use qpfast.py / cover13 / lpfast
(the certifier's own code), so it stays independent of the program that produced the certificate.

Exactness argument (why every verdict equals certcheck.py's verdict):

 (A) QP_hrep_int(a, d) == hw.QP_hrep(a, d) as sorted lists of integer 4-tuples.
     Let L be a positive integer with L*x integral for every coordinate x of every point used (the lcm of the
     denominators of all a_i and a_i + d_i; the constants 0 and 1 are integral).  Put p' = L p for each point.
     For a triple (a, b, c): u' = L u, v' = L v, so n' = u' x v' = L^2 n (n is hw's normal), and
     d' = -n'.a' = L^3 d.  For every point p: d' + n'.p' = L^3 (d + n.p), so every sign test (n == 0, all >= 0,
     all <= 0) has the same outcome, and the orientation flip is applied in the same cases.  The original row
     (d, n1, n2, n3) is therefore a POSITIVE rational multiple of (d', L n1', L n2', L n3') (factor L^3).  hw then scales
     its row to integers and divides by the gcd of the absolute values; the result is the unique primitive integer
     vector on that open ray.  We divide our integer vector by its gcd, which gives the primitive vector on the same
     ray, hence the identical tuple.  The set union over the 64 point choices, the extra row (0,1,0,0) and the sort are
     the same, so the two lists are equal.  (gcd > 0 because n != 0 on every row kept.)
 (B) check_solution_int(...) == certcheck.check_solution(...).
     The length test and 0 <= mu <= 1 (i.e. 0 <= q <= DEN, as mu = q/DEN) are the same tests.  Every point z the
     reference forms is z_i = base_i - q_T/DEN + sum_v (q_v/DEN) u_i [- (q_mu/DEN) LAM T s_i], with base_i a box/cube
     coordinate (L base_i integral), u_i, s_i integers and LAM*T = 3983/10000.  With M = L*DEN*10000 the vector
     Z = M z is integral and computed exactly in integers.  The reference tests  c0 * w + n.z >= 0  with w = 1 (plain)
     or w = 1 - LAM = 3/10 (merged).  Multiplying by the positive integer 10*M gives the equivalent integer test
     c0 * (10 w) * M + 10 * (n.Z) >= 0, with 10 w in {10, 3}.  Multiplying by a positive number preserves the sign,
     so each point passes iff it passes in the reference, and the verdict is the same.
 Supplements (independent review, 2026-10-30):
 * L*p is integral because hw.rects only produces coordinates in {0, 1, a_i, a_i + d_i}; the cube vertices used as
   bases in constraints() are in {0, 1}^3.  (hw's own `int(F(x)*den)` normalisation is exact as well.)
 * The reference literally tests  Hc*(w.num*Lz) + (Hn.zi)*w.den >= 0  with zi = Lz*z, Lz > 0 the lcm of z's
   denominators; that is (Lz*w.den) * (c0*w + n.z) >= 0 with a positive factor, i.e. the test used in (B).
 * Non-integer q entries: the reference raises TypeError (Fraction(float, int)); here they are rejected (False).
   bool entries are ints in both.  Neither direction can turn a rejected leaf into an accepted one.
 * Integrality is enforced by explicit raises (not asserts), so `python -O` changes nothing.
 Nothing is rounded anywhere; Python ints are unbounded.
 Driver: a ledger verdict is reused only if the chunk file's sha256 AND the sha256 of (certcheck_fast.py, certcheck.py,
 hw.py) match; coverage is computed from the start boxes check_file actually tiled.

Usage: python certcheck_fast.py OUTDIR [--only 151,153 | --only-file F] [--tag T] [--max-files N] [--verdict-only]
                                  [--xcheck FRAC] [--seed S]
  Ledgers: each process appends to OUTDIR/certcheck_strict_done[.T].jsonl (use a distinct --tag per concurrent process);
  all certcheck_strict_done*.jsonl are read.  A verdict is reused only if its sha256 (of the exact bytes that were
  checked) equals the current file's sha256 AND its checker_sha equals sha256(certcheck_fast.py+certcheck.py+hw.py).
  --verdict-only: check nothing, just print the verdict from the ledgers (NOT COMPLETE if any file lacks a match).
  --xcheck FRAC: also run the ORIGINAL certcheck.check_leaf on a random fraction FRAC of leaves and record any
  disagreement (there should be none)."""
import sys, os, json, gzip, time, itertools, math, random, io, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fractions import Fraction as F
import certcheck as K          # reference: structure, tiling, constants (unchanged)
import hw                       # only hw.rects (point lists); hw.QP_hrep is NOT used here

LAMT_NUM, LAMT_DEN = (K.LAM * K.T).numerator, (K.LAM * K.T).denominator     # 3983 / 10000
W_MERGED = 1 - K.LAM                                                          # 3/10
assert (K.LAM * K.T) == F(LAMT_NUM, LAMT_DEN) and W_MERGED == F(3, 10)

def _lcm(a, b): return a * b // math.gcd(a, b)

def box_scale(a, d):
    L = 1
    for i in range(6):
        L = _lcm(L, a[i].denominator); L = _lcm(L, (a[i] + d[i]).denominator)
    return L

def hrep_points_int(pts):
    """pts: integer 3-vectors. Same algorithm as hw.hrep_points, in integers; returns primitive rows (c0,c1,c2,c3)."""
    out = set()
    for a, b, c in itertools.combinations(pts, 3):
        u0, u1, u2 = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        v0, v1, v2 = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        n0 = u1 * v2 - u2 * v1; n1 = u2 * v0 - u0 * v2; n2 = u0 * v1 - u1 * v0
        if n0 == 0 and n1 == 0 and n2 == 0: continue
        dd = -(n0 * a[0] + n1 * a[1] + n2 * a[2])
        vals = [dd + n0 * p[0] + n1 * p[1] + n2 * p[2] for p in pts]
        if all(x >= 0 for x in vals): pass
        elif all(x <= 0 for x in vals): n0, n1, n2, dd = -n0, -n1, -n2, -dd
        else: continue
        out.add((dd, n0, n1, n2))
    return out

def QP_hrep_int(a, d):
    """== hw.QP_hrep(a, d) (see (A)).  a, d: Fractions."""
    L = box_scale(a, d)
    ss = hw.rects(a, d)
    S = [[tuple((F(x) * L).numerator if (F(x) * L).denominator == 1 else None for x in p) for p in face] for face in ss]
    if any(x is None for face in S for p in face for x in p): raise ValueError('box scale not integral')
    H = set()
    for i0 in range(4):
        for i2 in range(4):
            for i4 in range(4):
                pts = [S[0][i0], S[1][i0], S[2][i2], S[3][i2], S[4][i4], S[5][i4]]
                for dd, n0, n1, n2 in hrep_points_int(pts):
                    row = (dd, L * n0, L * n1, L * n2)        # positive multiple of hw's (d, n) -- see (A)
                    g = math.gcd(math.gcd(abs(row[0]), abs(row[1])), math.gcd(abs(row[2]), abs(row[3])))
                    H.add(tuple(x // g for x in row))
    H.add((0, 1, 0, 0))
    return sorted(H)

def check_solution_int(cons, nv, nmu, q, H, s_vec, L):
    """== certcheck.check_solution(cons, nv, nmu, q, H, s_vec) (see (B)).  L: box scale (L*base integral)."""
    if len(q) != nv + nmu: return False
    if any(not isinstance(x, int) for x in q): return False   # reference raises TypeError on these (bool is an int in both)
    mus = q[nv:]
    if any(m < 0 or m > K.DEN for m in mus): return False
    M = L * K.DEN * LAMT_DEN
    cq = K.DEN * LAMT_DEN            # M / L
    cl = L * LAMT_DEN                # M / DEN
    Hc = [h[0] for h in H]; Hn = [(h[1], h[2], h[3]) for h in H]
    k = 0
    for base, lin, Ti, merged in cons:
        Bi = []
        for i in range(3):
            bL = base[i] * L
            if bL.denominator != 1: raise ValueError('base not integral at scale L')
            Bi.append(bL.numerator)
        Z = [Bi[i] * cq - q[Ti + i] * cl + sum(q[v] * u[i] for v, u in lin) * cl for i in range(3)]
        if merged:
            mq = mus[k]; k += 1
            Z = [Z[i] - mq * LAMT_NUM * L * s_vec[i] for i in range(3)]
            w10 = 3
        else:
            w10 = 10
        cM = w10 * M
        for c0, (a1, a2, a3) in zip(Hc, Hn):
            if c0 * cM + 10 * (a1 * Z[0] + a2 * Z[1] + a3 * Z[2]) <= 0: return False   # STRICT
    return True

def check_leaf_fast(lf):
    a = [F(x) for x in lf['a']]; d = [F(x) for x in lf['d']]
    H = QP_hrep_int(a, d); R = hw.rects(a, d); L = box_scale(a, d)
    if lf['k'] == 'p':
        c = lf['c']; cons, nv, nmu = K.constraints(R, c['st'][0], c['st'][1], tuple(c['st'][2]), False)
        return check_solution_int(cons, nv, nmu, c['q'], H, None, L)
    if lf['k'] == 'c':
        for sg, c in ((1, lf['plus']), (-1, lf['minus'])):
            cons, nv, nmu = K.constraints(R, c['st'][0], c['st'][1], tuple(c['st'][2]), True)
            if not check_solution_int(cons, nv, nmu, c['q'], H, K.chord_dir(lf['ax'], sg), L): return False
        return True
    return False

GFLAG = os.environ.get('GAMING_FLAG', r'C:\Users\hackathon\.compute-broker\GAMING')
PAUSE_ON_GAMING = os.environ.get('CHECK_PAUSE_ON_GAMING', '1') == '1'

def check_file(path, xfrac=0.0, rng=None):
    t0 = time.time(); nsb = nl = nx = 0; bad = []; disagree = []; tiled = []
    data = open(path, 'rb').read(); fsha = hashlib.sha256(data).hexdigest()   # verdict is about exactly these bytes
    for line in gzip.open(io.BytesIO(data), 'rt'):
        while PAUSE_ON_GAMING and os.path.exists(GFLAG): time.sleep(60)   # Lee's GAMING rule, also mid-file
        r = json.loads(line); sb = tuple(r['sb']); nsb += 1
        a0 = [F(sb[i], 10) for i in range(6)]; d0 = [F(1, 10)] * 6
        Lv = [(tuple(F(x) for x in lf['a']), tuple(F(x) for x in lf['d'])) for lf in r['leaves']]
        if not K.tiles(a0, d0, Lv): bad.append(['tiling', list(sb)]); continue
        tiled.append(list(sb))
        for lf in r['leaves']:
            nl += 1
            v = check_leaf_fast(lf)
            if not v: bad.append(['leaf', list(sb), lf['a'], lf['d']])
            if xfrac and rng.random() < xfrac:
                nx += 1
                if K.check_leaf(lf) != v: disagree.append([list(sb), lf['a'], lf['d'], v])
    return {'sbs': nsb, 'leaves': nl, 'bad': len(bad), 'first_bad': bad[:5], 'xchecked': nx,
            'disagree': len(disagree), 'first_disagree': disagree[:5], 'sec': round(time.time() - t0, 1), 'tiled_sbs': tiled, 'sha256': fsha}

if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    out = sys.argv[1]; args = sys.argv[2:]
    maxf = int(args[args.index('--max-files') + 1]) if '--max-files' in args else 10 ** 9
    only = set(int(x) for x in args[args.index('--only') + 1].split(',')) if '--only' in args else None
    if '--only-file' in args: only = set(int(x) for x in open(args[args.index('--only-file') + 1]).read().split(',') if x.strip())
    xfrac = float(args[args.index('--xcheck') + 1]) if '--xcheck' in args else 0.0
    seed = int(args[args.index('--seed') + 1]) if '--seed' in args else 12345
    import hashlib
    def sha(path):
        h = hashlib.sha256()
        with open(path, 'rb') as fh:
            for blk in iter(lambda: fh.read(1 << 20), b''): h.update(blk)
        return h.hexdigest()
    here = os.path.dirname(os.path.abspath(__file__))
    CHK = hashlib.sha256(b''.join(open(os.path.join(here, *p), 'rb').read() for p in
                         (('certcheck_strict.py',), ('certcheck.py',), ('..', 'hw.py')))).hexdigest()
    tag = args[args.index('--tag') + 1] if '--tag' in args else ''
    led = os.path.join(out, 'certcheck_strict_done%s.jsonl' % ('.' + tag if tag else '')); done = {}
    # one ledger per concurrent process (concurrent appends to one file are not safe on Windows); all are read
    for lname in sorted(f for f in os.listdir(out) if f.startswith('certcheck_strict_done') and f.endswith('.jsonl')):
        for line in open(os.path.join(out, lname)):
            try: v = json.loads(line)
            except ValueError: continue            # torn line from a crash: ignore (the file is re-checked)
            if v.get('checker_sha') == CHK: done[v['file']] = v      # later entries win (append order)
    files = sorted(f for f in os.listdir(out) if f.endswith('.jsonl.gz'))
    fsha = {fn: sha(os.path.join(out, fn)) for fn in files}
    n = 0
    for fn in files:
        if only is not None and int(fn[6:11]) not in only: continue
        st = os.stat(os.path.join(out, fn)); v = done.get(fn)
        if v and v.get('sha256') == fsha[fn]: continue
        if n >= maxf or '--verdict-only' in args: break
        rng = random.Random('%d:%s' % (seed, fn))
        res = check_file(os.path.join(out, fn), xfrac, rng); n += 1
        res.update(file=fn, size=st.st_size, mtime=int(st.st_mtime), checker_sha=CHK,
                   checked=time.strftime('%Y-%m-%dT%H:%M:%S'), checker='certcheck_strict')
        with open(led, 'a') as fh: fh.write(json.dumps(res) + "\n")
        done[fn] = res
        print(fn, 'sbs', res['sbs'], 'leaves', res['leaves'], 'bad', res['bad'], 'xchk', res['xchecked'],
              'disagree', res['disagree'], 'sec', res['sec'], flush=True)
    if only is None:
        need = set(c for c in itertools.product(range(5), range(5), range(5), range(10), range(10), range(10))
                   if any(K.in_UD([2 * c[i] + e[i] for i in range(6)]) for e in itertools.product((0, 1), repeat=6)))
        # coverage from the start boxes check_file actually tiled, in the files whose current sha was verified
        cur = [done[f] for f in files if f in done and done[f].get('sha256') == fsha[f]]
        seen = set(); dup = 0
        for v in cur:
            for sb in v['tiled_sbs']:
                sb = tuple(sb)
                if sb in seen: dup += 1
                seen.add(sb)
        print('files', len(files), 'verified', len(cur), 'bad files', sum(1 for v in cur if v['bad']),
              'disagreements', sum(v['disagree'] for v in cur), 'leaves verified', sum(v['leaves'] for v in cur),
              'xchecked', sum(v['xchecked'] for v in cur))
        print('start boxes needed', len(need), 'present', len(seen & need), 'missing', len(need - seen),
              'extra', len(seen - need), 'dup', dup)
        ok = (not (need - seen) and not (seen - need) and not dup and len(cur) == len(files)
              and not any(v['bad'] or v['disagree'] for v in cur))
        print('COMPLETE AND VALID' if ok else 'NOT (YET) COMPLETE/VALID')

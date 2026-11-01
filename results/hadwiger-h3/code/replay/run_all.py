"""Driver: replay-check certificate chunk files.

usage: python run_all.py OUTDIR CERTDIR [CERTDIR ...] [--workers 8] [--sample FRACTION] [--seed S]
  Every chunk_*.jsonl.gz in each CERTDIR is checked (fullcert, and later fullcert_ext).
  Per chunk a JSON result goes to OUTDIR/<certdir-name>__<chunk>.json (resumable).
  --sample f: check only a stratified random fraction f of the start boxes in each chunk
  (tiling is still checked for the sampled start boxes; all chord leaves of sampled boxes).
Then: python summarize.py OUTDIR
Honours C:/Users/hackathon/.compute-broker/GAMING: while it exists only worker 0 runs.
"""
import gzip, json, os, sys, time, glob, random, argparse
from fractions import Fraction as Fr
from multiprocessing import Pool, current_process
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hcheck as H

GAMING = 'C:/Users/hackathon/.compute-broker/GAMING'
BATCH = 1500


def to_int640(s):
    x = Fr(s) * H.G
    assert x.denominator == 1, ('corner not on 1/640 grid', s)
    return int(x)


def check_record(r, ci_name):
    out = {'sb': r['sb'], 'n': 0, 'plain': 0, 'chord': 0, 'fast_ok': 0, 'exact_fallback_ok': 0,
           'fail': [], 'slots': {}, 'tiling': None, 'pol': r.get('pol')}
    leaves = r['leaves']
    out['n'] = len(leaves)
    LO = np.array([[to_int640(x) for x in L['a']] for L in leaves], dtype=np.int64)
    DD = np.array([[to_int640(x) for x in L['d']] for L in leaves], dtype=np.int64)
    sb = np.array(r['sb'], dtype=np.int64) * 64
    ok, msg = H.tiling_ok(sb, sb + 64, LO, LO + DD)
    out['tiling'] = msg if not ok else 'ok'
    if not ok:
        out['fail'].append({'leaf': None, 'why': 'tiling: ' + msg})
    plain = [i for i, L in enumerate(leaves) if L['k'] == 'p']
    chord = [i for i, L in enumerate(leaves) if L['k'] == 'c']
    other = [i for i, L in enumerate(leaves) if L['k'] not in 'pc']
    for i in other:
        out['fail'].append({'leaf': i, 'why': 'unknown kind %r' % leaves[i]['k']})
    out['plain'] = len(plain); out['chord'] = len(chord)
    for b in range(0, len(plain), BATCH):
        idx = plain[b:b + BATCH]
        Q = np.array([leaves[i]['c']['q'][:39] for i in idx], dtype=np.int64).reshape(-1, 13, 3)
        fok = H.fast_plain_batch(LO[idx], LO[idx] + DD[idx], Q)
        out['fast_ok'] += int(fok.sum())
        for j in np.nonzero(~fok)[0]:
            i = idx[j]
            e, why = H.exact_leaf(leaves[i])
            if e:
                out['exact_fallback_ok'] += 1
            else:
                out['fail'].append({'leaf': i, 'why': why, 'a': leaves[i]['a'], 'd': leaves[i]['d']})
    for i in chord:
        e, info = H.exact_leaf(leaves[i])
        if e:
            key = 'ax%d plus@%s minus@%s' % (leaves[i]['ax'], info['plus'] == leaves[i]['plus']['st'][1],
                                             info['minus'] == leaves[i]['minus']['st'][1])
            out['slots'][key] = out['slots'].get(key, 0) + 1
        else:
            out['fail'].append({'leaf': i, 'why': info, 'a': leaves[i]['a'], 'd': leaves[i]['d']})
    return out


def work(args):
    fn, outdir, sample, seed, strict = args
    H.STRICT = strict
    tag = os.path.basename(os.path.dirname(fn)) + '__' + os.path.basename(fn).replace('.jsonl.gz', '')
    of = os.path.join(outdir, tag + '.json')
    if os.path.exists(of):
        return tag, 'skip'
    wid = current_process()._identity[0] if current_process()._identity else 1
    rng = random.Random('%s|%s' % (seed, tag))
    t0 = time.time()
    res = {'file': fn, 'records': [], 'sample': sample, 'strict': strict}
    with gzip.open(fn, 'rt') as f:
        for line in f:
            while wid > 1 and os.path.exists(GAMING):
                time.sleep(30)
            r = json.loads(line)
            if sample < 1 and rng.random() >= sample:
                res['records'].append({'sb': r['sb'], 'skipped': True})
                continue
            res['records'].append(check_record(r, tag))
    res['sec'] = time.time() - t0
    tmp = of + '.tmp'
    with open(tmp, 'w') as g:
        json.dump(res, g)
    os.replace(tmp, of)
    nf = sum(len(x.get('fail', [])) for x in res['records'])
    return tag, 'done %.0fs fails=%d' % (res['sec'], nf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir'); ap.add_argument('certdirs', nargs='+')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--sample', type=float, default=1.0)
    ap.add_argument('--seed', default='replay')
    ap.add_argument('--only', default=None, help='comma list of chunk numbers')
    ap.add_argument('--strict', action='store_true', help='interior membership everywhere (see hcheck.STRICT)')
    a = ap.parse_args()
    assert a.workers <= 8
    os.makedirs(a.outdir, exist_ok=True)
    files = []
    for d in a.certdirs:
        files += sorted(glob.glob(os.path.join(d, 'chunk_*.jsonl.gz')))
    if a.only:
        keep = {int(x) for x in a.only.split(',')}
        files = [f for f in files if int(os.path.basename(f)[6:11]) in keep]
    # big files first for load balance
    files.sort(key=lambda f: -os.path.getsize(f))
    print('files', len(files), flush=True)
    with Pool(a.workers) as pool:
        for tag, st in pool.imap_unordered(work, [(f, a.outdir, a.sample, a.seed, a.strict) for f in files]):
            print(time.strftime('%H:%M:%S'), tag, st, flush=True)


if __name__ == '__main__':
    main()

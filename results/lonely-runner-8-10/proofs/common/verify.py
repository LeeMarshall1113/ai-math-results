"""One-command re-check of a proof package:   python proofs/common/verify.py proofs/k7

Trusted: gen_cnf.py (CNF generator), mass.py (integer inequality), the LRAT checkers
(cake_lpr, CakeML-verified; gate_lrat_check, ours), and the cited lemmas (PROOF.md).
NOT trusted: the manifest, the stored proofs, the SAT solver, every hash.

For the package's gate list it:
  1. checks lcm(2..k+1) * prod(gates) > B_k exactly (mass.py);
  2. for each gate P: regenerates G(k, L, P) with gen_cnf.py, obtains an LRAT proof
     (decompresses <art>/k{k}_L{L}_p{P}.lrat.xz, or re-runs CaDiCaL with the recorded options
     when the file is absent), and requires BOTH checkers to accept it against the
     regenerated CNF.  Hash mismatches are reported (reproducibility) but cannot make a
     check pass: acceptance comes only from the checkers.
Exit status 0 iff every step passes.

options: --art DIR (default <pkg>/artifacts)  --jobs J  --cake-jobs C  --only P [P ...]
         --no-cake (own checker only; the run is then reported as NOT fully verified)
env: GATE_LRAT_CHECK, CAKE_LPR (path; prefix 'wsl:' runs it in WSL with /mnt paths), CADICAL,
     CAKE_HEAP_MB (default 3000), WORK (temp dir; default <art>/tmp)
"""
import argparse
import hashlib
import json
import lzma
import os
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gen_cnf  # noqa: E402
import mass  # noqa: E402

EXE = '.exe' if os.name == 'nt' else ''
OWN = os.environ.get('GATE_LRAT_CHECK', str(HERE / 'checkers' / 'bin' / ('gate_lrat_check' + EXE)))
CAKE = os.environ.get('CAKE_LPR', str(HERE / 'checkers' / 'bin' / 'cake_lpr'))
CADICAL = os.environ.get('CADICAL', 'cadical')
HEAP = os.environ.get('CAKE_HEAP_MB', '3000')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 22), b''):
            h.update(blk)
    return h.hexdigest()


def wslpath(p):
    p = str(Path(p).resolve()).replace('\\', '/')
    return '/mnt/' + p[0].lower() + p[2:] if p[1] == ':' else p


def run_cake(cnf, lrat):
    args = ['--CML_HEAP_SIZE=' + HEAP, '--CML_STACK_SIZE=1000']
    if CAKE.startswith('wsl:'):
        cmd = ['wsl', '-d', os.environ.get('WSL_DISTRO', 'Ubuntu-24.04'), '--', CAKE[4:]] + args + [wslpath(cnf), wslpath(lrat)]
    else:
        cmd = [CAKE] + args + [str(cnf), str(lrat)]
    r = subprocess.run(cmd, capture_output=True)
    out = r.stdout.decode(errors='replace').replace('\0', '')
    err = r.stderr.decode(errors='replace').replace('\0', '')
    return ('s VERIFIED UNSAT' in out and r.returncode == 0), (out + err).strip()[-300:]


def check_gate(k, L, P, rec, art, work, use_cake, cake_lock):
    res = dict(P=P)
    t0 = time.perf_counter()
    stem = f'k{k}_L{L}_p{P}'
    cnf = work / (stem + '.cnf')
    lrat = work / (stem + '.lrat')
    b, meta = gen_cnf.generate(k, L, P)
    with open(cnf, 'w', newline='\n') as f:
        f.write(gen_cnf.dimacs(b, meta))
    res['cnf_sha256'] = sha(cnf)
    res['cnf_hash_matches_manifest'] = (res['cnf_sha256'] == rec.get('cnf_sha256'))
    xz = art / (stem + '.lrat.xz')
    raw = art / (stem + '.lrat')
    t = time.perf_counter()
    keep = False
    if raw.exists():
        lrat, keep = raw, True
        res['proof_source'] = 'stored-raw'
    elif xz.exists():
        with lzma.open(xz, 'rb') as fi, open(lrat, 'wb') as fo:
            shutil.copyfileobj(fi, fo, 1 << 22)
        res['proof_source'] = 'stored'
    else:
        opts = rec.get('solver_opts', '-q --lrat --no-binary').split()
        for extra in ([], ['--no-decompose'], ['--plain']):   # Windows CaDiCaL 3.0.1 crash workaround
            r = subprocess.run([CADICAL] + opts + extra + [str(cnf), str(lrat)], capture_output=True)
            if r.returncode in (10, 20):
                break
        res['solver_exit'] = r.returncode
        res['solver_opts_used'] = ' '.join(opts + extra)
        res['proof_source'] = 'regenerated'
    res['t_obtain'] = round(time.perf_counter() - t, 1)
    res['lrat_sha256'] = sha(lrat) if lrat.exists() else None
    res['lrat_hash_matches_manifest'] = res['lrat_sha256'] is not None and res['lrat_sha256'] == rec.get('lrat_sha256')
    t = time.perf_counter()
    r = subprocess.run([OWN, str(cnf), str(lrat)], capture_output=True, text=True)
    res['own_ok'] = r.returncode == 0
    res['own_msg'] = (r.stdout or r.stderr).strip()[-200:]
    res['t_own'] = round(time.perf_counter() - t, 1)
    if use_cake:
        with cake_lock:
            t = time.perf_counter()
            ok, msg = run_cake(cnf, lrat)
            res['t_cake'] = round(time.perf_counter() - t, 1)
        res['cake_ok'], res['cake_msg'] = ok, msg
    res['lrat_bytes'] = lrat.stat().st_size if lrat.exists() else 0
    for f in ((cnf,) if keep else (cnf, lrat)):
        if f.exists():
            f.unlink()
    res['t_total'] = round(time.perf_counter() - t0, 1)
    res['ok'] = res['own_ok'] and (res.get('cake_ok', False) if use_cake else True)
    return res


def verify_split(a, pkg, man, k, L, gates, todo, ok_mass, t0, work):
    """Split mode (k >= 8): per gate, the trusted level-1 split (level1.py, or its C++ port
    when LEVEL1_EXE is set) gives leaves (I, E); each leaf CNF G + units is regenerated, solved
    by CaDiCaL (untrusted) and its LRAT proof must be accepted by BOTH checkers."""
    import split
    results = []
    glog = pkg / 'verify_gates.jsonl'          # one line per checked gate, written as it finishes
    prior = {}
    if a.resume and glog.exists():
        for line in open(glog):
            r = json.loads(line)
            if r.get('ok') and r.get('cake_used'):
                prior[r['P']] = r
    for P in sorted(todo):
        if P in prior:
            results.append(prior[P])
            print('[2] p=%-4d ok=True (checked earlier, %s)' % (P, prior[P].get('when_utc')), flush=True)
            continue
        r = split.gate(k, L, P, a.jobs, not a.no_cake, work / 'split')
        r['ok'] = bool(r['gate'] and r['crashes'] == 0 and r['cake_rejects'] == 0)
        r['when_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        with open(glog, 'a') as f:
            f.write(json.dumps(r) + chr(10))
        results.append(r)
        print('[2] p=%-4d leaves=%-6d ok=%s  proofs %.0f MB  level1 %.0fs solve %.0fs own %.0fs cake %.0fs wall %.0fs'
              % (P, r['leaves'], r['ok'], r['bytes'] / 1e6, r['t_level1'], r['t_solve'], r['t_own'], r['t_cake'], r['t_wall']),
              flush=True)
    allok = ok_mass and all(r['ok'] for r in results) and (bool(a.only) or {r['P'] for r in results} >= set(gates))
    summary = dict(k=k, L=L, mode='split', gates=len(gates), mass_ok=ok_mass, all_ok=allok, cake_used=not a.no_cake,
                   wall_s=round(time.time() - t0, 1), proof_bytes=sum(r['bytes'] for r in results),
                   leaves=sum(r['leaves'] for r in results), steps=sum(r['steps'] for r in results),
                   cpu_solve_s=round(sum(r['t_solve'] for r in results), 1),
                   cpu_own_s=round(sum(r['t_own'] for r in results), 1),
                   cpu_cake_s=round(sum(r['t_cake'] for r in results), 1),
                   when_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    with open(pkg / 'verify_log.jsonl', 'a') as f:
        f.write(json.dumps(dict(summary=summary, results=results)) + '\n')
    print('[3]', json.dumps(summary))
    if allok and not a.no_cake and not a.only:
        print('VERIFIED: LRC(%d) follows from the cited lemmas + these %d gates (%d leaf certificates, both checkers).'
              % (k, len(gates), summary['leaves']))
    elif allok:
        print('PASSED (partial run)')
    else:
        print('FAILED')
    sys.exit(0 if allok else 1)


def main():
    if os.environ.get('SPLIT_SAMPLE'):   # sampled leaves are a measurement, never a proof (referee, 2026-10-22)
        sys.exit('refusing to verify: SPLIT_SAMPLE is set; unset it for a proof run')
    ap = argparse.ArgumentParser()
    ap.add_argument('pkg')
    ap.add_argument('--art')
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--cake-jobs', type=int, default=1)
    ap.add_argument('--only', type=int, nargs='+')
    ap.add_argument('--no-cake', action='store_true')
    ap.add_argument('--resume', action='store_true',
                    help='split mode: reuse gates already accepted by both checkers (verify_gates.jsonl)')
    a = ap.parse_args()
    pkg = Path(a.pkg).resolve()
    man = json.load(open(pkg / 'manifest.json'))
    k, L, gates = man['k'], man['L'], man['gates']
    art = Path(a.art or os.environ.get('ART', pkg / 'artifacts')).resolve()
    work = Path(os.environ.get('WORK', art / 'tmp'))
    work.mkdir(parents=True, exist_ok=True)
    recs = {r['P']: r for r in man.get('records', [])}
    t0 = time.time()
    ok_mass, info = mass.check(k, gates)
    print('[1] gate-set size: %d gates, ln prod = %.4f > needed %.4f : %s (exact integer test)'
          % (len(gates), info['ln_gates'], info['needed_from_gates'], ok_mass), flush=True)
    todo = a.only or gates
    if man.get('mode') == 'split':
        return verify_split(a, pkg, man, k, L, gates, todo, ok_mass, t0, work)
    sem = threading.Semaphore(a.cake_jobs)
    results = []
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(check_gate, k, L, P, recs.get(P, {}), art, work, not a.no_cake, sem)
                for P in sorted(todo, reverse=True)]
        for fu in futs:
            r = fu.result()
            results.append(r)
            print('[2] p=%-4d own=%s cake=%s  proof %s %.0f MB  times: obtain %.0fs own %.0fs cake %ss  hashes cnf=%s lrat=%s'
                  % (r['P'], r['own_ok'], r.get('cake_ok', 'skipped'), r['proof_source'], r['lrat_bytes'] / 1e6,
                     r['t_obtain'], r['t_own'], r.get('t_cake', '-'), r['cnf_hash_matches_manifest'],
                     r['lrat_hash_matches_manifest']), flush=True)
    allok = ok_mass and all(r['ok'] for r in results) and (bool(a.only) or set(r['P'] for r in results) >= set(gates))
    summary = dict(k=k, L=L, gates=len(gates), mass_ok=ok_mass, all_ok=allok, cake_used=not a.no_cake,
                   wall_s=round(time.time() - t0, 1), proof_bytes=sum(r['lrat_bytes'] for r in results),
                   cpu_own_s=round(sum(r['t_own'] for r in results), 1),
                   cpu_cake_s=round(sum(r.get('t_cake', 0) for r in results), 1),
                   when_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    with open(pkg / 'verify_log.jsonl', 'a') as f:
        f.write(json.dumps(dict(summary=summary, results=results)) + '\n')
    print('[3]', json.dumps(summary))
    if allok and not a.no_cake and not a.only:
        print('VERIFIED: LRC(%d) follows from the cited lemmas + these %d certificates (both checkers).' % (k, len(gates)))
    elif allok:
        print('PASSED (partial run: %s)' % ('own checker only' if a.no_cake else 'subset of gates'))
    else:
        print('FAILED')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()

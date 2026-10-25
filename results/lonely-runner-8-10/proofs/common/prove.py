"""Proof PRODUCER (untrusted).  For each prime: generate G(K,L,P) with gen_cnf.py, solve it with
CaDiCaL writing a text LRAT proof, check the proof with gate_lrat_check, xz-compress it into
the artifact directory, and append a JSON record.  A SAT answer is decoded and re-checked
against the definition (it certifies a NON-gate).  Nothing here is trusted: verify.py
re-generates every CNF and re-checks every proof.

usage: python prove.py K L --primes 31 37 ... --art DIR --out records.jsonl [--jobs 4]
env:   CADICAL (solver binary), GATE_LRAT_CHECK (checker binary)
"""
import argparse
import hashlib
import json
import lzma
import math
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gen_cnf  # noqa: E402

CADICAL = os.environ.get('CADICAL', 'cadical')
CHECK = os.environ.get('GATE_LRAT_CHECK', str(HERE / 'checkers' / 'bin' / 'gate_lrat_check.exe'))


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 22), b''):
            h.update(blk)
    return h.hexdigest()


def improper_def(x, K, L, P):
    N = L * P
    for i in range(K):
        g = N
        for j in range(K):
            if j != i:
                g = math.gcd(g, x[j])
        if g > 1:
            return False
    return not any(all((K + 1) * min(t * v % N, N - t * v % N) >= N for v in x) for t in range(N))


def solve(cnf, lrat, timeout):
    """CaDiCaL 3.0.1 on Windows sometimes crashes in 'decompose' (exit 127): retry."""
    for opts in ([], ['--no-decompose'], ['--plain']):
        cmd = [CADICAL, '-q'] + opts + (['--lrat', '--no-binary'] if lrat else []) + [str(cnf)]
        if lrat:
            cmd.append(str(lrat))
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return 'UNKNOWN', '', opts
        if r.returncode in (10, 20):
            st = 'SATISFIABLE' if r.returncode == 10 else 'UNSATISFIABLE'
            return st, r.stdout, opts
    return 'CRASH', '', opts


def one(K, L, P, art, timeout, keep_raw=False):
    art = Path(art)
    art.mkdir(parents=True, exist_ok=True)
    stem = art / f'k{K}_L{L}_p{P}'
    cnf, lrat = stem.with_suffix('.cnf'), stem.with_suffix('.lrat')
    b, meta = gen_cnf.generate(K, L, P)
    with open(cnf, 'w', newline='\n') as f:
        f.write(gen_cnf.dimacs(b, meta))
    rec = dict(K=K, L=L, P=P, vars=b.nv, clauses=len(b.clauses), cnf_sha256=sha(cnf))
    t = time.perf_counter()
    st, out, opts = solve(cnf, lrat, timeout)
    rec.update(status=st, solver='cadical-3.0.1', solver_opts=' '.join(['-q'] + opts + ['--lrat', '--no-binary']),
               t_solve=round(time.perf_counter() - t, 2))
    if st == 'SATISFIABLE':
        model = {int(v) for ln in out.splitlines() if ln.startswith('v ') for v in ln[2:].split() if int(v) > 0}
        D = [v for v in range(2, L * P // 2 + 1) if v % P]
        S = [v for i, v in enumerate(D) if meta['x_first'] + i in model]
        tup = [1] + S + [1] * (K - 1 - len(S))
        rec['improper_tuple'] = tup
        rec['tuple_checked'] = improper_def(tup, K, L, P)
        assert rec['tuple_checked'], tup
    elif st == 'UNSATISFIABLE':
        t = time.perf_counter()
        r = subprocess.run([CHECK, str(cnf), str(lrat)], capture_output=True, text=True)
        rec['t_check_own'] = round(time.perf_counter() - t, 2)
        assert r.returncode == 0, r.stderr
        steps, live = map(int, r.stdout.split())
        rec.update(lrat_steps=steps, lrat_bytes=lrat.stat().st_size, lrat_sha256=sha(lrat))
        xz = lrat.with_suffix('.lrat.xz')
        t = time.perf_counter()
        with open(lrat, 'rb') as fi, lzma.open(xz, 'wb', preset=6) as fo:
            shutil.copyfileobj(fi, fo, 1 << 22)
        rec.update(xz_bytes=xz.stat().st_size, xz_sha256=sha(xz), t_xz=round(time.perf_counter() - t, 2))
    if lrat.exists() and not keep_raw:
        lrat.unlink()
    cnf.unlink()
    return rec


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('K', type=int)
    ap.add_argument('L', type=int)
    ap.add_argument('--primes', type=int, nargs='+', required=True)
    ap.add_argument('--art', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--timeout', type=float, default=7200)
    a = ap.parse_args()
    done = set()
    if os.path.exists(a.out):
        done = {json.loads(l)['P'] for l in open(a.out) if l.strip()}
    todo = sorted((p for p in a.primes if p not in done), reverse=True)
    with ProcessPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(one, a.K, a.L, p, a.art, a.timeout): p for p in todo}
        for fu in as_completed(futs):
            try:
                rec = fu.result()
            except Exception as e:  # keep going; record the failure
                rec = dict(K=a.K, L=a.L, P=futs[fu], status='ERROR', error=repr(e))
            print(json.dumps(rec), flush=True)
            with open(a.out, 'a') as f:
                f.write(json.dumps(rec) + '\n')

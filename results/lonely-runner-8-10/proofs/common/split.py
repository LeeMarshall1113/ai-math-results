"""Split-mode certificates for one gate (used when one whole-gate proof is too large).

Producer AND checker in one: for gate P it builds G(K,L,P) (gen_cnf), the level-1 leaves
(level1.leaves, trusted), and for every leaf j the CNF  G + units(I_j) + units(-E_j)
(built here, trusted: plain concatenation).  Each leaf CNF is refuted by CaDiCaL (untrusted)
with a text LRAT proof, which is checked by gate_lrat_check and, unless --no-cake, by
cake_lpr; the proof is then deleted (proofs are regenerated deterministically on demand;
the SHA-256 of each is recorded).  A SAT leaf decodes to an improper tuple (non-gate).

usage: python split.py K L P [P ...] --out rec.jsonl [--jobs J] [--no-cake] [--work DIR]
"""
import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gen_cnf  # noqa: E402
import gen_cnf_composite  # noqa: E402
import level1  # noqa: E402


def generator(L):
    try:
        gen_cnf.prime_power_base(L)
        return gen_cnf
    except ValueError:
        return gen_cnf_composite
from verify import OWN, sha  # noqa: E402

NL = chr(10)
CADICAL = os.environ.get('CADICAL', 'cadical')


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


def leaf_cnf_text(b, I, E):
    units = [c - 1 for c in I] + [-(c - 1) for c in E]       # y_c is variable c-1
    head = 'p cnf %d %d\n' % (b.nv, len(b.clauses) + len(units))
    return head, ''.join('%d 0\n' % u for u in units)


def run_leaf(K, L, P, b, meta, body, j, I, E, work, use_cake, cake_lock):
    head, tail = leaf_cnf_text(b, I, E)
    cnf = work / f'k{K}p{P}_leaf{j}.cnf'
    lrat = cnf.with_suffix('.lrat')
    keep = os.environ.get('KEEP_LRAT_DIR')       # keep proofs on disk (cake_lpr check can run later)
    if keep:
        kd = Path(keep) / f'k{K}_p{P}'
        kd.mkdir(parents=True, exist_ok=True)
        lrat = kd / f'leaf{j}.lrat'
    with open(cnf, 'w', newline='\n') as f:
        f.write('c leaf %d of G(%d,%d,%d): I=%s E=%s\n' % (j, K, L, P, list(I), list(E)))
        f.write(head)
        f.write(body)
        f.write(tail)
    rec = dict(j=j, nI=len(I), nE=len(E))
    t = time.perf_counter()
    part = lrat.with_name(lrat.name + '.part')
    if keep and lrat.exists():                  # reuse a complete stored proof; the checkers below decide
        extra = []
        r = subprocess.CompletedProcess([], 20, '', '')
    else:
        for extra in ([], ['--no-decompose'], ['--plain']):
            r = subprocess.run([CADICAL, '-q', '--lrat', '--no-binary'] + extra + [str(cnf), str(part if keep else lrat)],
                               capture_output=True, text=True)
            if r.returncode in (10, 20):
                break
        if keep and r.returncode == 20:
            os.replace(part, lrat)
        elif keep and part.exists():
            part.unlink()
    rec['t_solve'] = round(time.perf_counter() - t, 2)
    rec['opts'] = ' '.join(extra)
    if r.returncode == 10:
        model = {int(v) for ln in r.stdout.splitlines() if ln.startswith('v ') for v in ln[2:].split() if int(v) > 0}
        if 'anchors' in meta:       # composite level: S includes the anchor; try paddings
            import itertools
            D = [v for v in range(1, L * P // 2 + 1) if v % P]
            S = [v for i, v in enumerate(D) if meta['x_first'] + i in model]
            tup = next((S + list(pd) for pd in itertools.combinations_with_replacement(S, K - len(S))
                        if improper_def(S + list(pd), K, L, P)), S)
        else:
            D = [v for v in range(2, L * P // 2 + 1) if v % P]
            S = [v for i, v in enumerate(D) if meta['x_first'] + i in model]
            tup = [1] + S + [1] * (K - 1 - len(S))
        rec.update(status='SAT', improper_tuple=tup,
                   tuple_checked=len(tup) == K and improper_def(tup, K, L, P))
    elif r.returncode == 20:
        rec['status'] = 'UNSAT'
        t = time.perf_counter()
        o = subprocess.run([OWN, str(cnf), str(lrat)], capture_output=True, text=True)
        rec['own_ok'] = o.returncode == 0
        rec['t_own'] = round(time.perf_counter() - t, 2)
        if rec['own_ok']:
            rec['steps'] = int(o.stdout.split()[0])
        rec['bytes'] = lrat.stat().st_size
        rec['lrat_sha256'] = sha(lrat)
        if use_cake:
            rec['files'] = (str(cnf), str(lrat))      # checked in one batch per gate
            return rec
    else:
        rec['status'] = 'CRASH rc=%d' % r.returncode
    for f in ((cnf,) if keep and rec['status'] == 'UNSAT' else (cnf, lrat)):
        if f.exists():
            f.unlink()
    return rec


def cake_batch(recs, work, tag):
    """Run cake_lpr on every (cnf, lrat) pair in ONE process launch (WSL start-up is ~1 s).
    Each pair is a separate cake_lpr run; a pair counts only if it prints 's VERIFIED UNSAT'."""
    from verify import CAKE, HEAP, wslpath
    todo = [r for r in recs if 'files' in r]
    if not todo:
        return 0.0
    wsl = CAKE.startswith('wsl:')
    exe = CAKE[4:] if wsl else CAKE
    conv = wslpath if wsl else str
    npar = max(1, int(os.environ.get('CAKE_PAR', '6')))
    chunks = [todo[i::npar] for i in range(npar) if todo[i::npar]]
    scripts = []
    for ci, chunk in enumerate(chunks):
        script = Path(work) / f'cake_{tag}_{ci}.sh'
        scripts.append(script)
        with open(script, 'w', newline=NL) as f:
            if wsl:
                # Windows host: copy each pair to a Linux temp dir first; cake_lpr's small reads over
                # WSL's /mnt/c (9p) ran at ~2 leaves/s.  Stale copies must never be checked in place
                # of a missing pair: delete first and require both copies to succeed.
                f.write('T=$(mktemp -d)' + NL)
                for r in chunk:
                    c, l = r['files']
                    f.write('rm -f $T/a.cnf $T/a.lrat; '
                            'if cat "%s" > $T/a.cnf && cat "%s" > $T/a.lrat && '
                            '%s --CML_HEAP_SIZE=%s --CML_STACK_SIZE=500 $T/a.cnf $T/a.lrat | grep -q "^s VERIFIED UNSAT$"; '
                            'then echo "OK %d"; else echo "FAIL %d"; fi' % (conv(c), conv(l), exe, HEAP, r['j'], r['j']) + NL)
                f.write('rm -rf $T' + NL)
            else:
                for r in chunk:
                    c, l = r['files']
                    f.write('if %s --CML_HEAP_SIZE=%s --CML_STACK_SIZE=500 "%s" "%s" | grep -q "^s VERIFIED UNSAT$"; '
                            'then echo "OK %d"; else echo "FAIL %d"; fi' % (exe, HEAP, c, l, r['j'], r['j']) + NL)
    t = time.perf_counter()
    procs = []
    for script in scripts:
        cmd = (['wsl', '-d', os.environ.get('WSL_DISTRO', 'Ubuntu-24.04'), '--', 'sh', conv(script)] if wsl
               else ['sh', str(script)])
        procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL))
    out = ''.join(pr.communicate()[0].decode(errors='replace').replace(chr(0), '') for pr in procs)
    dt = time.perf_counter() - t
    okset = {int(w.split()[1]) for w in out.splitlines() if w.startswith('OK ')}
    for r in todo:
        r['cake_ok'] = r['j'] in okset
        for fn in r.pop('files'):
            if os.environ.get('KEEP_LRAT_DIR') and fn.endswith('.lrat'):
                continue
            if os.path.exists(fn):
                os.unlink(fn)
    for script in scripts:
        script.unlink()
    return dt


def gate(K, L, P, jobs, use_cake, work, cake_jobs=1):
    t0 = time.perf_counter()
    G = generator(L)
    b, meta = G.generate(K, L, P)
    body = ''.join(' '.join(map(str, c)) + ' 0\n' for c in b.clauses)
    st = {}
    exe = os.environ.get('LEVEL1_EXE')          # C++ port (checkers/level1.cpp), same leaf list
    if exe:
        cached = Path(os.environ['LEAVES_DIR']) / f'k{K}_p{P}.leaves' if os.environ.get('LEAVES_DIR') else None
        if cached is not None and cached.exists():     # output of an earlier level1.exe run (sha recorded)
            out = cached.read_text()
            st['leaves_sha256'] = sha(cached)
        else:
            out = subprocess.run([exe, str(K), str(P)], capture_output=True, text=True, check=True).stdout
        lv = []
        for ln in out.splitlines():
            if ln.startswith('#'):
                n, d, nl = map(int, ln[1:].split())
                st.update(nodes=n, dead=d, leaves=nl, enumerator='level1.cpp')
            else:
                a, b2 = ln[1:].split('|')
                lv.append((tuple(map(int, a.split())), tuple(map(int, b2.strip()[1:].split()))))
        assert len(lv) == st['leaves']
    else:
        lv = level1.leaves(K, P, st)
        st['enumerator'] = 'level1.py'
    if os.environ.get('SPLIT_SAMPLE'):          # measurement only: every n-th leaf; never a proof
        n = int(os.environ['SPLIT_SAMPLE'])
        lv = lv[::n]
    t_l1 = time.perf_counter() - t0
    lock = threading.Semaphore(cake_jobs)
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(jobs) as ex:
        recs = list(ex.map(lambda a: run_leaf(K, L, P, b, meta, body, a[0], a[1][0], a[1][1], work, use_cake, lock),
                           enumerate(lv)))
    t_cake = cake_batch(recs, work, f'k{K}p{P}') if use_cake else 0.0
    unsat = [r for r in recs if r['status'] == 'UNSAT']
    sat = [r for r in recs if r['status'] == 'SAT']
    ok = (len(unsat) == len(lv) and all(r['own_ok'] for r in unsat)
          and (all(r.get('cake_ok') for r in unsat) if use_cake else True))
    h = hashlib.sha256(''.join(r.get('lrat_sha256') or '-' for r in recs).encode()).hexdigest()
    return dict(K=K, L=L, P=P, cnf_sha256=hashlib.sha256(G.dimacs(b, meta).encode()).hexdigest(),
                l1_nodes=st['nodes'], leaves_sha256=st.get('leaves_sha256'), l1_dead=st['dead'], l1_enum=st['enumerator'], leaves=len(lv), t_level1=round(t_l1, 1),
                gate=ok, nonsat=len(sat), improper=sat[0].get('improper_tuple') if sat else None,
                tuple_checked=sat[0].get('tuple_checked') if sat else None,
                crashes=sum(r['status'].startswith('CRASH') for r in recs),
                steps=sum(r.get('steps', 0) for r in unsat), bytes=sum(r.get('bytes', 0) for r in unsat),
                max_leaf_bytes=max((r.get('bytes', 0) for r in unsat), default=0),
                t_solve=round(sum(r['t_solve'] for r in recs), 1),
                t_own=round(sum(r.get('t_own', 0) for r in recs), 1),
                t_cake=round(t_cake, 1), cake_rejects=sum(1 for r in unsat if use_cake and not r.get('cake_ok')),
                cake_used=use_cake, proofs_sha256_of_list=h,
                t_wall=round(time.perf_counter() - t0, 1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('K', type=int)
    ap.add_argument('L', type=int)
    ap.add_argument('primes', type=int, nargs='+')
    ap.add_argument('--out', required=True)
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--cake-jobs', type=int, default=1)
    ap.add_argument('--no-cake', action='store_true')
    ap.add_argument('--work', default=str(HERE.parent / 'work'))
    a = ap.parse_args()
    done = set()
    if os.path.exists(a.out):
        done = {json.loads(l)['P'] for l in open(a.out) if l.strip()}
    for P in a.primes:
        if P in done:
            continue
        rec = gate(a.K, a.L, P, a.jobs, not a.no_cake, a.work, a.cake_jobs)
        print(json.dumps(rec), flush=True)
        with open(a.out, 'a') as f:
            f.write(json.dumps(rec) + '\n')

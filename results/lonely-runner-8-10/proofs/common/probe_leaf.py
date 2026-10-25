"""Measurement (untrusted): solve a few sampled leaves of a gate; report time, LRAT size, SAT/UNSAT.
usage: python probe_leaf.py K L P NSAMPLE   (never a proof: samples only)"""
import sys, subprocess, time, os, tempfile
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import split
K, L, P, N = map(int, sys.argv[1:5])
G = split.generator(L)
b, meta = G.generate(K, L, P)
body = ''.join(' '.join(map(str, c)) + ' 0\n' for c in b.clauses)
out = subprocess.run([str(HERE / 'checkers/bin/level1.exe'), str(K), str(P)], capture_output=True, text=True).stdout
lv = []
for ln in out.splitlines():
    if not ln.startswith('#'):
        a, e = ln[1:].split('|'); lv.append((tuple(map(int, a.split())), tuple(map(int, e.strip()[1:].split()))))
print('leaves', len(lv), flush=True)
work = Path(tempfile.mkdtemp(dir=str(HERE.parent / 'work'))) if (HERE.parent / 'work').exists() else Path(tempfile.mkdtemp())
step = max(1, len(lv) // N)
for j in range(0, len(lv), step)[:N] if False else list(range(0, len(lv), step))[:N]:
    r = split.run_leaf(K, L, P, b, meta, body, j, lv[j][0], lv[j][1], work, False, None)
    print({a: r.get(a) for a in ('j', 'status', 't_solve', 'steps', 'bytes', 't_own', 'improper_tuple', 'tuple_checked')}, flush=True)

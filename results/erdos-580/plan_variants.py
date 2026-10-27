"""Print run_certs job args with leaf-keep alternatives for trees whose rooted-core CNF is SAT.
Alternatives: every nonempty subset of parents, ordered by |F'| then by more marks.
usage: python plan_variants.py m i1 i2 ...   ->  m:i1@keep=..|keep=..  (one per line)"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lks import trees, core, core_keep
m = int(sys.argv[1])
for i in map(int, sys.argv[2:]):
    T = trees(m)[i]
    P = core(T, m)[2]
    c = []
    for r in range(1, len(P) + 1):
        for keep in itertools.combinations(P, r):
            cv, ce, mk = core_keep(T, m, keep)
            c.append((len(cv), -len(mk), "keep=" + ".".join(map(str, keep))))
    c.sort()
    print(f"{m}:{i}@" + "|".join(x[2] for x in c))

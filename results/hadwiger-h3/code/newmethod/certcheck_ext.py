"""certcheck_ext: coverage verdict for the repaired certificate fullcert/ + fullcert_ext/ (referee F1, 2026-10-31).

Leaf and tiling checks are NOT redone here: they are done by the unchanged certcheck_fast.py, run on each directory
(fullcert/ already done; fullcert_ext/ via  certcheck_fast.py fullcert_ext --only-file fullcert_ext/all_chunks.txt).
This script reads both directories' ledgers, accepts a file only if the ledger entry's sha256 equals the file's
current sha256 and its checker_sha equals sha256(certcheck_fast.py + certcheck.py + hw.py) (the same rule as
certcheck_fast.py), and then replaces certcheck_fast's U^D_20 coverage test by the symmetry-correct one:

  (1) fullcert/ tiles exactly the 39,380 U^D_20 start boxes (as before), fullcert_ext/ tiles exactly the boxes listed
      in dstar_ext_min.json, every start box once, no overlap between the two;
  (2) ORBIT COVERAGE: for every one of the 10^6 closed 1/10 boxes B of [0,1]^6 there is a cube symmetry g with g(B)
      a verified start box.  Each of the 48 symmetries acts on p as a signed coordinate permutation
      (p'_m = p_s(m) or 1 - p_s(m)), so it maps 1/10 grid boxes to 1/10 grid boxes; every p lies in some closed B,
      hence g(p) lies in a certified start box.  The action is derived here from the contact points themselves
      (independent of dstar_domain.py and of the certifier).

Usage: python certcheck_ext.py   (paths fixed: newmethod/fullcert, newmethod/fullcert_ext)"""
import os, sys, json, hashlib, itertools
from fractions import Fraction as F
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import certcheck as K

CHK = hashlib.sha256(b''.join(open(os.path.join(here, *p), 'rb').read() for p in
                     (('certcheck_fast.py',), ('certcheck.py',), ('..', 'hw.py')))).hexdigest()

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(1 << 20), b''): h.update(blk)
    return h.hexdigest()

def verified(out):
    done = {}
    for lname in sorted(f for f in os.listdir(out) if f.startswith('certcheck_fast_done') and f.endswith('.jsonl')):
        for line in open(os.path.join(out, lname)):
            try: v = json.loads(line)
            except ValueError: continue
            if v.get('checker_sha') == CHK: done[v['file']] = v
    files = sorted(f for f in os.listdir(out) if f.endswith('.jsonl.gz'))
    cur = []; unver = []
    for fn in files:
        v = done.get(fn)
        if v and v.get('sha256') == sha(os.path.join(out, fn)): cur.append(v)
        else: unver.append(fn)
    return files, cur, unver

def tiled(cur):
    seen = set(); dup = 0
    for v in cur:
        for sb in v['tiled_sbs']:
            sb = tuple(sb); dup += sb in seen; seen.add(sb)
    return seen, dup

def signed_perms():
    """the 48 symmetries as maps on p, derived from the labelled contact points of a generic p."""
    base = [F(v, 101) for v in (3, 7, 13, 19, 29, 37)]          # all distinct and distinct from 1 - x
    def contacts(p):
        p1, p2, p3, p4, p5, p6 = p
        return [((0, p1, p2), 0, 0), ((1, p1, p2), 0, 1), ((p3, 0, p4), 1, 0), ((p3, 1, p4), 1, 1),
                ((p5, p6, 0), 2, 0), ((p5, p6, 1), 2, 1)]
    out = set()
    for perm in itertools.permutations(range(3)):
        for fl in itertools.product((0, 1), repeat=3):
            face0 = {}
            for q, ax, sd in contacts(base):
                r = tuple((1 - q[perm[i]]) if fl[i] else q[perm[i]] for i in range(3))
                i = perm.index(ax); s2 = sd ^ fl[i]
                if r[i] != (1 if s2 else 0): raise AssertionError('contact point left its face')
                if s2 == 0: face0[i] = r
            newp = (face0[0][1], face0[0][2], face0[1][0], face0[1][2], face0[2][0], face0[2][1])
            m = []
            for y in newp:
                if y in base: m.append((base.index(y), 0))
                else: m.append((base.index(1 - y), 1))
            out.add(tuple(m))
    assert len(out) == 48
    return sorted(out)

def main():
    nm = here
    f0, cur0, un0 = verified(os.path.join(nm, 'fullcert'))
    f1, cur1, un1 = verified(os.path.join(nm, 'fullcert_ext'))
    S0, d0 = tiled(cur0); S1, d1 = tiled(cur1)
    NEED = set(c for c in itertools.product(range(5), range(5), range(5), range(10), range(10), range(10))
               if any(K.in_UD([2 * c[i] + e[i] for i in range(6)]) for e in itertools.product((0, 1), repeat=6)))
    EXT = set(tuple(c) for c in json.load(open(os.path.join(nm, 'dstar_ext_min.json')))['boxes'])
    bad = sum(1 for v in cur0 + cur1 if v['bad'] or v['disagree'])
    print('fullcert     files %d verified %d unverified %d leaves %d' % (len(f0), len(cur0), len(un0), sum(v['leaves'] for v in cur0)))
    print('fullcert_ext files %d verified %d unverified %d leaves %d' % (len(f1), len(cur1), len(un1), sum(v['leaves'] for v in cur1)))
    print('bad/disagreeing files', bad, '| unverified ext files', un1[:10])
    print('fullcert start boxes: need %d present %d missing %d extra %d dup %d' % (len(NEED), len(S0 & NEED), len(NEED - S0), len(S0 - NEED), d0))
    print('ext start boxes: listed %d present %d missing %d extra %d dup %d overlap-with-fullcert %d' % (
        len(EXT), len(S1 & EXT), len(EXT - S1), len(S1 - EXT), d1, len(S1 & S0)))
    S = S0 | S1; G = signed_perms(); nobox = []
    for c in itertools.product(range(10), repeat=6):
        if not any(tuple((9 - c[i]) if r else c[i] for i, r in g) in S for g in G): nobox.append(c)
    print('orbit coverage: 1/10 boxes of [0,1]^6 with no image in a verified start box:', len(nobox), nobox[:5])
    ok = (not un0 and not un1 and not bad and S0 == NEED and S1 == EXT and not d0 and not d1 and not (S0 & S1)
          and not nobox)
    print('COMPLETE AND VALID (symmetry-correct coverage)' if ok else 'NOT (YET) COMPLETE/VALID')

if __name__ == '__main__':
    sys.setrecursionlimit(10000)
    main()

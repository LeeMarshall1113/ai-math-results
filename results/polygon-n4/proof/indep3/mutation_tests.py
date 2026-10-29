"""Mutation tests for check3.py: every corrupted input below must be REJECTED.

  python mutation_tests.py

A mutant counts as rejected if the replay returns False or raises (malformed input).
Each test also runs the unmutated input first and requires it to be ACCEPTED, so a
rejection is attributable to the mutation.
"""
import copy, json, os, sys, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import poly as PL
import check3 as C
from bernstein import bern_nonneg

RESULTS = []


def run(fn):
    try:
        return bool(fn())
    except (C.Bad, KeyError, TypeError, ValueError, AssertionError, IndexError, AttributeError):
        return False


def case(name, good, bad):
    t = time.time()
    g = run(good) if good is not None else True
    b = run(bad)
    ok = g and not b
    RESULTS.append((name, ok))
    print('%-72s original %-8s mutant %-8s %s (%.1fs)' % (name[:72], 'accepted' if g else 'REJECTED',
          'ACCEPTED' if b else 'rejected', 'PASS' if ok else '** TEST FAILED **', time.time() - t), flush=True)


def eps_shift(N, num=1, den=10**6):
    """den*N - num, a positive multiple of N - num/den."""
    return PL.add(PL.scale(N, den), PL.const(-num))


def find(node, pred, path=()):
    """Paths to nodes of a proof tree satisfying pred."""
    out = []
    if not isinstance(node, dict):
        return out
    if pred(node):
        out.append(path)
    t = node.get('t')
    if 'p' in node and isinstance(node['p'], dict):
        out += find(node['p'], pred, path + ('p',))
    if t == 'factor':
        for i, kd in enumerate(node['odd']):
            out += find(kd['p'], pred, path + ('odd', i, 'p'))
    if t in ('split', 'splitat', 'simplex'):
        for i, kd in enumerate(node['kids']):
            out += find(kd, pred, path + ('kids', i))
    if t == 'blowup':
        for i, kd in enumerate(node['kids']):
            out += find(kd['p'], pred, path + ('kids', i, 'p'))
    return out


def at(node, path):
    for p in path:
        node = node[p]
    return node


def main():
    t0 = time.time()
    Phi = C.check_algebra()
    faces = C.all_faces()
    orb, norb = C.orbits(faces)
    easy, trees, tubes = C.load()

    # ---- 0. Bernstein routine on known polynomials
    t = PL.var(0)
    sq = PL.pw(PL.add(PL.scale(t, 2), PL.const(-1)), 2)            # (2t-1)^2 >= 0, zero at 1/2
    case('bern: (2t0-1)^2 - 1/100 (negative near 1/2)', lambda: bern_nonneg(sq),
         lambda: bern_nonneg(PL.add(PL.scale(sq, 100), PL.const(-1))))
    xy = PL.pw(PL.sub(PL.var(0), PL.var(1)), 2)
    # (t0-t1)^2 itself vanishes on a whole diagonal, which bisection alone never certifies
    # (that is what the blowup move is for), so the accepted original is (t0-t1)^2 + 1/100.
    case('bern: (t0-t1)^2 - 1/10^4 vs (t0-t1)^2 + 1/100', lambda: bern_nonneg(PL.add(PL.scale(xy, 100), PL.const(1))),
         lambda: bern_nonneg(PL.add(PL.scale(xy, 10**4), PL.const(-1))))

    # ---- 1. algebra
    Phi_bad = dict(Phi); k = max(Phi_bad); Phi_bad[k] += 1
    case('identity (*) with one coefficient of Phi changed', lambda: C.star_identity(Phi),
         lambda: C.star_identity(Phi_bad))
    Phi_sgn = dict(Phi); Phi_sgn[0] = -Phi_sgn[0]
    case('identity (*) with the constant term of Phi negated', lambda: C.star_identity(Phi),
         lambda: C.star_identity(Phi_sgn))

    # ---- 3. face polynomial cross-check
    import random
    f0 = sorted(easy)[0]
    N0 = C.face_poly(f0, Phi)
    Nb = dict(N0); kk = sorted(Nb)[len(Nb) // 2]; Nb[kk] += 1
    case('N_F with one coefficient off by 1 vs Faddeev-LeVerrier', lambda: (C.numeric_charpoly_check(f0, N0, Phi, random.Random(1), 3), True)[1],
         lambda: (C.numeric_charpoly_check(f0, Nb, Phi, random.Random(1), 3), True)[1])

    # ---- 4a. Bernstein faces and tight faces
    case('easy face %s: bern on N - 10^6 (huge shift)' % f0, lambda: bern_nonneg(N0),
         lambda: bern_nonneg(eps_shift(N0, 10**12)))
    # pick small tree faces (time the originals)
    timed = []
    for f, lst in sorted(trees.items()):
        N = C.face_poly(f, Phi)
        ts = time.time()
        ok = C.replay(N, lst[0][1])
        timed.append((time.time() - ts, f, lst[0][1], N, ok))
        if len(timed) >= 60:
            break
    timed.sort(key=lambda x: x[0])
    small = [x for x in timed if x[4]]
    _, fT, trT, NT, _ = small[0]
    case('tight face %s: plain bern on N' % fT, lambda: C.replay(NT, trT), lambda: bern_nonneg(NT))
    case('tight face %s: valid tree replayed against N - 10^-6' % fT, lambda: C.replay(NT, trT),
         lambda: C.replay(eps_shift(NT), trT))
    case('tight face %s: valid tree replayed against N - 10^-12' % fT, lambda: C.replay(NT, trT),
         lambda: C.replay(eps_shift(NT, 1, 10**12), trT))

    def mutate_first(ntype, mut, label, maxn=3):
        """Apply `mut` to the first node of type ntype in the smallest trees that have one."""
        done = 0
        for _, f, tr, N, _ in small:
            paths = find(tr, lambda n: n.get('t') == ntype)
            if not paths:
                continue
            tr2 = copy.deepcopy(tr)
            node = at(tr2, paths[0])
            desc = mut(node)
            if desc is None:
                continue
            case('%s %s: %s' % (label, f, desc), lambda: C.replay(N, tr), lambda: C.replay(N, tr2))
            done += 1
            if done >= maxn:
                return

    def flip_mono(n):
        n['at'] = 1 - n['at']; return 'mono direction flipped (at=%d)' % n['at']
    mutate_first('mono', flip_mono, 'tree')

    def mono_var(n):
        n['var'] = (n['var'] + 1) % 6; return 'mono variable changed to t%d' % n['var']
    mutate_first('mono', mono_var, 'tree', 2)

    def fac_sign(n):
        n['odd'][0]['s'] = -n['odd'][0]['s']; return 'sign s of odd factor %s.. flipped' % n['odd'][0]['f'][:20]
    mutate_first('factor', fac_sign, 'tree')

    def fac_coef(n):
        for kd in n['odd']:
            if len(kd['f']) > 6:
                kd['f'] = kd['f'] + ' + 1'
                return 'odd factor %s.. altered (+1)' % kd['f'][:20]
        return None
    mutate_first('factor', fac_coef, 'tree')

    def fac_drop(n):
        if not n['odd']:
            return None
        n['odd'] = n['odd'][1:]; return 'first odd factor dropped from the list'
    mutate_first('factor', fac_drop, 'tree')

    def blow_corner(n):
        j = sorted(n['corner'])[0]; n['corner'][j] = 1 - n['corner'][j]
        return 'blowup corner bit of t%s flipped' % j
    mutate_first('blowup', blow_corner, 'tree')

    def blow_dropkid(n):
        n['kids'] = n['kids'][1:]; return 'blowup cone dropped'
    mutate_first('blowup', blow_dropkid, 'tree', 1)

    def blow_S(n):
        if len(n['S']) < 3:
            return None
        j = n['S'][-1]; n['S'] = n['S'][:-1]; n['corner'].pop(str(j), None)
        n['kids'] = [kd for kd in n['kids'] if kd['m'] != j]
        return 'blowup variable t%d removed from S (cones over a smaller set)' % j
    mutate_first('blowup', blow_S, 'tree', 2)

    def splitat_r(n):
        n['r'] = [n['r'][1], n['r'][0]]; return 'splitat r inverted to %s (outside (0,1))' % n['r']
    mutate_first('splitat', splitat_r, 'tree', 1)

    def splitat_move(n):
        p, q = n['r']; n['r'] = [1, 64] if (p, q) != (1, 64) else [63, 64]
        return 'splitat point moved from %d/%d to %d/%d' % (p, q, n['r'][0], n['r'][1])
    mutate_first('splitat', splitat_move, 'tree', 2)

    def simplex_drop(n):
        n['kids'] = n['kids'][:4]; return 'simplex: one tetrahedron dropped'
    def simplex_corner(n):
        j = sorted(n['corner'])[0]; n['corner'][j] = 1 - n['corner'][j]
        return 'simplex corner bit of t%s flipped' % j
    # simplex nodes are rare: look in all trees, not just the small ones
    sx = None
    for f, lst in sorted(trees.items()):
        for fn, tr in lst:
            if find(tr, lambda n: n.get('t') == 'simplex'):
                sx = (f, tr); break
        if sx:
            break
    if sx:
        f, tr = sx
        N = C.face_poly(f, Phi)
        p = find(tr, lambda n: n.get('t') == 'simplex')[0]
        for mut in (simplex_drop, simplex_corner):
            tr2 = copy.deepcopy(tr); d = mut(at(tr2, p))
            case('tree %s: %s' % (f, d), lambda: C.replay(N, tr), lambda: C.replay(N, tr2))

    # ---- 4b. tube certificates (use the fastest one, tube_tubeC)
    fU = '**0*01*10**1'
    rec = [r for fn, r in tubes[fU]][0]
    NU = C.face_poly(fU, Phi)
    good = lambda: C.check_tube(NU, rec['tree'], rec['TG'])

    def kd_paths(node, pred, path=()):
        out = []
        if pred(node):
            out.append(path)
        if node.get('t') == 'split':
            for i, kd in enumerate(node['kids']):
                out += kd_paths(kd, pred, path + ('kids', i))
        return out
    case('tube %s: replayed against N - 10^-6' % fU, good, lambda: C.check_tube(eps_shift(NU), rec['tree'], rec['TG']))
    r2 = copy.deepcopy(rec)
    p = kd_paths(r2['tree'], lambda n: n.get('t') == 'reduce')[0]
    nd = at(r2['tree'], p); nd['S'] = [[l, 1 - v] for l, v in nd['S']]
    case('tube %s: pin value of a reduce leaf flipped' % fU, good, lambda: C.check_tube(NU, r2['tree'], r2['TG']))
    r3 = copy.deepcopy(rec)
    key = sorted(r3['TG'])[0]; del r3['TG'][key]
    case('tube %s: TG entry %s deleted' % (fU, key), good, lambda: C.check_tube(NU, r3['tree'], r3['TG']))
    r4 = copy.deepcopy(rec)
    ps = kd_paths(r4['tree'], lambda n: n.get('t') == 'reduce')
    for p in ps:
        at(r4['tree'], p).clear(); at(r4['tree'], p)['t'] = 'bern'
    case('tube %s: all %d reduce leaves replaced by plain bern' % (fU, len(ps)), good,
         lambda: C.check_tube(NU, r4['tree'], r4['TG']))
    r5 = copy.deepcopy(rec)
    p = kd_paths(r5['tree'], lambda n: n.get('t') == 'split')[1]
    at(r5['tree'], p)['kids'] = at(r5['tree'], p)['kids'][:1]
    case('tube %s: one k-d box dropped' % fU, good, lambda: C.check_tube(NU, r5['tree'], r5['TG']))
    r6 = copy.deepcopy(rec)
    sub = r6['TG'][sorted(r6['TG'])[0]]
    if isinstance(sub, dict) and sub.get('t') == 'tube':
        pp = kd_paths(sub['tree'], lambda n: n.get('t') == 'reduce')
        for q in pp:
            at(sub['tree'], q).clear(); at(sub['tree'], q)['t'] = 'bern'
        desc = 'nested TG tube: its %d reduce leaves replaced by bern' % len(pp)
    else:
        r6['TG'][sorted(r6['TG'])[0]] = {'t': 'bern'}
        desc = 'TG tree replaced by plain bern'
    case('tube %s: %s' % (fU, desc), good, lambda: C.check_tube(NU, r6['tree'], r6['TG']))
    # a tube with 'tree' leaves: replace them by bern
    for f, lst in sorted(tubes.items()):
        rr = lst[-1][1]
        ps = kd_paths(rr['tree'], lambda n: n.get('t') == 'tree')
        if ps and len(json.dumps(rr)) < 400000:
            N = C.face_poly(f, Phi)
            rr2 = copy.deepcopy(rr)
            for q in ps:
                at(rr2['tree'], q).clear(); at(rr2['tree'], q)['t'] = 'bern'
            case('tube %s: its %d tree leaves replaced by plain bern' % (f, len(ps)),
                 lambda: C.check_tube(N, rr['tree'], rr['TG']), lambda: C.check_tube(N, rr2['tree'], rr2['TG']))
            break

    # ---- 5. coverage
    ok_cov = lambda e, t, u: (C.coverage(orb, norb, e, t, u), True)[1]
    E = list(easy)
    case('coverage: one Bernstein face removed', lambda: ok_cov(E, trees, tubes), lambda: ok_cov(E[1:], trees, tubes))
    fT1 = sorted(trees)[0]
    case('coverage: a tree face also listed as Bernstein (overlap)', lambda: ok_cov(E, trees, tubes),
         lambda: ok_cov(E + [fT1], trees, tubes))
    other = next(g for g in faces if orb[g] == orb[E[0]] and g != E[0])
    case('coverage: Bernstein face replaced by another face of its orbit, plus the original', lambda: ok_cov(E, trees, tubes),
         lambda: ok_cov(E + [other], trees, tubes))
    case('coverage: a non-face (7 free coordinates) inserted', lambda: ok_cov(E, trees, tubes),
         lambda: ok_cov(E[1:] + ['*******00000'], trees, tubes))
    tb = dict(tubes); k1 = sorted(tb)[0]; del tb[k1]
    case('coverage: one tube face removed', lambda: ok_cov(E, trees, tubes), lambda: ok_cov(E, trees, tb))
    three = '***' + '*00' + '*00' + '*00'      # 6 free, but all 3 weights of row 0 free
    case('face rule: a face with 3 free weights in row 0 is not among the 38016 faces',
         lambda: True, lambda: three in set(faces))

    npass = sum(ok for _, ok in RESULTS)
    print('\nmutation tests: %d/%d passed (every mutant rejected, every original accepted)  %.0fs'
          % (npass, len(RESULTS), time.time() - t0))
    return 0 if npass == len(RESULTS) else 1


if __name__ == '__main__':
    sys.exit(main())

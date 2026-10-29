"""Independent replay of tube certificates (tube_*.jsonl), using verify.py's own
face polynomial (sympy), Bernstein check and proof-tree replay.

For a certificate {face, tree, TG}:
  * every TG[S] tree proves N|_{t_l = v, (l,v) in S} >= 0 on [0,1]^6 (verify.replay);
  * the k-d tree partitions [0,1]^6 (bisection at midpoints);
  * leaf 'bern': N on the box is >= 0 (verify.bern_nonneg on N o box map);
  * leaf 'reduce' S: for (l,v) in S in order, the box touches t_l = v, and the derivative
    of the current polynomial in t_l has the sign making it decrease towards t_l = v on the
    box (verify.bern_nonneg); then N(t) >= N(pi_S(t)) >= 0 by TG[S].
"""
import json, sys, glob
from fractions import Fraction
import sympy as sp
import flint
import verify


def check_tube(N, tree, TG, ctx):
    g = ctx.gens(); names = ctx.names()
    for key, sub in TG.items():
        S = json.loads(key)
        PG = N.subs({names[l]: v for l, v in S})
        if isinstance(sub, dict) and sub.get('t') == 'tube':
            if not check_tube(PG, sub['tree'], sub['TG'], ctx):
                print('TG tube fail', key); return False
        elif not verify.replay(PG, sub, ctx):
            print('TG fail', key); return False
    def walk(node, lo, hi):
        if node['t'] == 'split':
            l = node['var']; mid = (lo[l] + hi[l])/2
            a = walk(node['kids'][0], lo, hi[:l] + [mid] + hi[l+1:])
            b = walk(node['kids'][1], lo[:l] + [mid] + lo[l+1:], hi)
            return a and b and len(node['kids']) == 2
        args = [flint.fmpq(lo[i].numerator, lo[i].denominator) +
                flint.fmpq((hi[i]-lo[i]).numerator, (hi[i]-lo[i]).denominator)*g[i] for i in range(6)]
        P = N.compose(*args)
        if node['t'] == 'bern':
            return verify.bern_nonneg(P, ctx)
        if node['t'] == 'tree':
            return verify.replay(P, node['p'], ctx)
        if node['t'] == 'reduce':
            key = json.dumps(sorted([list(x) for x in node['S']]))
            if key not in TG:
                return False
            Q = P
            for l, v in node['S']:
                if (v == 0 and lo[l] != 0) or (v == 1 and hi[l] != 1):
                    return False
                D = Q.derivative(names[l])
                if not verify.bern_nonneg(D if v == 0 else -D, ctx):
                    return False
                Q = Q.subs({names[l]: 0 if v == 0 else 1})
            return True
        return False
    return walk(tree, [Fraction(0)]*6, [Fraction(1)]*6)

def main(files):
    Phi, abc = verify.check_identity()
    ctx = flint.fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(6)), 'lex')
    g = ctx.gens(); names = ctx.names()
    allok = True; done = set()
    for fn in files:
        for line in open(fn):
            r = json.loads(line)
            if not r.get('ok'):
                continue
            face = r['face']
            verify.validate_tube(r)
            Ps = verify.face_poly(face, Phi, abc)
            N = ctx.from_dict({tuple(m): flint.fmpq(int(sp.Rational(v).p), int(sp.Rational(v).q)) for m, v in Ps.terms()})
            ok = check_tube(N, r['tree'], r['TG'], ctx)
            print(face, 'OK' if ok else 'FAIL', flush=True)
            allok &= ok
            if ok: done.add(face)
    assert len(done) == 10 or len(files) < 5, 'expected 10 tube certificates'
    print('tube certificates verified:', len(done), 'all ok' if allok else 'SOME FAILED')
    if not allok:
        sys.exit(1)

if __name__ == '__main__':
    main(sys.argv[1:] or sorted(glob.glob('tube_*.jsonl')))

"""Run verify.py's independent checks on all Bernstein-only faces and on the first
NHARD blow-up trees (search output in hard[0-9]_*.jsonl / kd_*.jsonl)."""
import json, glob, sys, time, flint, sympy as sp
import verify
nhard = int(sys.argv[1])
Phi, abc = verify.check_identity(); print('identity OK', flush=True)
allf = verify.faces_all(); can = sorted(set(verify.canon(f) for f in allf))
print('faces', len(allf), 'canonical', len(can), flush=True)
assert can == sorted(l.strip() for l in open('faces_canon.txt') if l.strip())
assert len(can) == 836
ctx = flint.fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(6)), 'lex')
def poly(f):
    Ps = verify.face_poly(f, Phi, abc)
    return ctx.from_dict({tuple(m): flint.fmpq(int(sp.Rational(v).p), int(sp.Rational(v).q)) for m, v in Ps.terms()})
easy, trees, tubes = verify.load_certificates(can)
easy = sorted(easy)
bad = 0
for f in easy:
    if not verify.replay(poly(f), {'t': 'bern'}, ctx): bad += 1; print('FAIL easy', f, flush=True)
print('easy faces verified:', len(easy) - bad, 'of', len(easy), flush=True)
print('hard trees available:', len(trees), flush=True)
nb = 0; n = 0
for f, tr in list(trees.items())[:nhard]:
    t = time.time(); ok = verify.replay(poly(f), tr, ctx); n += 1; nb += (not ok)
    print('hard', f, ok, round(time.time() - t, 1), flush=True)
print('hard trees verified:', n - nb, 'of', n)
if bad or nb:
    sys.exit(1)

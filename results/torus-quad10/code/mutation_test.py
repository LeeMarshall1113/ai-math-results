"""Mutation tests: geomcheck must reject corrupted tori; encoders' clauses must be violated by corrupted chirotopes."""
import copy, sys
from geomcheck import check_torus
from soundness_test import T9, T11
import encode
bad = 0
def expect_reject(name, W):
    global bad
    r = check_torus(W['vertices'], W['faces'])
    print(f"{name}: {'rejected: ' + r['err'] if not r['valid'] else 'ACCEPTED (BUG)'}")
    bad += r['valid']
d = copy.deepcopy(T9); d['vertices'][0][2] += 1; expect_reject("warped face", d)
d = copy.deepcopy(T9)
for i in range(3, 6): d['vertices'][i][0] *= -1; d['vertices'][i][1] *= -1
expect_reject("reflected ring (intersecting)", d)
d = copy.deepcopy(T11); d['faces'][0] = [0, 6, 1, 2]; expect_reject("bow-tie face order", d)
d = copy.deepcopy(T11); d['vertices'][10] = [47, 3, 300]; expect_reject("moved vertex", d)
d = copy.deepcopy(T11); d['faces'][3] = d['faces'][3][::-1]; r = check_torus(d['vertices'], d['faces'])
print("reversed face cycle (must stay valid):", r['valid']); bad += not r['valid']
# chirotope mutation: flip one nonzero bracket of T11 used in an AD clause -> must violate the CNF
E, _ = encode.build(faces=T11['faces'])
chi = encode.chirotope(T11['vertices'])
if not encode.assignment(E, chi)[E.lit(E.sym, 1)]: chi = {k: -v for k, v in chi.items()}
T = tuple(sorted(E.sym)); chi2 = dict(chi); chi2[T] = 0
v = len(encode.violated(E, encode.assignment(E, chi2)))
print("zeroed AD bracket -> violated clauses:", v); bad += (v == 0)
f = tuple(sorted(T11['faces'][0])); chi3 = dict(chi); chi3[f] = 1
v = len(encode.violated(E, encode.assignment(E, chi3)))
print("non-planar face bracket -> violated clauses:", v); bad += (v == 0)
print("MUTATION TESTS PASSED" if not bad else "MUTATION TEST FAILURE"); sys.exit(1 if bad else 0)

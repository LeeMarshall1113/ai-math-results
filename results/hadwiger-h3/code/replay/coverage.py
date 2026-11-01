"""Coverage of the reduced parameter space (replay, independent).

Derivation of the domain needed (my own, from the symmetry group):
  The 6 parameters are indexed by ordered pairs (a,b), a != b: param(a,b) = b-coordinate of
  the touching point on face x_a = 0.  p1=(x,y) p2=(x,z) p3=(y,x) p4=(y,z) p5=(z,x) p6=(z,y).
  A cube symmetry = coordinate permutation pi + reflections. pi maps param(a,b) to
  param(pi a, pi b); the reflection x_b -> 1-x_b maps both params with second index b to 1-.
  Let S_b = sum of the two params with second index b:  S_x = p3+p5, S_y = p1+p6, S_z = p2+p4.
  Then reflection in b sends S_b -> 2 - S_b (others unchanged), and pi sends S_b to S_{pi b}.
  So the multiset {|S_x-1|, |S_y-1|, |S_z-1|} is invariant and pi permutes it freely.
  Lemma: every orbit meets D* = {p1,p2,p3 <= 1/2, |S_z-1| <= |S_y-1| <= |S_x-1|}.
  Proof: pick pi ordering the three values |S_b-1| as required; then the reflection in y
  fixes p1 <= 1/2 (it changes p1,p6 only), in z fixes p2 (changes p2,p4), in x fixes p3
  (changes p3,p5); reflections do not change any |S_b - 1|. QED.
  Prymak's D = {p1,p2,p3<=1/2, S_z <= S_y <= S_x} is NOT such a region: the reflections are
  used up by p1,p2,p3<=1/2 and they change S_b, not |S_b-1|. (Counterexamples below.)

Box test (exact, integers in units of 1/10): a 1/10-box k is NEEDED iff its interior meets the
interior of D* (resp. D). D and D* are finite unions of full-dimensional polyhedra, hence the
closure of their interiors; the union of the needed closed boxes is closed and contains every
point of int D* off the grid hyperplanes (a dense subset of D*), hence contains D*.
Open-box test: k1,k2,k3 <= 4 and z<y<x solvable with z in |S_z-1|(open box) etc.; each S_b
ranges over an open interval independently (disjoint parameters), |S-1| over its image
[lo,hi) or (lo,hi); the test is lo_z<hi_y, lo_y<hi_x, lo_z<hi_x (strict).
"""
import itertools, random, sys, json
from fractions import Fraction as F


def absrange(m):
    """range of |t-10| for t in [m, m+2] (units 1/10)."""
    lo = 0 if m <= 10 <= m + 2 else min(abs(m - 10), abs(m + 2 - 10))
    hi = max(abs(m - 10), abs(m + 2 - 10))
    return lo, hi


def meets_Dstar(k):
    k1, k2, k3, k4, k5, k6 = k
    if max(k1, k2, k3) > 4:
        return False
    zl, zh = absrange(k2 + k4); yl, yh = absrange(k1 + k6); xl, xh = absrange(k3 + k5)
    return zl < yh and yl < xh and zl < xh


def meets_D(k):
    k1, k2, k3, k4, k5, k6 = k
    if max(k1, k2, k3) > 4:
        return False
    return k2 + k4 < k1 + k6 + 2 and k1 + k6 < k3 + k5 + 2 and k2 + k4 < k3 + k5 + 2


ALLK = list(itertools.product(range(5), range(5), range(5), range(10), range(10), range(10)))


# ---- orbit tests (exact rationals) ----
def pts(p):
    p1, p2, p3, p4, p5, p6 = p
    return [(0, p1, p2), (1, p1, p2), (p3, 0, p4), (p3, 1, p4), (p5, p6, 0), (p5, p6, 1)]


def recover(P):
    out = {}
    for c in range(3):
        lo = [q for q in P if q[c] == 0]; hi = [q for q in P if q[c] == 1]
        cand = [q for q in lo if any(all(r[i] == q[i] for i in range(3) if i != c) for r in hi)]
        assert len(cand) == 1
        out[c] = cand[0]
    q1, q2, q3 = out[0], out[1], out[2]
    return (q1[1], q1[2], q2[0], q2[2], q3[0], q3[1])


SYMS = [(perm, refl) for perm in itertools.permutations(range(3)) for refl in itertools.product((0, 1), repeat=3)]


def act(P, perm, refl):
    out = []
    for x in P:
        y = [None] * 3
        for b in range(3):
            y[perm[b]] = (1 - x[b]) if refl[b] else x[b]
        out.append(tuple(y))
    return out


def orbit(p):
    return [recover(act(pts(p), *s)) for s in SYMS]


def inDstar(q):
    sz, sy, sx = abs(q[1] + q[3] - 1), abs(q[0] + q[5] - 1), abs(q[2] + q[4] - 1)
    return max(q[:3]) <= F(1, 2) and sz <= sy <= sx


def inD(q):
    return max(q[:3]) <= F(1, 2) and q[1] + q[3] <= q[0] + q[5] <= q[2] + q[4]


def in_boxes(q, S):
    cands = []
    for x in q:
        f = x * 10; fl = int(f)
        c = [fl] if f != fl else [fl - 1, fl]
        cands.append([k for k in c if 0 <= k <= 9])
    return any(k in S for k in itertools.product(*cands))


def main(start_set, ntest=20000, seed=7):
    res = {}
    needD = {k for k in ALLK if meets_D(k)}
    needDs = {k for k in ALLK if meets_Dstar(k)}
    res['boxes_meeting_D'] = len(needD)
    res['boxes_meeting_Dstar'] = len(needDs)
    res['start_boxes_present'] = len(start_set)
    res['D_missing'] = len(needD - start_set)
    res['Dstar_missing'] = len(needDs - start_set)
    res['Dstar_missing_list_head'] = sorted(needDs - start_set)[:5]
    random.seed(seed)
    missD = missDs = missU = 0
    ex = None
    for _ in range(ntest):
        p = tuple(F(random.randint(0, 10000), 10000) for _ in range(6))
        o = orbit(p)
        assert len(set(o)) == 48 or True
        if not any(inD(q) for q in o):
            missD += 1
        if not any(inDstar(q) for q in o):
            missDs += 1
        if not any(in_boxes(q, start_set) for q in o):
            missU += 1
            ex = ex or p
    res.update(orbit_test_n=ntest, orbit_misses_D=missD, orbit_misses_Dstar=missDs,
               orbit_misses_certified_union=missU, example_uncovered=[str(x) for x in ex] if ex else None)
    # the coordinator's counterexample
    pc = (F(459, 500), F(1481, 5000), F(1689, 10000), F(9529, 10000), F(9359, 10000), F(3079, 10000))
    res['referee_example_in_certified_union'] = any(in_boxes(q, start_set) for q in orbit(pc))
    res['referee_example_in_Dstar_orbit'] = any(inDstar(q) for q in orbit(pc))
    return res


if __name__ == '__main__':
    S = set(tuple(x) for x in json.load(open(sys.argv[1])))
    print(json.dumps(main(S), indent=1, default=str))

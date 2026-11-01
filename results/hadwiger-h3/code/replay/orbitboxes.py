"""Exact grid-orbit coverage: cube symmetries act on the 1/10-box grid of [0,1]^6 (each acts
as a signed coordinate permutation p'_m = p_s(m) or 1-p_s(m)). The certified region covers all
configurations up to symmetry iff every orbit of 1/10-boxes contains a certified start box.
usage: python orbitboxes.py startboxes.json [Dstar-only]"""
import itertools, json, sys
from fractions import Fraction as F
import coverage as C

def signed_perms():
    base = tuple(F(v, 97) for v in (3, 5, 7, 11, 13, 17))  # distinct, and distinct from 1-x
    out = []
    for s in C.SYMS:
        img = C.recover(C.act(C.pts(base), *s))
        m = []
        for y in img:
            if y in base: m.append((base.index(y), 0))
            else: m.append((base.index(1 - y), 1))
        out.append(m)
    return out

def main(S):
    G = signed_perms()
    def img(k, g): return tuple((9 - k[i]) if r else k[i] for i, r in g)
    seen = set(); missing_orbits = []
    for k in itertools.product(range(10), repeat=6):
        if k in seen: continue
        orb = {img(k, g) for g in G}
        seen |= orb
        if not (orb & S):
            missing_orbits.append(sorted(orb))
    reps_in_Dstar = sum(1 for o in missing_orbits if any(C.meets_Dstar(x) for x in o))
    reps = [min(x for x in o if C.meets_Dstar(x)) for o in missing_orbits if any(C.meets_Dstar(x) for x in o)]
    return dict(missing_Dstar_reps=reps, n_missing_orbits=len(missing_orbits), missing_reps_meeting_Dstar=reps_in_Dstar,
                example=missing_orbits[0][0] if missing_orbits else None,
                )

if __name__ == '__main__':
    S = set(tuple(x) for x in json.load(open(sys.argv[1])))
    print(json.dumps(main(S)))

"""check_explicit.py -- the two groups behind the exclusion of types 2, 4, 6, 8, written out by hand.

The certificates in ../certs use permutation generators found by GAP (classes 1132 and 1279 of the list of
conjugacy classes of subgroups of S_10).  This script builds the same two groups from the finite-field
descriptions given in the note (Section 5), checks their orders, their edge orbits and the Oliver-class witness
(normal p-subgroup N, Gamma = <g> N), writes four certificates in the format of ../certs, and replays them with
the standard-library checker replay_cert.py.  Nothing here uses GAP, nauty or a solver.

Vertices: 0 and 1 are the two hubs (the fixed edge e = {0,1}); the other eight vertices carry the field structure.
  Gamma_1 = AGL(1,8) x <(0 1)>          on W = F_8,        vertex of x is 2 + x.
  Gamma_2 = (F_2 x F_4 affine group) x <(0 1)>   on W = F_2 x F_4, vertex of (t,z) is 2 + 4t + z.

usage:  python note/check_explicit.py          (from the karp10 directory; about 20 seconds)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "replay"))
from replay_lib import closure, edge_orbits, check_G1, popcount, IDENT  # noqa: E402
import replay_cert  # noqa: E402


# ---- F_8 = F_2[x]/(x^3 + x + 1): elements 0..7 as bit vectors, addition = xor
def mul8(a, b):
    r = 0
    for i in range(3):
        if (b >> i) & 1:
            r ^= a << i
    for i in (4, 3):                       # reduce x^4, x^3
        if (r >> i) & 1:
            r ^= 0b1011 << (i - 3)
    return r


# ---- F_4 = F_2[w]/(w^2 + w + 1): elements 0,1,w=2,w+1=3
MUL4 = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]


def perm(f):
    """permutation of 0..9 from a function on vertices"""
    p = tuple(f(v) for v in range(10))
    assert sorted(p) == list(range(10))
    return p


HUBSWAP = perm(lambda v: {0: 1, 1: 0}.get(v, v))


def on_w8(h):                                  # h: F_8 -> F_8 acting on vertices 2..9
    return perm(lambda v: v if v < 2 else 2 + h(v - 2))


def on_w42(h):                                 # h: (t,z) -> (t',z') acting on vertices 2..9
    def f(v):
        if v < 2:
            return v
        t, z = divmod(v - 2, 4)
        t2, z2 = h(t, z)
        return 2 + 4 * t2 + z2
    return perm(f)


def group1():
    alpha = 2                                  # x, a primitive element of F_8
    assert len({alpha, mul8(alpha, alpha), mul8(mul8(alpha, alpha), alpha)}) == 3
    trans = [on_w8(lambda x, c=c: x ^ c) for c in (1, 2, 4)]
    g = on_w8(lambda x: mul8(alpha, x))
    gens = [trans[0], g, HUBSWAP]
    ngens = trans + [HUBSWAP]
    return dict(gens=gens, order=112, p=2, N_gens=ngens, g=g)


def group2():
    omega = 2
    trans = [on_w42(lambda t, z, s=s, w=w: (t ^ s, z ^ w)) for (s, w) in ((1, 0), (0, 1), (0, 2))]
    shear = [on_w42(lambda t, z, b=b: (t, z ^ (MUL4[b][1] if t else 0))) for b in (1, 2)]   # z -> z + b t
    g = on_w42(lambda t, z: (t, MUL4[omega][z]))                                            # z -> w z
    gens = trans + shear + [g, HUBSWAP]
    ngens = trans + shear + [HUBSWAP]
    return dict(gens=gens, order=192, p=2, N_gens=ngens, g=g)


def main():
    out = os.path.join(HERE, "explicit")
    os.makedirs(out, exist_ok=True)
    groups = {"GAMMA1_AGL18xC2": group1(), "GAMMA2_order192": group2()}
    for name, gd in groups.items():
        order, n = check_G1(gd["gens"], gd["order"], gd["p"], gd["N_gens"], gd["g"])
        orbs = edge_orbits(gd["gens"])
        sizes = sorted(popcount(o) for o in orbs)
        print(f"{name}: |Gamma| = {order}, |N| = {n} (normal {gd['p']}-subgroup), Gamma = <g> N, edge orbits {sizes}")
    assert sorted(popcount(o) for o in edge_orbits(groups["GAMMA1_AGL18xC2"]["gens"])) == [1, 16, 28]
    assert sorted(popcount(o) for o in edge_orbits(groups["GAMMA2_order192"]["gens"])) == [1, 12, 16, 16]
    # the hub edge {0,1} must be a singleton orbit of both groups (replay_cert checks this again)
    for name in groups:
        o = edge_orbits(groups[name]["gens"])
        assert any(popcount(m) == 1 and m == 1 for m in o), "edge {0,1} is not fixed"   # edge (0,1) has index 0

    def js(gd):
        return dict(gens=[list(p) for p in gd["gens"]], order=gd["order"], p=gd["p"],
                    N_gens=[list(p) for p in gd["N_gens"]], g=list(gd["g"]))

    G1, G2 = groups["GAMMA1_AGL18xC2"], groups["GAMMA2_order192"]
    certs = {
        "type2_explicit.json": dict(type=2, claim="type 2 excluded by Gamma_1 and Gamma_2 (explicit finite-field groups)",
                                    groups={"GAMMA1": js(G1), "GAMMA2": js(G2)},
                                    equations=[dict(group="GAMMA1", kind="lk", fixed_edge=[0, 1]),
                                               dict(group="GAMMA2", kind="lk", fixed_edge=[0, 1])]),
        "type8_explicit.json": dict(type=8, claim="type 8 excluded by Gamma_1 and Gamma_2 (explicit finite-field groups)",
                                    groups={"GAMMA1": js(G1), "GAMMA2": js(G2)},
                                    equations=[dict(group="GAMMA1", kind="lk", fixed_edge=[0, 1]),
                                               dict(group="GAMMA2", kind="lk", fixed_edge=[0, 1])]),
        "type4_explicit.json": dict(type=4, claim="type 4 excluded by Gamma_2 (explicit finite-field group)",
                                    groups={"GAMMA2": js(G2)},
                                    equations=[dict(group="GAMMA2", kind="lk", fixed_edge=[0, 1])]),
        "type6_explicit.json": dict(type=6, claim="type 6 excluded by Gamma_2 (explicit finite-field group)",
                                    groups={"GAMMA2": js(G2)},
                                    equations=[dict(group="GAMMA2", kind="lk", fixed_edge=[0, 1])]),
    }
    rc_all = 0
    for fn, cert in certs.items():
        path = os.path.join(out, fn)
        json.dump(cert, open(path, "w"))
        print(f"--- replay {fn}")
        rc = replay_cert.main(path, verbose=("type4" in fn))
        rc_all |= rc
    print("EXPLICIT GROUPS: ALL VALID" if rc_all == 0 else "EXPLICIT GROUPS: FAILED")
    return rc_all


if __name__ == "__main__":
    sys.exit(main())

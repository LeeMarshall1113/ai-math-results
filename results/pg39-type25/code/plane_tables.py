"""Plane-section tables for an involution g of PG(3,9): the sections of K with g-invariant planes are g|plane-invariant planar
(2,5)-sets (size 26 or 35); we list all of them (SAT enumeration in the standard plane with the standard planar involution) and
transport them to each invariant plane through a linear map A with A h0 A^-1 = g|plane.

Used for 2b (planes <l1,p2> and <l2,p1>, standard involution h0 = diag(1,1,-1)) and for 2d (baer3d.py does the same for the Baer case).
"""
import itertools
from pysat.solvers import Cadical195

from ff import GF
from pg import PG
from orbsat import OrbitInstance


def invariant_family(P2, M, frob):
    """all (2,5)-sets of size 26/35 of PG(2,9) invariant under the collineation (M, frob); {k: sorted tuples}"""
    perm = P2.perm_of_matrix(M, frob)
    inst = OrbitInstance(P2, [perm])
    out = {}
    for k in (26, 35):
        cnf = inst.build_cnf(k, (2, 5), amo=True)
        s = Cadical195(bootstrap_with=cnf.clauses)
        sols = []
        while s.solve():
            m = s.get_model()
            val = {abs(l): l > 0 for l in m}
            sols.append(tuple(sorted(p for i, o in enumerate(inst.orbs) if val[inst.z[i]] for p in o)))
            s.add_clause([(-inst.z[i] if val[inst.z[i]] else inst.z[i]) for i in range(len(inst.orbs))])
        s.delete()
        out[k] = sorted(sols)
    return out, inst


def linear_map_points(P2, P3, cols):
    """images of the points of PG(2,9) under x -> x1*cols[0] + x2*cols[1] + x3*cols[2] (vectors of GF(9)^4)"""
    F = P3.F
    mp = []
    for v in P2.points:
        w = (0, 0, 0, 0)
        for c, b in zip(v, cols):
            w = tuple(F.add[x][F.mul[c][y]] for x, y in zip(w, b))
        mp.append(P3.pidx[P3.normalize(w)])
    return mp


def allowed_rows(family, mp, fixed2, want_fixed):
    """transported sections whose fixed part (points of the standard plane in `fixed2`) maps onto `want_fixed`"""
    rows = []
    for k in (26, 35):
        for X in family[k]:
            if set(mp[p] for p in X if p in fixed2) == want_fixed:
                rows.append(tuple(sorted(mp[p] for p in X)))
    return sorted(rows)


def add_plane_tables(cnf, inst3, planes_rows):
    """planes_rows: list of (plane points, rows); selector clauses over the orbit variables inside the plane"""
    z = inst3.z
    for pts, rows in sorted(planes_rows, key=lambda t: tuple(sorted(t[0]))):
        ps = set(pts)
        orbs_in = sorted(set(inst3.oid[p] for p in pts))
        for o in orbs_in:
            assert all(p in ps for p in inst3.orbs[o])
        sels = []
        for r in rows:
            sv = cnf.new_var()
            sels.append(sv)
            rs = set(r)
            for o in orbs_in:
                cnf.add(-sv, z[o] if inst3.orbs[o][0] in rs else -z[o])
        cnf.add(*sels)

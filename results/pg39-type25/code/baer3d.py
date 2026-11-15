"""Class 2d (Baer involution) with plane-table constraints.

Every sigma-invariant plane pi of PG(3,9) is one of the 40 planes of PG(3,3) (Galois descent); K n pi is a sigma-invariant
planar (2,5)-set of size 26 or 35, i.e. one of the 3240 + 4212 sets Z of the standard plane (complete list: SAT enumeration,
LRAT certificate certs/baer_planar_complete_k*.lrat).  Z is invariant under PGL(3,3) (which commutes with sigma), so the
sections allowed in pi are the images of Z under any GF(3)-linear collineation standard plane -> pi.  With S = K n PG(3,3)
fixed (3 cases, baer_fixed.py), the allowed sections of pi are those whose fixed part equals S n pi.

CNF of a case = orbit CNF (class 2d, |K| = k) + S units + for each of the 40 planes a selector clause over its allowed sections.
"""
import json
import os
import itertools
from pysat.solvers import Cadical195

from ff import GF
from pg import PG
from orbsat import OrbitInstance

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def planar_Z(P2=None):
    """all sigma-invariant (2,5)-sets of size 26 / 35 of the standard plane; returns {k: [sorted tuples of P2 point indices]}"""
    F = GF(9)
    P2 = P2 or PG(3, 9, F)
    perm = P2.perm_of_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1)
    inst = OrbitInstance(P2, [perm])
    out = {}
    for k in (26, 35):
        cnf = inst.build_cnf(k, (2, 5), amo=True)
        s = Cadical195(bootstrap_with=cnf.clauses)
        sols = []
        while s.solve():
            m = s.get_model()
            val = {abs(l): l > 0 for l in m}
            X = tuple(sorted(p for i, o in enumerate(inst.orbs) if val[inst.z[i]] for p in o))
            sols.append(X)
            s.add_clause([(-inst.z[i] if val[inst.z[i]] else inst.z[i]) for i in range(len(inst.orbs))])
        s.delete()
        out[k] = sorted(sols)
    return out, inst, P2


def gf3_planes(P3):
    """the 40 planes of PG(3,3): (normal, sorted P3 indices of the 91 points, sorted indices of the 13 fixed points, GF(3)-basis)"""
    F = P3.F
    B = [i for i, v in enumerate(P3.points) if all(c in (0, 1, 2) for c in v)]
    out = []
    for a in P3.points:
        if not all(c in (0, 1, 2) for c in a):
            continue
        pts = [i for i, v in enumerate(P3.points) if P3._dot(a, v) == 0]
        fixed = [i for i in pts if i in set(B)]
        assert len(pts) == 91 and len(fixed) == 13
        # GF(3)-basis from fixed points
        basis, span = [], {(0, 0, 0, 0)}
        for i in fixed:
            v = P3.points[i]
            if v not in span:
                basis.append(v)
                new = set()
                for s in span:
                    for c in (0, 1, 2):
                        new.add(tuple(F.add[x][F.mul[c][y]] for x, y in zip(s, v)))
                span = new
            if len(basis) == 3:
                break
        out.append((a, pts, fixed, basis))
    return out


def transport(P2, P3, basis):
    mp = []
    F = P3.F
    for v in P2.points:
        w = (0, 0, 0, 0)
        for c, b in zip(v, basis):
            w = tuple(F.add[x][F.mul[c][y]] for x, y in zip(w, b))
        mp.append(P3.pidx[P3.normalize(w)])
    return mp


def plane_tables(Z, inst2, P2, P3, S):
    """for each plane: sorted list of allowed sections (sorted tuples of P3 indices)"""
    fixed2 = set(p for o in inst2.orbs if len(o) == 1 for p in o)
    Sset = set(S)
    tables = []
    for (a, pts, fixed, basis) in gf3_planes(P3):
        mp = transport(P2, P3, basis)
        assert sorted(mp) == sorted(pts)
        want = Sset & set(fixed)
        rows = []
        for k in (26, 35):
            for X in Z[k]:
                img = tuple(sorted(mp[p] for p in X))
                if set(mp[p] for p in X if p in fixed2) == want:
                    rows.append(img)
        tables.append((a, sorted(rows)))
    return tables


def add_tables_full(cnf, inst3, P3, tables):
    z = inst3.z
    planes = {a: pts for (a, pts, fixed, basis) in gf3_planes(P3)}
    for a, rows in sorted(tables, key=lambda t: tuple(sorted(planes[t[0]]))):
        pts = planes[a]
        orbs_in = sorted(set(inst3.oid[p] for p in pts))
        for o in orbs_in:
            assert all(p in set(pts) for p in inst3.orbs[o])
        sels = []
        for r in rows:
            sv = cnf.new_var()
            sels.append(sv)
            rs = set(r)
            for o in orbs_in:
                cnf.add(-sv, z[o] if inst3.orbs[o][0] in rs else -z[o])
        cnf.add(*sels)

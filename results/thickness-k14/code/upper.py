"""Upper bound theta(K14) <= 4, exactly: Seamons' K13 3-layer witness plus a 14th point whose
13 edges form layer 4 (a star: edges sharing an endpoint never cross).  Exact integer check.
Also re-checks the authors' k14_t4*.json witnesses with the referee's own checker.
The drawing is rebuilt in work/ and compared with certs/witnesses/k14_t4_referee.json."""
import json, os, glob, itertools
from controls import check_drawing, orient, accept

HERE = os.path.dirname(os.path.abspath(__file__))
WIT = os.path.join(HERE, "..", "certs", "witnesses")
W = json.load(open(os.path.join(WIT, "k13_seamons.json")))
P = [tuple(p) for p in W["points"]]
layer = {(min(u, v), max(u, v)): L for u, v, L in W["edges"]}
new = (1000003, 777781)
for a, b in itertools.combinations(range(13), 2):
    assert orient(P + [new], a, b, 13) != 0
P14 = P + [new]
for v in range(13):
    layer[(v, 13)] = 4
bad = check_drawing(P14, layer)
print("K14 = Seamons K13 + star point", new, ": layers", sorted(set(layer.values())), " crossings within layers:", bad)
assert bad == 0
out = {"points": P14, "edges": [[u, v, L] for (u, v), L in sorted(layer.items())]}
os.makedirs(os.path.join(HERE, "..", "work"), exist_ok=True)
json.dump(out, open(os.path.join(HERE, "..", "work", "k14_t4_referee.json"), "w"))
pub = json.load(open(os.path.join(WIT, "k14_t4_referee.json")))
print("rebuilt drawing equals certs/witnesses/k14_t4_referee.json:", pub == json.loads(json.dumps(out)))
assert pub == json.loads(json.dumps(out))
for f in sorted(glob.glob(os.path.join(WIT, "k14_t4_S*.json")) + [os.path.join(WIT, "k14_t4.json")]):
    w = json.load(open(f))
    Q = [tuple(p) for p in w["points"]]
    lay = {(min(u, v), max(u, v)): L for u, v, L in w["edges"]}
    bad, acc = check_drawing(Q, lay), accept(Q, lay, 4)
    print(os.path.basename(f), "exact crossings within layers:", bad, " accepted by referee encoding t=4:", acc)
    assert bad == 0 and acc and len(Q) == 14 and set(lay.values()) <= {1, 2, 3, 4}
acc = accept(P14, layer, 4)
print("referee K14 drawing accepted by referee encoding t=4:", acc)
assert acc
print("UPPER BOUND OK: exact 4-layer straight-line drawings of K14")

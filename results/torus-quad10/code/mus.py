"""Group-MUS of the core encoding GP+PL+AD: which GP relations / faces / adjacency conditions
are really needed?  Output feeds the human-readable proof attempt."""
import itertools, sys
from pysat.solvers import Solver
import encode

E = encode.Enc()
groups = []   # (name, clause list)
def grab(name, fn):
    n0 = len(E.cls); fn(); groups.append((name, E.cls[n0:])); del E.cls[n0:]
base_cls = list(E.cls); E.cls.clear()   # chirotope variable axioms (always on)
N, F = encode.N, encode.F
# GP per relation
prod = {}
def P(A, B):
    key = (A, B)
    if key in prod: return prod[key]
    pp, pn = E.new(), E.new()
    ap, an, bp, bn = E.lit(A, 1), E.lit(A, -1), E.lit(B, 1), E.lit(B, -1)
    base_cls.extend([[-ap, -bp, pp], [-an, -bn, pp], [-ap, -bn, pn], [-an, -bp, pn],
                     [-pp, ap, an], [-pp, bp, bn], [-pp, -ap, bp], [-pp, -an, bn],
                     [-pn, ap, an], [-pn, bp, bn], [-pn, -ap, bn], [-pn, -an, bp]])
    prod[key] = (pp, pn); return pp, pn
for s in itertools.combinations(range(N), 2):
    rest = [x for x in range(N) if x not in s]
    for y in itertools.combinations(rest, 4):
        y1, y2, y3, y4 = y
        t1 = P(s + (y1, y2), s + (y3, y4)); a = P(s + (y1, y3), s + (y2, y4)); t2 = (a[1], a[0])
        t3 = P(s + (y1, y4), s + (y2, y3))
        cl = []
        for (tp, tn) in (t1, t2, t3):
            cl.append([-tp, t1[1], t2[1], t3[1]]); cl.append([-tn, t1[0], t2[0], t3[0]])
        groups.append((("GP", s, y), cl))
for f in F:
    grab(("PL", f), lambda f=f: E.forbid(E.nz(f)))
for i, f in enumerate(F):
    for j, g in enumerate(F):
        s = set(f) & set(g)
        if i == j or len(s) != 2: continue
        a, b = [u for u in f if u in s]; c = [u for u in f if u not in s][0]
        x, y = [u for u in g if u not in s]
        def ad(a=a, b=b, c=c, x=x, y=y):
            E.forbid(E.eq((a, b, c, x), 0)); E.forbid(E.eq((a, b, c, y), 0))
            for sg in (1, -1): E.forbid(E.eq((a, b, c, x), sg), E.notv((a, b, c, y), sg))
        grab(("AD", f, g), ad)

S = Solver(name="cd19" if "cd19" in Solver.__init__.__doc__ else "cd15", bootstrap_with=base_cls)
sel = []
for k, (name, cl) in enumerate(groups):
    v = E.new(); sel.append(v)
    for c in cl: S.add_clause(c + [-v])
assert not S.solve(assumptions=sel)
core = set(S.get_core())
for _ in range(5):
    S.solve(assumptions=[v for v in sel if v in core]); core = set(S.get_core())
cur = [v for v in sel if v in core]
print("initial core groups:", len(cur), flush=True)
i = 0
while i < len(cur):
    trial = cur[:i] + cur[i+1:]
    if not S.solve(assumptions=trial):
        c2 = set(S.get_core()); cur = [v for v in trial if v in c2]
    else:
        i += 1
idx = {v: k for k, v in enumerate(sel)}
names = [groups[idx[v]][0] for v in cur]
from collections import Counter
print("MUS size:", len(names), Counter(n[0] for n in names))
for n in names: print(n)

# ---- how hard is the MUS by hand?  look for a shallow case-split refutation by unit propagation
mus_cls = list(base_cls)
for v in cur: mus_cls += groups[idx[v]][1]
from pysat.solvers import Solver as S2
T = S2(name="m22", bootstrap_with=mus_cls)
chivars = sorted({abs(l) for c in mus_cls for l in c if abs(l) <= 3 * 210})
def refuted(assm):
    st, _ = T.propagate(assumptions=assm); return not st
print("refuted by UP alone:", refuted([]))
def depth(assm, d):
    if refuted(assm): return True
    if d == 0: return False
    for v in chivars:
        if v in assm or -v in assm: continue
        if depth(assm + [v], d - 1) and depth(assm + [-v], d - 1):
            return True
    return False
for d in range(1, 4):
    ok = depth([], d); print("refutable with case-split depth", d, ":", ok, flush=True)
    if ok: break

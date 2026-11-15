"""Independent check: is any (2,5)-set of PG(2,9) invariant under a collineation of order 3?
Own GF(9) arithmetic, own geometry; stdlib only (python-sat used only as a second, independent solver if present).
GF(9)=GF(3)[i]/(i^2+1); element a+3b <-> a+b*i.  One thread, small memory."""
import itertools, random, sys, time
sys.stdout.reconfigure(encoding="utf-8")
random.seed(20261111)
def P(*a): print(*a, flush=True)

# ---------------- GF(9) ----------------
def mk():
    add=[[0]*9 for _ in range(9)]; mul=[[0]*9 for _ in range(9)]
    for x in range(9):
        for y in range(9):
            a,b=x%3,x//3; c,d=y%3,y//3
            add[x][y]=(a+c)%3+3*((b+d)%3)
            mul[x][y]=(a*c-b*d)%3+3*((a*d+b*c)%3)
    return add,mul
ADD,MUL=mk()
NEG=[next(y for y in range(9) if ADD[x][y]==0) for x in range(9)]
INV=[0]+[next(y for y in range(9) if MUL[x][y]==1) for x in range(1,9)]
FROB=[MUL[MUL[x][x]][x] for x in range(9)]       # x -> x^3
assert all(FROB[x]==(x%3)+3*((3-x//3)%3) for x in range(9))  # conj: a+bi -> a-bi
def sub(x,y): return ADD[x][NEG[y]]
def eorder(x):
    k,y=1,x
    while y!=1: y=MUL[y][x]; k+=1
    return k
assert max(eorder(x) for x in range(1,9))==8
GEN=next(x for x in range(1,9) if eorder(x)==8)

# ---------------- plane ----------------
def norm(v):
    for c in v:
        if c: s=INV[c]; return tuple(MUL[s][t] for t in v)
    return None
PTS=sorted({norm(v) for v in itertools.product(range(9),repeat=3) if any(v)})
PIDX={p:i for i,p in enumerate(PTS)}
def dot(u,v): return ADD[ADD[MUL[u[0]][v[0]]][MUL[u[1]][v[1]]]][MUL[u[2]][v[2]]]
LINES=[frozenset(i for i,p in enumerate(PTS) if dot(l,p)==0) for l in PTS]   # line with coordinate vector l
LIDX={L:i for i,L in enumerate(LINES)}
assert len(PTS)==91 and len(LINES)==91 and len(set(LINES))==91
assert all(len(L)==10 for L in LINES)
assert all(sum(i in L for L in LINES)==10 for i in range(91))
assert all(len(LINES[a]&LINES[b])==1 for a in range(91) for b in range(a))
P("plane: 91 points, 91 lines, 10 points/line, 10 lines/point, two lines meet in 1 point: OK")

# ---------------- collineations (M, f): x -> M * frob^f(x) ----------------
def mmul(A,B):
    return tuple(tuple(ADD[ADD[MUL[A[i][0]][B[0][j]]][MUL[A[i][1]][B[1][j]]]][MUL[A[i][2]][B[2][j]]] for j in range(3)) for i in range(3))
def mfrob(A): return tuple(tuple(FROB[x] for x in r) for r in A)
I3=((1,0,0),(0,1,0),(0,0,1))
def cmul(a,b):
    (A,f),(B,g)=a,b
    return (mmul(A,mfrob(B) if f else B),(f+g)%2)
def cpow(a,k):
    r=(I3,0)
    for _ in range(k): r=cmul(r,a)
    return r
def mvec(M,v): return tuple(ADD[ADD[MUL[M[i][0]][v[0]]][MUL[M[i][1]][v[1]]]][MUL[M[i][2]][v[2]]] for i in range(3))
def pperm(c):
    M,f=c
    return tuple(PIDX[norm(mvec(M,tuple(FROB[x] for x in p) if f else p))] for p in PTS)
def lperm(pp): return tuple(LIDX[frozenset(pp[i] for i in L)] for L in LINES)
def pord(pp):
    idn=tuple(range(91)); q=pp; k=1
    while q!=idn: q=tuple(pp[x] for x in q); k+=1
    return k
def det(M):
    t=lambda a,b,c:MUL[MUL[a][b]][c]
    s=0
    for (i,j,k),sg in (((0,1,2),1),((1,2,0),1),((2,0,1),1),((0,2,1),-1),((2,1,0),-1),((1,0,2),-1)):
        x=t(M[0][i],M[1][j],M[2][k]); s=ADD[s][x if sg>0 else NEG[x]]
    return s
def rank(M):
    A=[list(r) for r in M]; rk=0
    for c in range(3):
        pr=next((r for r in range(rk,3) if A[r][c]),None)
        if pr is None: continue
        A[rk],A[pr]=A[pr],A[rk]; iv=INV[A[rk][c]]; A[rk]=[MUL[iv][x] for x in A[rk]]
        for r in range(3):
            if r!=rk and A[r][c]:
                m=A[r][c]; A[r]=[sub(A[r][j],MUL[m][A[rk][j]]) for j in range(3)]
        rk+=1
    return rk
def rand_gl():
    while True:
        M=tuple(tuple(random.randrange(9) for _ in range(3)) for _ in range(3))
        if det(M): return M

# ---------------- step 1: classes of elements of order 3 ----------------
P("\n== Step 1: elements of order 3 in PGammaL(3,9) ==")
P("Theory: |PGL(3,9)|=%d, |PGammaL(3,9)|=%d."%(728*720*648//8,728*720*648//4))
P("  (g,f) of order 3 => f=0 (f is a homomorphism to Z/2, field automorphism has order 2).")
P("  In PGL: g^3=lambda*I; cubing is bijective on GF(9)* (gcd(3,8)=1), so g/mu has cube I; x^3-1=(x-1)^3 in char 3,")
P("  so g/mu is unipotent: Jordan type (2,1) [elation] or (3) [regular unipotent], and nothing else.")
P("  h u h^-1 = lambda u' with u,u' unipotent forces lambda=1, so PGL-classes of order-3 elements = GL-classes of unipotents.")
P("  Frobenius preserves Jordan type (rank of g-I), so PGammaL does not merge the two classes.")
E=((1,0,1),(0,1,0),(0,0,1))     # elation: I+E13
R=((1,1,0),(0,1,1),(0,0,1))     # regular unipotent: I+E12+E23
def comm_dim_count(g):
    rows=[]
    for i in range(3):
        for j in range(3):
            row=[0]*9
            for k in range(3):
                row[3*k+j]=ADD[row[3*k+j]][g[i][k]]      # (gX)_ij = sum_k g_ik X_kj
                row[3*i+k]=sub(row[3*i+k],g[k][j])       # (Xg)_ij = sum_k X_ik g_kj
            rows.append(row)
    A=rows; piv=[]; rk=0
    for c in range(9):
        pr=next((r for r in range(rk,len(A)) if A[r][c]),None)
        if pr is None: continue
        A[rk],A[pr]=A[pr],A[rk]; iv=INV[A[rk][c]]; A[rk]=[MUL[iv][x] for x in A[rk]]
        for r in range(len(A)):
            if r!=rk and A[r][c]:
                m=A[r][c]; A[r]=[sub(A[r][j],MUL[m][A[rk][j]]) for j in range(9)]
        piv.append(c); rk+=1
    free=[c for c in range(9) if c not in piv]; basis=[]
    for fc in free:
        v=[0]*9; v[fc]=1
        for r,pc in enumerate(piv): v[pc]=NEG[A[r][fc]]
        basis.append(v)
    cnt=0
    for co in itertools.product(range(9),repeat=len(basis)):
        X=[0]*9
        for c,b in zip(co,basis):
            for t in range(9): X[t]=ADD[X[t]][MUL[c][b[t]]]
        if det((tuple(X[0:3]),tuple(X[3:6]),tuple(X[6:9]))): cnt+=1
    return len(basis),cnt
GLORD=728*720*648
tot=0
for name,g in (("elation I+E13",E),("regular unipotent I+E12+E23",R)):
    assert cpow((g,0),3)[0]==I3 and g!=I3
    d,c=comm_dim_count(g); cs=GLORD//c; tot+=cs
    P("  %s: rank(g-I)=%d, commutant dim %d, |C_GL|=%d, class size=%d"%(name,rank(tuple(tuple(sub(g[i][j],1 if i==j else 0) for j in range(3)) for i in range(3))),d,c,cs))
P("  total = %d ; Steinberg count of nontrivial unipotents q^6-1 = %d -> %s"%(tot,9**6-1,"COMPLETE" if tot==9**6-1 else "MISMATCH"))
seen={}; allf0=True; n3=0
for _ in range(3000):
    c=(rand_gl(),random.randrange(2))
    o=pord(pperm(c))
    if o%3: continue
    h=c
    while True:   # reduce to exact order 3
        oo=pord(pperm(h))
        if oo==3: break
        h=cpow(h,oo//3) if oo%3==0 else h
        if pord(pperm(h))==oo: break
    if pord(pperm(h))!=3: continue
    n3+=1; M,f=h
    if f: allf0=False
    lam=cpow(h,3)[0][0][0]
    mu=next(m for m in range(1,9) if MUL[MUL[m][m]][m]==lam)
    mi=INV[mu]; U=tuple(tuple(MUL[mi][x] for x in r) for r in M)
    rk=rank(tuple(tuple(sub(U[i][j],1 if i==j else 0) for j in range(3)) for i in range(3)))
    seen[rk]=seen.get(rk,0)+1
P("  random PGammaL sample: %d order-3 elements found, all with f=0: %s, rank(g/mu - I) histogram %s (expected only ranks 1 and 2)"%(n3,allf0,seen))

# ---------------- solver ----------------
def orbits_of(pp):
    seen=[False]*91; orb=[]
    for s in range(91):
        if seen[s]: continue
        o=[]; x=s
        while not seen[x]: seen[x]=True; o.append(x); x=pp[x]
        orb.append(o)
    return orb
def build(orbs):
    oid=[0]*91
    for k,o in enumerate(orbs):
        for x in o: oid[x]=k
    cons=set()
    for L in LINES:
        cnt={}
        for x in L: cnt[oid[x]]=cnt.get(oid[x],0)+1
        cons.add(tuple(sorted(cnt.items())))
    return [list(c) for c in cons]
def dpll(orbs,cons,maxsol=10**9,tlimit=600):
    n=len(orbs); sz=[len(o) for o in orbs]; sols=[]; t0=time.time(); st={'nodes':0,'timeout':False,'capped':False}
    def prop(val):
        ch=True
        while ch:
            ch=False
            for c in cons:
                lo=0; un=0
                for o,m in c:
                    if val[o]==1: lo+=m
                    elif val[o] is None: un+=m
                hi=lo+un
                if not ((lo<=2<=hi) or (lo<=5<=hi)): return False
                for o,m in c:
                    if val[o] is None:
                        c1=(lo+m<=2<=hi) or (lo+m<=5<=hi)
                        c0=(lo<=2<=hi-m) or (lo<=5<=hi-m)
                        if not c1 and not c0: return False
                        if c1!=c0:
                            val[o]=1 if c1 else 0; ch=True
                            if val[o]==1: lo+=m
                            un-=m; hi=lo+un
        return True
    def rec(val):
        st['nodes']+=1
        if st['nodes']%500==0 and time.time()-t0>tlimit: st['timeout']=True
        if st['timeout'] or st['capped']: return
        if not prop(val): return
        free=[o for o in range(n) if val[o] is None]
        if not free:
            sols.append([o for o in range(n) if val[o]==1])
            if len(sols)>=maxsol: st['capped']=True
            return
        best=None
        for c in cons:
            fr=[o for o,m in c if val[o] is None]
            if fr and (best is None or len(fr)<len(best)): best=fr
        o=max(best,key=lambda x:sz[x])
        for b in (1,0):
            v=val[:]; v[o]=b; rec(v)
    rec([None]*n)
    return sols,st
def sat_enum(orbs,cons,limit=10**9):
    """independent cross-check: CNF with one clause per forbidden assignment of a line's orbit variables"""
    try: from pysat.solvers import Cadical153 as S
    except Exception:
        try: from pysat.solvers import Glucose4 as S
        except Exception: return None
    s=S()
    for c in cons:
        vs=[o for o,m in c]
        for bits in itertools.product((0,1),repeat=len(vs)):
            cnt=sum(m for (o,m),b in zip(c,bits) if b)
            if cnt not in (2,5): s.add_clause([-(v+1) if b else (v+1) for v,b in zip(vs,bits)])
    n=0
    while n<limit and s.solve():
        mod=s.get_model(); n+=1
        s.add_clause([-(v+1) if mod[v]>0 else (v+1) for v in range(len(orbs))])
    s.delete(); return n
def run(name,c,maxsol=10**9,tlimit=600):
    pp=pperm(c); o=pord(pp); lp=lperm(pp)
    po=orbits_of(pp); lo=orbits_of(lp)
    P("\n-- %s: order %d, %d point orbits (sizes %s), %d line orbits"%(name,o,len(po),sorted({len(x) for x in po}),len(lo)))
    cons=build(po); t=time.time()
    sols,st=dpll(po,cons,maxsol,tlimit)
    exh=not(st['timeout'] or st['capped'])
    P("   DPLL: %d invariant (2,5)-sets found, nodes %d, %.1fs, exhaustive=%s"%(len(sols),st['nodes'],time.time()-t,exh))
    sizes={}
    for s in sols:
        K=set(x for k in s for x in po[k])
        assert all(len(L&K) in (2,5) for L in LINES)
        sizes[len(K)]=sizes.get(len(K),0)+1
    P("   sizes of found sets (each verified directly against all 91 lines): %s"%sizes)
    if exh:
        n=sat_enum(po,cons)
        P("   SAT cross-check (pysat, clause per forbidden line assignment): %s solutions %s"%(n,"(AGREES)" if n==len(sols) else "(DISAGREES or unavailable)"))
    return sols,exh

# ---------------- step 2: order 3 ----------------
P("\n== Step 2: exhaustive search, order-3 representatives ==")
res3={}
for name,g in (("elation I+E13",E),("regular unipotent I+E12+E23",R)):
    c=(g,0); assert pord(pperm(c))==3
    sols,exh=run(name,c); res3[name]=(len(sols),exh)
P("   Why exhaustive: DPLL branches on each free orbit variable both ways; it prunes only when some line's interval [assigned, assigned+free] misses both 2 and 5,")
P("   and forces a variable only when the other value is infeasible; at leaves all counts are exact. A (2,5)-set invariant under <g> is a union of <g>-orbits, i.e. one of the 2^k unions searched.")

# ---------------- step 3: controls ----------------
P("\n== Step 3: controls (same code) ==")
inv=(((1,0,0),(0,1,0),(0,0,NEG[1])),0)
run("control order 2: homology diag(1,1,-1)",inv,maxsol=20000,tlimit=300)
run("control order 2: Frobenius (Baer involution)",(I3,1),maxsol=20000,tlimit=300)
def find_order(p,avoid=None,tries=400):
    for _ in range(tries):
        c=(rand_gl(),0); o=pord(pperm(c))
        if o%p==0:
            h=cpow(c,o//p)
            if avoid is None or sorted(len(x) for x in orbits_of(pperm(h)))!=avoid: return h
    return None
c13=find_order(13); run("control order 13",c13,tlimit=300)
c5=find_order(5); run("control order 5",c5,tlimit=300)
c5b=find_order(5,avoid=sorted(len(x) for x in orbits_of(pperm(c5))))
if c5b: run("control order 5 (a different cycle type)",c5b,tlimit=300)
else: P("   (all order-5 elements sampled have the same cycle type as the first control)")

# ---------------- step 4: conic groups ----------------
P("\n== Step 4: conic groups ==")
def sym2(a,b,c,d):
    m=MUL; ad=ADD
    return ((m[a][a],ad[m[a][b]][m[a][b]],m[b][b]),
            (m[a][c],ad[m[a][d]][m[b][c]],m[b][d]),
            (m[c][c],ad[m[c][d]][m[c][d]],m[d][d]))
def gen_group(gens):
    ps=[pperm((g,0)) for g in gens]; idn=tuple(range(91)); G={idn}; fr=[idn]
    while fr:
        nf=[]
        for x in fr:
            for g in ps:
                y=tuple(g[i] for i in x)
                if y not in G: G.add(y); nf.append(y)
        fr=nf
    return G,ps
M1=sym2(1,1,0,1); Mi=sym2(1,3,0,1); W=sym2(0,1,NEG[1],0)
dd=lambda a:sym2(a,0,0,INV[a])
dg=sym2(GEN,0,0,1)
groups=[("A4=PSL(2,3) via Sym^2 (PΩ(3,3), order 12)",[M1,W]),
        ("S4=PGL(2,3) (SO(3,3), order 24)",[M1,W,sym2(1,0,0,NEG[1])]),
        ("A6=PSL(2,9)=PΩ(3,9) (order 360)",[M1,Mi,W,dd(GEN)]),
        ("PGL(2,9)=SO(3,9) (order 720)",[M1,Mi,W,dd(GEN),dg])]
conic=set(PIDX[norm((MUL[s][s],MUL[s][t],MUL[t][t]))] for s,t in [(1,0)]+[(x,1) for x in range(9)])
P("   conic xz=y^2: %d points (expect 10); line intersection numbers with it: %s"%(len(conic),sorted({len(L&conic) for L in LINES})))
for name,gens in groups:
    G,ps=gen_group(gens)
    seen=[False]*91; orbs=[]
    for s in range(91):
        if seen[s]: continue
        o=[s]; seen[s]=True; k=0
        while k<len(o):
            for g in ps:
                y=g[o[k]]
                if not seen[y]: seen[y]=True; o.append(y)
            k+=1
        orbs.append(o)
    has3=any(pord(g)%3==0 for g in G)
    P("\n-- %s: |G|=%d, point-orbit sizes %s, has elements of order 3: %s"%(name,len(G),sorted(len(o) for o in orbs),has3))
    k=len(orbs); best=None; hits=[]; spectra=[]
    for mask in range(1,2**k-1):
        K=set()
        for j in range(k):
            if mask>>j&1: K.update(orbs[j])
        cnt={}
        for L in LINES:
            c=len(L&K); cnt[c]=cnt.get(c,0)+1
        bad=sum(v for c,v in cnt.items() if c not in (2,5))
        if bad==0: hits.append(len(K))
        if len(K) in (26,35): spectra.append((len(K),bad,dict(sorted(cnt.items())),[len(orbs[j]) for j in range(k) if mask>>j&1]))
        if best is None or bad<best[0]: best=(bad,len(K),dict(sorted(cnt.items())),[len(orbs[j]) for j in range(k) if mask>>j&1])
    P("   unions of orbits checked: %d ; (2,5)-sets among them: %s"%(2**k-2,hits if hits else "none"))
    P("   nearest miss over all proper unions: %d lines with count outside {2,5}; size %d, spectrum %s, orbit sizes %s"%best)
    for s in spectra: P("   union of size %d: bad lines %d, spectrum %s, orbits %s"%s)
    if not spectra: P("   no union of orbits has size 26 or 35")
P("\n== SUMMARY ==")
P("order-3 classes searched (count of invariant (2,5)-sets, exhaustive?): %s"%res3)

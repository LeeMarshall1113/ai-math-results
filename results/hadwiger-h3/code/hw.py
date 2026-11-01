"""Exact re-implementation of Prymak's L(P,tau) (arXiv 2112.10698) with float LP margin + exact rational check."""
import itertools, sys, csv, random
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog

vers=[[0,0,0],[0,0,1],[0,1,0],[1,0,0],[0,1,1],[1,0,1],[1,1,0],[1,1,1]]
edges=[[[0,0,0],[1,0,0]],[[1,0,0],[0,1,0]],[[0,1,0],[1,0,0]],[[0,0,0],[0,1,0]],[[0,0,1],[1,0,0]],[[1,0,1],[0,1,0]],[[0,1,1],[1,0,0]],[[0,0,1],[0,1,0]],[[0,0,0],[0,0,1]],[[1,0,0],[0,0,1]],[[1,1,0],[0,0,1]],[[0,1,0],[0,0,1]]]
v_edges=[[0,3],[3,6],[2,6],[0,2],[1,5],[5,7],[4,7],[1,4],[0,1],[3,5],[6,7],[2,4]]

def hrep_points(pts):
    """exact H-rep (c0,c1,c2,c3): c0+c.x>=0 for conv(pts), full-dim assumed; deduped."""
    out=set()
    for a,b,c in itertools.combinations(pts,3):
        u=[b[i]-a[i] for i in range(3)]; v=[c[i]-a[i] for i in range(3)]
        n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
        if n==[0,0,0]: continue
        d=-sum(n[i]*a[i] for i in range(3))
        vals=[d+sum(n[i]*p[i] for i in range(3)) for p in pts]
        if all(x>=0 for x in vals): pass
        elif all(x<=0 for x in vals): n=[-x for x in n]; d=-d
        else: continue
        # normalise to integers
        from math import gcd
        den=1
        for x in n+[d]:
            den=den*F(x).denominator//gcd(den,F(x).denominator)
        ints=[int(F(x)*den) for x in [d]+n]
        g=0
        for x in ints: g=gcd(g,abs(x))
        out.add(tuple(x//g for x in ints))
    return out

def rects(a,d):
    ss=[[[0,a[0],a[1]],[0,a[0]+d[0],a[1]],[0,a[0],a[1]+d[1]],[0,a[0]+d[0],a[1]+d[1]]],
        [[1,a[0],a[1]],[1,a[0]+d[0],a[1]],[1,a[0],a[1]+d[1]],[1,a[0]+d[0],a[1]+d[1]]],
        [[a[2],0,a[3]],[a[2]+d[2],0,a[3]],[a[2],0,a[3]+d[3]],[a[2]+d[2],0,a[3]+d[3]]],
        [[a[2],1,a[3]],[a[2]+d[2],1,a[3]],[a[2],1,a[3]+d[3]],[a[2]+d[2],1,a[3]+d[3]]],
        [[a[4],a[5],0],[a[4]+d[4],a[5],0],[a[4],a[5]+d[5],0],[a[4]+d[4],a[5]+d[5],0]],
        [[a[4],a[5],1],[a[4]+d[4],a[5],1],[a[4],a[5]+d[5],1],[a[4]+d[4],a[5]+d[5],1]]]
    return ss

def QP_hrep(a,d):
    ss=rects(a,d)
    H=set()
    for i0 in range(4):
        for i2 in range(4):
            for i4 in range(4):
                pts=[ss[0][i0],ss[1][i0],ss[2][i2],ss[3][i2],ss[4][i4],ss[5][i4]]
                H|=hrep_points([[F(x) for x in p] for p in pts])
    H.add((0,1,0,0))
    return sorted(H)

def build(a,d,tau):
    """return list of constraints as (coef over 60 vars (dict var->F), const F) meaning const+sum>=0"""
    H=QP_hrep(a,d); ss=rects(a,d); cons=[]
    def add_tr(x,k):
        for h in H:
            co={}; c=F(h[0])
            for i in range(3):
                c+=h[i+1]*F(x[i]); co[k+i]=co.get(k+i,0)+F(h[i+1])
            cons.append((co,c))
    def add_seg(ei,k,kk):
        v,u=edges[ei]
        for h in H:
            co={}; c=F(h[0])
            for i in range(3):
                c+=h[i+1]*F(v[i]); co[k+i]=co.get(k+i,0)+F(h[i+1])
                co[kk]=co.get(kk,0)+F(h[i+1])*u[i]
            cons.append((co,c))
    for v in range(8): add_tr(vers[v],3*v)
    for f in range(6):
        for t in range(4): add_tr(ss[f][t],24+3*f)
    l=42
    for i in range(12):
        if i not in tau:
            add_seg(i,3*v_edges[i][0],l); add_seg(i,3*v_edges[i][1],l); l+=1
        else:
            add_seg(i,3*v_edges[i][0],l); add_seg(i,3*v_edges[i][1],l+1)
            f=tau.index(i)
            add_seg(i,24+3*f,l); add_seg(i,24+3*f,l+1); l+=2
    return cons,len(H)

def margin_solve(cons,nv=60):
    """max s s.t. (const+co.w)/||co|| >= s (only constraints with nonzero co); returns (s,w)"""
    A=[];b=[]
    for co,c in cons:
        row=np.zeros(nv+1); nrm=np.sqrt(sum(float(v)**2 for v in co.values()))
        if nrm==0:
            if c<0: return -1,None
            continue
        for k,v in co.items(): row[k]=-float(v)/nrm
        row[nv]=1.0
        A.append(row); b.append(float(c)/nrm)
    obj=np.zeros(nv+1); obj[nv]=-1
    bounds=[(None,None)]*nv+[(None,1.0)]
    r=linprog(obj,A_ub=np.array(A),b_ub=np.array(b),bounds=bounds,method='highs-ds',options={'presolve':True})
    if r.status!=0: return None,None
    return r.x[nv],r.x[:nv]

def exact_check(cons,w,den=10**9):
    wq=[F(round(x*den),den) for x in w]
    for co,c in cons:
        if c+sum(v*wq[k] for k,v in co.items())<0: return False,wq
    return True,wq

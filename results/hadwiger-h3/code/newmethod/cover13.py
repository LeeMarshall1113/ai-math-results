"""Search for 13-translate covers of E u V_p by O_p at a single p (float LP exploration).
Structure: 8 vertex translates a(v); 5 face translates b_f (f != f0) each covering q_f and the middle of edge tau(f) (or no edge);
q_f0 is covered by the vertex translate a(v0), v0 a vertex of face f0.  Edges not in tau's image are covered by their two
endpoint translates.  Returns best uniform margin (>0 means closed-O_p cover with slack, hence open cover by int O_p too)."""
import itertools, sys, os
_u=r'C:\Users\hackathon\AppData\Roaming\Python\Python314\site-packages'
if os.path.isdir(_u) and _u not in sys.path: sys.path.append(_u)
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

CUBE=[tuple(v) for v in itertools.product((0,1),repeat=3)]
EDGES=[(a,b) for a,b in itertools.combinations(range(8),2) if sum(abs(x-y) for x,y in zip(CUBE[a],CUBE[b]))==1]
FACES=[(i,j) for i in range(3) for j in (0,1)]
def face_vertices(f):
    i,j=f; return [k for k in range(8) if CUBE[k][i]==j]
def face_edges(f):
    fv=set(face_vertices(f)); return [e for e,(a,b) in enumerate(EDGES) if a in fv and b in fv]

def Ap(p):
    p1,p2,p3,p4,p5,p6=p
    return [(0,p1,p2),(1,p1,p2),(p3,0,p4),(p3,1,p4),(p5,p6,0),(p5,p6,1)]

def hrep(verts):
    h=ConvexHull(np.array(verts,float)); E=h.equations
    # inside: n.x + c <= 0 ; normalise rows
    E=E/np.linalg.norm(E[:,:3],axis=1)[:,None]
    return E

def solve(p,f0,v0,tau,H=None,Hm=None,facepts=None):
    """Hm: optional H-rep for the merged translate a(v0) (a different body, used up to translation)."""
    """tau: dict face-index -> edge-index or None, for faces != f0. Returns margin."""
    Q=Ap(p)
    if H is None: H=hrep(Q)
    nv=0; idx={}
    def var(name,k=1):
        nonlocal nv; idx[name]=nv; nv+=k
    for v in range(8): var(('a',v),3)
    for f in range(6):
        if f!=f0: var(('b',f),3)
    for e in range(12):
        var(('t',e));
        if e in tau.values(): var(('s',e))
    rows=[];rhs=[]
    # constraint: point x(w) in translate T + O  <=>  n.(x - T) + c <= -s  for each facet
    def contain(const_pt, lin_pt, T):
        # point = const_pt + sum lin_pt[var]*dir ; lin_pt: list of (varidx, dirvec)
        HH=Hm if (Hm is not None and T==idx[('a',v0)]) else H
        for n in HH:
            r=np.zeros(nv_total+1); b=-n[3]-n[:3]@np.array(const_pt,float)
            for k,dvec in lin_pt: r[k]+=n[:3]@np.array(dvec,float)
            for i in range(3): r[T+i]-=n[i]
            r[nv_total]=1.0
            rows.append(r); rhs.append(b)
    nv_total=nv
    for v in range(8): contain(CUBE[v],[],idx[('a',v)])
    for f in range(6):
        T=idx[('a',v0)] if f==f0 else idx[('b',f)]
        for pt in (facepts[f] if facepts is not None else [Q[f]]): contain(pt,[],T)
    tauinv={e:f for f,e in tau.items() if e is not None}
    for e,(a,b) in enumerate(EDGES):
        u=tuple(CUBE[b][i]-CUBE[a][i] for i in range(3))
        if e not in tauinv:
            contain(CUBE[a],[(idx[('t',e)],u)],idx[('a',a)])
            contain(CUBE[a],[(idx[('t',e)],u)],idx[('a',b)])
        else:
            B=idx[('b',tauinv[e])]
            contain(CUBE[a],[(idx[('t',e)],u)],idx[('a',a)])
            contain(CUBE[a],[(idx[('s',e)],u)],idx[('a',b)])
            contain(CUBE[a],[(idx[('t',e)],u)],B)
            contain(CUBE[a],[(idx[('s',e)],u)],B)
    c=np.zeros(nv_total+1); c[-1]=-1
    r=linprog(c,A_ub=np.array(rows),b_ub=np.array(rhs),bounds=[(None,None)]*nv_total+[(None,1)],method='highs')
    return r.x[-1] if r.status==0 else None

def best13(p,verbose=False):
    H=hrep(Ap(p)); best=(-9,None)
    for f0 in range(6):
        for v0 in face_vertices(FACES[f0]):
            others=[f for f in range(6) if f!=f0]
            choices=[face_edges(FACES[f])+[None] for f in others]
            for combo in itertools.product(*choices):
                es=[e for e in combo if e is not None]
                if len(set(es))<len(es): continue
                tau=dict(zip(others,combo))
                m=solve(p,f0,v0,tau,H)
                if m is not None and m>best[0]:
                    best=(m,(f0,v0,tau))
                    if verbose: print(" ",m,best[1],flush=True)
    return best

if __name__=="__main__":
    import random
    random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 0)
    for t in range(int(sys.argv[2]) if len(sys.argv)>2 else 3):
        p=[random.random() for _ in range(6)]
        b=best13(p)
        print(" ".join("%.3f"%x for x in p),"best13 margin %.5f"%b[0],b[1],flush=True)

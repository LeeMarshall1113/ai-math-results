"""Referee check (2026-10-31): Prymak's domain D (and U^D_20, and the 39,380 certified 1/10 start boxes) is NOT a
fundamental domain for the cube-symmetry action on configurations p in [0,1]^6.
Exact rational arithmetic. Runs in ~1 min, 1 thread, <100 MB.  python symmetry_gap.py"""
import itertools, random
from fractions import Fraction as Fr
def pts(p):
    p1,p2,p3,p4,p5,p6=p
    return [(0,p1,p2),(1,p1,p2),(p3,0,p4),(p3,1,p4),(p5,p6,0),(p5,p6,1)]
def topar(P):  # read p back from the point lying on each face x_i = 0
    d={}
    for q in P:
        for i in range(3):
            if q[i]==0: d[i]=q
    return (d[0][1],d[0][2],d[1][0],d[1][2],d[2][0],d[2][1])
G=[(perm,fl) for perm in itertools.permutations(range(3)) for fl in itertools.product((0,1),repeat=3)]
def app(g,q):
    perm,fl=g; r=[q[perm[i]] for i in range(3)]
    return tuple(1-r[i] if fl[i] else r[i] for i in range(3))
def in_UD(k):
    k1,k2,k3,k4,k5,k6=k; return k2+k4<=k1+k6+1<=k3+k5+2
need=set(c for c in itertools.product(range(5),range(5),range(5),range(10),range(10),range(10))
   if any(in_UD([2*c[i]+e[i] for i in range(6)]) for e in itertools.product((0,1),repeat=6)))
def in_start(p):   # p in the closed union of the certified 1/10 start boxes (same set certcheck requires)
    cand=[[c for c in range((4 if i<3 else 9)+1) if Fr(c,10)<=x<=Fr(c+1,10)] for i,x in enumerate(p)]
    return any(c in need for c in itertools.product(*cand))
def inD(x): return all(y<=Fr(1,2) for y in x[:3]) and x[1]+x[3]<=x[0]+x[5]<=x[2]+x[4]
def inDstar(x):  # a correct fundamental domain: |s_z-1| <= |s_y-1| <= |s_x-1|
    return all(y<=Fr(1,2) for y in x[:3]) and abs(x[1]+x[3]-1)<=abs(x[0]+x[5]-1)<=abs(x[2]+x[4]-1)
print('start boxes', len(need))
p=[Fr(459,500),Fr(1481,5000),Fr(1689,10000),Fr(9529,10000),Fr(9359,10000),Fr(3079,10000)]
ims=[topar([app(g,q) for q in pts(p)]) for g in G]
print('example p =',[str(x) for x in p])
print(' images in D:',sum(inD(x) for x in ims),' in certified region:',sum(in_start(x) for x in ims),' in D*:',sum(inDstar(x) for x in ims))
random.seed(7); N=5000; bD=bS=bDs=0
for _ in range(N):
    q=[Fr(random.randint(1,9999),10000) for _ in range(6)]
    ims=[topar([app(g,x) for x in pts(q)]) for g in G]
    bD+=not any(inD(x) for x in ims); bS+=not any(in_start(x) for x in ims); bDs+=not any(inDstar(x) for x in ims)
print('random p (N=%d): no image in D %d, no image in certified region %d, no image in D* %d'%(N,bD,bS,bDs))

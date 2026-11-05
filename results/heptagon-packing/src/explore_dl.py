import numpy as np, sys
sys.path.insert(0, 'src')
from hept import *
d = kk_double_lattice()
print('a,b', d['a'], d['b'])
print('heptagon area', AREA, 'half-length parallelogram area', d['par_area'])
delta = AREA/(2*d['par_area']); print('delta_DL', delta)
closed = (2/97)*(-111+492*U-356*U**2); print('closed form', closed)
acrit = AREA/(2*delta); print('acrit per triangle', acrit)
# centers: lattice L = <t1,t2>, upright copies at L, reflected copies at 2p2 + L
t1, t2 = d['t1'], d['t2']
pts = []; refl = []
for i in range(-4, 5):
    for j in range(-4, 5):
        pts.append(i*t1+j*t2); refl.append(False)
        pts.append(2*d['p2']+i*t1+j*t2); refl.append(True)
pts = np.array(pts); refl = np.array(refl)
# overlap check & contacts
H = [verts(p, 0, r) for p, r in zip(pts, refl)]
o = np.argmin(np.linalg.norm(pts, axis=1))
ds = []
for j in range(len(pts)):
    if j == o: continue
    dist = np.linalg.norm(pts[j]-pts[o])
    if dist < 2.2:
        ds.append((round(dist, 5), round(sep(H[o], H[j]), 8), refl[j]))
print('neighbours of origin (dist, sep, refl):'); [print('  ', x) for x in sorted(ds)]
from scipy.spatial import Delaunay
tri = Delaunay(pts)
c = np.linalg.norm(pts[tri.simplices].mean(1), axis=1) < 3
ars = sorted(set(round(abs(tri_area(*pts[s])), 9) for s in tri.simplices[c]))
crs = sorted(set(round(circumradius(*pts[s]), 6) for s in tri.simplices[c]))
print('Delaunay triangle areas near origin', ars)
print('circumradii', crs)
for s in tri.simplices[c][:4]:
    P = pts[s]; L = [np.linalg.norm(P[(k+1)%3]-P[k]) for k in range(3)]
    ang = [np.degrees(np.arccos(np.dot(P[(k+1)%3]-P[k], P[(k+2)%3]-P[k])/np.linalg.norm(P[(k+1)%3]-P[k])/np.linalg.norm(P[(k+2)%3]-P[k]))) for k in range(3)]
    print(' tri edges', np.round(L, 5), 'angles', np.round(ang, 3), 'refl', refl[s])

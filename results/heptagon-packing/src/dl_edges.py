"""Pair configurations (r, a, b) of the oriented (ccw) edges of the two DL triangle types."""
import numpy as np, sys
sys.path.insert(0, 'src')
from hept import *
P7 = 2*np.pi/7
d = kk_double_lattice()
H = {'0': (np.zeros(2), 0.0), 't1': (d['t1'], 0.0), '2p2': (2*d['p2'], np.pi/7), '2p6': (2*d['p6'], np.pi/7)}
def cfg(i, j):
    (ci, ti), (cj, tj) = H[i], H[j]; v = cj-ci; phi = np.arctan2(v[1], v[0])
    return np.array([np.linalg.norm(v), np.mod(ti-phi, P7), np.mod(tj-phi, P7)])
def R(e): return np.array([e[0], np.mod(e[2]+np.pi/7, P7), np.mod(e[1]+np.pi/7, P7)])
for name, tri in (('s2', ('0', 't1', '2p2')), ('s6', ('0', '2p6', 't1'))):
    print(name)
    for k in range(3):
        i, j = tri[k], tri[(k+1) % 3]; e = cfg(i, j)
        print(f'  {i:>3}->{j:<3} r={e[0]:.6f} a={e[1]:.6f} b={e[2]:.6f}  self-reverse={np.allclose(R(e), e)}  R(e)=({R(e)[1]:.6f},{R(e)[2]:.6f})')

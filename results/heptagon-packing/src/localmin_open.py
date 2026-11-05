"""Float local minima of b started from open-box centres (canonical frame)."""
import numpy as np, sys, time
OPEN = sys.argv[1]; NS = int(sys.argv[2])
sys.argv = [sys.argv[0], '3', '2', '0', '0', 'nonobtuse', '0.3', '0', '0.05']
sys.path.insert(0, 'src')
exec(open('src/lp_transfer3.py').read().split("if __name__")[0])
from scipy.optimize import minimize
c = np.load('runs/lp4b_c.npy')
S = np.load(OPEN); Cs = S.mean(2)
def extra(x):  # canonical: c0c1 longest, nonobtuse
    return np.array([x[0]**2-(x[1]**2+x[2]**2), x[0]**2-((x[0]-x[1])**2+x[2]**2)])
rng = np.random.default_rng(0); idx = rng.choice(len(Cs), NS, replace=False); res = []
for i in idx:
    r = minimize(lambda x: bval(x, c), Cs[i], method='SLSQP',
                 constraints=[{'type': 'ineq', 'fun': lambda x: np.concatenate([bcons(x), extra(x)])}],
                 options={'maxiter': 300, 'ftol': 1e-14})
    if np.concatenate([bcons(r.x), extra(r.x)]).min() > -1e-9: res.append((bval(r.x, c)-ACRIT, r.x))
res.sort(key=lambda t: t[0])
print('open', len(S), 'feasible minima', len(res), 'b-acrit:', [f'{a:.2e}' for a, _ in res[:12]])
for a, x in res[:3]: print(f'{a:.3e}', np.round(x, 5))

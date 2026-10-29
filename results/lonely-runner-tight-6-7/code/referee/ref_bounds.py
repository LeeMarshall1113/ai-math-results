"""Second implementation: bookkeeping in exact integer arithmetic (plus mpmath logs for display).

For a primitive v with ML(v) <= 1/(k+1), not a dilate of any known tight T:
  (1) Lemma T (strict MSS, m = C(k+1,2)+1):  S = sum v_i < sum_J v_J <= m^(k-1).
  (2) prod v <= (S/k)^k  (AM-GM).
  (3) for each T, F_T := first nonzero I^T_j(v), j = 2..k, where
        I^T_j(v) = e_j(w) e_1(T^2)^j - e_1(w)^j e_j(T^2),  w = v^2.
      Both terms are >= 0, so |I_j| <= max(e_j(w) E1^j, e1(w)^j Ej) <= sum of them.
      e_1(w) <= S^2 and e_j(w) <= C(k,j) (e_1(w)/k)^j (Maclaurin).
  (4) every gate p with p | prod v divides rad(prod v) <= prod v; every other gate divides
      some F_T.  Hence  prod_{gates} p <= prod v * prod_T |F_T|.
Contradiction iff prod_{gates} p > Bound.  We check this with integers (Fractions).
"""
import json, math, sys, itertools, os
from fractions import Fraction as Fr
import mpmath as mp

mp.mp.dps = 50
KNOWN = {6: [(1, 2, 3, 4, 5, 6)],
         7: [(1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]}


def esym(xs, j):
    return sum(math.prod(c) for c in itertools.combinations(xs, j))


def bound(k, S, sharp=False):
    """Upper bound (Fraction) on prod v * prod_T max_j |I^T_j| given sum v <= S."""
    prod_bd = Fr(S, k) ** k
    tot = prod_bd
    per = []
    for T in KNOWN[k]:
        w = [a * a for a in T]
        E = {j: esym(w, j) for j in range(1, k + 1)}
        best = 0
        for j in range(2, k + 1):
            t1 = math.comb(k, j) * Fr(S * S, k) ** j * E[1] ** j
            t2 = Fr(S * S) ** j * E[j]
            ub = max(t1, t2) if sharp else t1 + t2
            best = max(best, ub)
        per.append(best)
        tot *= best
    return tot, prod_bd, per


def main(k, gates_file):
    m = math.comb(k + 1, 2) + 1
    S = m ** (k - 1)
    rows = [json.loads(l) for l in open(gates_file)]
    gates = sorted(r['p'] for r in rows if r['gate'])
    G = math.prod(gates)
    out = {}
    for label, SS, sharp in (('standard (S <= m^(k-1), |a-b| <= a+b)', S, False),
                             ('sharpened (S <= m^(k-1) - (2^k-k-1), |a-b| <= max)',
                              S - (2 ** k - k - 1), True)):
        B, P, per = bound(k, SS, sharp)
        lg = mp.log(G)
        lb = mp.log(B.numerator) - mp.log(B.denominator)
        out[label] = dict(log_gate_mass=mp.nstr(lg, 12), log_needed=mp.nstr(lb, 12),
                          margin_nats=mp.nstr(lg - lb, 12), exact_holds=(G * B.denominator > B.numerator),
                          log_prod=mp.nstr(mp.log(P.numerator) - mp.log(P.denominator), 10),
                          log_per_T=[mp.nstr(mp.log(x.numerator) - mp.log(x.denominator), 10) for x in per])
        # smallest prefix of gates that already suffices
        acc = 1
        for p in gates:
            acc *= p
            if acc * B.denominator > B.numerator:
                out[label]['prefix_ends_at'] = p
                break
        else:
            out[label]['prefix_ends_at'] = None
    print(json.dumps(dict(k=k, m=m, S0=S, n_gates=len(gates),
                          nongates=sorted(r['p'] for r in rows if not r['gate']),
                          pmax=max(r['p'] for r in rows), results=out), indent=1))


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2])

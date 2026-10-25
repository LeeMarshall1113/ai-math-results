"""Exact check that a gate set is large enough (TRUSTED, ~30 lines of integer arithmetic).

Claim checked:   lcm(2..k+1) * prod(gates)  >  B_k := ( C(k+1,2)^(k-1) / k )^k ,
done exactly as   (lcm * prod)  *  k^k   >   C(k+1,2)^(k(k-1))   in integers.
Also checks: every gate is a prime > k+1 and the gates are distinct (so they are coprime
to lcm(2..k+1) and to each other, and lcm * prod divides the product of every
counterexample).  Floats below are for display only.

usage: python mass.py k p1 p2 ...      (exit 0 iff the claim holds)
"""
import math
import sys


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def check(k, gates):
    assert len(set(gates)) == len(gates), 'repeated gate'
    assert all(is_prime(p) and p > k + 1 for p in gates), 'gate not a prime > k+1'
    lcm = 1
    for j in range(2, k + 2):
        lcm = lcm * j // math.gcd(lcm, j)
    prod = math.prod(gates)
    lhs = lcm * prod * k ** k
    rhs = math.comb(k + 1, 2) ** (k * (k - 1))
    ok = lhs > rhs
    lnB = k * ((k - 1) * math.log(math.comb(k + 1, 2)) - math.log(k))
    info = dict(k=k, n_gates=len(gates), lcm=lcm, ln_lcm=math.log(lcm),
                ln_gates=sum(math.log(p) for p in gates), ln_B=lnB,
                needed_from_gates=lnB - math.log(lcm), holds=ok)
    return ok, info


if __name__ == '__main__':
    k = int(sys.argv[1])
    gates = [int(a) for a in sys.argv[2:]]
    ok, info = check(k, gates)
    print('k=%(k)d: %(n_gates)d gates, ln(prod gates) = %(ln_gates).4f, need > ln B_k - ln lcm = '
          '%(ln_B).4f - %(ln_lcm).4f = %(needed_from_gates).4f : %(holds)s (exact integer test)' % info)
    sys.exit(0 if ok else 1)

"""Monte-Carlo cost model for the full level-1 run on Modal.

  python modal_cost.py [--speedup S] [--lrat-factor F] [--check-factor C] [--cutoff H] [--mem-gib 4]

Per-cube solve times (no proof, local CPU, seconds) measured so far:  A: 352, 365, 704, 1957;  C: 601, 2272, and one right-censored at 3600.
A lognormal is fitted to each class by maximum likelihood (censoring handled), the container time of a cube is
   solve * lrat_factor + cake_check_factor * solve        (cake_lpr time / solve time measured in the speed check),
divided by --speedup (config gain); with --cutoff H a cube that is not finished after H hours is restarted with a new seed (the abandoned
time is paid) .  Price: 1 core 0.0473 $/h + mem_gib * 0.008 $/h (modal billing rates, 2026-11-08).
"""
import math
import random
import sys

CPU_RATE, MEM_RATE = 0.0473, 0.008


def arg(name, default):
    return float(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def fit(obs, cens):
    """lognormal MLE with right-censored observations (crude grid search)"""
    import math
    best = None
    xs = [math.log(x) for x in obs]
    cs = [math.log(x) for x in cens]
    for mu in [m / 20.0 for m in range(100, 260)]:
        for sg in [s / 20.0 for s in range(6, 80)]:
            ll = 0.0
            for x in xs:
                ll += -math.log(sg) - 0.5 * ((x - mu) / sg) ** 2
            for c in cs:
                z = (c - mu) / (sg * math.sqrt(2))
                sf = 0.5 * math.erfc(z)
                ll += math.log(max(sf, 1e-300))
            if best is None or ll > best[0]:
                best = (ll, mu, sg)
    return best[1], best[2]


def main():
    speedup = arg("--speedup", 1.0)
    lrat_f = arg("--lrat-factor", 1.4)
    check_f = arg("--check-factor", 0.4)
    cutoff = arg("--cutoff", 0.0)
    mem = arg("--mem-gib", 4.0)
    rate = CPU_RATE + mem * MEM_RATE
    A = fit([352, 365, 704, 1957], [])
    C = fit([601, 2272], [3600])
    print("lognormal fits (seconds): A mu=%.2f sigma=%.2f (median %.0f s, mean %.0f s) ; C mu=%.2f sigma=%.2f (median %.0f s, mean %.0f s)" % (
        A[0], A[1], math.exp(A[0]), math.exp(A[0] + A[1] ** 2 / 2), C[0], C[1], math.exp(C[0]), math.exp(C[0] + C[1] ** 2 / 2)))
    rnd = random.Random(1)
    tot = []
    for _ in range(4000):
        h = 0.0
        for n, (mu, sg) in ((52, A), (119, C)):
            for _ in range(n):
                paid = 0.0
                while True:
                    t = math.exp(rnd.gauss(mu, sg)) / speedup
                    cont = (t * (lrat_f + check_f)) / 3600.0
                    if cutoff and t * lrat_f / 3600.0 > cutoff:
                        paid += cutoff
                        continue
                    paid += cont
                    break
                h += paid
        tot.append(h)
    tot.sort()
    q = lambda f: tot[int(f * (len(tot) - 1))]
    print("container hours: median %.0f, mean %.0f, p90 %.0f, p99 %.0f" % (q(0.5), sum(tot) / len(tot), q(0.9), q(0.99)))
    print("Modal cost at %.4f $/container-hour: median $%.1f, mean $%.1f, p90 $%.1f, p99 $%.1f" % (rate, q(0.5) * rate, sum(tot) / len(tot) * rate, q(0.9) * rate, q(0.99) * rate))


if __name__ == "__main__":
    main()

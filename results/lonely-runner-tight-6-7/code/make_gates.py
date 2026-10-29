"""Builds the published gate lists gates/gates_k6.json and gates/gates_k7.json from the
outputs of both implementations, after checking that they agree prime by prime.

Usage (after both check.sh runs):  python make_gates.py
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ok = True
    for k in (6, 7):
        a = {r['p']: r for r in json.load(open(os.path.join(HERE, 'original', 'out', f'k{k}.json')))}
        b = {r['p']: r for r in (json.loads(l) for l in open(os.path.join(HERE, 'referee', 'out', f'k{k}.jsonl')) if l.strip())}
        if sorted(a) != sorted(b):
            print(f'k={k}: prime sets differ'); ok = False; continue
        rows = []
        for p in sorted(a):
            ra, rb = a[p], b[p]
            na = ra['nontight'] if k == 6 else ra['nonexempt']
            agree = ra['gate'] == rb['gate'] and ra['level1_orbits'] == rb['level1'] and na == rb['nonexempt']
            if not agree:
                print(f'k={k} p={p}: DISAGREE {ra} {rb}'); ok = False
            rows.append(dict(p=p, gate=bool(ra['gate']), level1_orbits=ra['level1_orbits'],
                             nonexempt_orbits=na,
                             impl1_surviving_lifts=ra['lifted_survivors'],
                             impl2_surviving_orbits=rb['survivors']))
        out = dict(k=k, pmin=min(a), pmax=max(a), n_primes=len(rows),
                   n_gates=sum(r['gate'] for r in rows),
                   nongates=[r['p'] for r in rows if not r['gate']],
                   note=('gate = no non-exempt level-one orbit has a closed-improper, '
                         'non-Lemma-G lift mod (k+1)p; impl1 counts surviving lifts up to its '
                         'symmetry reduction, impl2 counts orbits with at least one surviving lift'),
                   primes=rows)
        os.makedirs(os.path.join(HERE, 'gates'), exist_ok=True)
        with open(os.path.join(HERE, 'gates', f'gates_k{k}.json'), 'w', newline='\n') as f:
            json.dump(out, f, indent=1)
        print(f'k={k}: {len(rows)} primes {min(a)}..{max(a)}, {out["n_gates"]} gates, '
              f'non-gates {out["nongates"]}')
    print('GATE LISTS', 'WRITTEN (implementations agree)' if ok else 'NOT CONSISTENT')
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)

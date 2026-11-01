"""Aggregate replay results and run the coverage check.
usage: python summarize.py OUTDIR  -> prints and writes OUTDIR/SUMMARY.json"""
import json, glob, os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage, orbitboxes

od = sys.argv[1]
tot = collections.Counter(); slots = collections.Counter(); fails = []
sbs = collections.Counter(); pol = collections.Counter(); files = 0; sec = 0.0
skipped = 0
for f in sorted(glob.glob(os.path.join(od, '*.json'))):
    if f.endswith('SUMMARY.json'):
        continue
    d = json.load(open(f)); files += 1; sec += d.get('sec', 0)
    tag = os.path.basename(f)[:-5]
    for k, r in enumerate(d['records']):
        sbs[tuple(r['sb'])] += 1
        if r.get('skipped'):
            skipped += 1
            continue
        for key in ('n', 'plain', 'chord', 'fast_ok', 'exact_fallback_ok'):
            tot[key] += r[key]
        tot['startboxes_checked'] += 1
        tot['tiling_ok'] += (r['tiling'] == 'ok')
        pol[r.get('pol')] += 1
        slots.update(r['slots'])
        for x in r['fail']:
            fails.append(dict(x, file=tag, record=k, sb=r['sb']))
S = set(sbs)
dups = [k for k, c in sbs.items() if c > 1]
cov = coverage.main(S, ntest=20000)
cov['grid_orbit_test'] = orbitboxes.main(S)
summary = dict(files=files, cpu_sec=round(sec), totals=dict(tot), skipped_startboxes=skipped,
               chord_slot_stats=dict(slots), policies=dict(pol), duplicate_startboxes=len(dups),
               dup_examples=dups[:5], n_fail=len(fails), fails=fails[:200], coverage=cov)
json.dump(summary, open(os.path.join(od, 'SUMMARY.json'), 'w'), indent=1, default=str)
print(json.dumps({k: v for k, v in summary.items() if k != 'fails'}, indent=1, default=str))
print('first fails:', json.dumps(fails[:10], default=str))

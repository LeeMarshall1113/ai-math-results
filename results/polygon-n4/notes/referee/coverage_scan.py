import json, glob, collections
types=collections.Counter(); ats=collections.Counter(); vs=collections.Counter(); vars_=collections.Counter()
def walk(n):
    if isinstance(n,dict):
        t=n.get('t'); types[t]+=1
        if t=='mono': ats[repr(n['at'])]+=1; vars_[n['var']]+=1
        if t=='reduce':
            for l,v in n['S']: vs[repr(v)]+=1; vars_[l]+=1
        if t in('split','splitat'): vars_[n['var']]+=1
        if 'TG' in n:
            for k,s in n['TG'].items():
                for l,v in json.loads(k): vs['TG'+repr(v)]+=1
        for k,v in n.items():
            if k!='t': walk(v)
    elif isinstance(n,list):
        for x in n: walk(x)
easy=set(); trees=set(); tubes=set()
for l in open('result_0_836.jsonl'):
    r=json.loads(l)
    if r['ok'] is True: easy.add(r['face'])
for fn in sorted(glob.glob('hard[0-9]_*.jsonl')+glob.glob('kd_*.jsonl')+glob.glob('iso*_*.jsonl')):
    for l in open(fn):
        r=json.loads(l)
        if r.get('ok') and r.get('tree'): trees.add(r['face']); walk(r['tree'])
for fn in ['tube_t158.jsonl','tube_tubeA.jsonl','tube_tubeB.jsonl','tube_tubeC.jsonl','tube_tubeD.jsonl']:
    for l in open(fn):
        r=json.loads(l)
        if r.get('ok'): tubes.add(r['face']); walk(r['tree']); walk(r['TG']); walk({'TG':r['TG']})
can=set(x.strip() for x in open('faces_canon.txt') if x.strip())
print(len(easy),len(trees),len(tubes),len(easy|trees|tubes), len(can), (easy|trees|tubes)==can, len(easy&trees),len(trees&tubes),len(easy&tubes))
print(can-(easy|trees|tubes))
print(types); print(ats); print(vs); print(vars_)

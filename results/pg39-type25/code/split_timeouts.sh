#!/bin/bash
# Level-2 split of every Modal TIMEOUT cube that has no certs/split_<ID>.json yet (4 python+PySAT processes in parallel), then the child jobs are added to modal_jobs.json.
cd "$(dirname "$0")"
IDS=$(/c/Python314/python - <<'PY'
import json, os
out = []
for l in open("../results/modal_ledger.jsonl"):
    r = json.loads(l)
    if r["status"] == "TIMEOUT" and r["id"][0] in "AC" and r["id"][1].isdigit() and len(r["id"]) == 6 and not os.path.exists("../certs/split_%s.json" % r["id"]):
        out.append(r["id"])
print(" ".join(out))
PY
)
echo "to split: $IDS"
for id in $IDS; do
  while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 5; done
  /c/Python314/python -u split_level2.py $id 8 300 > ../logs/split_$id.log 2>&1 &
done
wait
/c/Python314/python modal_children_jobs.py $(echo $IDS | tr ' ' ',') 7200

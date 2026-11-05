#!/bin/bash
cd /c/Users/hackathon/Documents/GitHub/erdos-triage/heptagon-packing
export PYTHONPATH=/c/Users/hackathon/AppData/Roaming/Python/Python314/site-packages
( /c/Python314/python src/local_replay.py obt 30; /c/Python314/python src/local_replay.py obt 13 ) &
( /c/Python314/python src/local_replay.py main 5; /c/Python314/python src/local_replay.py main 39c ) &
wait

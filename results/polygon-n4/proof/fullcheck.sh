set -e
/c/Python314/python -u vsome.py 500
/c/Python314/python -u verify_tube.py tube_t158.jsonl tube_tubeA.jsonl tube_tubeB.jsonl tube_tubeC.jsonl tube_tubeD.jsonl

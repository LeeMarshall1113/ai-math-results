"""Monolithic cover CNF of the level-2 children (speed check): cube CNF + one selector per child + (OR selectors).
  python bench_cover.py BASE.cnf CHILDREN.json OUT.cnf
"""
import json
import sys

from bench_make import cover_cnf


def main():
    base, ch, out = sys.argv[1:4]
    kids = json.load(open(ch))["children"]
    cover_cnf(base, kids, out)
    print("cover cnf with %d selectors written" % len(kids))


if __name__ == "__main__":
    main()

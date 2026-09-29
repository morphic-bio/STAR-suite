#!/usr/bin/env python3
"""Compare SAM bodies without discarding duplicate records or changing fields."""
import argparse
from collections import Counter
import json
from pathlib import Path


def compare(a, b):
    aa, bb = a.splitlines(), b.splitlines()
    if aa == bb:
        status = "identical"
    elif Counter(aa) == Counter(bb):
        status = "order_only"
    else:
        status = "different"
    return {"status": status, "records_a": len(aa), "records_b": len(bb)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("a", type=Path)
    parser.add_argument("b", type=Path)
    parser.add_argument("--require-order", action="store_true")
    args = parser.parse_args()
    result = compare(args.a.read_bytes(), args.b.read_bytes())
    print(json.dumps(result, sort_keys=True))
    return int(result["status"] == "different" or
               (args.require_order and result["status"] != "identical"))


if __name__ == "__main__":
    raise SystemExit(main())

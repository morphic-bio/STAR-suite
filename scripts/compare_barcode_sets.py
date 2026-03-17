#!/usr/bin/env python3
"""Compare two barcode sets and report overlap statistics."""

import argparse
import gzip
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def load_barcodes(path: Path, strip_suffix: bool, translation_map=None):
    barcodes = set()
    with open_text(path) as fh:
        for line in fh:
            bc = line.strip()
            if not bc:
                continue
            if strip_suffix and "-" in bc:
                base, suffix = bc.rsplit("-", 1)
                if suffix.isdigit():
                    bc = base
            if translation_map:
                bc = translation_map.get(bc, bc)
            barcodes.add(bc)
    return barcodes


def load_translation(path: Path, direction: str):
    mapping = {}
    with open_text(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            if direction == "left-to-right":
                mapping[parts[0]] = parts[1]
            else:
                mapping[parts[1]] = parts[0]
    return mapping


def main():
    parser = argparse.ArgumentParser(description="Compare two barcode sets.")
    parser.add_argument("file_a", help="First barcode file")
    parser.add_argument("file_b", help="Second barcode file")
    parser.add_argument("--label-a", default="A")
    parser.add_argument("--label-b", default="B")
    parser.add_argument("--strip-suffix", action="store_true")
    parser.add_argument("--barcode-translation", default="")
    parser.add_argument("--translation-direction", default="left-to-right")
    args = parser.parse_args()

    trans_map = None
    if args.barcode_translation:
        trans_map = load_translation(Path(args.barcode_translation), args.translation_direction)

    a = load_barcodes(Path(args.file_a), args.strip_suffix, trans_map)
    b = load_barcodes(Path(args.file_b), args.strip_suffix, trans_map)

    shared = a & b
    jaccard = len(shared) / len(a | b) if (a | b) else 0.0

    print(f"{args.label_a}: {len(a)}")
    print(f"{args.label_b}: {len(b)}")
    print(f"shared: {len(shared)}")
    print(f"{args.label_a}_only: {len(a - b)}")
    print(f"{args.label_b}_only: {len(b - a)}")
    print(f"jaccard: {jaccard:.6f}")


if __name__ == "__main__":
    main()

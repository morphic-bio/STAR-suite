#!/usr/bin/env python3
import argparse
import gzip
import os
import sys


def open_maybe_gzip(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def normalize_qname(raw):
    name = raw.strip().split()[0]
    if name.startswith("@"):
        name = name[1:]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def load_readnames(path):
    names = set()
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            names.add(normalize_qname(line))
    return names


def filter_fastq(in_path, out_path, names):
    kept = 0
    total = 0
    with open_maybe_gzip(in_path, "rt") as fin, open_maybe_gzip(out_path, "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                raise SystemExit(f"Truncated FASTQ at {in_path}, record {total + 1}")
            total += 1
            if normalize_qname(header) in names:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)
                kept += 1
    return kept, total


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTQ files by read name list (one per line)."
    )
    parser.add_argument("--readnames", required=True, help="Path to read name list")
    parser.add_argument("--r1", required=True, help="R1 FASTQ (gz or plain)")
    parser.add_argument("--r2", required=True, help="R2 FASTQ (gz or plain)")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    args = parser.parse_args()

    names = load_readnames(args.readnames)
    if not names:
        raise SystemExit("Read name list is empty")

    os.makedirs(args.out_dir, exist_ok=True)
    r1_out = os.path.join(args.out_dir, os.path.basename(args.r1))
    r2_out = os.path.join(args.out_dir, os.path.basename(args.r2))

    kept_r1, total_r1 = filter_fastq(args.r1, r1_out, names)
    kept_r2, total_r2 = filter_fastq(args.r2, r2_out, names)

    print(f"R1: kept {kept_r1}/{total_r1} reads -> {r1_out}")
    print(f"R2: kept {kept_r2}/{total_r2} reads -> {r2_out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def iter_fastq(path: Path):
    with gzip.open(path, 'rt') as handle:
        while True:
            name = handle.readline()
            if not name:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise SystemExit(f'Truncated FASTQ record in {path}')
            yield name, seq, plus, qual


def main() -> int:
    parser = argparse.ArgumentParser(description='Copy the first N paired FASTQ records from gzipped R1/R2 inputs')
    parser.add_argument('--r1', required=True)
    parser.add_argument('--r2', required=True)
    parser.add_argument('--out-r1', required=True)
    parser.add_argument('--out-r2', required=True)
    parser.add_argument('--read-limit', required=True, type=int)
    args = parser.parse_args()

    if args.read_limit <= 0:
        raise SystemExit('--read-limit must be > 0')

    out_r1 = Path(args.out_r1)
    out_r2 = Path(args.out_r2)
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    out_r2.parent.mkdir(parents=True, exist_ok=True)

    copied = 0
    with gzip.open(out_r1, 'wt') as h1, gzip.open(out_r2, 'wt') as h2:
        for rec1, rec2 in zip(iter_fastq(Path(args.r1)), iter_fastq(Path(args.r2))):
            h1.writelines(rec1)
            h2.writelines(rec2)
            copied += 1
            if copied >= args.read_limit:
                break

    print(copied)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

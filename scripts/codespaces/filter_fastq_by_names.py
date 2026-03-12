#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def normalize_name(raw: str) -> str:
    name = raw.strip().split()[0]
    if name.startswith('@'):
        name = name[1:]
    if name.endswith('/1') or name.endswith('/2'):
        name = name[:-2]
    return name


def load_names(path: Path) -> set[str]:
    names = set()
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            names.add(normalize_name(line))
    return names


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
    parser = argparse.ArgumentParser(description='Filter gzipped FASTQ records by normalized read name')
    parser.add_argument('--names', required=True)
    parser.add_argument('--r1', required=True)
    parser.add_argument('--r2', required=True)
    parser.add_argument('--out-r1', required=True)
    parser.add_argument('--out-r2', required=True)
    args = parser.parse_args()

    wanted = load_names(Path(args.names))
    out_r1 = Path(args.out_r1)
    out_r2 = Path(args.out_r2)
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    out_r2.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    with gzip.open(out_r1, 'wt') as h1, gzip.open(out_r2, 'wt') as h2:
        for rec1, rec2 in zip(iter_fastq(Path(args.r1)), iter_fastq(Path(args.r2))):
            if normalize_name(rec1[0]) not in wanted:
                continue
            h1.writelines(rec1)
            h2.writelines(rec2)
            kept += 1

    print(kept)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

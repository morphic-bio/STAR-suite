#!/usr/bin/env python3
import argparse
import gzip
import shutil
import sys
from pathlib import Path


def fail(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)


def open_maybe_gzip(path: Path, mode: str):
    if "b" not in mode and "t" not in mode:
        mode += "t"
    with path.open("rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, mode)
    return open(path, mode)


def derive_split_outputs(inpath: Path, split_dir: Path) -> tuple[Path, Path]:
    base = inpath.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if base.endswith(suffix):
            stem = base[: -len(suffix)]
            return split_dir / f"{stem}_Y{suffix}", split_dir / f"{stem}_noY{suffix}"
    return split_dir / f"{base}_Y", split_dir / f"{base}_noY"


def aggregate(split_paths: list[Path], outpath: Path) -> None:
    with gzip.open(outpath, "wt", compresslevel=6) as out_handle:
        for path in split_paths:
            if not path.exists():
                fail(f"Missing split FASTQ: {path}")
            with open_maybe_gzip(path, "rt") as in_handle:
                shutil.copyfileobj(in_handle, out_handle)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Aggregate remove_y_reads outputs in input FASTQ order"
    )
    parser.add_argument("--split-dir", required=True, help="Directory containing remove_y_reads outputs")
    parser.add_argument("--output-y", required=True, help="Aggregated Y FASTQ output path (.gz)")
    parser.add_argument("--output-noy", required=True, help="Aggregated noY FASTQ output path (.gz)")
    parser.add_argument("fastqs", nargs="+", help="Original FASTQ inputs in desired aggregation order")
    args = parser.parse_args()

    split_dir = Path(args.split_dir)
    if not split_dir.is_dir():
        fail(f"Split directory not found: {split_dir}")

    y_paths: list[Path] = []
    noy_paths: list[Path] = []
    for fastq in args.fastqs:
        fastq_path = Path(fastq)
        if not fastq_path.exists():
            fail(f"FASTQ input not found: {fastq_path}")
        y_path, noy_path = derive_split_outputs(fastq_path, split_dir)
        y_paths.append(y_path)
        noy_paths.append(noy_path)

    aggregate(y_paths, Path(args.output_y))
    aggregate(noy_paths, Path(args.output_noy))
    return 0


if __name__ == "__main__":
    sys.exit(main())

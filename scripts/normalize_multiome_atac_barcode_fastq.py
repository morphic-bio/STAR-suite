#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from pathlib import Path


TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract a 16 bp 10x Multiome ATAC barcode from a longer i5/barcode "
            "FASTQ read and write a normalized barcode FASTQ for Chromap."
        )
    )
    parser.add_argument("--input", required=True, help="Input barcode FASTQ(.gz)")
    parser.add_argument("--output", required=True, help="Output normalized FASTQ(.gz)")
    parser.add_argument(
        "--start",
        type=int,
        default=9,
        help="1-based start of the barcode window in the input read (default: 9)",
    )
    parser.add_argument(
        "--length",
        type=int,
        default=16,
        help="Barcode length to extract (default: 16)",
    )
    parser.add_argument(
        "--reverse-complement",
        action="store_true",
        help="Reverse-complement the extracted barcode and reverse the quality string.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite output if present")
    return parser.parse_args()


def open_read(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def open_write(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt")
    return open(path, "w", encoding="utf-8")


def rc(seq: str) -> str:
    return seq.translate(TRANS)[::-1].upper()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    if args.start < 1:
        raise SystemExit("--start must be 1-based and positive")
    if args.length < 1:
        raise SystemExit("--length must be positive")
    if output_path.exists() and not args.force:
        raise SystemExit(f"Output exists; pass --force to overwrite: {output_path}")

    offset = args.start - 1
    end = offset + args.length
    n_reads = 0
    min_input_len: int | None = None

    with open_read(input_path) as src, open_write(output_path) as dst:
        while True:
            name = src.readline()
            if not name:
                break
            seq = src.readline()
            plus = src.readline()
            qual = src.readline()
            if not qual:
                raise SystemExit(f"Truncated FASTQ record after {n_reads} reads in {input_path}")
            seq = seq.rstrip("\n")
            qual = qual.rstrip("\n")
            if len(seq) < end or len(qual) < end:
                raise SystemExit(
                    f"Read {n_reads + 1} is too short for requested window "
                    f"{args.start}:{end}: sequence={len(seq)} quality={len(qual)}"
                )
            barcode = seq[offset:end]
            barcode_qual = qual[offset:end]
            if args.reverse_complement:
                barcode = rc(barcode)
                barcode_qual = barcode_qual[::-1]
            min_input_len = len(seq) if min_input_len is None else min(min_input_len, len(seq))
            dst.write(name)
            dst.write(barcode + "\n")
            dst.write(plus)
            dst.write(barcode_qual + "\n")
            n_reads += 1

    print(f"Wrote {output_path}")
    print(f"reads={n_reads}")
    print(f"barcode_start_1based={args.start}")
    print(f"barcode_length={args.length}")
    print(f"reverse_complement={str(args.reverse_complement).lower()}")
    print(f"min_input_read_length={min_input_len if min_input_len is not None else 0}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path


def fail(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def open_maybe_gzip(path: Path):
    with path.open("rb") as handle:
        magic = handle.read(2)
    return gzip.open if magic == b"\x1f\x8b" else open


def normalize_header(raw: str, ignore_comments: bool) -> str:
    header = raw.rstrip("\n")
    if not header.startswith("@"):
        fail(f"Invalid FASTQ header: {header!r}")
    if ignore_comments:
        return header.split()[0]
    return header


def read_record(handle):
    header = handle.readline()
    if not header:
        return None
    seq = handle.readline()
    plus = handle.readline()
    qual = handle.readline()
    if not qual:
        fail("Truncated FASTQ record")
    return header, seq, plus, qual


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare FASTQ record streams with optional header-comment normalization"
    )
    parser.add_argument("--actual", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--ignore-header-comments", action="store_true")
    parser.add_argument("--label", default="FASTQ")
    args = parser.parse_args()

    actual = Path(args.actual)
    reference = Path(args.reference)
    if not actual.exists():
        fail(f"Missing actual FASTQ: {actual}")
    if not reference.exists():
        fail(f"Missing reference FASTQ: {reference}")

    actual_open = open_maybe_gzip(actual)
    reference_open = open_maybe_gzip(reference)
    with actual_open(actual, "rt") as act, reference_open(reference, "rt") as ref:
        idx = 0
        while True:
            act_rec = read_record(act)
            ref_rec = read_record(ref)
            if act_rec is None and ref_rec is None:
                break
            if act_rec is None or ref_rec is None:
                fail(f"{args.label}: FASTQ length mismatch at record {idx + 1}")
            idx += 1
            act_header, act_seq, _act_plus, act_qual = act_rec
            ref_header, ref_seq, _ref_plus, ref_qual = ref_rec
            if normalize_header(act_header, args.ignore_header_comments) != normalize_header(
                ref_header, args.ignore_header_comments
            ):
                fail(f"{args.label}: header mismatch at record {idx}")
            if act_seq != ref_seq:
                fail(f"{args.label}: sequence mismatch at record {idx}")
            if act_qual != ref_qual:
                fail(f"{args.label}: quality mismatch at record {idx}")

    print(f"PASS: {args.label} matched {idx} records")
    return 0


if __name__ == "__main__":
    sys.exit(main())

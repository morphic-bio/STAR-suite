#!/usr/bin/env python3
import argparse
import gzip
import io
import sys
from pathlib import Path


def fail(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)


def normalize_qname(raw: str) -> str:
    name = raw.split()[0]
    if name.startswith("@"):
        name = name[1:]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    elif len(name) > 2 and name[-2:] in {".1", ".2"}:
        name = name[:-2]
    return name


def open_maybe_gzip(path: Path):
    with path.open("rb") as handle:
        magic = handle.read(2)
    return gzip.open if magic == b"\x1f\x8b" else open


class DeterministicFastqWriter:
    def __init__(self, path: Path):
        self.path = path
        self._raw = path.open("wb")
        self._gz = gzip.GzipFile(filename="", mode="wb", fileobj=self._raw, compresslevel=6, mtime=0)
        self._txt = io.TextIOWrapper(self._gz, encoding="ascii", newline="")

    def writelines(self, lines):
        self._txt.writelines(lines)

    def close(self):
        self._txt.flush()
        self._txt.close()
        self._gz.close()
        self._raw.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


def read_fastq_record(handle, path: Path):
    record = [handle.readline() for _ in range(4)]
    if not record[0]:
        return None
    if any(line == "" for line in record[1:]):
        fail(f"Truncated FASTQ record in {path}")
    if not record[0].startswith("@"):
        fail(f"Unexpected FASTQ header in {path}: {record[0]!r}")
    return record


def write_manifest(path: Path, input_r1: Path, input_r2: Path, sample_a: str, sample_b: str,
                   pairs_per_sample: int, total_pairs: int):
    with path.open("w", encoding="ascii", newline="") as handle:
        handle.write("sample\tinput_r1\tinput_r2\tstart_pair_0based\tend_pair_0based\tpairs_written\n")
        handle.write(f"{sample_a}\t{input_r1}\t{input_r2}\t0\t{pairs_per_sample - 1}\t{pairs_per_sample}\n")
        handle.write(
            f"{sample_b}\t{input_r1}\t{input_r2}\t{pairs_per_sample}\t{2 * pairs_per_sample - 1}\t{pairs_per_sample}\n"
        )
        handle.write(f"summary\t{input_r1}\t{input_r2}\t0\t{total_pairs - 1}\t{total_pairs}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Split one PE FASTQ pair into two consecutive pseudo-samples")
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ(.gz)")
    parser.add_argument("--r2", required=True, help="Input R2 FASTQ(.gz)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--sample-a", required=True, help="Sample name for first consecutive chunk")
    parser.add_argument("--sample-b", required=True, help="Sample name for second consecutive chunk")
    parser.add_argument("--pairs-per-sample", required=True, type=int, help="Pairs to write into each output sample")
    args = parser.parse_args()

    if args.pairs_per_sample <= 0:
        fail("--pairs-per-sample must be > 0")

    r1 = Path(args.r1)
    r2 = Path(args.r2)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not r1.is_file():
        fail(f"Missing R1 FASTQ: {r1}")
    if not r2.is_file():
        fail(f"Missing R2 FASTQ: {r2}")

    out_a_r1 = outdir / f"{args.sample_a}_R1_001.fastq.gz"
    out_a_r2 = outdir / f"{args.sample_a}_R2_001.fastq.gz"
    out_b_r1 = outdir / f"{args.sample_b}_R1_001.fastq.gz"
    out_b_r2 = outdir / f"{args.sample_b}_R2_001.fastq.gz"

    opener_r1 = open_maybe_gzip(r1)
    opener_r2 = open_maybe_gzip(r2)
    total_needed = args.pairs_per_sample * 2
    total_seen = 0

    with opener_r1(r1, "rt") as in_r1, opener_r2(r2, "rt") as in_r2, \
        DeterministicFastqWriter(out_a_r1) as a_r1, DeterministicFastqWriter(out_a_r2) as a_r2, \
        DeterministicFastqWriter(out_b_r1) as b_r1, DeterministicFastqWriter(out_b_r2) as b_r2:
        while total_seen < total_needed:
            rec_r1 = read_fastq_record(in_r1, r1)
            rec_r2 = read_fastq_record(in_r2, r2)
            if rec_r1 is None or rec_r2 is None:
                fail(f"Not enough input pairs in {r1} / {r2}; needed {total_needed}, found {total_seen}")
            q1 = normalize_qname(rec_r1[0].rstrip())
            q2 = normalize_qname(rec_r2[0].rstrip())
            if q1 != q2:
                fail(f"Mate qname mismatch at pair {total_seen}: {q1!r} != {q2!r}")
            if total_seen < args.pairs_per_sample:
                a_r1.writelines(rec_r1)
                a_r2.writelines(rec_r2)
            else:
                b_r1.writelines(rec_r1)
                b_r2.writelines(rec_r2)
            total_seen += 1

    write_manifest(outdir / "split_manifest.tsv", r1, r2, args.sample_a, args.sample_b, args.pairs_per_sample, total_seen)
    print(f"Wrote {args.pairs_per_sample} pairs to {out_a_r1.name} / {out_a_r2.name}")
    print(f"Wrote {args.pairs_per_sample} pairs to {out_b_r1.name} / {out_b_r2.name}")
    print(f"Total pairs consumed: {total_seen}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

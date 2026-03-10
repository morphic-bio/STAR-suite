#!/usr/bin/env python3
import argparse
import gzip
import subprocess
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


def normalize_qname(raw: str) -> str:
    name = raw.split()[0]
    if name.startswith("@"):
        name = name[1:]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    elif len(name) > 2 and name[-2:] in {".1", ".2"}:
        name = name[:-2]
    return name


def load_y_names(bam_path: Path) -> set[str]:
    proc = subprocess.run(
        ["samtools", "view", str(bam_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    names = set()
    for line in proc.stdout.splitlines():
        if not line:
            continue
        names.add(normalize_qname(line.split("\t", 1)[0]))
    return names


def derive_outputs(inpath: Path, outdir: Path) -> tuple[Path, Path]:
    base = inpath.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if base.endswith(suffix):
            stem = base[: -len(suffix)]
            return outdir / f"{stem}_Y{suffix}", outdir / f"{stem}_noY{suffix}"
    return outdir / f"{base}_Y", outdir / f"{base}_noY"


def split_fastq(path: Path, outdir: Path, y_names: set[str]) -> None:
    y_path, noy_path = derive_outputs(path, outdir)
    with open_maybe_gzip(path, "rt") as inp, y_path.open("wt") as y_out, noy_path.open("wt") as noy_out:
        total = 0
        y_count = 0
        noy_count = 0
        while True:
            header = inp.readline()
            if not header:
                break
            seq = inp.readline()
            plus = inp.readline()
            qual = inp.readline()
            if not qual:
                fail(f"Truncated FASTQ record in {path}")
            target = y_out if normalize_qname(header.rstrip()) in y_names else noy_out
            target.write(header)
            target.write(seq)
            target.write(plus)
            target.write(qual)
            total += 1
            if target is y_out:
                y_count += 1
            else:
                noy_count += 1
    print(f"{path}: total={total} Y={y_count} noY={noy_count}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Reference FASTQ splitter using normalized BAM qnames")
    parser.add_argument("--ybam", required=True, help="Y-only BAM file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("fastqs", nargs="+", help="FASTQ files to split")
    args = parser.parse_args()

    ybam = Path(args.ybam)
    outdir = Path(args.outdir)
    if not ybam.is_file():
        fail(f"Missing Y BAM: {ybam}")
    outdir.mkdir(parents=True, exist_ok=True)

    y_names = load_y_names(ybam)
    print(f"Loaded {len(y_names)} normalized Y read names")
    if not y_names:
        fail("No Y read names found in BAM")

    for fastq in args.fastqs:
        split_fastq(Path(fastq), outdir, y_names)
    return 0


if __name__ == "__main__":
    sys.exit(main())

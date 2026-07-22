#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path


def fastq_pairs(source: Path) -> list[tuple[Path, Path]]:
    pairs = []
    for r1 in sorted(source.glob("*_R1_001.fastq.gz")):
        r2 = Path(str(r1).replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
        if not r2.is_file():
            raise SystemExit(f"missing mate for {r1}")
        pairs.append((r1, r2))
    if not pairs:
        raise SystemExit(f"no paired FASTQs under {source}")
    return pairs


def copy_records(source: Path, destination: Path, count: int) -> int:
    destination.parent.mkdir(parents=True, exist_ok=True)
    copied = 0
    with gzip.open(source, "rt") as reader, gzip.open(destination, "wt") as writer:
        while copied < count:
            record = [reader.readline() for _ in range(4)]
            if not record[0]:
                break
            if any(line == "" for line in record):
                raise SystemExit(f"truncated FASTQ record in {source}")
            writer.writelines(record)
            copied += 1
    return copied


def stage_library(source: Path, destination: Path, total_reads: int) -> dict:
    pairs = fastq_pairs(source)
    quotient, remainder = divmod(total_reads, len(pairs))
    rows = []
    copied_total = 0
    for index, (r1, r2) in enumerate(pairs):
        target = quotient + (1 if index < remainder else 0)
        copied_r1 = copy_records(r1, destination / r1.name, target)
        copied_r2 = copy_records(r2, destination / r2.name, target)
        if copied_r1 != copied_r2:
            raise SystemExit(f"mate count mismatch for {r1.name}: {copied_r1} != {copied_r2}")
        rows.append({"r1": r1.name, "r2": r2.name, "reads": copied_r1})
        copied_total += copied_r1
    return {"source": str(source), "reads": copied_total, "lanes": rows}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a deterministic UCSF GEX+guide fixture for scheduler tutorials."
    )
    parser.add_argument("--dataset-root", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--reads-per-library", type=int, default=1_000_000)
    args = parser.parse_args()
    if args.reads_per_library <= 0:
        raise SystemExit("--reads-per-library must be positive")

    source_sample = Path(args.dataset_root) / args.sample
    output_sample = Path(args.output_root) / args.sample
    manifest = {
        "schema": "star_suite.ucsf_slurm_tutorial_fixture/v1",
        "sample": args.sample,
        "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "method": "deterministic_prefix_per_lane",
        "reads_per_library": args.reads_per_library,
        "libraries": {
            "GEX": stage_library(
                source_sample / "GEX",
                output_sample / "GEX",
                args.reads_per_library,
            ),
            "guides": stage_library(
                source_sample / "guides",
                output_sample / "guides",
                args.reads_per_library,
            ),
        },
    }
    output_sample.mkdir(parents=True, exist_ok=True)
    (output_sample / "FIXTURE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

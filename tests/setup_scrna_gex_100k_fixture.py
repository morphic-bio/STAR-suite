#!/usr/bin/env python3
"""Create a deterministic, paired 100K GEX fixture without reading whole lanes."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(stream):
    lines = [stream.readline() for _ in range(4)]
    if not lines[0].startswith(b"@") or not lines[2].startswith(b"+"):
        raise ValueError("missing or malformed FASTQ record")
    if len(lines[1].rstrip()) != len(lines[3].rstrip()):
        raise ValueError("FASTQ sequence/quality length mismatch")
    return lines


def name(lines):
    return lines[0].split()[0].removesuffix(b"/1").removesuffix(b"/2")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", type=Path, action="append", required=True)
    parser.add_argument("--r2", type=Path, action="append", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if len(args.r1) != len(args.r2):
        parser.error("provide the same number of R1 and R2 lanes in matching order")
    args.outdir.mkdir(parents=True, exist_ok=False)
    manifest = {"schema": "scrna-gex-100k-v1", "read_pairs": 100000, "lanes": []}
    for index, (r1, r2) in enumerate(zip(args.r1, args.r2)):
        pairs = 100000 // len(args.r1) + (index < 100000 % len(args.r1))
        paths = [args.outdir / f"L{index + 1:03d}_R{mate}.fastq.gz" for mate in (1, 2)]
        with gzip.open(r1, "rb") as a, gzip.open(r2, "rb") as b:
            with paths[0].open("wb") as out1, paths[1].open("wb") as out2:
                with gzip.GzipFile(filename="", fileobj=out1, mode="wb", mtime=0) as x:
                    with gzip.GzipFile(filename="", fileobj=out2, mode="wb", mtime=0) as y:
                        for _ in range(pairs):
                            left, right = record(a), record(b)
                            if name(left) != name(right):
                                raise ValueError(f"mismatched read pair in {r1} / {r2}")
                            x.write(b"".join(left))
                            y.write(b"".join(right))
        manifest["lanes"].append({
            "read_pairs": pairs, "source_r1": str(r1), "source_r2": str(r2),
            "r1": paths[0].name, "r2": paths[1].name,
            "r1_sha256": sha256(paths[0]), "r2_sha256": sha256(paths[1]),
        })
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared exactly 100000 matched read pairs: {args.outdir}")


if __name__ == "__main__":
    main()

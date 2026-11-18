#!/usr/bin/env python3
"""Create the 100K read-pair fixture of one STAR-perturb add-on regression.

    python3 tests/setup_perturb_addons_100k_fixture.py --assay splitseq \
        --dataset GSE352893 --run SRR41469160 \
        --cdna SRR41469160_1.fastq.gz --barcode SRR41469160_2.fastq.gz \
        --outdir /storage/downsampled_100K/splitseq_gse352893

Writes the first 100,000 read pairs of a public library as cDNA.fastq.gz and
barcode.fastq.gz with a manifest.json that holds their checksums. The gzip
header carries no name or time, so the same source gives the same bytes.
tests/run_perturb_addons_100k_regression.py reads the fixture. Only the head
of the source files is read.
"""
import argparse
import gzip
import json
from pathlib import Path

from setup_scrna_gex_100k_fixture import name, record, sha256

READ_PAIRS = 100000


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assay", required=True, choices=("splitseq", "pipseq"))
    parser.add_argument("--dataset", required=True, help="public accession of the dataset, for example GSE352893")
    parser.add_argument("--run", required=True, help="public accession of the sequencing run")
    parser.add_argument("--cdna", type=Path, required=True, help="gzipped FASTQ of the cDNA read")
    parser.add_argument("--barcode", type=Path, required=True, help="gzipped FASTQ of the barcode read")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    paths = [args.outdir / "cDNA.fastq.gz", args.outdir / "barcode.fastq.gz"]
    with gzip.open(args.cdna, "rb") as a, gzip.open(args.barcode, "rb") as b:
        with paths[0].open("wb") as out1, paths[1].open("wb") as out2:
            with gzip.GzipFile(filename="", fileobj=out1, mode="wb", mtime=0) as x:
                with gzip.GzipFile(filename="", fileobj=out2, mode="wb", mtime=0) as y:
                    for _ in range(READ_PAIRS):
                        left, right = record(a), record(b)
                        if name(left) != name(right):
                            raise ValueError(f"mismatched read pair in {args.cdna} / {args.barcode}")
                        x.write(b"".join(left))
                        y.write(b"".join(right))
    manifest = {
        "schema": "perturb-addons-100k-v1", "assay": args.assay, "dataset": args.dataset, "run": args.run,
        "read_pairs": READ_PAIRS, "selection": "the first 100,000 read pairs of the source files",
        "source_cdna": str(args.cdna), "source_barcode": str(args.barcode),
        "cdna": paths[0].name, "barcode": paths[1].name,
        "cdna_sha256": sha256(paths[0]), "barcode_sha256": sha256(paths[1]),
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {READ_PAIRS} matched read pairs: {args.outdir}")


if __name__ == "__main__":
    main()

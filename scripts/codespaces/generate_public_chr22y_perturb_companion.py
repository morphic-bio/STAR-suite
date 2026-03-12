#!/usr/bin/env python3
import argparse
import csv
import gzip
import itertools
import os
from pathlib import Path

GUIDES = [
    ("guide_a", "GuideA", "ACGTACGTACGTACGTACGT"),
    ("guide_b", "GuideB", "TTGCAATGCAATGCAATGCA"),
]


def iter_fastq(path: Path):
    with gzip.open(path, "rt") as handle:
        while True:
            name = handle.readline().rstrip("\n")
            if not name:
                break
            seq = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            qual = handle.readline().rstrip("\n")
            yield name, seq, plus, qual


def enumerate_pairs(fastq_dir: Path):
    r1s = sorted(fastq_dir.glob("*_R1_001.fastq.gz"))
    if not r1s:
        raise SystemExit(f"No R1 FASTQs found under {fastq_dir}")
    pairs = []
    for r1 in r1s:
        r2 = Path(str(r1).replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
        if not r2.is_file():
            raise SystemExit(f"Missing R2 mate for {r1}")
        pairs.append((r1, r2))
    return pairs


def write_feature_ref(path: Path):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["id", "name", "sequence", "feature_type"])
        for guide_id, name, seq in GUIDES:
            writer.writerow([guide_id, name, seq, "CRISPR Guide Capture"])


def write_config(path: Path, ref_root: Path, feature_ref: Path, gex_dir: Path, guide_dir: Path):
    path.write_text(
        "[gene-expression]\n"
        f"reference,{ref_root}\n"
        "create-bam,false\n"
        "chemistry,auto\n\n"
        "[feature]\n"
        f"reference,{feature_ref}\n\n"
        "[libraries]\n"
        "fastq_id,fastqs,feature_types\n"
        f"public_chr22y,{gex_dir},Gene Expression\n"
        f"public_chr22y,{guide_dir},CRISPR Guide Capture\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a perturb companion library from the public chr22+chrY GEX fixture")
    parser.add_argument("--gex-fastq-dir", required=True)
    parser.add_argument("--ref-root", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--read-limit", type=int, default=4000)
    args = parser.parse_args()

    gex_fastq_dir = Path(args.gex_fastq_dir).resolve()
    ref_root = Path(args.ref_root).resolve()
    outdir = Path(args.outdir).resolve()
    gex_out = outdir / "gex" / "public_chr22y_perturb"
    guide_out = outdir / "guides" / "public_chr22y_perturb"
    outdir.mkdir(parents=True, exist_ok=True)
    gex_out.mkdir(parents=True, exist_ok=True)
    guide_out.mkdir(parents=True, exist_ok=True)

    pairs = enumerate_pairs(gex_fastq_dir)
    renamed_pairs = []
    for idx, (r1, r2) in enumerate(pairs, start=1):
        out_r1 = gex_out / f"public_chr22y_perturb_S1_L{idx:03d}_R1_001.fastq.gz"
        out_r2 = gex_out / f"public_chr22y_perturb_S1_L{idx:03d}_R2_001.fastq.gz"
        for dst, src in ((out_r1, r1), (out_r2, r2)):
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(src, dst)
        renamed_pairs.append((out_r1, out_r2))

    whitelist = set()
    feature_ref = outdir / "feature_ref.csv"
    whitelist_path = outdir / "whitelist_TRU.txt"
    config_path = outdir / "config.csv"
    manifest_path = outdir / "MANIFEST.txt"

    read_count = 0
    source_r1 = renamed_pairs[0][0]
    guide_r1 = guide_out / "public_chr22y_guides_S1_L001_R1_001.fastq.gz"
    guide_r2 = guide_out / "public_chr22y_guides_S1_L001_R2_001.fastq.gz"

    with gzip.open(guide_r1, "wt") as out_r1, gzip.open(guide_r2, "wt") as out_r2:
        for (name, seq, plus, qual), guide in zip(iter_fastq(source_r1), itertools.cycle(GUIDES)):
            cb = seq[:16]
            whitelist.add(cb)
            out_r1.write(f"{name}\n{seq}\n{plus}\n{qual}\n")
            guide_seq = guide[2]
            out_r2.write(f"{name.replace('/1', '/2')}\n{guide_seq}\n+\n{'I' * len(guide_seq)}\n")
            read_count += 1
            if read_count >= args.read_limit:
                break

    write_feature_ref(feature_ref)
    whitelist_path.write_text("\n".join(sorted(whitelist)) + "\n", encoding="utf-8")
    write_config(config_path, ref_root, feature_ref, gex_out, guide_out)

    manifest_path.write_text(
        "Public chr22+chrY perturb companion\n"
        f"Generated reads: {read_count}\n"
        f"GEX dir: {gex_out}\n"
        f"Guide dir: {guide_out}\n"
        f"Whitelist: {whitelist_path}\n"
        f"Feature ref: {feature_ref}\n"
        f"CR config: {config_path}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

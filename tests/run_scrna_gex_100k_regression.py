#!/usr/bin/env python3
"""Fixture-backed, non-Flex/non-perturb GEX regression against reviewed goldens."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

from compare_flex_hash_screen_mex import read_mex
from setup_scrna_gex_100k_fixture import sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GOLDEN = ROOT / "tests/fixtures/scrna_gex_100k_golden.json"
QC_KEYS = ("Number of input reads", "Uniquely mapped reads number",
           "Number of reads mapped to multiple loci", "Number of reads mapped to too many loci")


def digest_rows(rows):
    digest = hashlib.sha256()
    for row in sorted(rows):
        digest.update((json.dumps(row, separators=(",", ":")) + "\n").encode())
    return digest.hexdigest()


def matrix_summary(path):
    features, barcodes, counts = read_mex(path)
    if any(value < 0 for value in counts.values()):
        raise AssertionError(f"negative UMI count: {path}")
    return {
        "genes": len(features), "barcodes": len(barcodes), "nnz": len(counts),
        "umis": sum(counts.values()),
        "features_sha256": digest_rows(features), "barcodes_sha256": digest_rows(barcodes),
        "counts_sha256": digest_rows((gene, cb, count) for (gene, cb), count in counts.items()),
    }, counts


def run(command, log):
    log.with_suffix(".command.json").write_text(json.dumps(command, indent=2) + "\n")
    with log.open("w") as stream:
        subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)


def verify_report(actual, expected, profiles):
    if actual["inputs"] != expected["inputs"]:
        raise AssertionError("fixture/reference differs from the pinned golden")
    for profile in profiles:
        if actual["profiles"][profile] != expected["profiles"][profile]:
            raise AssertionError(f"{profile}: raw/filtered matrix, cell set, or mapping QC drift")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, default=os.environ.get("STAR_BIN", ROOT / "core/legacy/source/STAR"))
    parser.add_argument("--fixture", type=Path, default=os.environ.get("SCRNA_GEX_FIXTURE", "/storage/downsampled_100K/pbmc10k_gex"))
    parser.add_argument("--genome", type=Path, default=os.environ.get("GENOME_DIR", "/storage/autoindex_110_44/bulk_index"))
    parser.add_argument("--whitelist", type=Path, default=os.environ.get("CB_WHITELIST", "/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt"))
    parser.add_argument("--outdir", type=Path, default=os.environ.get("SCRNA_GEX_OUTDIR"))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--golden", type=Path, default=DEFAULT_GOLDEN)
    parser.add_argument("--profiles", nargs="+", choices=("vanilla", "modern", "modern_bam"),
                        default=["vanilla", "modern", "modern_bam"])
    parser.add_argument("--record-golden", type=Path,
                        help="explicit bootstrap only; never reports a regression pass")
    parser.add_argument("--verify-report", type=Path,
                        help="check a completed report without rerunning STAR")
    args = parser.parse_args()
    if args.verify_report:
        verify_report(json.loads(args.verify_report.read_text()),
                      json.loads(args.golden.read_text()), args.profiles)
        print(f"PASS: completed 100K report matches golden: {args.verify_report}")
        return
    if args.threads < 1:
        parser.error("--threads must be positive")
    if not args.star.is_file() or not os.access(args.star, os.X_OK):
        parser.error(f"STAR executable not found: {args.star}")
    if args.record_golden is not None and args.record_golden.exists():
        parser.error("refusing to overwrite an existing golden")
    if args.outdir is None:
        import tempfile
        args.outdir = Path(tempfile.mkdtemp(prefix="star-scrna-gex-100k-"))
    else:
        args.outdir.mkdir(parents=True, exist_ok=False)
    args.star = args.star.resolve()
    fixture = json.loads((args.fixture / "manifest.json").read_text())
    if fixture["schema"] != "scrna-gex-100k-v1" or fixture["read_pairs"] != 100000:
        raise ValueError("requires a prepared 100K paired GEX fixture")
    if sum(lane["read_pairs"] for lane in fixture["lanes"]) != 100000:
        raise ValueError("lane totals do not equal 100K")
    for lane in fixture["lanes"]:
        for mate in ("r1", "r2"):
            if sha256(args.fixture / lane[mate]) != lane[f"{mate}_sha256"]:
                raise ValueError(f"fixture checksum mismatch: {lane[mate]}")
    identity = {
        "lanes": [{k: v for k, v in lane.items() if not k.startswith("source_")}
                  for lane in fixture["lanes"]],
        "whitelist_sha256": sha256(args.whitelist),
        # Annotation/index metadata plus exact output digests guard reference
        # drift without rehashing a 30-GB index on every smoke invocation.
        "reference_metadata": {name: sha256(args.genome / name) for name in
                               ("genomeParameters.txt", "chrNameLength.txt", "geneInfo.tab", "exonInfo.tab")},
    }
    expected = None
    if args.record_golden is None:
        expected = json.loads(args.golden.read_text())
        if expected["inputs"] != identity:
            raise ValueError("fixture/reference differs from the pinned golden")
    report = {
        "schema": "scrna-gex-100k-regression-v1", "inputs": identity,
        "binary_sha256": sha256(args.star),
        "version": subprocess.check_output([str(args.star), "--version"], text=True).strip(),
        "source_revision": subprocess.check_output([str(args.star), "--source-revision"], text=True).strip(),
        "profiles": {},
    }
    common = [str(args.star), "--genomeDir", str(args.genome),
              "--runThreadN", str(args.threads), "--readFilesIn",
              ",".join(str(args.fixture / lane["r2"]) for lane in fixture["lanes"]),
              ",".join(str(args.fixture / lane["r1"]) for lane in fixture["lanes"]),
              "--readFilesCommand", "zcat", "--soloType", "CB_UMI_Simple",
              "--soloCBstart", "1", "--soloCBlen", "16", "--soloUMIstart", "17",
              "--soloUMIlen", "12", "--soloBarcodeReadLength", "0",
              "--soloCBwhitelist", str(args.whitelist), "--soloStrand", "Forward",
              "--soloFeatures", "Gene", "GeneFull", "--clipAdapterType", "CellRanger4",
              "--clip3pPolyG", "yes"]
    profiles = {
        "vanilla": ["--outSAMtype", "None"],
        "modern": ["--defaultCoreScrna", "Yes", "--outSAMtype", "None"],
        "modern_bam": ["--defaultCoreScrna", "Yes", "--outSAMtype", "BAM", "SortedByCoordinate",
                       "--outSAMattributes", "NH", "HI", "AS", "nM", "CB", "UB",
                       "--limitBAMsortRAM", "1000000000"],
    }
    for profile, flags in profiles.items():
        if profile not in args.profiles:
            continue
        out = args.outdir / profile
        out.mkdir()
        print(f"Running plain GEX 100K: {profile}", flush=True)
        run(common + flags + ["--outFileNamePrefix", str(out) + "/"], out / "run.log")
        qc = {}
        for line in (out / "Log.final.out").read_text().splitlines():
            parts = line.split("|")
            if len(parts) == 2 and parts[0].strip() in QC_KEYS:
                qc[parts[0].strip()] = int(parts[1].strip())
        if qc.get("Number of input reads") != 100000 or qc.get("Uniquely mapped reads number", 0) < 10000:
            raise AssertionError(f"GEX input/mapping sanity failed: {qc}")
        result = {"qc": qc, "matrices": {}}
        for feature in ("Gene", "GeneFull"):
            raw_stats, raw = matrix_summary(out / f"Solo.out/{feature}/raw")
            if raw_stats["umis"] < 1000:
                raise AssertionError(f"empty or implausible raw {feature} matrix")
            filtered_dir = out / f"Solo.out/{feature}/filtered"
            if filtered_dir.is_dir():
                filtered_stats, filtered = matrix_summary(filtered_dir)
            else:
                with (out / f"Solo.out/{feature}/Summary.csv").open() as stream:
                    summary = dict(csv.reader(stream))
                if int(summary["Estimated Number of Cells"]) != 0:
                    raise AssertionError("missing filtered matrix with nonzero called cells")
                filtered_stats, filtered = {"not_exported": "zero_called_cells", "barcodes": 0, "umis": 0}, {}
            if any(raw.get(key) != value for key, value in filtered.items()):
                raise AssertionError("filtered counts are not a subset of raw counts")
            result["matrices"][feature] = {"raw": raw_stats, "filtered": filtered_stats}
        if result["matrices"]["GeneFull"]["raw"]["umis"] < result["matrices"]["Gene"]["raw"]["umis"]:
            raise AssertionError("GeneFull lost exonic molecules")
        if profile == "modern_bam":
            run(["samtools", "quickcheck", str(out / "Aligned.sortedByCoord.out.bam")], out / "bam-check.log")
        if profile == "modern_bam" and "modern" in report["profiles"] and result != report["profiles"]["modern"]:
            raise AssertionError("BAM output changed GEX matrices or cell calls")
        report["profiles"][profile] = result
        (args.outdir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        if expected is not None and result != expected["profiles"][profile]:
            raise AssertionError(f"{profile}: raw/filtered matrix, cell set, or mapping QC drift; inspect report.json")
    if args.record_golden is not None:
        with args.record_golden.open("x") as stream:
            stream.write(json.dumps(report, indent=2) + "\n")
        print(f"BASELINE RECORDED (not a regression pass): {args.record_golden}")
    else:
        verify_report(report, expected, args.profiles)
        (args.outdir / "PASS").write_text("Selected plain GEX profiles match the pinned baseline.\n")
        print(f"PASS: plain GEX 100K; report: {args.outdir / 'report.json'}")


if __name__ == "__main__":
    main()

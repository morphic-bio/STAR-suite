#!/usr/bin/env python3
"""Pinned 100K-read regression of the experimental STAR-perturb add-ons on public reads.

    python3 tests/run_perturb_addons_100k_regression.py [--profiles splitseq pipseq]

Two profiles, each one STAR run on the first 100,000 read pairs of a public
library with the barcode lists packaged in share/star-suite/barcode_lists:

  splitseq  --soloSplitseqMode Fixed, v3 short-linker layout, 48_wells lists;
            v1.13.0 counting policy (rescue off), GEO GSE352893, run SRR41469160
  pipseq    --soloPipseqMode Illumina, GeneFull_Ex50pAS, the read-derived tier
            lists; GEO GSE295097, run SRR33231853

The fixtures are made by tests/setup_perturb_addons_100k_fixture.py and are
not tracked. The golden file pins, per profile, the mapping counts of
Log.final.out, the raw and filtered matrices (sizes, totals and digests of
features, barcodes and counts), every row of the mode's metrics table
(non-integer numbers to nine significant digits) and, for the PIP mode, the
digest of Pipseq.groups.tsv. Any drift fails.

--record-golden PATH writes a new golden and never reports a pass.
"""
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
DEFAULT_GOLDEN = ROOT / "tests/fixtures/perturb_addons_100k_golden.json"
QC_KEYS = ("Number of input reads", "Uniquely mapped reads number",
           "Number of reads mapped to multiple loci", "Number of reads mapped to too many loci")
PROFILES = {
    "splitseq": {"fixture": "splitseq_gse352893", "dataset": "GSE352893", "run": "SRR41469160", "feature": "GeneFull",
                 "lists": ["splitseq_v3/48_wells/round1.tsv", "splitseq_v3/48_wells/round2.tsv", "splitseq_v3/48_wells/round3.tsv"],
                 "metrics": "Solo.out/Splitseq.metrics.tsv"},
    "pipseq": {"fixture": "pipseq_gse295097", "dataset": "GSE295097", "run": "SRR33231853", "feature": "GeneFull_Ex50pAS",
               "lists": ["pipseq_readderived/readderived_tier%d.tsv" % i for i in (1, 2, 3, 4)],
               "metrics": "Solo.out/GeneFull_Ex50pAS/Pipseq.metrics.tsv"},
}


def digest_rows(rows):
    digest = hashlib.sha256()
    for row in sorted(rows):
        digest.update((json.dumps(row, separators=(",", ":")) + "\n").encode())
    return digest.hexdigest()


def matrix_summary(path):
    features, barcodes, counts = read_mex(path)
    if any(value < 0 for value in counts.values()):
        raise AssertionError(f"negative count: {path}")
    return {
        "genes": len(features), "barcodes": len(barcodes), "nnz": len(counts),
        "counts": sum(counts.values()),
        "features_sha256": digest_rows(features), "barcodes_sha256": digest_rows(barcodes),
        "counts_sha256": digest_rows((gene, cb, count) for (gene, cb), count in counts.items()),
    }, counts


def normalized(value):
    """Metric values as written, except non-integer numbers, which are kept to nine significant digits."""
    if value.lstrip("-").isdigit():
        return value
    try:
        return format(float(value), ".9g")
    except ValueError:
        return value


def command(profile, star, genome, threads, fixture, lists, out):
    reads = ["--readFilesIn", str(fixture / "cDNA.fastq.gz"), str(fixture / "barcode.fastq.gz"), "--readFilesCommand", "zcat"]
    common = [str(star), "--genomeDir", str(genome), "--runThreadN", str(threads)] + reads
    tail = ["--soloStrand", "Forward", "--outSAMtype", "None", "--outFileNamePrefix", str(out) + "/"]
    if profile == "splitseq":
        return common + ["--soloSplitseqMode", "Fixed", "--soloSplitseqCBposition", "50_8", "30_8", "10_8",
                         "--soloSplitseqUMIposition", "0_10", "--soloSplitseqWhitelist", *lists,
                         "--soloCrMultimapRescue", "no"] + tail
    return common + ["--soloPipseqMode", "Illumina", "--soloPipseqWhitelist", *lists,
                     "--soloFeatures", "GeneFull_Ex50pAS"] + tail


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--star", type=Path, default=os.environ.get("STAR_BIN", ROOT / "core/legacy/source/STAR"))
    parser.add_argument("--fixtures", type=Path, default=os.environ.get("PERTURB_ADDONS_FIXTURES", "/storage/downsampled_100K"),
                        help="folder that holds splitseq_gse352893/ and pipseq_gse295097/")
    parser.add_argument("--genome", type=Path, default=os.environ.get("GENOME_DIR", "/storage/autoindex_110_44/bulk_index"))
    parser.add_argument("--lists", type=Path, default=os.environ.get("PERTURB_ADDONS_LISTS", ROOT / "share/star-suite/barcode_lists"),
                        help="the packaged barcode_lists folder")
    parser.add_argument("--outdir", type=Path, default=os.environ.get("PERTURB_ADDONS_OUTDIR"))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--golden", type=Path, default=DEFAULT_GOLDEN)
    parser.add_argument("--profiles", nargs="+", choices=tuple(PROFILES),
                        default=os.environ.get("PERTURB_ADDONS_PROFILES", "splitseq pipseq").split())
    parser.add_argument("--record-golden", type=Path, help="explicit bootstrap only; never reports a regression pass")
    args = parser.parse_args()
    for profile in args.profiles:
        if profile not in PROFILES:
            parser.error(f"unknown profile: {profile}")
    if args.threads < 1:
        parser.error("--threads must be positive")
    args.star = Path(args.star)
    if not args.star.is_file() or not os.access(args.star, os.X_OK):
        parser.error(f"STAR executable not found: {args.star}")
    if args.record_golden is not None and args.record_golden.exists():
        parser.error("refusing to overwrite an existing golden")
    if args.outdir is None:
        import tempfile
        args.outdir = Path(tempfile.mkdtemp(prefix="star-perturb-addons-100k-"))
    else:
        args.outdir = Path(args.outdir)
        args.outdir.mkdir(parents=True, exist_ok=False)
    args.star = args.star.resolve()
    genome = Path(args.genome)
    reference = {name: sha256(genome / name) for name in
                 ("genomeParameters.txt", "chrNameLength.txt", "geneInfo.tab", "exonInfo.tab")}
    expected = None if args.record_golden is not None else json.loads(Path(args.golden).read_text())
    report = {
        "schema": "perturb-addons-100k-regression-v1",
        "binary_sha256": sha256(args.star),
        "version": subprocess.check_output([str(args.star), "--version"], text=True).strip(),
        "source_revision": subprocess.check_output([str(args.star), "--source-revision"], text=True).strip(),
        "reference_metadata": reference,
        "profiles": {},
    }
    if expected is not None and expected["reference_metadata"] != reference:
        raise ValueError("the genome index differs from the one of the pinned golden")
    for profile in args.profiles:
        spec = PROFILES[profile]
        fixture = Path(args.fixtures) / spec["fixture"]
        manifest = json.loads((fixture / "manifest.json").read_text())
        if (manifest["schema"], manifest["assay"], manifest["run"], manifest["read_pairs"]) != \
                ("perturb-addons-100k-v1", profile, spec["run"], 100000):
            raise ValueError(f"{fixture} is not the prepared 100K fixture of {spec['run']}")
        for mate in ("cdna", "barcode"):
            if sha256(fixture / manifest[mate]) != manifest[f"{mate}_sha256"]:
                raise ValueError(f"fixture checksum mismatch: {fixture / manifest[mate]}")
        lists = [Path(args.lists) / name for name in spec["lists"]]
        inputs = {"dataset": spec["dataset"], "run": spec["run"], "read_pairs": manifest["read_pairs"],
                  "cdna_sha256": manifest["cdna_sha256"], "barcode_sha256": manifest["barcode_sha256"],
                  "lists_sha256": {name: sha256(path) for name, path in zip(spec["lists"], lists)}}
        if expected is not None and expected["profiles"][profile]["inputs"] != inputs:
            raise ValueError(f"{profile}: fixture or barcode lists differ from the pinned golden")
        out = args.outdir / profile
        out.mkdir()
        argv = command(profile, args.star, genome, args.threads, fixture, [str(p) for p in lists], out)
        (out / "run.command.json").write_text(json.dumps(argv, indent=2) + "\n")
        print(f"Running {profile} 100K", flush=True)
        with (out / "run.log").open("w") as stream:
            subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT, check=True)
        qc = {}
        for line in (out / "Log.final.out").read_text().splitlines():
            parts = line.split("|")
            if len(parts) == 2 and parts[0].strip() in QC_KEYS:
                qc[parts[0].strip()] = int(parts[1].strip())
        if qc.get("Number of input reads") != 100000 or qc.get("Uniquely mapped reads number", 0) < 10000:
            raise AssertionError(f"{profile}: input/mapping sanity failed: {qc}")
        solo = out / "Solo.out" / spec["feature"]
        raw_stats, raw = matrix_summary(solo / "raw")
        if raw_stats["counts"] < 1000:
            raise AssertionError(f"{profile}: empty or implausible raw matrix")
        filtered_stats, filtered = matrix_summary(solo / "filtered")
        if any(raw.get(key) != value for key, value in filtered.items()):
            raise AssertionError(f"{profile}: filtered counts are not a subset of raw counts")
        with (out / spec["metrics"]).open() as stream:
            metrics = {row[0]: normalized(row[1]) for row in csv.reader(stream, delimiter="\t") if len(row) == 2}
        if int(metrics["reads_total"]) != 100000:
            raise AssertionError(f"{profile}: the metrics table does not account for 100,000 reads")
        result = {"inputs": inputs, "qc": qc, "matrices": {"raw": raw_stats, "filtered": filtered_stats}, "metrics": metrics}
        if profile == "pipseq":
            result["groups_sha256"] = sha256(solo / "Pipseq.groups.tsv")
        report["profiles"][profile] = result
        (args.outdir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        if expected is not None and result != expected["profiles"][profile]:
            raise AssertionError(f"{profile}: matrix, cell set, metrics or mapping counts drifted; inspect {args.outdir / 'report.json'}")
    if args.record_golden is not None:
        with args.record_golden.open("x") as stream:
            stream.write(json.dumps(report, indent=2) + "\n")
        print(f"BASELINE RECORDED (not a regression pass): {args.record_golden}")
    else:
        (args.outdir / "PASS").write_text("Selected add-on profiles match the pinned baseline.\n")
        print(f"PASS: perturb add-ons 100K ({', '.join(args.profiles)}); report: {args.outdir / 'report.json'}")


if __name__ == "__main__":
    main()

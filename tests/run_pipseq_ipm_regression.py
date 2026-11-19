#!/usr/bin/env python3
"""Run the barcode-selected PIP calibration fixture against frozen ledger and MEX expectations."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

from compare_flex_hash_screen_mex import read_mex
from run_perturb_addons_100k_regression import PROFILES, ROOT, command
from setup_scrna_gex_100k_fixture import sha256


GOLDEN = ROOT / "tests/fixtures/pipseq_ipm_golden.json"
COUNT_FIELDS = ("reads", "occupied_bins", "unique_bi_imi", "estimated_transcripts")


def read_groups(path, expected=False):
    fields = COUNT_FIELDS + (("transcripts_ipm_1",) if expected else ())
    groups = {}
    with path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            key = (row["gene_id"], row["barcode"])
            if key in groups:
                raise ValueError(f"duplicate group {key} in {path}")
            groups[key] = tuple(int(row[f]) for f in fields)
    return groups


def verify(solo, expected_file, golden):
    expected = read_groups(expected_file, expected=True)
    groups = read_groups(solo / "Pipseq.groups.tsv")
    if expected.keys() != groups.keys():
        raise ValueError(f"ledger keys differ: {len(expected.keys() - groups.keys())} missing, "
                         f"{len(groups.keys() - expected.keys())} unexpected")
    reference = golden["expected"]
    if len(groups) != reference["groups"]:
        raise ValueError("unexpected number of cell/gene groups")
    sensitive = 0
    for key, (reads, bins, fragments, count, floor_count) in expected.items():
        if groups[key] != (reads, bins, fragments, count):
            raise ValueError(f"ledger mismatch at {key}: got {groups[key]}, "
                             f"expected {(reads, bins, fragments, count)}")
        sensitive += count != floor_count
    if not sensitive or sensitive != reference["groups_whose_count_differs_from_ipm_1"]:
        raise ValueError(f"unexpected number of IPM-sensitive groups: {sensitive}")

    raw = solo / "raw"
    if sha256(raw / "features.tsv") != golden["features_sha256"]:
        raise ValueError("unexpected raw matrix feature axis")
    _, barcodes, counts = read_mex(raw)
    if set(barcodes) != {bc for _, bc in expected} or len(barcodes) != reference["barcodes"]:
        raise ValueError("unexpected raw matrix barcode axis")
    if counts.keys() != expected.keys():
        raise ValueError("raw matrix coordinates differ from expected cell/gene keys")
    for key, count in counts.items():
        if count != expected[key][3]:
            raise ValueError(f"raw matrix mismatch at {key}: got {count}, expected {expected[key][3]}")

    with (solo / "Pipseq.metrics.tsv").open() as stream:
        metrics = dict(csv.reader(stream, delimiter="\t"))
    required = {"ipm_source": "estimated", "ipm_floor_applied": "0",
                "ipm_calibration_groups": str(reference["calibration_groups_5_to_32_bins"]),
                "reads_total": str(golden["read_pairs"]), "reads_counted": str(reference["counted_reads"]),
                "estimated_transcripts": str(reference["estimated_transcripts"])}
    for name, value in required.items():
        if metrics[name] != value:
            raise ValueError(f"{name}: got {metrics[name]}, expected {value}")
    for name in ("ipm_raw", "ipm_effective"):
        if not math.isclose(float(metrics[name]), reference["ipm_raw"], rel_tol=1e-12):
            raise ValueError(f"{name}: got {metrics[name]}, expected {reference['ipm_raw']}")
    return {"groups": len(groups), "matrix_entries": len(counts), "ipm_sensitive_groups": sensitive,
            "ipm_raw": metrics["ipm_raw"], "ipm_effective": metrics["ipm_effective"],
            **{k: metrics[k] for k in required}}


def check_files(directory, hashes):
    for name, digest in hashes.items():
        if sha256(directory / name) != digest:
            raise ValueError(f"fixture or reference checksum mismatch: {directory / name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, default=os.environ.get("STAR_BIN", ROOT / "core/legacy/source/STAR"))
    parser.add_argument("--fixture", type=Path, default=os.environ.get("PIPSEQ_IPM_FIXTURE", "/storage/pipseq_v5/GSE295097/ipm_fixture_20cb_v1"))
    parser.add_argument("--genome", type=Path, default=os.environ.get("PIPSEQ_IPM_GENOME", "/storage/references_public/pipseq_gencode_v44/index"))
    parser.add_argument("--lists", type=Path, default=os.environ.get("PERTURB_ADDONS_LISTS", ROOT / "share/star-suite/barcode_lists"))
    parser.add_argument("--outdir", type=Path, default=os.environ.get("OUT_ROOT"))
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    star = args.star.resolve()
    if not star.is_file() or not os.access(star, os.X_OK):
        parser.error(f"missing STAR executable: {star}")
    golden = json.loads(GOLDEN.read_text())
    check_files(args.fixture, golden["fixture_sha256"])
    check_files(args.genome, golden["index_sha256"])
    check_files(args.lists, golden["lists_sha256"])
    if args.outdir is None:
        out = Path(tempfile.mkdtemp(prefix="star-pipseq-ipm-"))
    else:
        out = args.outdir.resolve()
        out.mkdir(parents=True, exist_ok=False)
    lists = [str(args.lists.resolve() / name) for name in PROFILES["pipseq"]["lists"]]
    argv = command("pipseq", star, args.genome.resolve(), args.threads, args.fixture.resolve(), lists, out)
    argv += ["--soloPipseqIPM", "0", "--soloPipseqTrim", "Kit", "--soloCellFilter", "None"]
    (out / "run.command.json").write_text(json.dumps(argv, indent=2) + "\n")
    with (out / "run.log").open("w") as log:
        subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=True)
    solo = out / "Solo.out" / PROFILES["pipseq"]["feature"]
    report = verify(solo, args.fixture / "expected.tsv", golden)
    report["source_revision"] = subprocess.check_output([str(star), "--source-revision"], text=True).strip()
    report["binary_sha256"] = sha256(star)
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "PASS").write_text("The ledger, calibrated counts and raw matrix match the frozen expectations.\n")
    print(f"PASS: {report['groups']} groups, {report['ipm_sensitive_groups']} IPM-sensitive ({out})")


if __name__ == "__main__":
    main()

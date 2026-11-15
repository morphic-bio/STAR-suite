#!/usr/bin/env python3
"""Additional PIP integration cases using an existing synthetic smoke fixture."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--star", required=True, type=Path)
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--case", choices=["all", "no-cells"], default="all")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    fixture = args.fixture.resolve()
    commands = []
    base = [args.star.resolve(), "--genomeDir", fixture / "index",
            "--soloPipseqMode", "Illumina", "--soloPipseqWhitelist",
            *[fixture / ("tier%d.txt" % i) for i in range(4)],
            "--outSAMtype", "None", "--soloStrand", "Forward", "--alignIntronMax", "1000",
            "--outFilterMatchNminOverLread", "0.5", "--outFilterScoreMinOverLread", "0.5"]

    def run(label, r2, r1, extra, success=True):
        target = args.out / label
        target.mkdir()
        argv = [str(v) for v in base + ["--readFilesIn", r2, r1,
                "--outFileNamePrefix", str(target.resolve()) + "/"] + extra]
        commands.append({"label": label, "argv": argv})
        (args.out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        with (target / "stdout.log").open("w") as log:
            result = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
        if success != (result.returncode == 0):
            raise AssertionError((target / "stdout.log").read_text()[-4000:])
        return target / "Solo.out/GeneFull"

    original_r2, original_r1 = fixture / "R2.fastq", fixture / "R1.fastq"
    # Explicitly reject every detected cell to exercise empty filtered output
    # with a nonempty raw matrix (different from having no counted reads).
    no_cells = run("no_cells", original_r2, original_r1,
                   ["--soloCellFilter", "CellRanger2.2", "1", "0.99", "0.01", "--runThreadN", "1"])
    rows = [x for x in (no_cells / "filtered/matrix.mtx").read_text().splitlines() if not x.startswith("%")]
    assert rows == ["2 0 0"], rows
    assert (no_cells / "filtered/barcodes.tsv").read_text() == ""
    if args.case == "no-cells":
        (args.out / "RESULTS.json").write_text(json.dumps({"status": "PASS", "case": "no-cells"}, indent=2) + "\n")
        print("PASS nonempty raw matrix with zero called cells")
        return
    # Exercise the shared EmptyDrops adapter on estimated-transcript counts.
    ed = run("emptydrops", original_r2, original_r1, ["--soloCellFilter", "EmptyDrops_CR", "--runThreadN", "1"])
    assert (ed / "filtered/matrix.mtx").exists()
    assert (ed / "filtered/barcodes.tsv").read_text().splitlines() == ["ACGTACGTGACTGCCTAGCATGCATGCA"]
    # Shallow input with no assigned barcode should yield valid empty outputs.
    empty_r1, empty_r2 = args.out / "invalid_R1.fastq", args.out / "invalid_R2.fastq"
    sequence = original_r2.read_text().splitlines()[1]
    empty_r1.write_text("@invalid\n" + "N" * 42 + "\n+\n" + "I" * 42 + "\n")
    empty_r2.write_text("@invalid\n" + sequence + "\n+\n" + "I" * len(sequence) + "\n")
    empty = run("empty", empty_r2, empty_r1, ["--soloCellFilter", "TopCells", "1", "--runThreadN", "1"])
    for kind in ("raw", "filtered"):
        rows = [x for x in (empty / kind / "matrix.mtx").read_text().splitlines() if not x.startswith("%")]
        assert rows == ["2 0 0"], rows
        assert (empty / kind / "barcodes.tsv").read_text() == ""
    # Cross-chunk/thread deduplication: duplicate input far exceeds a small
    # input chunk; all copies must retain the original fragment cardinalities.
    repeated_r1, repeated_r2 = args.out / "repeated_R1.fastq", args.out / "repeated_R2.fastq"
    for source, dest in ((original_r1, repeated_r1), (original_r2, repeated_r2)):
        data = source.read_text()
        with dest.open("w") as out:
            for _ in range(500): out.write(data)
    threaded = run("threaded", repeated_r2, repeated_r1,
                   ["--runThreadN", "4", "--limitIObufferSize", "240000", "1000000",
                    "--soloPipseqIPM", "1.35", "--soloCellFilter", "None"])
    ledger = list(csv.DictReader((threaded / "Pipseq.groups.tsv").open(), delimiter="\t"))
    assert {r["gene_id"]: int(r["unique_bi_imi"]) for r in ledger} == {"PLUS": 22, "MINUS": 20}, ledger
    assert {r["gene_id"]: int(r["estimated_transcripts"]) for r in ledger} == {"PLUS": 16, "MINUS": 10}, ledger
    metrics = dict(csv.reader((threaded / "Pipseq.metrics.tsv").open(), delimiter="\t"))
    assert int(metrics["reads_counted"]) == 45000, metrics
    assert int(metrics["reads_total"]) == 46500, metrics
    # All 64 bins but no eligible IPM groups: fail clearly instead of choosing
    # an undocumented fallback. An explicit IPM permits the same data.
    body = "ACGTACGTATGGACTGCGAGCTAGCATCGAGTGCATGCA"
    saturated_r1, saturated_r2 = args.out / "saturated_R1.fastq", args.out / "saturated_R2.fastq"
    r1, r2 = [], []
    for bin_index in range(64):
        bi = "".join("ACGT"[(bin_index >> shift) & 3] for shift in (4, 2, 0))
        r1.append("@sat%d\n%s\n+\n%s\n" % (bin_index, body + bi, "I" * 42))
        r2.append("@sat%d\n%s\n+\n%s\n" % (bin_index, sequence, "I" * len(sequence)))
    saturated_r1.write_text("".join(r1)); saturated_r2.write_text("".join(r2))
    run("missing_calibration", saturated_r2, saturated_r1, ["--runThreadN", "1"], success=False)
    assert "No PIP cell/gene groups with 5-32 bins" in (args.out / "missing_calibration/stdout.log").read_text()
    supplied = run("saturated_supplied", saturated_r2, saturated_r1,
                   ["--runThreadN", "1", "--soloPipseqIPM", "2"])
    ledger = list(csv.DictReader((supplied / "Pipseq.groups.tsv").open(), delimiter="\t"))
    assert len(ledger) == 1 and int(ledger[0]["estimated_transcripts"]) == 64, ledger
    receipt = {"status": "PASS", "star_sha256": hashlib.sha256(args.star.read_bytes()).hexdigest(),
               "cases": ["EmptyDrops_CR", "empty_matrix", "cross_thread_dedup", "missing_calibration", "saturation"]}
    (args.out / "RESULTS.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS PIP EmptyDrops, empty outputs, threaded merging and calibration edges")


if __name__ == "__main__":
    main()

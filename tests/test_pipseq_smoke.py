#!/usr/bin/env python3
"""Synthetic end-to-end PIP RNA acceptance fixture; no external data required."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    commands = []

    def run(label, argv, success=True):
        commands.append({"label": label, "argv": [str(x) for x in argv]})
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        with (out / (label + ".stdout")).open("w") as log:
            result = subprocess.run([str(x) for x in argv], stdout=log, stderr=subprocess.STDOUT)
        if success and result.returncode:
            raise AssertionError(label + " failed: " + (out / (label + ".stdout")).read_text()[-4000:])
        if not success and result.returncode == 0:
            raise AssertionError(label + " unexpectedly succeeded")
        return (out / (label + ".stdout")).read_text()

    rng = random.Random(7319)
    genome = "".join(rng.choice("ACGT") for _ in range(8000))
    fasta = out / "genome.fa"
    fasta.write_text(">chrSynthetic\n" + genome + "\n")
    gtf = out / "genes.gtf"
    gtf.write_text(
        'chrSynthetic\tfixture\texon\t501\t3500\t.\t+\t.\tgene_id "PLUS"; transcript_id "PLUS.1"; gene_name "PLUS";\n'
        'chrSynthetic\tfixture\texon\t4501\t7500\t.\t-\t.\tgene_id "MINUS"; transcript_id "MINUS.1"; gene_name "MINUS";\n')
    index = out / "index"
    index.mkdir()
    run("index", [args.star, "--runMode", "genomeGenerate", "--genomeDir", index,
                  "--genomeFastaFiles", fasta, "--sjdbGTFfile", gtf, "--sjdbOverhang", "49",
                  "--genomeSAindexNbases", "3", "--genomeChrBinNbits", "12",
                  "--runThreadN", "1", "--outFileNamePrefix", str(out / "index_")])
    tiers = ["ACGTACGT", "GACTGC", "CTAGCA", "TGCATGCA"]
    lists = []
    for i, sequence in enumerate(tiers):
        path = out / ("tier%d.txt" % i)
        path.write_text(sequence + "\n")
        lists.append(path)
    body = "ACGTACGTATGGACTGCGAGCTAGCATCGAGTGCATGCA"
    cell = "".join(tiers)
    records = []
    def bi(n):
        return "".join("ACGT"[(n >> s) & 3] for s in (4, 2, 0))
    def add(sequence, bin_index=0, phase=0, barcode_body=body, suffix=None):
        records.append(("r%d" % len(records), sequence,
                        "G" * phase + barcode_body + (bi(bin_index) if suffix is None else suffix) + "TTTT"))
    # PLUS: 11 bins, two starts per bin, three read copies per fragment.
    for b in range(11):
        for start in (700, 900):
            for copy in range(3):
                add(genome[start:start + 50], b, (b + copy) % 4)
    # Same aligned start with different lengths/clipping remains one IMI.
    add(genome[700:755], 0)
    add("N" * 5 + genome[700:750], 0)
    # MINUS: ten bins and two 5' endpoints, duplicate lengths sharing an end.
    complement = str.maketrans("ACGT", "TGCA")
    rc = lambda s: s.translate(complement)[::-1]
    for b in range(10):
        for end in (5100, 5300):
            add(rc(genome[end - 50:end]), b, b % 4)
    add(rc(genome[5045:5100]), 0)
    changed = list(body)
    for start in (0, 11, 20, 31):
        changed[start] = "A" if changed[start] != "A" else "C"
    add(genome[700:750], 0, barcode_body="".join(changed))
    add(genome[700:750], suffix="NNN") # retain barcode, reject invalid BI
    add(genome[700:750], barcode_body="N" * 39)
    add("N" * 50) # unmapped but valid barcode/BI
    r2, r1 = out / "R2.fastq", out / "R1.fastq"
    for path, column in ((r2, 1), (r1, 2)):
        path.write_text("".join("@%s\n%s\n+\n%s\n" % (row[0], row[column], "I" * len(row[column])) for row in records))

    mode = [args.star, "--genomeDir", index,
            "--soloPipseqMode", "Illumina", "--soloPipseqWhitelist", *lists,
            "--outSAMtype", "None", "--soloStrand", "Forward",
            "--alignIntronMax", "1000", "--outFilterMatchNminOverLread", "0.5",
            "--outFilterScoreMinOverLread", "0.5"]
    # These cDNA reads carry no technical sequence: the mode's trimming is off
    # here and exercised separately below.
    common = mode + ["--readFilesIn", r2, r1, "--soloPipseqTrim", "None"]

    def matrix(run_dir, feature="GeneFull"):
        raw = run_dir / "Solo.out" / feature / "raw"
        assert (raw / "barcodes.tsv").read_text().splitlines() == [cell]
        features = [x.split("\t")[0] for x in (raw / "features.tsv").read_text().splitlines()]
        lines = [x for x in (raw / "matrix.mtx").read_text().splitlines() if not x.startswith("%")]
        assert lines[0] == "2 1 2"
        return {features[int(g) - 1]: int(n) for g, c, n in (x.split() for x in lines[1:])}

    receipts = {}
    for label, threads, ipm, filtering in (("fixed", 1, "1.35", ["None"]),
                                           ("automatic", 2, "0", ["TopCells", "1"])):
        target = out / label
        target.mkdir()
        run(label, common + ["--runThreadN", str(threads), "--soloPipseqIPM", ipm,
                            "--soloCellFilter", *filtering, "--outFileNamePrefix", str(target) + "/"])
        counts = matrix(target)
        ledger = list(csv.DictReader((target / "Solo.out/GeneFull/Pipseq.groups.tsv").open(), delimiter="\t"))
        assert {r["gene_id"]: int(r["unique_bi_imi"]) for r in ledger} == {"PLUS": 22, "MINUS": 20}, ledger
        assert {r["gene_id"]: int(r["occupied_bins"]) for r in ledger} == {"PLUS": 11, "MINUS": 10}
        metrics = dict(csv.reader((target / "Solo.out/GeneFull/Pipseq.metrics.tsv").open(), delimiter="\t"))
        assert metrics["barcode_corrected"] == "1", metrics
        assert metrics["barcode_correction"] == "exact_then_one_substitution_per_tier", metrics
        assert metrics["cdna_trim"] == "none", metrics
        exclusive = ("barcode_exact", "barcode_corrected", "barcode_invalid", "barcode_nonmatch",
                     "barcode_phase_tie", "barcode_ambiguous")
        assert sum(int(metrics[key]) for key in exclusive) == int(metrics["reads_total"]), metrics
        assert (target / "Solo.out/GeneFull/Summary.csv").read_text().startswith("Number of Reads,"), label
        assert metrics["bi_invalid_with_valid_barcode"] == "1", metrics
        assert metrics["reads_unmapped"] == "1", metrics
        if label == "fixed":
            assert counts == {"PLUS": 16, "MINUS": 10}, counts
        else:
            assert metrics["ipm_calibration_groups"] == "2"
            assert counts == {"PLUS": 11, "MINUS": 10}, counts
            filtered = target / "Solo.out/GeneFull/filtered"
            assert (filtered / "matrix.mtx").read_bytes() == (target / "Solo.out/GeneFull/raw/matrix.mtx").read_bytes()
        receipts[label] = {"counts": counts, "metrics": metrics}
    # The mode's own cDNA trimming (the default): every read gets the leading
    # base; some get the template-switch oligo, whole or as its last 14 bases;
    # some run through poly-A into further sequence. Trimmed, they are the reads
    # above, so the ledger and the matrix must not change.
    tso = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
    kit_r2 = out / "kit_R2.fastq"
    with kit_r2.open("w") as handle:
        for i, row in enumerate(records):
            sequence = "T" + (tso, tso[-14:], "")[i % 3] + row[1] + ("A" * 20 + "GCTTAGGCATCG" if i % 2 == 0 else "")
            handle.write("@%s\n%s\n+\n%s\n" % (row[0], sequence, "I" * len(sequence)))
    target = out / "kit_trim"
    target.mkdir()
    run("kit_trim", mode + ["--readFilesIn", kit_r2, r1, "--runThreadN", "1", "--soloPipseqIPM", "1.35",
                            "--soloCellFilter", "None", "--outFileNamePrefix", str(target) + "/"])
    for name in ("Pipseq.groups.tsv", "raw/matrix.mtx", "raw/barcodes.tsv"):
        assert (target / "Solo.out/GeneFull" / name).read_bytes() == (out / "fixed/Solo.out/GeneFull" / name).read_bytes(), name
    kit_metrics = dict(csv.reader((target / "Solo.out/GeneFull/Pipseq.metrics.tsv").open(), delimiter="\t"))
    assert kit_metrics["cdna_trim"] == "kit", kit_metrics
    assert "[PIP effective] clipAdapterType = PIPseq" in (target / "Log.out").read_text()
    for key in ("reads_counted", "reads_unmapped", "reads_total"):
        assert kit_metrics[key] == receipts["fixed"]["metrics"][key], (key, kit_metrics)
    # Untrimmed, the same reads lose alignments: the trimming is doing the work.
    target = out / "kit_untrimmed"
    target.mkdir()
    run("kit_untrimmed", mode + ["--readFilesIn", kit_r2, r1, "--soloPipseqTrim", "None", "--runThreadN", "1",
                                 "--soloPipseqIPM", "1.35", "--soloCellFilter", "None",
                                 "--outFileNamePrefix", str(target) + "/"])
    untrimmed = dict(csv.reader((target / "Solo.out/GeneFull/Pipseq.metrics.tsv").open(), delimiter="\t"))
    assert int(untrimmed["reads_counted"]) < int(kit_metrics["reads_counted"]), untrimmed
    receipts["kit_trim"] = {"metrics": kit_metrics, "reads_counted_untrimmed": untrimmed["reads_counted"]}
    # BAM output: CR and UR as read, CB after correction, written with each
    # alignment; counting is unchanged.
    at = mode.index("--outSAMtype")
    mode_bam = mode[:at] + mode[at + 2:]
    for kind, name in (("Unsorted", "Aligned.out.bam"), ("SortedByCoordinate", "Aligned.sortedByCoord.out.bam")):
        target = out / ("bam_" + kind.lower())
        target.mkdir()
        run("bam_" + kind.lower(), mode_bam + ["--readFilesIn", kit_r2, r1, "--runThreadN", "2", "--soloPipseqIPM", "1.35",
                                        "--soloCellFilter", "None", "--outSAMtype", "BAM", kind,
                                        "--outSAMattributes", "NH", "HI", "AS", "nM", "CR", "CB", "UR", "GX", "GN", "sS",
                                        "--outSAMunmapped", "Within", "--outFileNamePrefix", str(target) + "/"])
        for item in ("Pipseq.groups.tsv", "raw/matrix.mtx"):
            assert (target / "Solo.out/GeneFull" / item).read_bytes() == (out / "fixed/Solo.out/GeneFull" / item).read_bytes(), (kind, item)
        if shutil.which("samtools"):
            view = subprocess.run(["samtools", "view", str(target / name)], check=True, capture_output=True, text=True).stdout
            rows = [dict(f.split(":", 2)[::2] for f in line.split("\t")[11:]) | {"name": line.split("\t")[0], "flag": int(line.split("\t")[1])}
                    for line in view.splitlines()]
            assert len(rows) == len(records), (kind, len(rows))
            barcodes = {row[0]: row[2] for row in records}
            for row in rows:
                read = barcodes[row["name"]]
                assert row["sS"] == read, row
                phase = len(read) - len(read.lstrip("G")) if read.startswith("G") and not read.startswith("GACT") else 0
                if "N" * 39 in read:
                    assert row["CR"] == "-" and row["CB"] == "-" and row["UR"] == "-", row
                    continue
                assert row["CB"] == cell and len(row["CR"]) == 28 and row["UR"] == read[phase + 39:phase + 42], row
                if not row["flag"] & 4 and row["NH"] == "1":
                    assert row["GX"] in ("PLUS", "MINUS"), row
            assert sum(1 for row in rows if row["CR"] != row["CB"] and row["CB"] == cell) == 1  # the corrected read
            receipts["bam_" + kind.lower()] = {"records": len(rows)}
    ub = run("ub_undefined", mode_bam + ["--readFilesIn", kit_r2, r1, "--outSAMtype", "BAM", "Unsorted",
                                     "--outSAMattributes", "NH", "CB", "UB", "--outFileNamePrefix", str(out / "ub_")], success=False)
    assert "not defined in PIP mode" in ub
    clip_conflict = run("clip_conflict", mode + ["--readFilesIn", kit_r2, r1, "--clipAdapterType", "Hamming",
                                                 "--outFileNamePrefix", str(out / "clip_conflict_")], success=False)
    assert "conflicts with PIP mode" in clip_conflict
    # Exon-priority gene assignment is accepted and counts the same fully exonic reads.
    target = out / "ex50pas"
    target.mkdir()
    run("ex50pas", common + ["--soloFeatures", "GeneFull_Ex50pAS", "--runThreadN", "1", "--soloPipseqIPM", "1.35",
                             "--soloCellFilter", "None", "--outFileNamePrefix", str(target) + "/"])
    assert matrix(target, "GeneFull_Ex50pAS") == {"PLUS": 16, "MINUS": 10}
    conflict = run("conflict", common + ["--soloUMIdedup", "1MM_CR", "--outFileNamePrefix", str(out / "conflict_")], success=False)
    assert "conflicts with PIP mode" in conflict
    bad_ipm = run("bad_ipm", common + ["--soloPipseqIPM", "0.5", "--outFileNamePrefix", str(out / "bad_ipm_")], success=False)
    assert "finite value >= 1" in bad_ipm
    # A mode-off fixture exercises existing simple-barcode counting on this binary.
    ordinary_r1 = out / "ordinary_R1.fastq"
    ordinary_r1.write_text("".join("@%s\n%s\n+\n%s\n" % (r[0], "ACGTACGTACGTACGT" + "ACGTACGTAC", "I" * 26) for r in records))
    ordinary_wl = out / "ordinary_wl.txt"
    ordinary_wl.write_text("ACGTACGTACGTACGT\n")
    target = out / "ordinary"
    target.mkdir()
    run("ordinary", [args.star, "--genomeDir", index, "--readFilesIn", r2, ordinary_r1,
                     "--soloType", "CB_UMI_Simple", "--soloCBwhitelist", ordinary_wl,
                     "--soloFeatures", "GeneFull", "--soloUMIdedup", "Exact", "--soloCellFilter", "None",
                     "--outSAMtype", "None", "--runThreadN", "1", "--outFileNamePrefix", str(target) + "/"])
    lines = (target / "Solo.out/GeneFull/raw/matrix.mtx").read_text().splitlines()
    entries = [line for line in lines if line and not line.startswith("%")][1:]
    assert sorted(int(line.split()[2]) for line in entries) == [1, 1], entries
    assert not (target / "Solo.out/GeneFull/Pipseq.metrics.tsv").exists()
    receipts["status"] = "PASS"
    receipts["star_sha256"] = hashlib.sha256(args.star.read_bytes()).hexdigest()
    (out / "RESULTS.json").write_text(json.dumps(receipts, indent=2) + "\n")
    print("PASS PIP RNA end-to-end fixture, cDNA trimming and mode-off counting")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Synthetic split-pool FASTQ -> standard UMI counts -> default EmptyDrops."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--resume-receipt", type=Path,
                        help="reuse completed cases from a retained validation manifest with the same STAR hash")
    args = parser.parse_args()
    star, out = args.star.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    commands, results = [], {}
    previous, prior_commands = None, {}
    if args.resume_receipt:
        receipt = json.loads(args.resume_receipt.read_text())
        assert receipt["star_sha256"] == hashlib.sha256(star.read_bytes()).hexdigest(), "STAR binary changed"
        previous = Path(receipt["artifacts"]) / "smoke"
        prior_commands = {item["label"]: item for item in json.loads((previous / "commands.json").read_text())}

    def run(label, argv, error=None):
        commands.append({"label": label, "argv": [str(x) for x in argv]})
        prior = prior_commands.get(label)
        if prior and prior.get("exit_code") == 0 and label in ("index", "default", "threads4", "exact"):
            normalized = [s.replace(str(previous), str(out)) for s in prior["argv"]]
            assert normalized == commands[-1]["argv"], "reused command arguments changed: " + label
            target = out / label
            target.rmdir()  # the caller created an empty destination
            shutil.copytree(previous / label, target)
            commands[-1].update(exit_code=0, reused_from=str(previous / label))
            (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
            return
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        with (out / (label + ".stdout")).open("w") as log:
            proc = subprocess.run(commands[-1]["argv"], stdout=log, stderr=subprocess.STDOUT)
        commands[-1]["exit_code"] = proc.returncode
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        text = (out / (label + ".stdout")).read_text()
        if error:
            assert proc.returncode and error in text, (label, proc.returncode, text[-4000:])
        else:
            assert proc.returncode == 0, (label, text[-4000:])

    rng = random.Random(2703)
    genome = "".join(rng.choice("ACGT") for _ in range(8000))
    (out / "genome.fa").write_text(">chrSynthetic\n" + genome + "\n")
    (out / "genes.gtf").write_text(
        'chrSynthetic\tfixture\texon\t501\t3500\t.\t+\t.\tgene_id "PLUS"; transcript_id "PLUS.1"; gene_name "PLUS";\n'
        'chrSynthetic\tfixture\texon\t4501\t7500\t.\t-\t.\tgene_id "MINUS"; transcript_id "MINUS.1"; gene_name "MINUS";\n')
    index = out / "index"
    index.mkdir()
    run("index", [star, "--runMode", "genomeGenerate", "--genomeDir", index,
                  "--genomeFastaFiles", out / "genome.fa", "--sjdbGTFfile", out / "genes.gtf",
                  "--sjdbOverhang", "49", "--genomeSAindexNbases", "3", "--genomeChrBinNbits", "12",
                  "--runThreadN", "1", "--outFileNamePrefix", str(out / "index_")])
    # Same final 16 bases for both cells: truncation to a 32-bit CB would merge them.
    lists = [out / ("round%d.tsv" % i) for i in (1, 2, 3)]
    lists[0].write_text("ACGTACGT TGCATGCA\nTGCATGCA TGCATGCA\nGCTAGCTA GCTAGCTA\n")
    lists[1].write_text("GACTGACT\n")
    lists[2].write_text("CTAGCTAG\nCTAGCTAA\n")
    cell_a, cell_b = "TGCATGCAGACTGACTCTAGCTAG", "GCTAGCTAGACTGACTCTAGCTAG"
    u1, u1_error, u2 = "ACGTACGTAC", "ACGTACGTAT", "GCTAGCTAGC"
    plus = genome[700:750]
    minus = genome[5000:5050].translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def barcode(first="ACGTACGT", second="GACTGACT", third="CTAGCTAG", umi=u1):
        bases = list("A" * 86)
        for offset, sequence in ((0, umi), (10, third), (48, second), (78, first)):
            bases[offset:offset + len(sequence)] = sequence
        return "".join(bases)

    records = [(plus, barcode(), cell_a, u1)] * 4
    records += [(plus, barcode(umi=u1_error), cell_a, u1_error),
                (plus, barcode(umi=u2), cell_a, u2),
                (plus, barcode(first="TGCATGCA"), cell_a, u1),  # replacement alias
                (plus, barcode(first="CCGTACGT", second="AACTGACT", third="ATAGCTAG"), cell_a, u1),
                (minus, barcode(umi=u2), cell_a, u2),
                (plus, barcode(first="GCTAGCTA"), cell_b, u1)]
    # Include a cell with only corrected evidence: no whole-tuple exact anchor is required.
    cell_c = "GCTAGCTAGACTGACTCTAGCTAA"
    records += [(plus, barcode(first="ACTAGCTA", third="CTAGCTAA"), cell_c, u2)]
    invalid = [(plus, barcode(third="CTAGCTAC"), None, None),  # equally close to two tier entries
               (plus, barcode(first="NCGTACGT"), None, None),
               (plus, barcode(second="AAAAAAAA"), None, None),
               (plus, barcode()[:85], None, None),
               (plus, barcode(umi="N" + u1[1:]), None, None),
               (plus, barcode(umi="AAAAAAAAAA"), None, None)]
    records += invalid + [("N" * 50, barcode(), None, None)]
    # Multiple input chunks exercise independent per-thread decoding and sparse key remapping.
    records *= 1500

    def fastqs(prefix, rows, flat=False):
        rna, cb = out / (prefix + "_RNA.fastq"), out / (prefix + "_CB.fastq")
        with rna.open("w") as a, cb.open("w") as b:
            for i, (seq, bar, canonical, umi) in enumerate(rows):
                if flat:
                    bar = canonical + umi
                for stream, sequence in ((a, seq), (b, bar)):
                    stream.write("@r%d\n%s\n+\n%s\n" % (i, sequence, "I" * len(sequence)))
        return rna, cb

    rna, cb = fastqs("split", records)
    if previous:
        for path in (rna, cb, *lists, out / "genome.fa", out / "genes.gtf"):
            assert hashlib.sha256(path.read_bytes()).digest() == hashlib.sha256((previous / path.name).read_bytes()).digest(), "reused fixture input changed"
    base = [star, "--genomeDir", index, "--outSAMtype", "None", "--alignIntronMax", "1000",
            "--soloStrand", "Forward", "--limitIObufferSize", "1000000", "1000000"]
    split = ["--soloSplitseqMode", "Fixed", "--soloSplitseqCBposition", "78_8", "48_8", "10_8",
             "--soloSplitseqUMIposition", "0_10", "--soloSplitseqWhitelist", *lists]
    common = base + ["--readFilesIn", rna, cb] + split

    def matrix(directory):
        barcodes = (directory / "barcodes.tsv").read_text().splitlines()
        genes = [s.split("\t")[0] for s in (directory / "features.tsv").read_text().splitlines()]
        lines = [s for s in (directory / "matrix.mtx").read_text().splitlines() if not s.startswith("%")]
        dims = [int(s) for s in lines[0].split()]
        assert dims == [len(genes), len(barcodes), len(lines) - 1], dims
        return {(barcodes[int(c)-1], genes[int(g)-1]): int(n)
                for g, c, n in (s.split() for s in lines[1:])}

    expected = {(cell_a, "PLUS"): 2, (cell_a, "MINUS"): 1, (cell_b, "PLUS"): 1, (cell_c, "PLUS"): 1}
    for label, extra in (("default", ["--runThreadN", "1"]),
                         ("threads4", ["--runThreadN", "4"]),
                         ("exact", ["--runThreadN", "2", "--soloUMIdedup", "Exact", "--soloCellFilter", "None"]),
                         # A deterministic test-only threshold forces no calls. EmptyDrops defaults
                         # are exercised above; its bootstrap deliberately supersedes legacy knee settings.
                         ("zero_called", ["--runThreadN", "1", "--soloCellFilter", "CellRanger2.2", "1", "0.99", "0.01"])):
        target = out / label
        target.mkdir()
        run(label, common + extra + ["--outFileNamePrefix", str(target) + "/"])
        actual = matrix(target / "Solo.out/GeneFull/raw")
        want = dict(expected)
        if label == "exact":
            want[(cell_a, "PLUS")] = 3
        assert actual == want, (label, actual, want)
        metrics = dict(csv.reader((target / "Solo.out/Splitseq.metrics.tsv").open(), delimiter="\t"))
        for key in ("read_too_short", "barcode_non_acgt", "barcode_no_match", "barcode_ambiguous", "umi_non_acgt", "umi_homopolymer"):
            assert metrics[key] == "1500", (label, key, metrics)
        assert metrics["barcode_corrected"] == "3000", metrics
        exclusive = ("barcode_exact", "barcode_corrected", "read_too_short", "barcode_non_acgt",
                     "barcode_no_match", "barcode_ambiguous", "umi_non_acgt", "umi_homopolymer")
        assert sum(int(metrics[key]) for key in exclusive) == int(metrics["reads_total"]), metrics
        assert metrics["round1_corrected"] == "3000" and metrics["round2_corrected"] == "1500" and metrics["round3_corrected"] == "1500"
        filtered = target / "Solo.out/GeneFull/filtered"
        if label == "exact":
            assert not filtered.exists()
        else:
            called = matrix(filtered)
            assert set(called).issubset(actual)
            if label == "zero_called":
                assert called == {} and (filtered / "barcodes.tsv").read_text() == ""
        if label in ("default", "threads4"):
            assert "[SPLIT effective] soloCellFilter = EmptyDrops_CR" in (target / "Log.out").read_text()
        assert not (target / "Aligned.out.cb_ub.bin").exists()
        results[label] = {"counts": {c + "/" + g: n for (c, g), n in actual.items()}, "metrics": metrics}
    assert results["default"] == results["threads4"]

    # Mode-off reference: the same accepted records, flattened to canonical CB+UMI.
    flat_rna, flat_cb = fastqs("flat", [r for r in records if r[2]], flat=True)
    target = out / "ordinary"
    target.mkdir()
    run("ordinary", base + ["--readFilesIn", flat_rna, flat_cb, "--soloType", "CB_UMI_Simple",
                            "--soloCBwhitelist", "None", "--soloCBlen", "24", "--soloUMIstart", "25", "--soloUMIlen", "10",
                            "--soloFeatures", "GeneFull", "--soloCellFilter", "None", "--runThreadN", "1",
                            "--outFileNamePrefix", str(target) + "/"])
    assert matrix(target / "Solo.out/GeneFull/raw") == expected
    assert not (target / "Solo.out/Splitseq.metrics.tsv").exists()

    empty_rna, empty_cb = fastqs("invalid", [invalid[1]])
    target = out / "empty"
    target.mkdir()
    run("empty", base + ["--readFilesIn", empty_rna, empty_cb] + split +
        ["--runThreadN", "1", "--outFileNamePrefix", str(target) + "/"])
    assert matrix(target / "Solo.out/GeneFull/raw") == {}
    assert matrix(target / "Solo.out/GeneFull/filtered") == {}

    for label, extra, error in (
        ("bad_cb_mode", ["--soloCBmatchWLtype", "1MM"], "conflicts with split-seq mode"),
        ("both_modes", ["--soloPipseqMode", "Illumina"], "mutually exclusive"),
        ("bad_layout", ["--soloSplitseqUMIposition", "10_10"], "must not overlap"),
        ("bam_unsupported", ["--outSAMtype", "BAM", "Unsorted"], "requires --outSAMtype None"),
        ("orphan", ["--soloSplitseqMode", "None"], "require --soloSplitseqMode Fixed"),
    ):
        # STAR rejects repeated parameter assignments; replace the existing argument when overriding.
        argv = list(common)
        name = extra[0]
        if name in argv:
            start = argv.index(name)
            end = start + 1
            while end < len(argv) and not str(argv[end]).startswith("--"):
                end += 1
            del argv[start:end]
        run(label, argv + extra + ["--outFileNamePrefix", str(out / (label + "_"))], error=error)
    results["status"] = "PASS"
    results["star_sha256"] = hashlib.sha256(star.read_bytes()).hexdigest()
    results["input_sha256"] = {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [rna, cb, *lists, out / "genome.fa", out / "genes.gtf"]}
    (out / "RESULTS.json").write_text(json.dumps(results, indent=2) + "\n")
    print("PASS split-seq FASTQ decoding, sparse cell identity, UMI collapse, EmptyDrops and mode-off regression")


if __name__ == "__main__":
    main()

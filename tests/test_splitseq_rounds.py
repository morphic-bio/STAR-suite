#!/usr/bin/env python3
"""--soloSplitseqMode Fixed with two and with four barcode rounds, FASTQ to matrix.

    python3 tests/test_splitseq_rounds.py --star core/legacy/source/STAR --out /tmp/out

The mode takes two to four rounds. tests/test_splitseq_smoke.py runs three.
This test runs the two ends of the range on a small generated genome: a
two-round layout and a four-round layout whose segments add up to the largest
cell barcode the mode allows (31 bases), each with exact reads, a correction
in every round, a replacement alias, reads that must be rejected, one and
four threads, and the default cell caller. A four-round layout of 32 bases
must be refused. No external data is read.
"""
import argparse
import csv
from pathlib import Path
import random
import subprocess


def matrix(directory):
    barcodes = (directory / "barcodes.tsv").read_text().splitlines()
    genes = [s.split("\t")[0] for s in (directory / "features.tsv").read_text().splitlines()]
    lines = [s for s in (directory / "matrix.mtx").read_text().splitlines() if not s.startswith("%")]
    dims = [int(s) for s in lines[0].split()]
    assert dims == [len(genes), len(barcodes), len(lines) - 1], dims
    return {(barcodes[int(c) - 1], genes[int(g) - 1]): int(n) for g, c, n in (s.split() for s in lines[1:])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    star, out = args.star.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)

    def run(label, argv, error=None):
        log = out / (label + ".stdout")
        with log.open("w") as stream:
            proc = subprocess.run([str(x) for x in argv], stdout=stream, stderr=subprocess.STDOUT)
        text = log.read_text()
        if error:
            assert proc.returncode and error in text, (label, proc.returncode, text[-3000:])
        else:
            assert proc.returncode == 0, (label, text[-3000:])

    rng = random.Random(9143)
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
    plus = genome[700:750]
    minus = genome[5000:5050].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    u1, u2, u3 = "ACGTACGTAC", "GCTAGCTAGC", "TTGACCATGA"
    base = [star, "--genomeDir", index, "--outSAMtype", "None", "--alignIntronMax", "1000",
            "--soloStrand", "Forward", "--soloSplitseqMode", "Fixed", "--soloSplitseqUMIposition", "0_10"]

    # Segment values. "alias" is listed in round 1 with a replacement: it counts as the first value.
    a1, alias, b1 = "ACGTACGT", "TGCATGCA", "GCTAGCTA"
    a2, b2 = "GACTGACT", "CAGTCAGT"
    a3, a4 = "CTAGCTAG", "TCAGGTC"            # round 4 has seven bases: 8 + 8 + 8 + 7 = 31
    lists = {name: out / (name + ".tsv") for name in ("round1", "round2", "round3", "round4")}
    lists["round1"].write_text("%s\n%s %s\n%s\n" % (a1, alias, a1, b1))
    lists["round2"].write_text("%s\n%s\n" % (a2, b2))
    lists["round3"].write_text(a3 + "\n")
    lists["round4"].write_text(a4 + "\n")

    def sub(word, position=0):
        """word with one substitution"""
        return word[:position] + ("A" if word[position] != "A" else "C") + word[position + 1:]

    def fastqs(prefix, rows):
        rna, cb = out / (prefix + "_RNA.fastq"), out / (prefix + "_CB.fastq")
        with rna.open("w") as x, cb.open("w") as y:
            for i, (sequence, barcode) in enumerate(rows):
                for stream, s in ((x, sequence), (y, barcode)):
                    stream.write("@r%d\n%s\n+\n%s\n" % (i, s, "I" * len(s)))
        return rna, cb

    def check(label, layout, used, rows, expected, corrected, rejected):
        """Run a layout with one and four threads and with the default caller."""
        rna, cb = fastqs(label, rows * 400)   # several input chunks
        common = base + ["--readFilesIn", rna, cb, "--soloSplitseqCBposition", *layout,
                         "--soloSplitseqWhitelist", *[lists[name] for name in used]]
        results = []
        # two runs without a caller, then one with the mode's default caller (no --soloCellFilter)
        for threads, default_caller in ((1, False), (4, False), (2, True)):
            target = out / ("%s_t%d%s" % (label, threads, "_caller" if default_caller else ""))
            target.mkdir()
            run(target.name, common + ([] if default_caller else ["--soloCellFilter", "None"])
                + ["--runThreadN", str(threads), "--outFileNamePrefix", str(target) + "/"])
            actual = matrix(target / "Solo.out/GeneFull/raw")
            assert actual == expected, (target.name, actual, expected)
            metrics = dict(csv.reader((target / "Solo.out/Splitseq.metrics.tsv").open(), delimiter="\t"))
            assert int(metrics["reads_total"]) == len(rows) * 400, metrics
            exclusive = ("barcode_exact", "barcode_corrected", "read_too_short", "barcode_non_acgt",
                         "barcode_no_match", "barcode_ambiguous", "umi_non_acgt", "umi_homopolymer")
            assert sum(int(metrics[key]) for key in exclusive) == int(metrics["reads_total"]), metrics
            assert int(metrics["barcode_corrected"]) == corrected["reads"] * 400, metrics
            for n, count in enumerate(corrected["per_round"], 1):
                assert int(metrics["round%d_corrected" % n]) == count * 400, (n, metrics)
            for n in range(len(layout) + 1, 5):   # the table has four round rows; unused rounds stay at zero
                assert metrics["round%d_corrected" % n] == "0", (n, metrics)
            for key, count in rejected.items():
                assert int(metrics[key]) == count * 400, (key, metrics)
            filtered = target / "Solo.out/GeneFull/filtered"
            if default_caller:
                assert filtered.exists() and "[SPLIT effective] soloCellFilter = EmptyDrops_CR" in (target / "Log.out").read_text()
                assert set(matrix(filtered)).issubset(actual)
            else:
                assert not filtered.exists()
            results.append((actual, {k: v for k, v in metrics.items()}))
        assert results[0] == results[1], label + ": one and four threads differ"

    # ---------------------------------------------------------------- two rounds
    # barcode read: UMI 0-9, round 2 at 10-17, round 1 at 20-27
    def two(first=a1, second=a2, umi=u1):
        bases = list("T" * 30)
        for offset, s in ((0, umi), (10, second), (20, first)):
            bases[offset:offset + len(s)] = s
        return "".join(bases)

    cell_a, cell_b = a1 + a2, b1 + b2
    rows = [(plus, two())] * 3                          # one molecule, three reads
    rows += [(plus, two(umi=u2)),
             (plus, two(first=alias, umi=u3)),          # alias: the same cell, a third molecule
             (plus, two(first=sub(a1), second=sub(a2))),  # one substitution in each round: molecule u1 again
             (minus, two(umi=u2)),
             (plus, two(first=b1, second=b2)),
             (plus, two(second="AAAAAAAA")),            # no match in round 2
             (plus, two(first="NN" + a1[2:])),          # two non-ACGT bases
             (plus, two()[:27])]                        # too short for round 1
    check("two_rounds", ["20_8", "10_8"], ["round1", "round2"], rows,
          {(cell_a, "PLUS"): 3, (cell_a, "MINUS"): 1, (cell_b, "PLUS"): 1},
          {"reads": 1, "per_round": [1, 1]},
          {"barcode_no_match": 1, "barcode_non_acgt": 1, "read_too_short": 1})
    assert len(cell_a) == 16

    # --------------------------------------------------------------- four rounds
    # barcode read: UMI 0-9, round 4 at 10-16 (seven bases), round 3 at 20-27, round 2 at 30-37, round 1 at 40-47
    def four(first=a1, second=a2, third=a3, fourth=a4, umi=u1):
        bases = list("T" * 50)
        for offset, s in ((0, umi), (10, fourth), (20, third), (30, second), (40, first)):
            bases[offset:offset + len(s)] = s
        return "".join(bases)

    cell_a, cell_b = a1 + a2 + a3 + a4, b1 + b2 + a3 + a4
    rows = [(plus, four())] * 3
    rows += [(plus, four(umi=u2)),
             (plus, four(first=alias, umi=u3)),
             (plus, four(first=sub(a1), second=sub(a2), third=sub(a3), fourth=sub(a4))),  # a correction in every round
             (plus, four(fourth=sub(a4, 6), umi=u2)),   # last base of the last round: molecule u2 again
             (minus, four(umi=u2)),
             (plus, four(first=b1, second=b2)),
             (plus, four(third="AAAAAAAA")),            # no match in round 3
             (plus, four(fourth="NN" + a4[2:])),        # two non-ACGT bases in round 4
             (plus, four()[:47])]                       # too short for round 1
    layout = ["40_8", "30_8", "20_8", "10_7"]
    check("four_rounds", layout, ["round1", "round2", "round3", "round4"], rows,
          {(cell_a, "PLUS"): 3, (cell_a, "MINUS"): 1, (cell_b, "PLUS"): 1},
          {"reads": 2, "per_round": [1, 1, 1, 2]},
          {"barcode_no_match": 1, "barcode_non_acgt": 1, "read_too_short": 1})
    assert len(cell_a) == 31

    # four rounds of eight bases are 32 bases of cell barcode: refused
    lists8 = out / "round4_eight.tsv"
    lists8.write_text(a4 + "A\n")
    rna, cb = fastqs("too_long", [(plus, four())])
    run("too_long", base + ["--readFilesIn", rna, cb, "--soloSplitseqCBposition", "40_8", "30_8", "20_8", "10_8",
                            "--soloSplitseqWhitelist", lists["round1"], lists["round2"], lists["round3"], lists8,
                            "--runThreadN", "1", "--outFileNamePrefix", str(out / "too_long_")], error="31")
    # five rounds are refused as well
    run("five_rounds", base + ["--readFilesIn", rna, cb, "--soloSplitseqCBposition", "44_4", "40_4", "30_8", "20_8", "10_7",
                               "--soloSplitseqWhitelist", lists["round1"], lists["round1"], lists["round2"], lists["round3"], lists["round4"],
                               "--runThreadN", "1", "--outFileNamePrefix", str(out / "five_rounds_")], error="2-4")
    print("PASS split-seq with two and four rounds: counts, corrections per round, aliases, rejections, threads, caller, limits")


if __name__ == "__main__":
    main()

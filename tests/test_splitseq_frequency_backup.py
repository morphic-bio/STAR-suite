#!/usr/bin/env python3
"""Check split-seq frequency tie resolution across the millionth FASTQ pair and two lanes."""
import argparse
import csv
from pathlib import Path
import random
import subprocess

from test_splitseq_rounds import matrix


def fastq(name, seq):
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"


def write_lane(path, records, pad=0):
    mates = [path.with_suffix(ext) for ext in (".r1", ".r2")]
    with mates[0].open("w") as cdna, mates[1].open("w") as barcode:
        for name, sequence, tag in records:
            cdna.write(fastq(name, sequence))
            barcode.write(fastq(name, tag))
        for fh, seq in [(cdna, "N"), (barcode, "AC")]:
            record = fastq("pad", seq)
            blocks, tail = divmod(pad, 10000)
            for _ in range(blocks):
                fh.write(record * 10000)
            fh.write(record * tail)
    return mates


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--star", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    star, out = args.star.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    index = out / "index"
    index.mkdir()
    genome = "".join(random.Random(9143).choices("ACGT", k=6000))
    (out / "genome.fa").write_text(f">chrTest\n{genome}\n")
    (out / "genes.gtf").write_text('chrTest\tfixture\texon\t501\t3500\t.\t+\t.\tgene_id "PLUS"; transcript_id "PLUS.1";\n')
    subprocess.run([star, "--runMode", "genomeGenerate", "--genomeDir", index,
                    "--genomeFastaFiles", out / "genome.fa", "--sjdbGTFfile", out / "genes.gtf",
                    "--sjdbOverhang", "49", "--genomeSAindexNbases", "3", "--genomeChrBinNbits", "12",
                    "--outFileNamePrefix", str(index) + "/"], check=True)
    seconds = ["CGA", "GAC", "TAG", "TTC"]
    wl1, wl2 = out / "round1", out / "round2"
    wl1.write_text("AAA\nAAT\n")
    wl2.write_text("\n".join(seconds) + "\n")
    plus = genome[700:750]
    reads = [(word, plus, "ACGTAAC" + word) for word in seconds]
    def train(words):
        return [(f"exact{i}", "N", "AAAA" + word) for i, word in enumerate(words)]
    reads += train(["AAACGA", "AAAGAC", "AAAGAC", "AATGAC", "AAATAG", "AATTAG", "AAATTC"])
    lane1 = write_lane(out / "lane1", reads, pad=999998 - len(reads))
    # First two entries of lane 2 are in the prefix. The third is outside it.
    later = train(["AATCGA", "AAACGA", "AATCGA", "AATGAC", "AATGAC", "AAATAG", "AATTTC"])
    later.append(("counted_exact", plus, "TGCAAAACGA"))
    lane2 = write_lane(out / "lane2", later)
    common = [star, "--genomeDir", index, "--outSAMtype", "None", "--soloStrand", "Forward",
              "--readFilesIn", *[f"{a},{b}" for a, b in zip(lane1, lane2)], "--alignIntronMax", "1000",
              "--soloSplitseqMode", "FixedEdit2", "--soloSplitseqCBposition", "4_3", "7_3",
              "--soloSplitseqWhitelist", wl1, wl2, "--soloSplitseqUMIposition", "0_4",
              "--soloUMIdedup", "Exact", "--soloCellFilter", "None", "--limitIObufferSize", "1000000", "1000000"]
    # Full/prefix counts for the two candidate wells. A different full-input
    # leader or full-input tie still admits the prefix winner.
    support = dict(zip(seconds, [("3/2", "2/1"), ("2/2", "3/1"), ("2/1", "1/1"), ("1/1", "1/0")]))
    for label, threads, backup in [("one", 1, True), ("four", 4, True), ("off", 2, False)]:
        dest = out / label
        dest.mkdir()
        extra = [] if backup else ["--soloSplitseqFrequencyBackup", "no"]
        subprocess.run(common + ["--outFileNamePrefix", str(dest) + "/", "--runThreadN", str(threads)] + extra,
                       check=True)
        assert (dest / "Log.final.out").exists()
        with (dest / "Splitseq.ties.tsv").open() as f:
            trace = list(csv.DictReader(f, delimiter="\t"))
        assert sorted(r["read_name"] for r in trace) == sorted(seconds)
        for r in trace:
            key = r["read_name"]
            a, b = support[key]
            accepted = backup and key in ("CGA", "GAC", "TTC")
            assert r["candidate_exact_counts"] == (f"AAA{key}={a},AAT{key}={b}" if backup else "-"), r
            assert r["assignment_method"] == ("prefix_frequency" if accepted else "unresolved"), r
            assert r["assigned_barcode"] == ("AAA" + key if accepted else "-"), r
            assert float(r["best_posterior"]) < 0.95
            if accepted and key == "GAC":
                assert float(r["assigned_posterior"]) < float(r["best_posterior"])
        expected = {("AAACGA", "PLUS"): 2, ("AAAGAC", "PLUS"): 1, ("AAATTC", "PLUS"): 1}
        assert matrix(dest / "Solo.out/GeneFull/raw") == (expected if backup else {("AAACGA", "PLUS"): 1})
    print("PASS split-seq barcode decisions, first-million boundary and UMI counts across lanes and threads")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Check the released PIP RNA gene policies with isolated, hand-labeled reads."""
import argparse
import csv
import random
import subprocess
from pathlib import Path
from typing import NamedTuple

from compare_flex_hash_screen_mex import read_mex


NO_GENE = "reads_no_unique_gene"
MULTIMAPPED = "reads_multimapped"
COUNTED = "reads_counted"
CELL = "ACGTACGTGACTGCCTAGCATGCATGCA"
BARCODE_READ = "ACGTACGTATGGACTGCGAGCTAGCATCGAGTGCATGCA" + "ACT"


class Case(NamedTuple):
    name: str
    reads: tuple
    gene_full: str
    exon_priority: str
    fragments: int = 1


def write_reference(out):
    rng = random.Random(1742)
    genome = list("".join(rng.choices("ACGT", k=10000)))
    genome[6200:6202] = "GT"
    genome[6398:6400] = "AG"
    genome[9100:9220] = genome[8400:8520]
    genome = "".join(genome)
    (out / "genome.fa").write_text(">chrPip\n" + genome + "\n")

    # All coordinates in this test are zero-based, half-open. Each gene has
    # one plus-strand transcript. FOCUS and HOST span every exon-threshold read.
    exons = {
        "OUTER": [(400, 600), (2200, 2400)], "INNER": [(950, 1120)],
        "HOST": [(2600, 2700), (4200, 4300)], "FOCUS": [(3300, 3600), (3850, 3950)],
        "END": [(5000, 5100), (5310, 5320)], "SPLICE": [(6000, 6200), (6400, 6600)],
        "SHARED_A": [(7500, 7800)], "SHARED_B": [(7500, 7800)],
        "COPY_A": [(8400, 8520)], "COPY_B": [(9100, 9220)],
    }
    with (out / "genes.gtf").open("w") as gtf:
        for gene, blocks in exons.items():
            for start, end in blocks:
                gtf.write(f'chrPip\tfixture\texon\t{start + 1}\t{end}\t.\t+\t.\t'
                          f'gene_id "{gene}"; transcript_id "{gene}.1"; gene_name "{gene}";\n')
    return genome


def cases(genome):
    def read(start):
        return genome[start:start + 100]

    # The three junction reads occupy one BI. Two share a 5' coordinate;
    # the third starts one base earlier, so the group has two fragments.
    junction = genome[6150:6200] + genome[6400:6450]
    shifted = genome[6149:6200] + genome[6400:6449]
    antisense = read(5000).translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return [
        Case("nested_exon", (read(1000),), NO_GENE, "INNER"),
        Case("outer_intron", (read(1500),), "OUTER", "OUTER"),
        Case("exonic_49", (read(3551),), NO_GENE, NO_GENE),
        Case("exonic_50", (read(3550),), NO_GENE, NO_GENE),
        Case("exonic_51", (read(3549),), NO_GENE, "FOCUS"),
        # Each read crosses the END gene boundary. At most ten bases are
        # exonic; 19, 20 or 21 bases overlap its gene body.
        Case("gene_body_19", (read(5301),), "END", NO_GENE),
        Case("gene_body_20", (read(5300),), "END", NO_GENE),
        Case("gene_body_21", (read(5299),), "END", NO_GENE),
        Case("antisense", (antisense,), NO_GENE, NO_GENE),
        Case("splice", (junction, junction, shifted), "SPLICE", "SPLICE", fragments=2),
        Case("shared_exon", (read(7550),), NO_GENE, NO_GENE),
        Case("two_genomic_placements", (read(8410),), MULTIMAPPED, MULTIMAPPED),
    ]


def run_star(command, target):
    with (target / "stdout.log").open("w") as log:
        result = subprocess.run(command + ["--outFileNamePrefix", str(target) + "/"],
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"STAR exited {result.returncode} in {target}: "
                           + (target / "stdout.log").read_text()[-2000:])


def check_counts(solo, case, outcome):
    with (solo / "Pipseq.metrics.tsv").open() as stream:
        metrics = dict(csv.reader(stream, delimiter="\t"))
    with (solo / "Pipseq.groups.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    _, barcodes, matrix = read_mex(solo / "raw")

    accepted = outcome not in (NO_GENE, MULTIMAPPED)
    reason = COUNTED if accepted else outcome
    n = len(case.reads)
    assert int(metrics["reads_total"]) == int(metrics["barcode_exact"]) == n, solo
    for key in (COUNTED, NO_GENE, MULTIMAPPED, "reads_rejected_barcode", "reads_rejected_bi", "reads_unmapped"):
        assert int(metrics[key]) == (n if key == reason else 0), (solo, key, metrics[key], reason)

    expected_rows = []
    expected_matrix = {}
    if accepted:
        expected_rows.append({"barcode": CELL, "gene_id": outcome, "reads": str(n),
                              "occupied_bins": "1", "unique_bi_imi": str(case.fragments),
                              "estimated_transcripts": "1"})
        expected_matrix[(outcome, CELL)] = 1
    assert rows == expected_rows, (solo, rows, expected_rows)
    assert barcodes == ([CELL] if accepted else []), (solo, barcodes)
    assert matrix == expected_matrix, (solo, matrix, expected_matrix)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    star = args.star.resolve()
    genome = write_reference(out)
    index = out / "index"
    index.mkdir()
    run_star([star, "--runMode", "genomeGenerate", "--genomeDir", index, "--genomeFastaFiles", out / "genome.fa",
              "--sjdbGTFfile", out / "genes.gtf", "--sjdbOverhang", "99", "--genomeSAindexNbases", "3",
              "--genomeChrBinNbits", "12", "--runThreadN", "1"], index)

    tiers = []
    for number, sequence in enumerate(("ACGTACGT", "GACTGC", "CTAGCA", "TGCATGCA"), 1):
        path = out / f"tier{number}.tsv"
        path.write_text(sequence + "\n")
        tiers.append(path)
    base = [star, "--genomeDir", index, "--soloPipseqMode", "Illumina", "--soloPipseqWhitelist", *tiers,
            "--soloStrand", "Forward", "--soloPipseqTrim", "None", "--soloCellFilter", "None",
            "--outSAMtype", "None", "--runThreadN", "1", "--alignEndsType", "EndToEnd", "--outFilterMismatchNmax", "0"]
    probes = cases(genome)
    for case in probes:
        source = out / case.name
        source.mkdir()
        cdna, barcode = source / "cDNA.fastq", source / "barcode.fastq"
        for path, reads in ((cdna, case.reads), (barcode, [BARCODE_READ] * len(case.reads))):
            with path.open("w") as fastq:
                for i, sequence in enumerate(reads):
                    fastq.write(f"@{case.name}_{i}\n{sequence}\n+\n{'I' * len(sequence)}\n")
        for feature, outcome in (("GeneFull", case.gene_full), ("GeneFull_Ex50pAS", case.exon_priority)):
            target = source / feature
            target.mkdir()
            # Omit --soloFeatures for the bare CLI. This also pins its default.
            flags = [] if feature == "GeneFull" else ["--soloFeatures", feature]
            run_star(base + ["--readFilesIn", cdna, barcode] + flags, target)
            check_counts(target / "Solo.out" / feature, case, outcome)
    print(f"PASS: {len(probes)} PIP gene-assignment cases under both feature policies")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Analytical GEX-only counts for every supported UMI deduplication method."""
import argparse
import gzip
import os
from pathlib import Path
import random
import subprocess
import tempfile

from compare_flex_hash_screen_mex import read_mex


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", default=os.environ.get("STAR_BIN", "core/legacy/source/STAR"))
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()
    star = str(Path(args.star).resolve())
    out = args.outdir or Path(tempfile.mkdtemp(prefix="star-scrna-gex-counts-"))
    if args.outdir:
        out.mkdir(parents=True, exist_ok=False)
    index = out / "index"
    index.mkdir()
    rng = random.Random(20260927)
    genome = "".join(rng.choice("ACGT") for _ in range(6000))
    (out / "ref.fa").write_text(">chrTest\n" + genome + "\n")
    with (out / "ref.gtf").open("w") as stream:
        for gene, start, end in (("A", 1001, 1200), ("A", 1601, 1800), ("B", 3501, 3800)):
            stream.write(f'chrTest\ttest\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{gene}.1"; gene_name "{gene}";\n')
    cells = ["ACGTACGTACGTACGT", "TGCATGCATGCATGCA", "AACCGGTTAACCGGTT", "TTGGCCAATTGGCCAA"]
    ambient = "ATGCCGTAGCTAACGT"
    (out / "whitelist.txt").write_text("\n".join(cells + [ambient]) + "\n")
    # A: 7 exonic reads -> 3 exact UMIs -> 2 corrected UMIs; 2 intronic
    # reads -> 1 additional GeneFull UMI. B: 3 reads -> 1 UMI.
    molecules = [(1010, "ACGTACGTACGT", 4), (1010, "TCGTACGTACGT", 1),
                 (1010, "TTAGCGATCGGA", 2), (1300, "GGCTAACTGTAC", 2),
                 (3510, "GTAGACCACTTA", 3)]
    records = []
    for cb in cells:
        for start, umi, copies in molecules:
            records.extend([(cb + umi, genome[start:start + 75] + "G" * 40)] * copies)
    records.append((ambient + "GATCGTACCTAG", genome[1010:1085] + "G" * 40))
    with gzip.open(out / "R1.gz", "wt") as r1, gzip.open(out / "R2.gz", "wt") as r2:
        for i, (barcode_read, cdna) in enumerate(records):
            for stream, seq in ((r1, barcode_read), (r2, cdna)):
                stream.write(f"@fixture_{i}\n{seq}\n+\n{'I' * len(seq)}\n")

    def execute(command, path):
        with path.open("w") as stream:
            subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)

    execute([star, "--runMode", "genomeGenerate", "--genomeDir", str(index),
             "--genomeFastaFiles", str(out / "ref.fa"), "--sjdbGTFfile", str(out / "ref.gtf"),
             "--genomeSAindexNbases", "3", "--genomeChrBinNbits", "10",
             "--sjdbOverhang", "74", "--runThreadN", "1"], out / "index.log")
    modes = {"default": (2, 3, 1), "1MM_CR": (2, 3, 1),
             "1MM_Directional": (2, 3, 1), "1MM_Directional_UMItools": (2, 3, 1),
             "1MM_All": (2, 3, 1), "Exact": (3, 4, 1), "NoDedup": (7, 9, 3),
             "MultiGeneUMI_CR": (2, 3, 1)}
    for mode, (gene_a, full_a, gene_b) in modes.items():
        run_dir = out / mode
        run_dir.mkdir()
        command = [star, "--genomeDir", str(index), "--runThreadN", "2",
                   "--readFilesIn", str(out / "R2.gz"), str(out / "R1.gz"),
                   "--readFilesCommand", "zcat", "--outSAMtype", "None",
                   "--outFileNamePrefix", str(run_dir) + "/", "--soloType", "CB_UMI_Simple",
                   "--soloCBwhitelist", str(out / "whitelist.txt"), "--soloCBlen", "16",
                   "--soloUMIlen", "12", "--soloUMIstart", "17", "--soloStrand", "Forward",
                   "--soloFeatures", "Gene", "GeneFull", "--soloCellFilter", "CellRanger2.2", "4", "0.99", "1",
                   "--clipAdapterType", "CellRanger4", "--clip3pPolyG", "yes",
                   "--outFilterScoreMinOverLread", "0.8"]
        if mode == "MultiGeneUMI_CR":
            command += ["--soloUMIdedup", "1MM_CR", "--soloUMIfiltering", mode]
        elif mode != "default":
            command += ["--soloUMIdedup", mode]
        execute(command, run_dir / "run.log")
        for feature, a_count in (("Gene", gene_a), ("GeneFull", full_a)):
            for surface in ("raw", "filtered"):
                _, barcodes, observed = read_mex(run_dir / f"Solo.out/{feature}/{surface}")
                observed = {(gene, cb.split("-")[0]): n for (gene, cb), n in observed.items()}
                expected = {(gene, cb): n for cb in cells for gene, n in (("A", a_count), ("B", gene_b))}
                if surface == "raw":
                    expected[("A", ambient)] = 1
                else:
                    assert {cb.split("-")[0] for cb in barcodes} == set(cells)
                assert observed == expected, f"{mode}/{feature}/{surface}: {observed} != {expected}"
        print(f"PASS: {mode}, exact Gene/GeneFull and filtered-cell counts", flush=True)
    (out / "PASS").write_text("All UMI methods preserve analytical GEX counts.\n")
    print(f"Results: {out}")


if __name__ == "__main__":
    main()

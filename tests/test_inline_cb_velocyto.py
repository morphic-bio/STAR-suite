#!/usr/bin/env python3
"""Synthetic no-BAM Velocyto regression with inline correction and native OCM."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile

from compare_flex_hash_screen_mex import read_mex


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star", type=Path, default=Path(os.environ.get("STAR_BIN", "core/legacy/source/STAR")))
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--index", type=Path, help="Reuse this fixture's previously generated index")
    parser.add_argument("--expect-old-failure", action="store_true",
                        help="Diagnostic only: verify the old inline-storage failure")
    args = parser.parse_args()
    star = args.star.resolve()
    out = args.outdir or Path(tempfile.mkdtemp(prefix="star-inline-velo-"))
    if args.outdir:
        out.mkdir(parents=True, exist_ok=False)
    index = args.index or out / "index"
    if not args.index:
        index.mkdir()
    rng = random.Random(20260928)
    genome = "".join(rng.choice("ACGT") for _ in range(6000))
    (out / "ref.fa").write_text(">chrTest\n" + genome + "\n")
    with (out / "ref.gtf").open("w") as stream:
        for start, end in ((1001, 1200), (1601, 1800)):
            stream.write(f'chrTest\ttest\texon\t{start}\t{end}\t.\t+\t.\tgene_id "A"; transcript_id "A.1"; gene_name "A";\n')
    cells = ["AAACCCTGTAAGCGCG", "AAACCCGCAACTAGAC", "AAACCATTCACCTGGG", "AAACCAAAGCATTGAT"]
    (out / "whitelist.txt").write_text("\n".join(cells) + "\n")
    (out / "config.csv").write_text("[samples]\nsample_id,ocm_barcode_ids\n" +
                                     "".join(f"sample{i},OB{i}\n" for i in range(1, 5)))
    # One junction-spanning and one intronic molecule per barcode, with PCR copies.
    molecules = [(genome[1170:1200] + genome[1600:1645], "ACGTACGTACGT", 3),
                 (genome[1300:1375], "TTAGCGATCGGA", 2)]
    with gzip.open(out / "R1.gz", "wt") as r1, gzip.open(out / "R2.gz", "wt") as r2:
        read_id = 0
        for cell in cells:
            for cdna, umi, copies in molecules:
                for _ in range(copies):
                    for stream, seq in ((r1, cell + umi), (r2, cdna)):
                        stream.write(f"@fixture_{read_id}\n{seq}\n+\n{'I' * len(seq)}\n")
                    read_id += 1

    records = []

    def execute(command, path):
        with path.open("w") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        records.append({"argv": command, "exit": result.returncode})
        (out / "commands.json").write_text(json.dumps(records, indent=2) + "\n")
        return result.returncode

    command = [str(star), "--runMode", "genomeGenerate", "--genomeDir", str(index),
               "--outFileNamePrefix", str(index) + "/",
               "--genomeFastaFiles", str(out / "ref.fa"), "--sjdbGTFfile", str(out / "ref.gtf"),
               "--genomeSAindexNbases", "3", "--genomeChrBinNbits", "10",
               "--sjdbOverhang", "74", "--runThreadN", "1"]
    if not args.index:
        assert execute(command, out / "index.log") == 0, out / "index.log"
    for mode in ("plain", "inline", "ocm"):
        run = out / mode
        run.mkdir()
        command = [str(star), "--genomeDir", str(index), "--runThreadN", "2",
                   "--readFilesIn", str(out / "R2.gz"), str(out / "R1.gz"),
                   "--readFilesCommand", "zcat", "--outSAMtype", "None",
                   "--outFileNamePrefix", str(run) + "/", "--soloType", "CB_UMI_Simple",
                   "--soloCBwhitelist", str(out / "whitelist.txt"), "--soloCBlen", "16",
                   "--soloUMIstart", "17", "--soloUMIlen", "12", "--soloBarcodeReadLength", "0",
                   "--soloStrand", "Forward", "--soloFeatures", "GeneFull", "Velocyto",
                   "--soloCellFilter", *(["None"] if mode == "ocm" else ["TopCells", "4"]),
                   "--soloCrMultimapRescue", "yes", "--soloUMIdedup", "1MM_CR",
                   "--soloUMIfiltering", "MultiGeneUMI_CR", "--soloCBmatchWLtype", "Exact",
                   "--soloInlineHashMode", "no", "--soloInlineCBCorrection", "no" if mode == "plain" else "yes"]
        if mode == "ocm":
            command += ["--ocmMultiEnable", "yes", "--ocmMultiConfig", str(out / "config.csv"),
                        "--ocmMultiBarcodeMode", "flex", "--soloCrGexFeature", "genefull"]
        rc = execute(command, run / "run.log")
        if mode == "inline" and args.expect_old_failure:
            assert rc == 111, (rc, run / "run.log")
            assert "gene-like source did not populate per-read CB/UMI storage" in (run / "run.log").read_text()
            (out / "EXPECTED_OLD_FAILURE").write_text("Inline correction suppresses required readInfo.\n")
            print("Confirmed old storage failure after plain Velocyto control passed")
            return
        assert rc == 0, (rc, run / "run.log")
        _, _, counts = read_mex(run / "Solo.out/GeneFull/raw")
        assert {(gene, cb[:16]): value for (gene, cb), value in counts.items()} == {
            ("A", cell): 2 for cell in cells}, counts
        velo = run / "Solo.out/Velocyto/raw"
        barcodes = (velo / "barcodes.tsv").read_text().splitlines()
        assert {cb[:16] for cb in barcodes} == set(cells), barcodes
        assert all(len(cb.split("-")[0]) == (24 if mode == "ocm" else 16) for cb in barcodes), barcodes
        for layer, expected in (("spliced", 1), ("unspliced", 1), ("ambiguous", 0)):
            features, layer_barcodes, counts = read_mex(velo, matrix_name=f"{layer}.mtx")
            assert len(features) == 1 and layer_barcodes == barcodes
            assert {(gene, cb[:16]): value for (gene, cb), value in counts.items() if value} == {
                ("A", cell): expected for cell in cells if expected}, (mode, layer, counts)
        assert not list(run.glob("*.bam")), "Test must exercise the no-BAM path"
        print(f"PASS: {mode}, exact GeneFull and spliced/unspliced/ambiguous counts", flush=True)
    (out / "PASS.json").write_text(json.dumps({"binary_sha256": hashlib.sha256(star.read_bytes()).hexdigest(),
                                               "modes": ["plain", "inline", "ocm"]}, indent=2) + "\n")


if __name__ == "__main__":
    main()

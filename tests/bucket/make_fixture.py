#!/usr/bin/env python3
"""Create a deterministic paired multi-tag FASTQ fixture for CB bucket tests."""

import argparse
from pathlib import Path


BARCODES = (
    "ACGTACGTACGTACGT", "TGCATGCATGCATGCA", "GATTACAGATTACAGA", "CTAGGCTACTAGGCTA",
    "AACCGGTTAACCGGTT", "TTGGCCAATTGGCCAA", "AGCTAGCTAGCTAGCT", "TCGATCGATCGATCGA",
    "CATGCATGCATGCATG", "GTACGTACGTACGTAC", "ATGCATGCATGCATGC", "CGTACGTACGTACGTA",
    "AAGCTTGCAAGCTTGC", "TTCGAACGTTCGAACG", "GGCATACGGGCATACG", "CCGTATGCCCGTATGC",
)

UMIS = (
    "ACGTACGTACGT", "TGCATGCATGCA", "GATTACAGATTA", "CTAGGCTACTAG",
    "AACCGGTTAACC", "TTGGCCAATTGG", "AGCTAGCTAGCT", "TCGATCGATCGA",
)

TAGS = (
    ("sample_a", "AACCGGTT"),
    ("sample_b", "TTGGCCAA"),
    ("sample_c", "AGCTAGCT"),
    ("sample_d", "TCGATCGA"),
)

PROBES = (
    ("probe_1", "ACGTTGCAACGATCGTACGTTAGCTAGGCTAA", 1),
    ("probe_2", "GCTAGCTAGGATCCGATCGTACCTAGGCTAAC", 2),
    ("probe_3", "TTACCGGATGCTAACCGTATCGGATACCGTTA", 3),
    ("probe_4", "CGATTCGGAACCTTAGCGATACCTTGGCATAC", 1),
    ("probe_5", "ATCGGCTAACGTTCGATGGCATTCGATACGGT", 2),
    ("probe_6", "GGTACCATCGATGCTAGCTTACCGATGGCATA", 0),
)


def write_fastq_record(handle, name, sequence):
    handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    (outdir / "barcodes.tsv").write_text(
        "".join(f"{index}\t{sequence}\n" for index, sequence in enumerate(BARCODES)),
        encoding="ascii",
    )
    (outdir / "tags.tsv").write_text(
        "".join(f"{index}\t{label}\t{sequence}\n"
                for index, (label, sequence) in enumerate(TAGS, 1)),
        encoding="ascii",
    )
    (outdir / "probes.tsv").write_text(
        "".join(f"{index}\t{name}\t{sequence}\t{region}\n"
                for index, (name, sequence, region) in enumerate(PROBES, 1)),
        encoding="ascii",
    )

    truth_lines = ["read_id\tcb_idx\tumi\tgene_idx\ttag_idx\tregion\n"]
    read_id = 0
    with (outdir / "R1.fastq").open("w", encoding="ascii", newline="\n") as r1, \
            (outdir / "R2.fastq").open("w", encoding="ascii", newline="\n") as r2:
        for cb_idx, barcode in enumerate(BARCODES):
            for gene_idx, (_, probe, region) in enumerate(PROBES, 1):
                tag_idx = ((cb_idx + gene_idx - 1) % len(TAGS)) + 1
                tag = TAGS[tag_idx - 1][1]
                umi_indices = [(cb_idx * 3 + gene_idx) % len(UMIS)]
                if (cb_idx + gene_idx) % 3 == 0:
                    umi_indices.append((umi_indices[0] + 3) % len(UMIS))
                for umi_slot, umi_index in enumerate(umi_indices):
                    umi = UMIS[umi_index]
                    copies = 1 + ((cb_idx + gene_idx + umi_slot) % 3)
                    for copy in range(copies):
                        name = f"bucket_{read_id:06d}/1"
                        write_fastq_record(r1, name, barcode + umi)
                        # Sample tag starts at zero-based offset 68, matching the
                        # production Flex fixture layout.
                        r2_sequence = probe + ("ACGT" * 9)[:36] + tag + "GATTACA"
                        write_fastq_record(r2, name.replace("/1", "/2"), r2_sequence)
                        truth_lines.append(
                            f"{read_id}\t{cb_idx}\t{umi}\t{gene_idx}\t{tag_idx}\t{region}\n"
                        )
                        read_id += 1

    (outdir / "truth.tsv").write_text("".join(truth_lines), encoding="ascii")
    (outdir / "README.txt").write_text(
        "R1 is CB16+UMI12. R2 is probe32+padding36+tag8+tail7. "
        "truth.tsv is generator provenance; reference_counter.py derives its result from FASTQ.\n",
        encoding="ascii",
    )


if __name__ == "__main__":
    main()

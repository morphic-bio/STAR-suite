#!/usr/bin/env python3
"""Generate synthetic FASTQ + whitelist + feature_ref for NXT/TRU auto-detection tests.

Creates:
  fixtures/whitelist_tru.txt       — 50 barcodes in TRU namespace (1-column)
  fixtures/whitelist_nxt.txt       — same 50 barcodes in NXT namespace (1-column)
  fixtures/feature_ref.csv         — 2 CRISPR guides with anchored pattern
  fixtures/feature_ref_lineage.csv — 2 lineage barcodes with different anchor
  fixtures/reads_tru/              — 1000 reads with TRU barcodes (CRISPR guides)
  fixtures/reads_nxt/              — 1000 reads with NXT barcodes (CRISPR guides)
  fixtures/reads_mixed/            — 500 TRU + 500 NXT
  fixtures/reads_tiny/             — 5 reads (drain finalization test)
  fixtures/reads_lib_a_nxt/        — 1000 reads: NXT barcodes + CRISPR guides
  fixtures/reads_lib_b_tru/        — 1000 reads: TRU barcodes + lineage barcodes
"""

import os, random, gzip

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
N_BARCODES = 50
BARCODE_LEN = 16
UMI_LEN = 12
N_READS = 1000
N_TINY = 5

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

GUIDE_ANCHOR = "TTCCAGCTTAGCTCTTAAAC"
GUIDES = [
    ("Guide_A", "GGGTGGTGCCCATCCTGGTC"),
    ("Guide_B", "TATTTCCTGGTTCGCCGGCC"),
]

LINEAGE_ANCHOR = "AAGCAGTGGTATCAACGCAG"
LINEAGE_BARCODES = [
    ("Lineage_1", "ACGTACGTACGTACGTACGT"),
    ("Lineage_2", "TGCATGCATGCATGCATGCA"),
]

random.seed(42)


def random_seq(length):
    return "".join(random.choice(BASES) for _ in range(length))


def translate_nxt(barcode):
    """Complement positions 7-8 (0-indexed)."""
    bc = list(barcode)
    bc[7] = bc[7].translate(COMPLEMENT)
    bc[8] = bc[8].translate(COMPLEMENT)
    return "".join(bc)


def qual_string(length):
    return "I" * length


def write_fastq_gz(path, records):
    """Write list of (name, seq, qual) tuples to gzipped FASTQ."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        for name, seq, qual in records:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def generate_read_pair(read_idx, barcode, anchor, feature_seq):
    name = f"read_{read_idx}"
    umi = random_seq(UMI_LEN)
    r1_seq = barcode + umi
    r1_qual = qual_string(len(r1_seq))
    r2_seq = anchor + feature_seq
    r2_qual = qual_string(len(r2_seq))
    return (name, r1_seq, r1_qual), (name, r2_seq, r2_qual)


def make_reads(barcodes_tru, barcodes_nxt, mode, n_reads,
               anchor=GUIDE_ANCHOR, features=GUIDES):
    r1_records = []
    r2_records = []
    for i in range(n_reads):
        _, feature_seq = features[i % len(features)]
        if mode == "tru":
            bc = random.choice(barcodes_tru)
        elif mode == "nxt":
            bc = random.choice(barcodes_nxt)
        elif mode == "mixed":
            bc = random.choice(barcodes_tru if i % 2 == 0 else barcodes_nxt)
        else:
            raise ValueError(f"Unknown mode: {mode}")
        r1, r2 = generate_read_pair(i, bc, anchor, feature_seq)
        r1_records.append(r1)
        r2_records.append(r2)
    return r1_records, r2_records


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    barcodes_tru = [random_seq(BARCODE_LEN) for _ in range(N_BARCODES)]
    # Ensure no self-complementary pairs at positions 7-8
    for i, bc in enumerate(barcodes_tru):
        while bc[7] == bc[7].translate(COMPLEMENT) and bc[8] == bc[8].translate(COMPLEMENT):
            barcodes_tru[i] = random_seq(BARCODE_LEN)
            bc = barcodes_tru[i]

    barcodes_nxt = [translate_nxt(bc) for bc in barcodes_tru]

    # Whitelists
    wl_tru_path = os.path.join(OUTDIR, "whitelist_tru.txt")
    wl_nxt_path = os.path.join(OUTDIR, "whitelist_nxt.txt")
    with open(wl_tru_path, "w") as f:
        for bc in barcodes_tru:
            f.write(bc + "\n")
    with open(wl_nxt_path, "w") as f:
        for bc in barcodes_nxt:
            f.write(bc + "\n")

    # Feature reference CSV
    ref_path = os.path.join(OUTDIR, "feature_ref.csv")
    with open(ref_path, "w") as f:
        f.write("id,name,read,pattern,sequence,feature_type,target_gene_id,target_gene_name\n")
        for gid, (gname, gseq) in enumerate(GUIDES, 1):
            f.write(f"{gname},{gname},R2,{GUIDE_ANCHOR}(BC),{gseq},"
                    f"CRISPR Guide Capture,target_{gid},{gname}\n")

    # Lineage barcode feature reference CSV
    lineage_ref_path = os.path.join(OUTDIR, "feature_ref_lineage.csv")
    with open(lineage_ref_path, "w") as f:
        f.write("id,name,read,pattern,sequence,feature_type,target_gene_id,target_gene_name\n")
        for lid, (lname, lseq) in enumerate(LINEAGE_BARCODES, 1):
            f.write(f"{lname},{lname},R2,{LINEAGE_ANCHOR}(BC),{lseq},"
                    f"Custom,lineage_{lid},{lname}\n")

    # Generate single-library read sets
    for mode, n in [("tru", N_READS), ("nxt", N_READS), ("mixed", N_READS), ("tru", N_TINY)]:
        dirname = f"reads_{mode}" if n != N_TINY else "reads_tiny"
        r1_records, r2_records = make_reads(barcodes_tru, barcodes_nxt, mode, n)
        read_dir = os.path.join(OUTDIR, dirname)
        write_fastq_gz(os.path.join(read_dir, "sample1_S1_L001_R1_001.fastq.gz"), r1_records)
        write_fastq_gz(os.path.join(read_dir, "sample1_S1_L001_R2_001.fastq.gz"), r2_records)

    # Multi-library read sets:
    # Library A = CRISPR guides with NXT barcodes (simulates 3HT CRISPR)
    r1a, r2a = make_reads(barcodes_tru, barcodes_nxt, "nxt", N_READS,
                          anchor=GUIDE_ANCHOR, features=GUIDES)
    lib_a_dir = os.path.join(OUTDIR, "reads_lib_a_nxt")
    write_fastq_gz(os.path.join(lib_a_dir, "libA_S1_L001_R1_001.fastq.gz"), r1a)
    write_fastq_gz(os.path.join(lib_a_dir, "libA_S1_L001_R2_001.fastq.gz"), r2a)

    # Library B = lineage barcodes with TRU barcodes
    r1b, r2b = make_reads(barcodes_tru, barcodes_nxt, "tru", N_READS,
                          anchor=LINEAGE_ANCHOR, features=LINEAGE_BARCODES)
    lib_b_dir = os.path.join(OUTDIR, "reads_lib_b_tru")
    write_fastq_gz(os.path.join(lib_b_dir, "libB_S1_L001_R1_001.fastq.gz"), r1b)
    write_fastq_gz(os.path.join(lib_b_dir, "libB_S1_L001_R2_001.fastq.gz"), r2b)

    print(f"Generated synthetic fixtures in {OUTDIR}")
    print(f"  Whitelists: {wl_tru_path}, {wl_nxt_path}")
    print(f"  Feature refs: {ref_path}, {lineage_ref_path}")
    print(f"  Single-library: reads_tru/ reads_nxt/ reads_mixed/ reads_tiny/")
    print(f"  Multi-library: reads_lib_a_nxt/ (CRISPR+NXT) reads_lib_b_tru/ (lineage+TRU)")
    print(f"  Barcodes: {N_BARCODES}, Reads: {N_READS} (+ {N_TINY} tiny)")


if __name__ == "__main__":
    main()

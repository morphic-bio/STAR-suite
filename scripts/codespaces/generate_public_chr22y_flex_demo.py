#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path

SAMPLE_TAGS = [("PB01", "AAAACCCC"), ("PB02", "TTTTGGGG")]
READ_LENGTH = 100
PROBE_LENGTH = 50
TAG_OFFSET = 68


def iter_fastq(path: Path):
    with gzip.open(path, "rt") as handle:
        while True:
            name = handle.readline().rstrip("\n")
            if not name:
                break
            seq = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            qual = handle.readline().rstrip("\n")
            yield name, seq, plus, qual


def load_fasta(path: Path):
    seqs = {}
    current = None
    chunks = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(table)[::-1]


def parse_gene_id(attr: str) -> str:
    for field in attr.split(";"):
        field = field.strip()
        if field.startswith("gene_id"):
            parts = field.split('"')
            if len(parts) >= 2:
                return parts[1]
    return ""


def choose_probe_sequences(fasta_path: Path, gtf_path: Path):
    genome = load_fasta(fasta_path)
    wanted = ["chr22", "chrY", "22", "Y"]
    chosen = []
    seen = set()
    open_handle = gzip.open(gtf_path, "rt") if str(gtf_path).endswith(".gz") else gtf_path.open()
    with open_handle as handle:
        for chrom in wanted:
            handle.seek(0)
            for raw in handle:
                if not raw or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[0] != chrom or fields[2] != "exon":
                    continue
                gene_id = parse_gene_id(fields[8])
                if not gene_id or gene_id in seen:
                    continue
                start0 = int(fields[3]) - 1
                end0 = int(fields[4])
                seq = genome.get(chrom, "")[start0:end0]
                if len(seq) < PROBE_LENGTH:
                    continue
                seq = seq[:PROBE_LENGTH]
                if fields[6] == "-":
                    seq = reverse_complement(seq)
                chosen.append((gene_id, seq))
                seen.add(gene_id)
                break
    if len(chosen) < 2:
        raise SystemExit(f"Could not find one chr22 and one chrY exon of length >= {PROBE_LENGTH} in {gtf_path}")
    return chosen[:2]


def enumerate_pairs(fastq_dir: Path):
    r1s = sorted(fastq_dir.glob("*_R1_001.fastq.gz"))
    if not r1s:
        raise SystemExit(f"No R1 FASTQs found under {fastq_dir}")
    pairs = []
    for r1 in r1s:
        r2 = Path(str(r1).replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
        if not r2.is_file():
            raise SystemExit(f"Missing R2 mate for {r1}")
        pairs.append((r1, r2))
    return pairs


def write_probe_set(path: Path, probes):
    with path.open("w", newline="") as handle:
        handle.write("#probe_set_file_format=2.0\n")
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["gene_id", "probe_seq", "probe_id", "included", "region"])
        for idx, (gene_id, probe_seq) in enumerate(probes, start=1):
            writer.writerow([gene_id, probe_seq, f"{gene_id}|demo|{idx}", "TRUE", "spliced"])


def write_config(path: Path, ref_root: Path, probe_set: Path, fastq_dir: Path, catalog: Path):
    path.write_text(
        "[gene-expression]\n"
        f"reference,{ref_root}\n"
        f"probe-set,{probe_set}\n"
        "create-bam,false\n\n"
        "[star]\n"
        f"sample-probe-catalog,{catalog}\n"
        f"sample-probe-offset,{TAG_OFFSET}\n\n"
        "[libraries]\n"
        "fastq_id,fastqs,feature_types\n"
        f"public_chr22y_flex,{fastq_dir},Gene Expression\n\n"
        "[samples]\n"
        "sample_id,probe_barcode_ids,description\n"
        "PB01,PB01,Demo sample 1\n"
        "PB02,PB02,Demo sample 2\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a chr22+chrY Flex microfixture from the public derived GEX surface")
    parser.add_argument("--gex-fastq-dir", required=True)
    parser.add_argument("--ref-root", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--read-limit", type=int, default=4000)
    args = parser.parse_args()

    gex_fastq_dir = Path(args.gex_fastq_dir).resolve()
    ref_root = Path(args.ref_root).resolve()
    outdir = Path(args.outdir).resolve()
    fastq_out = outdir / "gex" / "public_chr22y_flex"
    outdir.mkdir(parents=True, exist_ok=True)
    fastq_out.mkdir(parents=True, exist_ok=True)

    fasta = ref_root / "fasta" / "genome.fa"
    gtf = ref_root / "genes" / "genes.gtf"
    probes = choose_probe_sequences(fasta, gtf)

    probe_set = outdir / "probe_set.csv"
    whitelist = outdir / "whitelist_TRU.txt"
    catalog = outdir / "sample_probe_catalog.tsv"
    config = outdir / "config.csv"
    manifest = outdir / "MANIFEST.txt"

    write_probe_set(probe_set, probes)
    catalog.write_text("PB01\tAAAACCCC\tPB01\nPB02\tTTTTGGGG\tPB02\n", encoding="utf-8")

    source_r1 = enumerate_pairs(gex_fastq_dir)[0][0]
    out_r1 = fastq_out / "public_chr22y_flex_S1_L001_R1_001.fastq.gz"
    out_r2 = fastq_out / "public_chr22y_flex_S1_L001_R2_001.fastq.gz"

    whitelist_values = set()
    read_count = 0
    with gzip.open(out_r1, "wt") as handle_r1, gzip.open(out_r2, "wt") as handle_r2:
        for idx, (name, seq, plus, qual) in enumerate(iter_fastq(source_r1)):
            cb = seq[:16]
            whitelist_values.add(cb)
            sample_id, sample_tag = SAMPLE_TAGS[idx % len(SAMPLE_TAGS)]
            gene_id, probe_seq = probes[idx % len(probes)]
            filler_before = "A" * max(0, TAG_OFFSET - len(probe_seq))
            tail_len = READ_LENGTH - len(probe_seq) - len(filler_before) - len(sample_tag)
            tail = "T" * max(0, tail_len)
            r2_seq = (probe_seq + filler_before + sample_tag + tail)[:READ_LENGTH]
            handle_r1.write(f"{name}\n{seq}\n{plus}\n{qual}\n")
            handle_r2.write(f"{name.replace('/1', '/2')}\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")
            read_count += 1
            if read_count >= args.read_limit:
                break

    whitelist.write_text("\n".join(sorted(whitelist_values)) + "\n", encoding="utf-8")
    write_config(config, ref_root, probe_set, fastq_out, catalog)
    manifest.write_text(
        "Public chr22+chrY Flex microfixture\n"
        f"Generated reads: {read_count}\n"
        f"FASTQ dir: {fastq_out}\n"
        f"Probe set: {probe_set}\n"
        f"Sample probe catalog: {catalog}\n"
        f"Whitelist: {whitelist}\n"
        f"CR config: {config}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

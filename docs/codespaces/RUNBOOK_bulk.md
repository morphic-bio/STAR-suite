# Bulk RNA-seq Codespaces Runbook

Step-by-step walkthrough for the bulk RNA-seq demo in GitHub Codespaces.
Copy-paste each block in order. Every command should succeed before moving on.

## Prerequisites

You need a running Codespace with the STAR-suite repository.
The Dockerfile installs all required tools (gcc, samtools, pigz, fastq-dump).

## Step 0: Build STAR (one-time, ~2 min)

```bash
make -C core/legacy/source -j$(nproc) STAR
```

Verify:

```bash
core/legacy/source/STAR --version
```

Expected: `1.1.0`

For upstream STAR provenance, run:

```bash
core/legacy/source/STAR --upstream-version
```

---

## Step 1: Download reference and build index (~1 min)

The demo uses a small human reference (chr22 + chrY only) that fits in Codespaces.

```bash
bash scripts/codespaces/fetch_public_chr22y_reference.sh
```

```bash
bash scripts/codespaces/build_public_chr22y_index.sh --threads 4
```

Verify:

```bash
ls .codespaces-demo/indices/public_human_chr22y_star/Genome
```

Expected: file exists, ~106 MB.

---

## Step 2: Download public FASTQ fixture (~10 sec)

Downloads 25,000 paired-end reads from SRA (SRR4422207, human esophagus).

```bash
bash scripts/codespaces/run_bulk_public_demo.sh \
  --dry-run \
  --genome-dir .codespaces-demo/indices/public_human_chr22y_star
```

This downloads the FASTQs and writes `RUN_COMMAND.sh` without running STAR.

Verify:

```bash
ls .codespaces-demo/runs/bulk_public_demo/RUN_COMMAND.sh
cat .codespaces-demo/runs/bulk_public_demo/RUN_MANIFEST.txt
```

---

## Step 3: Run STAR (~2 min)

```bash
bash scripts/codespaces/run_bulk_public_demo.sh \
  --run \
  --genome-dir .codespaces-demo/indices/public_human_chr22y_star
```

Or equivalently, execute the saved command:

```bash
bash .codespaces-demo/runs/bulk_public_demo/RUN_COMMAND.sh
```

---

## Step 4: Verify outputs

```bash
bash scripts/codespaces/verify_bulk_public_demo_outputs.sh \
  --outdir .codespaces-demo/runs/bulk_public_demo
```

Expected output:

```
Bulk public demo outputs look good.

Verified files:
  .../RUN_COMMAND.sh
  .../Aligned.sortedByCoord.out.bam
  .../ReadsPerGene.out.tab

Quick summary:
  BAM reads: ~3000
  ReadsPerGene rows: ~2050
```

---

## Step 5: Inspect results

View the gene count table:

```bash
head -n 12 .codespaces-demo/runs/bulk_public_demo/ReadsPerGene.out.tab
```

The first four rows are summary statistics:
- `N_unmapped` — reads that did not map
- `N_multimapping` — reads mapping to multiple loci
- `N_noFeature` — mapped reads not overlapping any gene
- `N_ambiguous` — reads overlapping multiple genes

The remaining rows show per-gene counts (columns: unstranded, sense, antisense).

Check the top expressed genes:

```bash
awk '$2 > 0 && $1 !~ /^N_/' .codespaces-demo/runs/bulk_public_demo/ReadsPerGene.out.tab \
  | sort -k2 -nr | head -10
```

**Note:** Most reads (23K of 25K) are unmapped because this is a full-genome
dataset aligned to a chr22+chrY-only reference. This is expected for the demo.
With a full genome index, the mapping rate would be ~85-95%.

Inspect the BAM:

```bash
samtools flagstat .codespaces-demo/runs/bulk_public_demo/Aligned.sortedByCoord.out.bam
```

---

## What STAR did

The demo ran this pipeline in a single STAR invocation:

1. **Read decompression** — gzipped FASTQs via `zcat`
2. **Adapter trimming** — Illumina adapters via `--trimCutadapt Yes`
3. **Alignment** — splice-aware alignment to the genome
4. **BAM output** — coordinate-sorted BAM (`Aligned.sortedByCoord.out.bam`)
5. **Gene counting** — per-gene read counts (`ReadsPerGene.out.tab`)

No external tools (Trim Galore, featureCounts, etc.) were needed.

---

## Using your own data

To switch to real data, change two things:

1. **Build a full genome index** using your FASTA + GTF:
   ```bash
   core/legacy/source/STAR --runMode genomeGenerate \
     --genomeDir /path/to/my_index \
     --genomeFastaFiles /path/to/genome.fa \
     --sjdbGTFfile /path/to/genes.gtf \
     --runThreadN $(nproc)
   ```

2. **Point to your FASTQs and index**:
   ```bash
   core/legacy/source/STAR \
     --runMode alignReads \
     --runThreadN $(nproc) \
     --genomeDir /path/to/my_index \
     --readFilesIn sample_R2.fastq.gz sample_R1.fastq.gz \
     --readFilesCommand zcat \
     --trimCutadapt Yes \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outFileNamePrefix /path/to/output/
   ```

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `STAR binary not found` | Run `make -C core/legacy/source -j$(nproc) STAR` |
| `fastq-dump: command not found` | Install SRA toolkit: `apt-get install sra-toolkit` |
| `pigz: command not found` | Install pigz: `apt-get install pigz` |
| Most reads unmapped | Expected with chr22+chrY index. Use full genome for real data. |
| `EXDEV: cross-device link` on Log.out | Harmless warning; STAR still runs correctly. |
| Step 3 hangs | Check `--threads` (default 4). Codespaces may have limited cores. |

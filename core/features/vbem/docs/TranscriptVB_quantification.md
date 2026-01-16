# TranscriptVB: Salmon-Style Transcript Quantification for Bulk RNA-seq

## Overview

TranscriptVB provides integrated Salmon-style transcript quantification directly within STAR alignment. This feature performs **Variational Bayes (VB) or Expectation-Maximization (EM)** quantification using equivalence classes built inline during alignment.

> **Note**: This feature is designed for **Bulk RNA-seq** workflows, not single-cell/FLEX applications. For single-cell quantification, use STARsolo.

## Key Features

- **Inline EC Building**: Equivalence classes are constructed during alignment, eliminating the need for a separate quantification step
- **Salmon-Compatible Output**: Produces `quant.sf` files in the same format as Salmon
- **VB/EM Algorithms**: Uses the same Variational Bayes algorithm as Salmon (with EM as an option)
- **GC Bias Collection**: Optional fragment GC observation for bias-aware quantification
- **Multi-threaded**: EC tables are built per-thread and merged efficiently

## Usage

### Basic Usage

```bash
STAR --runMode alignReads \
     --genomeDir /path/to/genome \
     --readFilesIn reads_1.fastq.gz reads_2.fastq.gz \
     --readFilesCommand zcat \
     --quantMode TranscriptVB \
     --outFileNamePrefix output_
```

### With GC Bias Collection

```bash
STAR --runMode alignReads \
     --genomeDir /path/to/genome \
     --readFilesIn reads_1.fastq.gz reads_2.fastq.gz \
     --readFilesCommand zcat \
     --quantMode TranscriptVB \
     --quantVBgcBias 1 \
     --outFileNamePrefix output_
```

### Using EM Instead of VB

```bash
STAR --runMode alignReads \
     --genomeDir /path/to/genome \
     --readFilesIn reads_1.fastq.gz reads_2.fastq.gz \
     --readFilesCommand zcat \
     --quantMode TranscriptVB \
     --quantVBem 1 \
     --outFileNamePrefix output_
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--quantMode TranscriptVB` | - | Enable TranscriptVB quantification |
| `--quantVBgcBias` | 0 | Enable GC bias collection (0/1) |
| `--quantVBem` | 0 | Use EM instead of VB (0/1) |
| `--quantVBLibType` | A | Library format for compatibility (A=auto-detect, IU/ISF/ISR/U fixed) |
| `--quantVBAutoDetectWindow` | 1000 | Number of reads used for auto-detect (when `--quantVBLibType A`) |
| `--quantVBprior` | 0.01 | Dirichlet prior for VB algorithm |
| `--quantVBgenes` | 1 | Output gene-level quantification (0/1) |

### Library Format Auto-Detect Warnings

When `--quantVBLibType A` is used, STAR logs format votes and may emit warnings in these cases:
- Winner is weak (<85% of votes).
- Winner is outward/same-strand (OSF/OSR/MSF/MSR) or labeled UNKNOWN.
- UNKNOWN votes exceed 15% of total.

These warnings often indicate mate order issues (R1/R2 swapped) or unexpected library prep. If you see them, verify mate order (for these pipelines STAR expects R2 then R1) or set `--quantVBLibType` explicitly.

## Output Files

| File | Description |
|------|-------------|
| `*quant.sf` | Salmon-compatible transcript-level quantification results |
| `*quant.genes.sf` | Gene-level quantification results (aggregated from transcripts) |
| `*Aligned.out.sam/bam` | Standard STAR alignment output |
| `*Log.out` | Log file with quantification statistics |

### quant.sf Format

```
Name    Length    EffectiveLength    TPM    NumReads
ENST00000456328    1657    1458.000    0.000000    0.000
ENST00000450305    632    433.000    12.345678    1.234
...
```

### quant.genes.sf Format

When `--quantMode TranscriptVB` is used, STAR also produces `quant.genes.sf` containing gene-level expression estimates. This file aggregates transcript-level results to genes.

| Column | Description |
|--------|-------------|
| Name | Gene ID (from GTF gene_id attribute) |
| Length | TPM-weighted average of transcript lengths |
| EffectiveLength | TPM-weighted average of transcript effective lengths |
| TPM | Sum of transcript TPMs for this gene |
| NumReads | Sum of transcript counts for this gene |

**Example:**
```
Name	Length	EffectiveLength	TPM	NumReads
ENSG00000268674	510	295.805	0	0
ENSG00000271254	3081.77	2865.89	9.91087	491.292
ENSG00000275063	351	150.012	0	0
```

**Disabling Gene Output:**

To disable gene-level output:

```bash
--quantVBgenes 0
```

**Notes:**
- Length/EffectiveLength use TPM-weighted averaging for expressed genes
- Genes with zero expression use unweighted averages
- The output format matches Salmon's `--geneMap` output (trailing zeros trimmed)

## Algorithm Details

### Equivalence Class Construction

TranscriptVB builds equivalence classes inline during alignment:

1. Each read's transcriptomic alignments are collected via `quantTranscriptome()`
2. Alignments are grouped by transcript ID set (the EC signature)
3. ECs with identical signatures are merged, incrementing the count
4. Thread-local EC tables are merged after mapping completes

### Variational Bayes Quantification

The VB algorithm iterates:

1. **E-step**: Compute expected transcript assignments using `exp(digamma(α))` weights
2. **M-step**: Update Dirichlet parameters `α = prior + expected_counts`
3. **Convergence**: Stop when relative change < 10⁻⁸ or max iterations (200)

### Effective Length Computation

Effective lengths are computed using fragment length distribution (FLD) adjustment:

```
effective_length = max(1, length - mean_fragment_length + 1)
```

When `--quantVBgcBias 1` is enabled, GC observations from properly-paired fragments are collected for potential bias correction.

## Integration Test Results

TranscriptVB was validated against Salmon using the same STAR alignment BAM file:

### Test Dataset
- **Source**: nf-core test data (GSE110004, yeast)
- **Reads**: SRR6357070 paired-end (50,000 read pairs)
- **Reference**: Yeast transcriptome (124 transcripts)

### Comparison: STAR-TranscriptVB vs Salmon (same BAM input)

| Metric | STAR-TranscriptVB | Salmon |
|--------|------------------|--------|
| Equivalence Classes | 100 | 101 |
| Total Reads Quantified | 30,059 | 29,653 |
| Expressed Transcripts | 89 | 90 |

### Correlation Statistics

| Correlation | All Transcripts | Jointly Expressed |
|-------------|-----------------|-------------------|
| **Spearman** | 0.968 | **0.997** |
| **Pearson** | 0.999 | **0.999** |

### Interpretation

- **0.997+ Spearman correlation** on jointly expressed transcripts demonstrates near-identical quantification
- **100 vs 101 ECs** (< 1% difference) shows equivalent equivalence class construction
- Small count differences arise from:
  - Minor alignment filtering differences
  - Effective length computation variations
  - Floating-point precision in VB convergence

## Requirements

- STAR genome index must include transcriptome annotations (GTF file)
- Paired-end reads recommended for GC bias collection
- Works with both SAM and BAM output modes

## Comparison with Other Methods

| Feature | TranscriptVB | Salmon | RSEM |
|---------|-------------|--------|------|
| Alignment | STAR (genomic) | Quasi-mapping | Bowtie2 |
| EC Building | Inline | Separate | Separate |
| Algorithm | VB/EM | VB/EM | EM |
| GC Bias | Collected | Full correction | No |
| Output Format | quant.sf | quant.sf | .isoforms.results |

## When to Use TranscriptVB

**Use TranscriptVB when:**
- You want transcript quantification integrated into STAR alignment
- You need Salmon-compatible output format
- You're doing bulk RNA-seq analysis
- You want to avoid running a separate quantification tool

**Use Salmon directly when:**
- You need full GC bias correction with transcript sequences
- You want sequence-bias or position-bias correction
- You prefer quasi-mapping over genomic alignment

## Troubleshooting

### No quant.sf output
- Ensure `--quantMode TranscriptVB` is specified
- Check that genome index includes transcriptome info (GTF at index generation)

### All zeros in output
- Verify reads map to annotated transcripts
- Check alignment rate in Log.final.out

### Low correlation with Salmon
- Ensure comparing against Salmon alignment-mode (`salmon quant -a`)
- Not Salmon quasi-mapping mode (different alignment approach)

## References

- Salmon: Patro et al. (2017) "Salmon provides fast and bias-aware quantification of transcript expression"
- STAR: Dobin et al. (2013) "STAR: ultrafast universal RNA-seq aligner"
- Variational Bayes for RNA-seq: Based on collapsed VB inference for Dirichlet-Multinomial models

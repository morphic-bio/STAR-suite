# tximport Algorithm Notes

**Reference version**: tximport 1.38.2 (Bioconductor 3.22)
**Default mode**: `countsFromAbundance = "lengthScaledTPM"` for Salmon-type input

## Overview

tximport summarizes transcript-level quantification to gene-level using a
tx2gene mapping. The `lengthScaledTPM` mode produces counts that are:
1. Proportional to TPM (abundance)
2. Scaled by effective length
3. Normalized to preserve total read count

## Input Columns (Salmon/STAR quant.sf)

| Column           | Description                          | Used for             |
|------------------|--------------------------------------|----------------------|
| `Name`           | Transcript ID                        | tx2gene lookup       |
| `Length`         | Transcript length (bp)               | Fallback only        |
| `EffectiveLength`| Effective length after bias model    | **Primary length**   |
| `TPM`            | Transcripts Per Million              | Abundance            |
| `NumReads`       | Estimated read count (fractional)    | Original counts      |

**For Salmon**: tximport uses `EffectiveLength` as the length column.

## Algorithm: Gene-Level Summarization

### Step 1: Transcript Filtering

```r
# Filter transcripts not in tx2gene
sub.idx <- txId %in% tx2gene$tx
abundanceMatTx <- abundanceMatTx[sub.idx, ]
countsMatTx <- countsMatTx[sub.idx, ]
lengthMatTx <- lengthMatTx[sub.idx, ]
```

**Behavior**: Transcripts in quant.sf but missing from tx2gene are silently dropped
(with a message: "transcripts missing from tx2gene: N").

### Step 2: Gene-Level Aggregation

```r
# Sum TPM and counts by gene
abundanceMat <- rowsum(abundanceMatTx, geneId)  # gene TPM = sum of transcript TPMs
countsMat <- rowsum(countsMatTx, geneId)        # gene counts = sum of transcript counts

# TPM-weighted average effective length
weightedLength <- rowsum(abundanceMatTx * lengthMatTx, geneId)
lengthMat <- weightedLength / abundanceMat
```

**Gene effective length** = Σ(TPM_tx × EffLen_tx) / Σ(TPM_tx) for transcripts in that gene

### Step 3: Handle Missing/Zero Lengths

When a gene has zero abundance (TPM = 0), the weighted average length is 0/0 = NaN.
tximport falls back to the unweighted average:

```r
aveLengthSamp <- rowMeans(lengthMatTx)           # avg length per transcript (across samples)
aveLengthSampGene <- tapply(aveLengthSamp, geneId, mean)  # avg of those per gene
lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGene)
```

### Step 4: lengthScaledTPM Counts

The `makeCountsFromAbundance()` function transforms counts:

```r
# For lengthScaledTPM:
countsSum <- colSums(countsMat)                  # total gene counts per sample
newCounts <- abundanceMat * rowMeans(lengthMat)  # TPM × average gene length
newSum <- colSums(newCounts)
countsMat <- t(t(newCounts) * (countsSum / newSum))  # scale to preserve total
```

## Mathematical Formula (Single Sample)

For gene g with transcripts T_g:

```
# Step 2: Gene TPM and length
gene_tpm[g] = Σ(tpm[t]) for t ∈ T_g
gene_efflen[g] = Σ(tpm[t] × efflen[t]) / gene_tpm[g]  (TPM-weighted avg)

# Step 4: lengthScaledTPM counts
raw_counts[g] = Σ(numreads[t]) for t ∈ T_g
length_scaled[g] = gene_tpm[g] × gene_efflen[g]
total_raw = Σ(raw_counts[g]) for all g
total_scaled = Σ(length_scaled[g]) for all g

final_counts[g] = length_scaled[g] × (total_raw / total_scaled)
```

**Key insight**: The final counts are proportional to `TPM × length`, normalized
to preserve the total read count. This corrects for length bias in TPM.

## Edge Cases

### Zero-TPM Genes
- Gene TPM = 0 → effective length uses fallback (unweighted average)
- Gene counts = 0 (unchanged)

### Transcripts Missing from tx2gene
- Silently filtered out
- Message printed: "transcripts missing from tx2gene: N"
- No error unless ALL transcripts are missing

### Zero/Invalid Effective Length
- If EffectiveLength = 0 for a transcript, it contributes 0 to weighted average
- Fallback handles the NaN case at gene level

## Output Columns

| Column           | Description                                    |
|------------------|------------------------------------------------|
| `Name`           | Gene ID (from tx2gene)                         |
| `Length`         | Same as EffectiveLength (tximport's length matrix) |
| `EffectiveLength`| TPM-weighted average effective length          |
| `TPM`            | Gene TPM (sum of transcript TPMs)              |
| `NumReads`       | lengthScaledTPM counts (see formula above)     |

**Note**: tximport's output `length` matrix contains TPM-weighted effective length,
identical to `EffectiveLength`. Both columns output the same value for compatibility.

## Gene Ordering

tximport outputs genes in **alphabetical order** by gene ID:
- R's `rowsum()` function sorts groups using locale-aware collation
- This is case-insensitive alphabetical order (e.g., "HRA1" < "snR18" < "YAL001C")

## Implementation Notes for C++

1. Use `EffectiveLength` as primary length (Salmon behavior)
2. Sum transcript TPM/counts by gene using a map
3. Compute TPM-weighted effective length per gene
4. Handle zero-TPM genes with unweighted fallback
5. Apply the count scaling formula at the end
6. Preserve gene ordering from tx2gene (first occurrence)

## Verification

To verify parity with R tximport:

```r
library(tximport)
tx2gene <- read.delim("tx2gene.tsv", header=FALSE, col.names=c("tx","gene"))
txi <- tximport("quant.sf", type="salmon", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM")
write.table(data.frame(Name=rownames(txi$counts),
                       Length=txi$length[,1],
                       EffectiveLength=txi$length[,1],
                       TPM=txi$abundance[,1],
                       NumReads=txi$counts[,1]),
            "quant.genes.tximport.sf", sep="\t", quote=FALSE, row.names=FALSE)
```

Note: tximport's output `length` matrix contains the TPM-weighted effective length,
not the raw transcript length. For the output file, we use this as both Length and
EffectiveLength (they're conceptually the same at gene level with this weighting).



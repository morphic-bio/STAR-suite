# TranscriptVB Testing Plan

## Overview

This document outlines the testing strategy for the TranscriptVB (Salmon-style) transcript quantification feature. Tests are organized by category with specific test cases, expected outcomes, and validation criteria.

> **Status**: Integration testing complete. Additional testing recommended before production use.

---

## Test Categories

1. [Basic Functionality](#1-basic-functionality)
2. [Edge Cases](#2-edge-cases)
3. [Numerical Stability](#3-numerical-stability)
4. [Library Types](#4-library-types)
5. [Performance & Scaling](#5-performance--scaling)
6. [Regression Testing](#6-regression-testing)
7. [Salmon Parity](#7-salmon-parity)

---

## 1. Basic Functionality

### Test 1.1: Basic PE Quantification
**Purpose**: Verify basic paired-end quantification works

```bash
STAR --runMode alignReads \
     --genomeDir $GENOME_DIR \
     --readFilesIn reads_1.fq.gz reads_2.fq.gz \
     --readFilesCommand zcat \
     --quantMode TranscriptVB \
     --runThreadN 4 \
     --outFileNamePrefix test_basic_
```

**Validation**:
- [ ] `quant.sf` file created
- [ ] Header line: `Name\tLength\tEffectiveLength\tTPM\tNumReads`
- [ ] TPM column sums to ~1,000,000
- [ ] NumReads > 0 for expressed transcripts
- [ ] Log shows "Quantification converged: yes"

### Test 1.2: GC Bias Collection
**Purpose**: Verify GC observation collection for PE reads

```bash
STAR --runMode alignReads \
     --genomeDir $GENOME_DIR \
     --readFilesIn reads_1.fq.gz reads_2.fq.gz \
     --readFilesCommand zcat \
     --quantMode TranscriptVB \
     --quantVBgcBias 1 \
     --runThreadN 4 \
     --outFileNamePrefix test_gc_
```

**Validation**:
- [ ] Log shows "GC bias: collected N fragment observations" where N > 0
- [ ] Log shows "GC bias: using FLD-adjusted effective lengths"
- [ ] Effective lengths differ from raw lengths

### Test 1.3: EM Mode
**Purpose**: Verify EM algorithm works as alternative to VB

```bash
STAR --runMode alignReads \
     --genomeDir $GENOME_DIR \
     --readFilesIn reads_1.fq.gz reads_2.fq.gz \
     --readFilesCommand zcat \
     --quantMode TranscriptVB \
     --quantVBem 1 \
     --runThreadN 4 \
     --outFileNamePrefix test_em_
```

**Validation**:
- [ ] Quantification completes successfully
- [ ] Results correlate highly with VB mode (r > 0.99)

### Test 1.4: Custom Prior
**Purpose**: Verify VB prior parameter affects results

```bash
# Low prior (sparse solution)
STAR --quantMode TranscriptVB --quantVBprior 0.001 ...

# High prior (smoothed solution)  
STAR --quantMode TranscriptVB --quantVBprior 1.0 ...
```

**Validation**:
- [ ] Both runs complete successfully
- [ ] Low prior produces more zeros
- [ ] High prior produces fewer zeros

---

## 2. Edge Cases

### Test 2.1: Single-End Reads
**Purpose**: Verify SE quantification (no GC bias collection expected)

```bash
STAR --runMode alignReads \
     --genomeDir $GENOME_DIR \
     --readFilesIn reads.fq.gz \
     --readFilesCommand zcat \
     --quantMode TranscriptVB \
     --outFileNamePrefix test_se_
```

**Validation**:
- [ ] Quantification completes
- [ ] No crash from PE-specific code paths
- [ ] Results are reasonable (compare with Salmon SE mode)

### Test 2.2: Very Short Transcripts
**Purpose**: Verify handling of transcripts shorter than mean fragment length

**Test Data**: Include transcripts < 100bp in reference

**Validation**:
- [ ] No division by zero errors
- [ ] Effective length >= 1 for all transcripts
- [ ] Short transcripts can still have counts

### Test 2.3: Zero-Count Transcripts
**Purpose**: Verify transcripts with no alignments handled correctly

**Validation**:
- [ ] TPM = 0.0 for unaligned transcripts
- [ ] NumReads = 0.0 for unaligned transcripts
- [ ] No NaN or Inf values in output

### Test 2.4: Highly Multi-Mapping Reads
**Purpose**: Verify handling of reads mapping to many transcripts

**Test Data**: Use gene families or repetitive elements

**Validation**:
- [ ] No memory explosion
- [ ] EC count remains reasonable
- [ ] Quantification converges

### Test 2.5: Empty Input
**Purpose**: Verify graceful handling of no alignable reads

```bash
# Use reads that won't map to reference
STAR --quantMode TranscriptVB --readFilesIn unmappable_reads.fq.gz ...
```

**Validation**:
- [ ] No crash
- [ ] All TPM = 0
- [ ] Warning or info message about low mapping rate

### Test 2.6: Single Transcript Reference
**Purpose**: Verify minimal reference works

**Test Data**: Reference with only 1-2 transcripts

**Validation**:
- [ ] Quantification works
- [ ] Results match expected (100% to single transcript if only one)

---

## 3. Numerical Stability

### Test 3.1: Underflow Protection
**Purpose**: Verify very low abundance transcripts don't cause underflow

**Test Data**: Large transcriptome with many zero/near-zero transcripts

**Validation**:
- [ ] No `-inf` or `nan` in output
- [ ] Very low TPM values are positive (not negative)
- [ ] Log-likelihood is finite

### Test 3.2: Convergence Behavior
**Purpose**: Verify VB converges within iteration limit

```bash
# Check different datasets
for dataset in small medium large; do
    STAR --quantMode TranscriptVB ... 2>&1 | grep "iterations"
done
```

**Validation**:
- [ ] All datasets converge (iterations < 200)
- [ ] Final log-likelihood is reported
- [ ] No "did not converge" warnings

### Test 3.3: Numerical Precision
**Purpose**: Verify deterministic results across runs

```bash
# Run twice with same input
STAR --quantMode TranscriptVB --runThreadN 1 ... -o run1_
STAR --quantMode TranscriptVB --runThreadN 1 ... -o run2_

diff run1_quant.sf run2_quant.sf
```

**Validation**:
- [ ] Identical results with single thread
- [ ] Results within floating-point tolerance with multiple threads

---

## 4. Library Types

### Test 4.1: Unstranded Library
**Purpose**: Verify unstranded PE library quantification

**Validation**:
- [ ] Reads assigned to both strands appropriately
- [ ] Antisense transcripts handled correctly

### Test 4.2: Stranded Library (FR-firststrand)
**Purpose**: Verify stranded library support

**Note**: Current implementation uses unstranded model. Test to document behavior.

**Validation**:
- [ ] Document whether strand info is used
- [ ] Compare with Salmon stranded mode if applicable

---

## 5. Performance & Scaling

### Test 5.1: Thread Scaling
**Purpose**: Verify multi-threaded EC merging works correctly

```bash
for threads in 1 2 4 8 16; do
    time STAR --quantMode TranscriptVB --runThreadN $threads ...
done
```

**Validation**:
- [ ] All thread counts produce same results (within tolerance)
- [ ] Runtime scales reasonably with threads
- [ ] No race conditions or crashes

### Test 5.2: Memory Usage
**Purpose**: Profile memory with large transcriptomes

**Test Data**: Human transcriptome (~200K transcripts)

```bash
/usr/bin/time -v STAR --quantMode TranscriptVB ... 2>&1 | grep "Maximum resident"
```

**Validation**:
- [ ] Memory usage is reasonable (< 2x base STAR memory)
- [ ] No memory leaks (stable across multiple runs)

### Test 5.3: Large Dataset Performance
**Purpose**: Benchmark on realistic dataset

**Test Data**: 50M PE reads, human transcriptome

**Validation**:
- [ ] Completes in reasonable time
- [ ] Compare runtime to `STAR + Salmon` two-step approach

---

## 6. Regression Testing

### Test 6.1: Golden Output Comparison
**Purpose**: Detect unintended changes in quantification

**Setup**:
1. Run TranscriptVB on reference dataset
2. Save `quant.sf` as golden reference
3. Future runs compare against golden

```bash
# Generate golden (one-time)
STAR --quantMode TranscriptVB ... -o golden_
cp golden_quant.sf tests/golden/yeast_pe_quant.sf

# Regression test
STAR --quantMode TranscriptVB ... -o test_
python3 compare_quant.py tests/golden/yeast_pe_quant.sf test_quant.sf
```

**Validation**:
- [ ] Pearson correlation > 0.9999 with golden
- [ ] Max absolute difference in TPM < 1.0
- [ ] Same number of expressed transcripts

### Test 6.2: Version Compatibility
**Purpose**: Ensure genome indices work across versions

**Validation**:
- [ ] Document minimum genome index version
- [ ] Test with indices from different STAR versions

---

## 7. Salmon Parity

### Test 7.1: Same-BAM Comparison (Primary Validation)
**Purpose**: Validate against Salmon using identical alignments

```bash
# Generate transcriptome BAM
STAR --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted ...

# Run Salmon on BAM
salmon quant -t transcriptome.fa -l A -a Aligned.toTranscriptome.out.bam -o salmon_out

# Run TranscriptVB
STAR --quantMode TranscriptVB ...

# Compare
python3 compare_salmon_star.py salmon_out/quant.sf star_quant.sf
```

**Validation Criteria**:
| Metric | Threshold |
|--------|-----------|
| Spearman correlation (all) | > 0.95 |
| Spearman correlation (expressed) | > 0.99 |
| Pearson correlation | > 0.99 |
| EC count difference | < 5% |

### Test 7.2: Large Dataset Parity
**Purpose**: Verify parity holds on larger, more complex data

**Test Data**: Human RNA-seq (e.g., ENCODE or GTEx sample)

**Validation**:
- [ ] Correlation > 0.95 across ~60K expressed transcripts
- [ ] No systematic biases by expression level
- [ ] Similar handling of multi-mappers

### Test 7.3: Parameter Equivalence
**Purpose**: Verify equivalent parameters produce equivalent results

| TranscriptVB | Salmon Equivalent |
|--------------|-------------------|
| `--quantVBprior 0.01` | Default VB prior |
| `--quantVBem 1` | `--useEM` |

---

## Test Data Sources

### Small Test Data (Quick Tests)
- **Location**: `/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/`
- **Species**: Yeast
- **Files**: `SRR6357070_1.fastq.gz`, `SRR6357070_2.fastq.gz`
- **Transcripts**: 124

### Medium Test Data (Integration Tests)
- **Source**: ENCODE or similar
- **Species**: Human chr22 subset
- **Reads**: 1-5M PE reads

### Large Test Data (Performance Tests)
- **Source**: Full human RNA-seq
- **Species**: Human (GRCh38)
- **Reads**: 50M+ PE reads
- **Transcripts**: ~200K

---

## Test Execution Scripts

### Quick Validation Script

```bash
#!/bin/bash
# quick_test.sh - Run basic validation tests

STAR_BIN=/path/to/STAR
GENOME_DIR=/path/to/genome
READS_1=/path/to/reads_1.fq.gz
READS_2=/path/to/reads_2.fq.gz
OUTDIR=/tmp/transcriptvb_test

mkdir -p $OUTDIR
cd $OUTDIR

echo "=== Test 1: Basic TranscriptVB ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --runThreadN 4 \
    --outFileNamePrefix basic_

if [ -f basic_quant.sf ]; then
    echo "✓ quant.sf created"
    TPM_SUM=$(tail -n +2 basic_quant.sf | awk '{sum+=$4} END {print sum}')
    echo "  TPM sum: $TPM_SUM (expected ~1000000)"
else
    echo "✗ quant.sf NOT created"
    exit 1
fi

echo "=== Test 2: GC Bias Mode ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --quantVBgcBias 1 \
    --runThreadN 4 \
    --outFileNamePrefix gc_

grep -q "GC bias: collected" gc_Log.out && echo "✓ GC observations collected" || echo "✗ No GC observations"

echo "=== Test 3: EM Mode ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --quantVBem 1 \
    --runThreadN 4 \
    --outFileNamePrefix em_

[ -f em_quant.sf ] && echo "✓ EM mode completed" || echo "✗ EM mode failed"

echo "=== All quick tests completed ==="
```

### Salmon Comparison Script

```python
#!/usr/bin/env python3
# compare_salmon_star.py - Compare TranscriptVB with Salmon output

import pandas as pd
import sys
from scipy.stats import spearmanr, pearsonr

def compare(salmon_file, star_file):
    salmon = pd.read_csv(salmon_file, sep='\t')
    star = pd.read_csv(star_file, sep='\t')
    
    merged = salmon.merge(star, on='Name', suffixes=('_salmon', '_star'))
    
    # All transcripts
    r_spear_all, _ = spearmanr(merged['NumReads_salmon'], merged['NumReads_star'])
    r_pear_all, _ = pearsonr(merged['NumReads_salmon'], merged['NumReads_star'])
    
    # Jointly expressed
    expressed = merged[(merged['NumReads_salmon'] > 0) & (merged['NumReads_star'] > 0)]
    r_spear_expr, _ = spearmanr(expressed['NumReads_salmon'], expressed['NumReads_star'])
    r_pear_expr, _ = pearsonr(expressed['NumReads_salmon'], expressed['NumReads_star'])
    
    print(f"=== Salmon vs TranscriptVB Comparison ===")
    print(f"Total transcripts: {len(merged)}")
    print(f"Expressed in both: {len(expressed)}")
    print(f"\nAll transcripts:")
    print(f"  Spearman: {r_spear_all:.4f}")
    print(f"  Pearson:  {r_pear_all:.4f}")
    print(f"\nJointly expressed:")
    print(f"  Spearman: {r_spear_expr:.4f}")
    print(f"  Pearson:  {r_pear_expr:.4f}")
    
    # Pass/fail
    if r_spear_expr > 0.99 and r_pear_expr > 0.99:
        print("\n✓ PASS: High correlation with Salmon")
        return 0
    else:
        print("\n✗ FAIL: Correlation below threshold")
        return 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} salmon_quant.sf star_quant.sf")
        sys.exit(1)
    sys.exit(compare(sys.argv[1], sys.argv[2]))
```

---

## Current Test Status

| Test Category | Status | Notes |
|---------------|--------|-------|
| Basic PE | ✅ Pass | Tested with yeast data |
| GC Bias Collection | ✅ Pass | 245 observations collected |
| EM Mode | 🔲 Not tested | |
| Single-End | 🔲 Not tested | |
| Large Dataset | 🔲 Not tested | |
| Salmon Parity | ✅ Pass | r=0.997 on same BAM |
| Thread Scaling | 🔲 Not tested | |
| Regression Suite | 🔲 Not created | |

---

## Test Environment

**Validated Configuration**:
- STAR version: 2.7.11b (STAR-Flex)
- Salmon version: (for comparison)
- Test OS: Linux (Ubuntu)
- Reference: Yeast transcriptome (124 transcripts)

---

## Appendix: Integration Test Results

### Yeast PE Dataset (GSE110004)

**Comparison: STAR-TranscriptVB vs Salmon (same BAM)**

| Metric | STAR-TranscriptVB | Salmon |
|--------|------------------|--------|
| Equivalence Classes | 100 | 101 |
| Total Reads | 30,059 | 29,653 |
| Expressed Transcripts | 89 | 90 |
| **Spearman (expressed)** | **0.997** | - |
| **Pearson (expressed)** | **0.999** | - |

**Conclusion**: Near-identical quantification when using same input alignments.


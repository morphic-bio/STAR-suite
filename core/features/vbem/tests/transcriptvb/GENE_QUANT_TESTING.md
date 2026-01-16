# Gene Quant Output Testing Guide

**Date:** December 29, 2025  
**Feature:** TranscriptVB Gene-Level Quantification (`quant.genes.sf`)

---

## Overview

This document describes the test suite for validating gene-level quantification output (`quant.genes.sf`) in TranscriptVB mode.

---

## Test Scripts

### 1. `quick_test.sh` - Basic Validation

**Location:** `tests/transcriptvb/quick_test.sh`

**Usage:**
```bash
./quick_test.sh [STAR_BIN] [GENOME_DIR] [READS_1] [READS_2]
```

**What it tests:**
- **Test 1b**: Gene-level output creation and format validation
  - Verifies `quant.genes.sf` is created by default
  - Checks header format: `Name\tLength\tEffectiveLength\tTPM\tNumReads`
  - Validates no scientific notation
  - Validates no NaN/Inf values
  - Checks decimal precision per column (Length: 3, EffectiveLength: 4, TPM: 6, NumReads: 3)
  - Compares gene TPM sum with transcript TPM sum

- **Test 5**: Gene quant disabled flag
  - Verifies `--quantVBgenes 0` prevents `quant.genes.sf` creation

**Example:**
```bash
cd tests/transcriptvb
./quick_test.sh
```

---

### 2. `spot_check_gene_format.sh` - Format Validation

**Location:** `tests/transcriptvb/spot_check_gene_format.sh`

**Usage:**
```bash
./spot_check_gene_format.sh [quant.genes.sf] [reference.genes.sf]
```

**What it checks:**
1. **Header Format**: Exact match with expected header
2. **No Scientific Notation**: Grep for `[eE][+-]` patterns
3. **No NaN/Inf**: Case-insensitive check for invalid values
4. **Decimal Precision**: Spot-check first 10 lines for per-column precision limits
5. **Zero Value Formatting**: Ensures zeros are formatted as `0` not `0.000000`
6. **Sample Output**: Displays first 5 data lines

**Example:**
```bash
cd tests/transcriptvb
./spot_check_gene_format.sh basic_quant.genes.sf
```

---

### 3. `salmon_parity_test.sh` - Parity Testing

**Location:** `tests/transcriptvb/salmon_parity_test.sh`

**Usage:**
```bash
./salmon_parity_test.sh [GENOME_DIR] [TRANSCRIPTOME_FASTA] [READS_1] [READS_2]
```

**What it tests:**
- **Step 4**: Transcript-level comparison (existing)
- **Step 5**: Gene-level comparison (NEW)
  - Compares STAR `quant.genes.sf` with Salmon `quant.genes.sf` (if available)
  - Requires Salmon to be run with `--geneMap` flag to generate gene-level output

**To generate Salmon gene-level output:**
```bash
salmon quant \
    -t transcriptome.fa \
    -l A \
    -a aligned.bam \
    --geneMap gene_map.txt \
    -o salmon_out
```

**Example:**
```bash
cd tests/transcriptvb
./salmon_parity_test.sh
```

---

## Expected Output Format

### Header
```
Name	Length	EffectiveLength	TPM	NumReads
```

### Example Rows
```
ENSG00000268674	510	295.805	0	0
ENSG00000271254	3081.77	2865.89	9.91087	491.292
ENSG00000275063	351	150.012	0	0
```

### Format Rules
- **Length**: Max 3 decimal places (e.g., `510`, `3081.77`)
- **EffectiveLength**: Max 4 decimal places (e.g., `295.805`, `2865.89`)
- **TPM**: Max 6 decimal places (e.g., `0`, `9.91087`)
- **NumReads**: Max 3 decimal places (e.g., `0`, `491.292`)
- **No scientific notation**: All values use decimal format
- **Trailing zeros trimmed**: `510` not `510.000`, `0` not `0.000000`

---

## Validation Checklist

When testing gene quant output, verify:

- [ ] `quant.genes.sf` exists after TranscriptVB run
- [ ] Header matches exactly: `Name	Length	EffectiveLength	TPM	NumReads`
- [ ] Line count = number of genes + 1 (header)
- [ ] No scientific notation (`[eE][+-]`)
- [ ] No NaN or Inf values
- [ ] Decimal precision within limits per column
- [ ] Zero values formatted as `0` (not `0.000000`)
- [ ] Gene TPM sum ≈ transcript TPM sum (within tolerance)
- [ ] Gene NumReads sum ≈ transcript NumReads sum (within tolerance)
- [ ] `--quantVBgenes 0` prevents file creation

---

## Running Tests

### Quick Validation
```bash
cd tests/transcriptvb
./quick_test.sh
```

### Format Spot-Check
```bash
cd tests/transcriptvb
./spot_check_gene_format.sh output_quant.genes.sf
```

### Full Parity Test (requires Salmon)
```bash
cd tests/transcriptvb
# First, generate Salmon gene-level output with --geneMap
./salmon_parity_test.sh
```

---

## Troubleshooting

### quant.genes.sf not created
- Check that `--quantMode TranscriptVB` is specified
- Verify `--quantVBgenes` is not set to `0`
- Check Log.out for errors

### Precision issues
- Verify `formatNum()` is called with correct precision per column
- Check that trailing zeros are trimmed correctly

### Scientific notation found
- Ensure `std::fixed` is used in `formatNum()`
- Verify precision limits are appropriate

### TPM sum mismatch
- Gene TPM sum should approximately equal transcript TPM sum
- Small differences (< 0.1%) are expected due to floating-point rounding

---

## Reference Files

- **Salmon reference**: `/mnt/pikachu/globus_mnt/JAX_RNAseq13_processed/Counts/CITED_HOXB8_WT1_GT25-03671_TTCGCAGTNNNNNNNNN-TCCGATCA_S35_L002_R1_001_val/quant.genes.sf`
- **Format specification**: See `docs/TranscriptVB_quantification.md`

---

## Related Documentation

- **Implementation Plan**: `plans/gene_quant_implementation_plan.md`
- **Design Decisions**: `plans/gene_quant_questions.md`
- **User Documentation**: `docs/TranscriptVB_quantification.md`

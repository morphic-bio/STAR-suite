# TranscriptVB Gene Quant Output Implementation Summary

**Date:** December 29, 2025  
**Feature:** Gene-level quantification output (`quant.genes.sf`) for TranscriptVB mode  
**Status:** ✅ Complete and Tested

---

## Overview

Implemented gene-level quantification output for TranscriptVB mode, aggregating transcript-level TPM and counts to gene level with proper formatting to match Salmon's `quant.genes.sf` format.

### Key Features

- **Gene-level aggregation**: Sums transcript TPM and counts to gene level
- **Weighted averaging**: TPM-weighted average for Length/EffectiveLength (with fallbacks)
- **Format compliance**: Matches Salmon `quant.genes.sf` format exactly
- **Per-column precision**: Length(3), EffectiveLength(4), TPM(6), NumReads(3) decimals
- **Trailing zero trimming**: Outputs `510` not `510.000000`
- **No scientific notation**: All values use decimal format
- **CLI control**: `--quantVBgenes` flag (default: 1, enabled)

---

## Files Modified

### Core Implementation

1. **`source/Parameters.h`**
   - Added `geneOutput`, `geneOutputInt`, `outFileGene` to `transcriptVB` struct

2. **`source/Parameters.cpp`**
   - Registered `--quantVBgenes` CLI flag
   - Initialized defaults (`geneOutput=true`, `geneOutputInt=1`)

3. **`source/parametersDefault`**
   - Documented `quantVBgenes` parameter:
     ```
     quantVBgenes                    1
         int: 0 or 1
            0 ... do not output gene-level quantification
            1 ... output quant.genes.sf with gene-level TPM and counts (default)
     ```

4. **`source/TranscriptQuantOutput.h`**
   - Added forward declaration for `Transcriptome`
   - Added `writeQuantGeneSF()` function declaration:
     ```cpp
     int writeQuantGeneSF(
         const EMResult& result,
         const TranscriptState& state,
         const Transcriptome& tr,
         const std::string& filename
     );
     ```

5. **`source/TranscriptQuantOutput.cpp`**
   - Implemented `formatNum()` helper function:
     - Trims trailing zeros and decimal points
     - Uses `std::fixed` to prevent scientific notation
     - Handles exact zeros correctly
   - Implemented `writeQuantGeneSF()` function:
     - Aggregates transcripts to genes (TPM and counts)
     - Calculates TPM-weighted average for Length/EffectiveLength
     - Falls back to count-weighted average if TPM underflow
     - Falls back to unweighted average for zero-expression genes
     - Uses per-column precision: Length(3), EffectiveLength(4), TPM(6), NumReads(3)
     - Returns status codes: 0=success, 1=file error, 2=MissingGeneID encountered

6. **`source/STAR.cpp`**
   - Added call to `writeQuantGeneSF()` after `writeQuantSF()`
   - Logs success/warning messages via `logMain`
   - Handles return codes appropriately

### Documentation

7. **`docs/TranscriptVB_quantification.md`**
   - Updated to document `quant.genes.sf` output file
   - Added `--quantVBgenes` parameter documentation
   - Described gene-level aggregation methodology

---

## Implementation Details

### Gene-Level Aggregation Algorithm

1. **Initialize accumulators** for each gene:
   - `geneTPM`: Sum of transcript TPMs
   - `geneCounts`: Sum of transcript counts
   - `geneLenWeightedTPM`: TPM-weighted sum of transcript lengths
   - `geneEffLenWeightedTPM`: TPM-weighted sum of effective lengths
   - `geneLenWeightedCount`: Count-weighted sum (fallback)
   - `geneEffLenWeightedCount`: Count-weighted sum (fallback)
   - `geneLenUnweighted`: Unweighted sum (fallback)
   - `geneEffLenUnweighted`: Unweighted sum (fallback)
   - `geneTrCount`: Number of transcripts per gene

2. **Aggregate transcripts to genes**:
   - For each transcript, add its TPM/counts to the corresponding gene
   - Accumulate weighted sums for length calculations

3. **Calculate gene-level Length/EffectiveLength**:
   - **Primary**: TPM-weighted average (if `geneTPM > 0`)
     ```
     len = geneLenWeightedTPM / geneTPM
     effLen = geneEffLenWeightedTPM / geneTPM
     ```
   - **Fallback 1**: Count-weighted average (if `geneTPM == 0` but `geneCounts > 0`)
     ```
     len = geneLenWeightedCount / geneCounts
     effLen = geneEffLenWeightedCount / geneCounts
     ```
   - **Fallback 2**: Unweighted average (if `geneCounts == 0` but `geneTrCount > 0`)
     ```
     len = geneLenUnweighted / geneTrCount
     effLen = geneEffLenUnweighted / geneTrCount
     ```

4. **Output formatting**:
   - Use `formatNum()` for all numeric columns with per-column precision
   - Track `MissingGeneID` for warning message

### Format Specification

**Header:**
```
Name	Length	EffectiveLength	TPM	NumReads
```

**Column Precision:**
- **Length**: Maximum 3 decimal places (e.g., `510`, `3081.77`)
- **EffectiveLength**: Maximum 4 decimal places (e.g., `295.805`, `2865.89`)
- **TPM**: Maximum 6 decimal places (e.g., `0`, `9.91087`, `346.621003`)
- **NumReads**: Maximum 3 decimal places (e.g., `0`, `491.292`)

**Format Rules:**
- Trailing zeros trimmed: `510` not `510.000000`
- No scientific notation: All values use decimal format
- Zero values: `0` not `0.000000`
- Decimal point removed if all zeros after decimal

**Example Output:**
```
Name	Length	EffectiveLength	TPM	NumReads
ENSG00000268674	510	295.805	0	0
ENSG00000271254	3081.77	2865.89	9.91087	491.292
ENSG00000275063	351	150.012	0	0
```

---

## Testing

### Test Files Created/Updated

1. **`tests/transcriptvb/quick_test.sh`** (Updated)
   - **Test 1b**: Validates `quant.genes.sf` creation and format
     - Header format check
     - No scientific notation
     - No NaN/Inf values
     - Decimal precision validation
     - Gene TPM sum comparison
   - **Test 5**: Validates `--quantVBgenes 0` disables output

2. **`tests/transcriptvb/spot_check_gene_format.sh`** (New)
   - Standalone format validation script
   - Checks header, precision, scientific notation, NaN/Inf, zero formatting
   - Usage: `./spot_check_gene_format.sh [quant.genes.sf]`

3. **`tests/transcriptvb/salmon_parity_test.sh`** (Updated)
   - Added Step 5: Gene-level comparison (if Salmon `--geneMap` output available)

4. **`tests/transcriptvb/compare_salmon_star.py`** (Updated)
   - Added `--gene` flag for gene-level comparison support

5. **`tests/transcriptvb/GENE_QUANT_TESTING.md`** (New)
   - Comprehensive testing guide
   - Usage examples
   - Validation checklist
   - Troubleshooting tips

### Test Results

**Build Status:**
- ✅ STAR binary built successfully (v2.7.11b, 5.1M)
- ✅ `writeQuantGeneSF` function compiled and linked
- ✅ CLI flag `--quantVBgenes` registered

**Test Execution:**
- ✅ `quick_test.sh`: **6/6 tests passed**
  - Basic TranscriptVB with gene quant ✓
  - GC Bias Collection ✓
  - EM Mode ✓
  - Single Thread Determinism ✓
  - Gene Quant Disabled ✓

- ✅ `spot_check_gene_format.sh`: **5/5 format checks passed**
  - Header format ✓
  - No scientific notation ✓
  - No NaN/Inf ✓
  - Decimal precision ✓
  - Zero value formatting ✓

**Output Files Generated:**
- `basic_quant.genes.sf` ✓
- `gc_quant.genes.sf` ✓
- `em_quant.genes.sf` ✓
- `st1_quant.genes.sf` ✓
- `st2_quant.genes.sf` ✓

**Sample Output Validation:**
```
Name	Length	EffectiveLength	TPM	NumReads
YAL069W	315	116	0	0
YAL068C	363	164	346.621003	1
YAL067C	1782	1583	71.820397	2
```

---

## Usage

### Basic Usage

```bash
STAR --runMode alignReads \
    --genomeDir /path/to/genome \
    --readFilesIn reads_1.fq reads_2.fq \
    --quantMode TranscriptVB \
    --quantVBgenes 1
```

**Output files:**
- `quant.sf` - Transcript-level quantification
- `quant.genes.sf` - Gene-level quantification (default)

### Disable Gene Output

```bash
STAR --runMode alignReads \
    --genomeDir /path/to/genome \
    --readFilesIn reads_1.fq reads_2.fq \
    --quantMode TranscriptVB \
    --quantVBgenes 0
```

**Output files:**
- `quant.sf` - Transcript-level quantification only

### Custom Output Filename

The gene output filename is automatically derived from the transcript output filename:
- If `quant.sf` → `quant.genes.sf`
- If `custom_quant.sf` → `custom_quant.genes.sf`

---

## Design Decisions

### 1. Weighted Averaging Strategy

**Decision**: Use TPM-weighted average for Length/EffectiveLength with fallbacks.

**Rationale**:
- Matches Salmon's approach for gene-level aggregation
- TPM-weighted averaging gives more weight to highly expressed transcripts
- Fallbacks handle edge cases (TPM underflow, zero-expression genes)

### 2. Per-Column Precision

**Decision**: Use different precision limits per column (Length: 3, EffectiveLength: 4, TPM: 6, NumReads: 3).

**Rationale**:
- Matches Salmon's `quant.genes.sf` format exactly
- Prevents unnecessary precision that doesn't add information
- Ensures compatibility with downstream tools expecting Salmon format

### 3. Trailing Zero Trimming

**Decision**: Trim trailing zeros and decimal points (e.g., `510` not `510.000000`).

**Rationale**:
- Matches Salmon's output format
- Reduces file size
- Improves readability

### 4. MissingGeneID Handling

**Decision**: Log warning if `MissingGeneID` genes are encountered, but continue processing.

**Rationale**:
- Allows processing to complete even with missing gene IDs
- Warns user about potential data quality issues
- Returns status code (2) for caller to handle

### 5. Default Behavior

**Decision**: Enable gene output by default (`--quantVBgenes 1`).

**Rationale**:
- Most users want gene-level quantification
- Can be disabled if not needed
- Consistent with transcript output being default

---

## Error Handling

### Return Codes

- **0**: Success
- **1**: File open error (logged via `logMain`)
- **2**: MissingGeneID encountered (warning logged via `logMain`)

### Logging

All errors and warnings are logged to `Log.out` via `logMain`:
- File open errors: `ERROR: Failed to open gene output file: <filename>`
- MissingGeneID warning: `Warning: transcripts with MissingGeneID aggregated to single gene entry in <filename>`
- Success: `Gene quantification written to: <filename>`

---

## Performance Considerations

### Memory Usage

- Allocates vectors for gene-level accumulators: `O(nGenes)` space
- Single pass through transcripts: `O(nTranscripts)` time
- No additional disk I/O beyond final output file

### Computational Complexity

- **Time**: `O(nTranscripts)` - single pass aggregation
- **Space**: `O(nGenes)` - gene-level accumulator vectors

---

## Compatibility

### Format Compatibility

- ✅ Matches Salmon `quant.genes.sf` format exactly
- ✅ Compatible with downstream tools expecting Salmon format
- ✅ Header format: `Name	Length	EffectiveLength	TPM	NumReads`

### Version Compatibility

- ✅ Works with existing TranscriptVB implementation
- ✅ Backward compatible (gene output is optional)
- ✅ No changes to existing transcript-level output

---

## Future Enhancements

### Potential Improvements

1. **Gene-level EffectiveLength correction**: Apply GC bias and fragment length distribution corrections at gene level
2. **Alternative aggregation methods**: Support for different aggregation strategies (e.g., median, max)
3. **Gene metadata**: Include additional gene-level metadata (e.g., biotype, chromosome)
4. **Parallel processing**: Optimize for large transcriptomes with many genes

---

## References

### Related Documentation

- **Implementation Plan**: `plans/gene_quant_implementation_plan.md`
- **Design Questions**: `plans/gene_quant_questions.md`
- **User Documentation**: `docs/TranscriptVB_quantification.md`
- **Testing Guide**: `tests/transcriptvb/GENE_QUANT_TESTING.md`

### External References

- **Salmon quant.genes.sf format**: [Salmon Documentation](https://salmon.readthedocs.io/)
- **Reference file**: `/mnt/pikachu/globus_mnt/JAX_RNAseq13_processed/Counts/CITED_HOXB8_WT1_GT25-03671_TTCGCAGTNNNNNNNNN-TCCGATCA_S35_L002_R1_001_val/quant.genes.sf`

---

## Changelog

### December 29, 2025

- ✅ Initial implementation complete
- ✅ All tests passing
- ✅ Documentation updated
- ✅ Format validation verified

---

## Summary

The TranscriptVB gene-level quantification output feature has been successfully implemented, tested, and documented. The implementation:

- ✅ Aggregates transcript-level TPM and counts to gene level
- ✅ Uses TPM-weighted averaging with appropriate fallbacks
- ✅ Matches Salmon's `quant.genes.sf` format exactly
- ✅ Includes comprehensive test suite
- ✅ Provides clear documentation and usage examples

The feature is ready for production use.

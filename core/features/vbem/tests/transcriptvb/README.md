# TranscriptVB Test Suite

Test suite for validating STAR-Flex's integrated transcript quantification feature (TranscriptVB).

## Overview

This test suite provides:
- **Quick validation tests** - Basic functionality checks (~2 minutes)
- **Salmon parity tests** - Comparison with Salmon on same data
- **Regression tests** - Compare against golden outputs
- **Edge case tests** - Handle unusual inputs gracefully

## Quick Start

```bash
# Run quick validation
./quick_test.sh

# Run Salmon parity test
./salmon_parity_test.sh

# Generate golden references (one-time)
./regression_test.sh generate

# Run regression tests
./regression_test.sh test
```

## Test Scripts

### `quick_test.sh`

Fast validation of core functionality.

**Tests**:
1. Basic VB quantification
2. GC bias collection
3. EM mode
4. Single-thread determinism

**Usage**:
```bash
./quick_test.sh [STAR_BIN] [GENOME_DIR] [READS_1] [READS_2]

# With defaults
./quick_test.sh
```

**Runtime**: ~2 minutes

---

### `salmon_parity_test.sh`

Validates TranscriptVB against Salmon using the same input alignments.

**Steps**:
1. Generate transcriptome BAM with STAR
2. Run Salmon on the BAM
3. Run TranscriptVB
4. Compare outputs

**Usage**:
```bash
./salmon_parity_test.sh [GENOME_DIR] [TRANSCRIPTOME_FASTA] [READS_1] [READS_2]
```

**Requirements**:
- Salmon installed and in PATH
- Transcriptome FASTA file
- Optional gene map file for gene-level parity (use `make_gene_map_from_gtf.sh`)

**Runtime**: ~5 minutes

---

### `regression_test.sh`

Compares current outputs against saved golden references.

**Commands**:
```bash
# Generate golden references (do once after validating implementation)
./regression_test.sh generate

# Test against golden
./regression_test.sh test

# Update golden after code changes (backs up old)
./regression_test.sh update
```

**Thresholds**:
- Spearman correlation: ≥ 0.9999
- Pearson correlation: ≥ 0.9999
- Max TPM difference: ≤ 1.0

**Runtime**: ~5 minutes

---

### `edge_case_tests.sh`

Tests handling of unusual/edge case inputs.

**Tests**:
- Empty input files
- Very few reads
- Unmappable reads (random sequences)
- Single-end reads
- Different prior values
- Multi-threading

**Usage**:
```bash
./edge_case_tests.sh [STAR_BIN] [GENOME_DIR]
```

**Runtime**: ~10 minutes

---

### `compare_salmon_star.py`

Python utility for detailed comparison of quant.sf files.

**Usage**:
```bash
python3 compare_salmon_star.py salmon_quant.sf star_quant.sf [--verbose]
```

**Output**:
- Transcript counts
- Spearman/Pearson correlations
- Top differences (with --verbose)
- Pass/fail status

---

## Test Data

### Default Test Data

Located at `/mnt/pikachu/test-datasets-rnaseq/`:

| File | Description |
|------|-------------|
| `testdata/GSE110004/SRR6357070_*.fastq.gz` | Yeast PE reads |
| `reference/transcriptome.fasta.gz` | Yeast transcriptome |
| `reference/genes.gtf.gz` | Yeast GTF |

### Golden References

After running `regression_test.sh generate`, golden files are saved to:
```
tests/transcriptvb/golden/
├── vb_quant.sf      # VB mode output
├── em_quant.sf      # EM mode output
├── gc_quant.sf      # GC bias mode output
└── metadata.txt     # Test metadata
```

---

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `STAR_BIN` | Path to STAR binary | `/mnt/pikachu/STAR-Flex/source/STAR` |
| `SALMON_BIN` | Path to Salmon binary | `salmon` |
| `GENOME_DIR` | STAR genome index | `/tmp/star_vb_test/star_new_index` |

---

## Running All Tests

```bash
# Full test suite
./quick_test.sh && \
./edge_case_tests.sh && \
./regression_test.sh test && \
./salmon_parity_test.sh
```

---

## Adding New Tests

1. Add test case to appropriate script
2. Update `TESTING_TranscriptVB.md` documentation
3. If adding golden outputs, run `regression_test.sh update`

---

## Troubleshooting

### "Genome index not found"
Build genome index first:
```bash
STAR --runMode genomeGenerate \
    --genomeDir $GENOME_DIR \
    --genomeFastaFiles transcriptome.fasta \
    --sjdbGTFfile genes.gtf
```

### "Salmon not found"
Install Salmon or set `SALMON_BIN`:
```bash
conda install -c bioconda salmon
# or
export SALMON_BIN=/path/to/salmon
```

### "Test failed but output looks correct"
May need to update golden references:
```bash
./regression_test.sh update
```

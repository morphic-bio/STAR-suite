# TranscriptVB Test Suite

Test suite for validating STAR-Flex's integrated transcript quantification feature (TranscriptVB).

## Overview

This test suite provides:
- **Quick validation tests** - Basic functionality checks (~2 minutes)
- **Salmon parity tests** - Comparison with Salmon on same data
- **Regression tests** - Compare against golden outputs
- **Edge case tests** - Handle unusual inputs gracefully
- **PE bulk auto-detect validation** - Real PE regression checks for `--quantVBLibType A`

## Quick Start

```bash
# Run quick validation
./quick_test.sh

# Run Salmon parity test
./salmon_parity_test.sh

# Validate archived bad PE artifacts fail and fixed artifacts pass
../run_pe_bulk_regression_replay.sh

# Run same-BAM chr22 Salmon parity smoke
../run_transcriptvb_chr22_parity_smoke.sh

# Run public human male bulk PE smoke (after preparing fixture)
../run_public_bulk_pe_smoke.sh

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

### `../run_transcriptvb_chr22_parity_smoke.sh`

Runs the cached chr22 real-read fixture through STAR TranscriptVB auto-detect,
then runs Salmon on the exact STAR transcriptome BAM with fixed `ISR` and
`-p 1`. Salmon is pinned to one quantification thread for this smoke because its
alignment-mode online fragment-length distribution is thread-order sensitive on
the small chr22 fixture. STAR defaults to 32 threads to match the retained chr22
parity artifact; this is a semantic parity smoke, not a thread-determinism test.

**Checks**:
- detected library format is `ISR`
- `NumReads` Pearson is at least `0.99998`
- TPM Pearson is at least `0.99995`
- total assigned reads match within `0.05`
- half-L1 read-equivalent movement is at most `20`
- at most `10` transcripts differ by more than one read-equivalent

**Usage**:
```bash
../run_transcriptvb_chr22_parity_smoke.sh
```

**Runtime**: ~1-2 minutes with the cached fixture/reference.

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

### `validate_pe_autodetect_output.py`

Validates PE TranscriptVB auto-detect logs and optional Salmon parity metrics.

**Checks**:
- detected format is not `UNKNOWN`
- `dropped_incompat` stays at zero
- `dropped_missing_mate_fields` stays at zero
- `dropped_unknown_obs_fmt` stays at zero
- optional Salmon parity thresholds

**Usage**:
```bash
python3 validate_pe_autodetect_output.py --outdir /path/to/star_out
python3 validate_pe_autodetect_output.py \
    --outdir /path/to/star_out \
    --salmon-quant /path/to/salmon/quant.sf \
    --min-spearman-all 0.95 \
    --min-spearman-expressed 0.99 \
    --min-pearson-expressed 0.99
```

---

## Test Data

### Default Test Data

Located at `/mnt/pikachu/test-datasets-rnaseq/`:

| File | Description |
|------|-------------|
| `testdata/GSE110004/SRR6357070_*.fastq.gz` | Yeast PE reads |
| `reference/transcriptome.fasta.gz` | Yeast transcriptome |
| `reference/genes.gtf.gz` | Yeast GTF |

### Public PE Bulk Fixture

Local fixture download script:

```bash
scripts/download_public_bulk_fixture.sh
```

Default pinned source:

| Field | Value |
|------|-------|
| GEO series | `GSE88509` |
| GEO sample | `GSM2344101` |
| SRA run | `SRR4422207` |
| Sample | human male esophagus muscularis mucosa |
| Layout | paired-end |

Default local location:

```bash
/tmp/public_bulk_fixture_SRR4422207/
```

For the public PE smoke, the default STAR `--readFilesIn` order is `R2,R1`
via `PUBLIC_BULK_STAR_MATE_ORDER=R2R1`, matching the TranscriptVB PE
documentation for mate-order-sensitive libraries.

For the pinned `SRR4422207` public smoke fixture, the expected auto-detected
format is `ISF`.

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

# Testing STAR-Flex

This document describes how to test the STAR-Flex pipeline.

## Test Scripts

| Script | Purpose |
|--------|---------|
| `tests/run_flex_multisample_test.sh` | Full flex pipeline with multi-sample output |
| `tests/run_baseline_test.sh` | Baseline regression (flex disabled) |
| `tests/run_flex_inline_test.sh` | Basic inline hash test |

## Gold Standard

Baseline MEX output files are bundled in `tests/gold_standard/`:

```
tests/gold_standard/
├── raw/                    # Raw MEX output
│   ├── barcodes.tsv
│   ├── features.tsv
│   └── matrix.mtx
├── per_sample/             # Per-sample filtered MEX
│   ├── BC004/Gene/filtered/
│   ├── BC006/Gene/filtered/
│   ├── BC007/Gene/filtered/
│   ├── BC008/Gene/filtered/
│   └── flexfilter_summary.tsv
└── README.md
```

These files were generated from a 100K downsampled SC2300771 dataset using the reference STAR-Flex implementation.

## Running the Multi-Sample Test

The primary test validates the full flex pipeline:

```bash
cd /path/to/STAR-Flex
./tests/run_flex_multisample_test.sh
```

### Prerequisites

The test requires:

1. **STAR-Flex binary**: Built at `source/STAR`
2. **Reference genome**: Flex-enabled STAR index at `/storage/flex_filtered_reference/star_index`
3. **Input data**: 100K downsampled FASTQ files at `/storage/downsampled_100K/SC2300771/`
4. **Whitelist files**: CB whitelist, sample whitelist, probe list, sample probes

### Expected Output

Successful run produces:

```
=== Testing STAR-Flex with multi-sample output ===
...
✓ Solo.out/Gene/raw/matrix.mtx exists
  Entries: 137649
...
✓ per_sample directory exists
  ✓ BC004: matrix.mtx exists (7174 entries)
    ✓ Entry count matches gold standard
  ✓ BC006: matrix.mtx exists (7629 entries)
    ✓ Entry count matches gold standard
...
=== Full MEX comparison ===
✓ barcodes.tsv matches gold standard
✓ features.tsv matches gold standard
✓ matrix.mtx matches gold standard
```

### Validation Checks

The test compares:

1. **Raw MEX**: Entry counts and file contents vs gold standard
2. **Per-sample MEX**: Entry counts for each sample tag
3. **File identity**: Byte-for-byte comparison of barcodes, features, matrix

## Running Baseline Regression

To verify flex changes don't break standard Solo:

```bash
./tests/run_baseline_test.sh
```

This runs STAR with `--flex no` (default) and validates standard output.

## Manual Testing

### Minimal Flex Test

```bash
source/STAR \
  --genomeDir /path/to/flex_reference \
  --readFilesIn R2.fastq R1.fastq \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /path/to/whitelist.txt \
  --flex yes \
  --soloFlexExpectedCellsPerTag 1000 \
  --outFileNamePrefix test_output/
```

### Verify Output

```bash
# Check raw MEX was generated
ls -la test_output/Solo.out/Gene/raw/

# Check matrix dimensions
head -3 test_output/Solo.out/Gene/raw/matrix.mtx

# Check flex log messages
grep -i flex test_output/Log.out
```

## Troubleshooting

### "inlineHash_ is empty"

The inline hash received no reads. Check:
- CB whitelist path is correct
- Read files are in correct order (R2, R1 for Flex)
- Barcode parameters match your chemistry

### Per-sample output missing

Sample detection failed. Check:
- `--soloSampleWhitelist` path exists
- `--soloSampleProbes` and `--soloSampleProbeOffset` are set
- Sample tags in whitelist match expected sequences

### Entry count mismatch

Gene indexing or filtering differs. Check:
- `--soloProbeList` matches reference build
- FlexFilter parameters (`--soloFlexExpectedCellsPerTag`)

## Expected Performance

For the 100K downsampled test set:

| Metric | Expected |
|--------|----------|
| Raw barcodes | ~20,000 |
| Raw entries | ~137,649 |
| Cells per sample | 200-400 |
| Runtime | ~2-5 minutes |

## Continuous Integration

For CI pipelines, the test can be run with exit code checking:

```bash
#!/bin/bash
set -e
./tests/run_flex_multisample_test.sh 2>&1 | tee test.log

# Check for failures
if grep -q "⚠\|✗" test.log; then
    echo "Test failures detected"
    exit 1
fi
```


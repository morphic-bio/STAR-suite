# run_flexfilter_mex - Standalone MEX FlexFilter Wrapper

**Purpose**: Read composite MEX files (CB+TAG barcodes) and run libflexfilter pipeline to produce per-sample split MEX outputs.

---

## Overview

This tool provides a standalone entry point to the libflexfilter pipeline, consuming the composite MEX output from STAR's inline-hash MEX writer and producing filtered, per-sample MEX directories.

**Pipeline stages:**
1. **Load composite MEX** - Reads matrix.mtx, barcodes.tsv (CB16+TAG8), features.tsv
2. **Tag extraction** - Parses composite barcodes to extract tag sequences
3. **Dominance** - (Optional) Filter to dominant tag per cell
4. **OrdMag** - Order-of-magnitude filtering per tag
5. **EmptyDrops** - Multinomial-based cell calling per tag
6. **Per-sample output** - Write filtered MEX for each tag/sample

---

## Building

From the repository root:

```bash
# Build via main Makefile (requires STAR built first)
cd source && make flexfilter

# Or directly from this directory
cd tools/flexfilter && make
```

**Output**: `tools/flexfilter/run_flexfilter_mex` binary

---

## Usage

### Basic Example

```bash
./run_flexfilter_mex \
  --mex-dir /storage/100K/SC2300771/alignment/InlineHashDedup \
  --total-expected 12000 \
  --output-prefix /storage/100K/SC2300771/flexfilter_output
```

### With Allowed Tags Filter

```bash
./run_flexfilter_mex \
  --mex-dir /path/to/composite_mex \
  --allowed-tags /path/to/allowed_tags.tsv \
  --total-expected 12000 \
  --output-prefix /path/to/output
```

### With Custom Parameters

```bash
./run_flexfilter_mex \
  --mex-dir /path/to/composite_mex \
  --total-expected 12000 \
  --output-prefix /path/to/output \
  --ordmag-nsamples 5000 \
  --ordmag-target-pct 0.99 \
  --ed-lower 100 \
  --ed-niters 10000 \
  --ed-fdr-threshold 0.001 \
  --total-partitions 115000 \
  --recovery-factor 0.606
```

---

## Required Inputs

### --mex-dir

Path to directory containing composite MEX files:
- `matrix.mtx` - Matrix Market sparse matrix (1-based indices)
- `barcodes.tsv` - Composite barcodes (CB16+TAG8 format, 24 chars)
- `features.tsv` - Gene IDs/names (tab-separated)

**Composite barcode format**:
```
AAACCCAAGAAACACTACGTACGT  (CB16 + TAG8 = 24 chars)
AAACCCAAGAAACACTTGCATGCA
AAACCCAAGAAACACTGGCCGGCC
```

### --total-expected

Total expected cells across all tags/samples.

This is used to allocate expected cells per tag (proportional to UMI counts).

**Example**: If you expect 12,000 total cells across all samples, use `--total-expected 12000`

### --output-prefix

Output directory prefix for per-sample MEX outputs.

**Output structure** (sample labels are taken verbatim from the whitelist when provided; otherwise auto-derived; whitelist order is preserved):
```
<output-prefix>/
  SampleA/
    Gene/
      filtered/
        matrix.mtx
        barcodes.tsv
        features.tsv
        OrdMag/
          filter_summary.json
          ordmag_result.tsv
        EmptyDrops/
          emptydrops_result.tsv
  SampleB/
    ...
  filter_summary.json  (aggregate)
```

---

## Optional Inputs

### --allowed-tags

Path to TSV file with allowed tag sequences (one per line).

**Format**:
```
ACGTACGT
TGCATGCA
GGCCGGCC
```

Only composite barcodes with these tag sequences will be processed.

### --cells-per-tag-json

Path to JSON file specifying expected cells per tag (overrides auto-allocation).

**Format** (TODO - not yet implemented):
```json
{
  "ACGTACGT": 4000,
  "TGCATGCA": 5000,
  "GGCCGGCC": 3000
}
```

---

## Parameters

### OrdMag Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ordmag-nsamples` | 2500 | Number of random samples for OrdMag |
| `--ordmag-ambient-umi` | 10 | Max ambient UMI count |
| `--ordmag-target-pct` | 0.99 | Target percentile for threshold |
| `--dominance-ratio` | 0 | Dominance filter before split (0 = disabled) |
| `--tag-weights-umi` | (off) | Use UMI-weighted expected cell allocation (instead of ncells_simple) |

### EmptyDrops Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ed-lower` | 10 | Lower UMI bound for EmptyDrops |
| `--ed-niters` | 10000 | Monte Carlo iterations |
| `--ed-fdr-threshold` | 0.001 | FDR threshold for cell calling |

### Occupancy Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--total-partitions` | 115000 | Total partitions (for occupancy) |
| `--recovery-factor` | 0.606 | Recovery factor (1/1.65) |
| `--occupancy-percentile` | 0.999 | Percentile for occupancy cutoff |
| `--low-umi-threshold` | 0 | Low UMI threshold (0=auto) |

### Debug/Testing Flags

| Flag | Description |
|------|-------------|
| `--disable-occupancy` | Disable occupancy filter (for testing) |
| `--debug-tag-log` | Enable detailed per-tag logging |
| `--debug-output-dir <path>` | Directory for debug outputs |

---

## Output Format

### Per-Sample MEX Directories

For each tag/sample, the tool creates:

```
<output-prefix>/<sample_label>/Gene/filtered/
  matrix.mtx          # Filtered sparse matrix (1-based)
  barcodes.tsv        # Filtered composite barcodes
  features.tsv        # Gene IDs/names (same as input)
```

### OrdMag Outputs

```
<output-prefix>/<sample_label>/Gene/filtered/OrdMag/
  filter_summary.json # OrdMag filtering summary
  ordmag_result.tsv   # Per-barcode OrdMag scores
```

### EmptyDrops Outputs

```
<output-prefix>/<sample_label>/Gene/filtered/EmptyDrops/
  emptydrops_result.tsv # Per-barcode EmptyDrops results
```

### Aggregate Summary

```
<output-prefix>/filter_summary.json
```

Contains aggregate statistics across all tags/samples.

---

## Validation

### Check Output Structure

```bash
ls -R <output-prefix>/
```

Expected: Per-sample directories with MEX files; directory names match whitelist labels (or auto-derived labels if no whitelist is supplied).

### Verify Barcode Counts

```bash
wc -l <output-prefix>/<sample_label>/Gene/filtered/barcodes.tsv
```

Should match the number of passing cells for that sample.

### Compare to Expected

```bash
cat <output-prefix>/filter_summary.json
```

Check total passing cells vs `--total-expected`

---

## Example Workflow

### 1. Run STAR with Inline-Hash MEX Writer

```bash
STAR ... (your params) \
  --outFileNamePrefix /storage/output/
```

**Output**: `/storage/output/InlineHashDedup_*.{mtx,tsv}`

### 2. Run FlexFilter

```bash
./run_flexfilter_mex \
  --mex-dir /storage/output \
  --total-expected 12000 \
  --output-prefix /storage/output/flexfilter
```

**Output**: Per-sample MEX in `/storage/output/flexfilter/<sample_label>/`

### 3. Load into Seurat/Scanpy

```r
# R/Seurat (sample label preserved)
library(Seurat)
sample_a <- Read10X("/storage/output/flexfilter/SampleA/Gene/filtered/")
seurat_obj <- CreateSeuratObject(sample_a, ...)
```

```python
# Python/Scanpy
import scanpy as sc
adata = sc.read_10x_mtx("/storage/output/flexfilter/SampleA/Gene/filtered/")
```

---

## Troubleshooting

### Error: "Could not open MEX directory"

Check that `--mex-dir` points to a directory containing `matrix.mtx`, `barcodes.tsv`, `features.tsv`.

### Error: "No composite barcodes found"

Check barcode format in `barcodes.tsv` - should be 24 characters (CB16 + TAG8).

### Warning: "No allowed tags filter"

If `--allowed-tags` not specified, all tags in the composite MEX will be processed.

### Low cell counts

- Check `--total-expected` is reasonable
- Try adjusting `--ordmag-target-pct` or `--ed-fdr-threshold`
- Use `--debug-tag-log` to see per-stage filtering

---

## Notes

- **No Solo involvement**: This tool is completely independent of STAR's Solo pipeline
- **No replayer involvement**: Reads MEX directly, no keys.bin processing
- **Deterministic**: Given same inputs/params, produces identical outputs
- **Composite barcodes**: Handles CB16+TAG8 format (24 chars total)
- **Tag extraction**: Automatically extracts last 8 chars as tag sequence

---

## Status

✅ **Implemented**: Basic CLI framework, FlexFilter integration  
✅ **Tested**: Compiles with libflex  
⏳ **TODO**: JSON parsing for `--cells-per-tag-json`  
⏳ **TODO**: Add synthetic test case  

**Date**: 2025-11-23  
**Version**: 1.0

# Feature Barcode Tools (process_features vendor)

`process_features` is the canonical feature-processing implementation in this
repo. Both FASTQ feature extraction and MEX feature calling are built from the
same `libprocess_features` core.

## Build

```bash
make process-features-tools
```

Binaries are produced in `core/features/process_features/`:
- `assignBarcodes`
- `demux_fastq`
- `demux_bam`
- `call_features`

Compatibility alias:
```bash
make feature-barcodes-tools
```
This target is a compatibility alias and builds the same
`core/features/process_features` tools.

## Quick Usage (assignBarcodes)

```bash
core/features/process_features/assignBarcodes \
  --whitelist /path/to/whitelist.txt \
  --featurelist /path/to/feature_ref.csv \
  --directory /path/to/output \
  /path/to/fastq_dir
```

## Quick Usage (call_features from MEX)

```bash
core/features/process_features/call_features --gmm \
  /path/to/mex_dir \
  /path/to/output_dir
```

## Feature Offset Detection

By default, `assignBarcodes` and the `pf_api` library automatically detect the optimal feature offset from the `pattern` column in your feature reference CSV.

### How It Works

1. **Pattern Column**: If your feature CSV has a `pattern` column with `(BC)` markers (e.g., `NNNN(BC)NNNN`), the offset is extracted as the position of `(BC)`.

2. **Auto-Detection**: At startup, the tool scans all feature offsets:
   - If all features share the same offset → uses it as global offset (fast path)
   - If multiple offsets detected (>5% heterogeneity) → stops with an error

3. **User Override**: You can explicitly control offset behavior:
   ```bash
   # Use specific global offset (skips auto-detection)
   --feature_constant_offset 26
   
   # Force per-feature offsets (slower for large feature sets)
   --force-individual-offsets
   ```

### Error: Multiple Offsets Detected

If your feature reference has heterogeneous offsets, you'll see:

```
ERROR: Multiple feature offsets detected in pattern column.
       Dominant offset: 26 (used by 9500 features)
       Other offsets detected (threshold: 5% of dominant):
         offset 30: 500 features (5.3%)

To proceed, choose one of:
  1. --force-individual-offsets   Use per-feature offsets (slower for large feature sets)
  2. --feature_constant_offset 26  Use dominant offset globally (faster)
```

**Resolution:**
- Use `--force-individual-offsets` if features genuinely have different offsets (e.g., mixed assay types)
- Use `--feature_constant_offset N` to apply the dominant offset globally (faster for large feature sets)

### Error: Conflicting Flags

You cannot specify both `--feature_constant_offset` and `--force-individual-offsets`:

```
Error: Cannot specify both --force-individual-offsets and --feature_constant_offset.
       Use --force-individual-offsets for per-feature offsets from pattern column,
       or --feature_constant_offset N for a single global offset.
```

### pf_api Usage

The library API (`pf_api.h`) follows the same behavior:

```c
pf_config *config = pf_config_create();

// Option 1: Use auto-detect (default, no setter call needed)
// Option 2: Explicit global offset
pf_config_set_feature_offset(config, 26);
// Option 3: Per-feature offsets from pattern
pf_config_set_use_feature_offset_array(config, 1);

pf_context *ctx = pf_init(config);
// ...
pf_error err = pf_process_fastq_dir(ctx, fastq_dir, output_dir, &stats);
if (err == PF_ERR_MULTI_OFFSET_DETECTED) {
    fprintf(stderr, "%s\n", pf_get_error(ctx));
}
if (err == PF_ERR_OFFSET_CONFLICT) {
    fprintf(stderr, "%s\n", pf_get_error(ctx));
}
```

## MEX Stub Conversion

`assignBarcodes` outputs `features.txt` and `barcodes.txt`. To create 10x-style TSVs:

```bash
tools/feature_barcodes/assignbarcodes_mex_stub.py \
  --assign-out /path/to/output \
  --feature-csv /path/to/feature_ref.csv
```

This writes `features.tsv` and `barcodes.tsv` alongside the existing `matrix.mtx`.

## Notes
- Canonical implementation path: `core/features/process_features`.
- `core/features/feature_barcodes` is retained only as a compatibility path.

---

## CRISPR Feature Calling (CR-Compat Mode)

When running STAR with `--pfMultiConfig` and CRISPR Guide Capture features, STAR automatically runs GMM-based feature calling after EmptyDrops filtering.

### Parameter: `--crMinUmi N`

**Default:** 3 (general STAR-suite default)

Controls the minimum UMI threshold for feature calling. Adjust based on assay type:

| Assay Type | Recommended Value | Rationale |
|------------|-------------------|-----------|
| **CRISPR Guide Capture (general)** | 3 (default) | Baseline setting across perturb/lineage contexts |
| **A375 CR-parity fixture** | 10 (override) | Fixture-specific parity target |
| **Lineage Barcodes** | 2-3 | Stable features with minimal noise - lower threshold captures more signal |
| **FLEX Probes** | 10 | Similar characteristics to CRISPR guides |

### Example

```bash
# CRISPR Guide Capture (general default)
STAR ... --pfMultiConfig config.csv

# A375 CR-parity fixture (explicit override)
STAR ... --pfMultiConfig config.csv --crMinUmi 10

# Lineage Barcodes (lower threshold)
STAR ... --pfMultiConfig config.csv --crMinUmi 3
```

### Output

CR-compat mode produces `outs/crispr_analysis/`:
- `protospacer_calls_per_cell.csv` - Per-cell feature assignments
- `protospacer_calls_summary.csv` - Calling statistics
- `protospacer_umi_thresholds.csv` - GMM-derived UMI thresholds

See `tests/crispr_feature_calling_comparison_report.md` for validation details.

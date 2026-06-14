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

## ADT / Protein MEX (`--output-mode adt_mex`)

For Multiomics Suite protein quantification, `assignBarcodes` can emit a gzipped
10x MEX directory directly (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`)
with `Antibody Capture` feature rows and provenance sidecars. See
`docs/RUNBOOK_PROCESS_FEATURES_ADT_MEX.md`.

In **pf-multi** configs, ADT/protein is another feature-library arm (same path as
gRNA or LARRY): declare a non-GEX `feature_types` value such as `Antibody Capture`,
`ADT`, or `Protein`, provide `star_feature_ref`, and the permits-based runner
calls `assignBarcodes` with ADT MEX mode for that library only. Multiomics Suite
consumes the resulting `protein.mex_dir`; STAR-suite does not orchestrate CLR or
downstream packaging.

---

## CRISPR Feature Calling (CR-Compat Mode)

When running STAR with `--pfMultiConfig` and CRISPR Guide Capture features, STAR automatically runs GMM-based feature calling after EmptyDrops filtering.

### Parameter: `--crMinUmi N`

**Default:** 3 (general STAR Suite default)

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

### Planned: Ambient-FDR Guide Calling

The current CR-compatible caller should remain available for Cell Ranger parity.
For perturb-seq QC, STAR-suite should also support an ambient-FDR guide caller
that uses raw guide counts in non-cell barcodes as the noise floor and emits
tunable q-values for the final called-cell set:

```text
call at FDR alpha = guide_qvalue <= alpha
```

This mode should run in addition to compatibility calling, not replace it. Raw
non-cell barcodes should be used only to estimate ambient guide rates; q-values
and calls should be stored only for filtered cells. The intended output is a
sparse q-value matrix plus per-cell calls at a default FDR, so downstream QC and
MuData builders can expose an FDR slider without re-running guide assignment.
The implementation contract is tracked in
`docs/HANDOFF_GUIDE_AMBIENT_FDR_CALLER_20260613.md`.

---

## Per-Library NXT/TRU Chemistry Override

10x 3' HT experiments may use mixed NXT/TRU chemistries across libraries
(e.g., mRNA GEX with TRU barcodes, gRNA CRISPR Guide Capture with NXT
barcodes). STAR supports per-library chemistry specification via the
`star_chemistry` column in the pfMultiConfig `[libraries]` section.

### pfMultiConfig `star_chemistry` Column

| Value | Behavior |
|-------|----------|
| `TRU` | Explicit TRU chemistry; auto-detection skipped for this library |
| `NXT` | Explicit NXT chemistry; auto-detection skipped for this library |
| `auto` | Auto-detect chemistry from reads for this library |
| *(empty or absent)* | Inherit from `--crChemistry` flag (default: `auto`) |

**Precedence**: `star_chemistry` column > `--crChemistry` flag > auto-detect.

### Example

```csv
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry
/path/to/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU
/path/to/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT
/path/to/LARRY,DE_30KO,Custom,Custom,TRU

[feature]
ref,/path/to/feature_ref.csv
```

The `star_` prefix avoids collision with any future 10x-native column.
Cell Ranger silently ignores unknown columns, so the same config file
works with both STAR and Cell Ranger.

### Merge Normalization

At merge time, all feature library barcodes are normalized to the output
namespace (TRU by default, controlled by `--crOutputChemistry`). This
ensures the merged MEX has a consistent barcode namespace regardless of
per-library chemistry differences.

---

## Per-Library Feature Reference and Library ID

When multiple feature libraries of different types are present (e.g., gRNA +
lineage barcodes), each library can specify its own feature reference CSV and
a stable output identifier.

### pfMultiConfig Columns

| Column | Description |
|--------|-------------|
| `star_feature_ref` | Per-library feature reference CSV path. When set, the global `[feature] ref` is not used for this library, and type-based filtering is skipped. |
| `star_library_id` | Stable provenance/output key for the library. Auto-generated as `{sample}_{feature_types}_{index}` when absent. Must be unique across libraries. |

**Feature ref precedence**: `star_feature_ref` (per-library) > `--crFeatureRef` (CLI flag) > `[feature] ref` (config global).

### Data-Driven Feature Routing

Any non-GEX `feature_types` value in the config is automatically routed to
`assignBarcodes`. This means custom types like "Lineage", "Custom", or
"CRISPRa Activation" work without requiring hardcoded support. Known 10x
types (CRISPR Guide Capture, Antibody Capture, Multiplexing Capture) map to
their canonical `featureRefType`; unknown types use their `feature_types`
verbatim.

### Example

```csv
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
/path/to/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,gex_de
/path/to/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,/path/to/ref_grna.csv,grna_de
/path/to/LARRY,DE_30KO,Custom,Custom,TRU,/path/to/ref_larry.csv,larry_de
/path/to/ADT,DE_30KO,Protein,Protein,TRU,/path/to/ref_protein.csv,adt_de

[feature]
ref,/path/to/ref_grna.csv
```

In this example:
- The gRNA library uses its own `star_feature_ref` (no filtering needed).
- The LARRY lineage library uses a separate reference CSV with LARRY barcodes.
- The ADT/protein library is routed like any other feature library; pf-multi
  enables `assignBarcodes --output-mode adt_mex` for `Antibody Capture`, `ADT`,
  or `Protein` rows and writes protein MEX sidecars under
  `cr_assign/<feature_types>/<star_library_id>/`.
- The global `[feature] ref` is a fallback for libraries without `star_feature_ref`.

### Validation

- `star_feature_ref` must point to an existing file (hard error at parse time).
- `star_library_id` must be unique across all libraries (hard error on duplicates).
- Missing trailing fields are padded with empty strings (backward compatible).

---

## Per-Library Max Hamming Distance

When running multiple feature libraries with very different sequence lengths
(e.g., 8-nt gRNA guides vs 40-nt lineage barcodes), a single global Hamming
distance is suboptimal. `star_max_hamming` lets each library specify its own
max Hamming distance independently.

### pfMultiConfig Column

| Column | Description |
|--------|-------------|
| `star_max_hamming` | Per-library max Hamming distance override. When set, overrides the global `--crAssignMaxHamming` for this library only. Empty or absent = use global. |

**Precedence**: `star_max_hamming` (per-library) > `--crAssignMaxHamming` (CLI flag).

### Example

```csv
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id,star_max_hamming
/path/to/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,gex_de,
/path/to/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,/path/to/ref_grna.csv,grna_de,1
/path/to/LARRY,DE_30KO,Custom,Custom,TRU,/path/to/ref_larry.csv,larry_de,5
```

In this example:
- The gRNA library (29 guides, 8 nt) uses Hamming=1, appropriate for short sequences.
- The LARRY library (245K barcodes, 40 nt) uses Hamming=5, enabling deeper fuzzy matching.
- The GEX row has no `star_max_hamming`; GEX does not use feature assignment.

### Guidance

| Sequence length | Recommended max Hamming |
|-----------------|-------------------------|
| 8–10 nt | 1 |
| 15–20 nt | 2 |
| 30–40 nt | 3–5 |

Higher Hamming distances on short sequences risk spurious matches (e.g.,
Hamming=5 on 8-nt sequences allows 62.5% of positions to mismatch). The
prehash memory budget also scales with Hamming distance and library size,
so per-library control avoids wasting memory on unnecessary prehash tiers.

---

## Optimized Benchmark Parameters (CR-Compat)

When running STAR benchmarks in CR-compat mode (GEX + features), use these
parameters for optimal throughput and parity. **Do not hardcode thread
counts**; use the dynamic interface so threads are auto-sized.

### Required Threading Parameters

```bash
--runThreadN 32                    # (or nproc)
--dynamicThreadInterface 1         # parallel phases: pf-preload overlaps Solo
--crAssignConsumerThreads -1       # AUTO-SIZE from runThreadN (DO NOT hardcode 4)
--crAssignSearchThreads 1          # 1 search thread per consumer (prevents oversubscription)
```

**Why these values matter:**

| Parameter | Bad value | Effect | Correct |
|-----------|-----------|--------|---------|
| `crAssignConsumerThreads` | `4` (hardcoded) | Only 5 threads active during feature assignment (~15% CPU) | `-1` (auto) |
| `crAssignSearchThreads` | `4` | Oversubscription: 31×4=124 threads, thrashing | `1` |
| `dynamicThreadInterface` | `0` (off) | Solo and pf-preload run sequentially, not overlapped | `1` |

### Full Reference Command

```bash
/usr/bin/time -v STAR \
  --runThreadN 32 \
  --dynamicThreadInterface 1 \
  --genomeDir /path/to/genome \
  --readFilesIn R2_files R1_files \
  --readFilesCommand zcat \
  --outSAMtype None \
  --clipAdapterType CellRanger4 \
  --clip3pPolyG yes \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist /path/to/3M-february-2018_TRU.txt \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand Forward \
  --soloFeatures GeneFull \
  --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes \
  --pfMultiConfig /path/to/config.csv \
  --crChemistry TRU \
  --crOutputChemistry TRU \
  --crWhitelist /path/to/3M-february-2018_TRU.txt \
  --crMinUmi 3 \
  --crAssignMaxHamming 1 \
  --crAssignFeatureOffset 0 \
  --crAssignLimitSearch -1 \
  --crAssignMinCounts 0 \
  --crAssignMaxBarcodeMismatches 5 \
  --crAssignFeatureN 0 \
  --crAssignBarcodeN 1 \
  --crAssignConsumerThreads -1 \
  --crAssignSearchThreads 1 \
  --crFeatureRef /path/to/feature_ref.csv
```

### Notes

- **No BAM for benchmarks**: Use `--outSAMtype None` unless BAM output is
  needed. Saves significant I/O and disk.
- **GeneFull only**: Skip `Gene` and `Velocyto` unless specifically required.
  Each adds a full pass over the read array.
- **Poly-G trimming**: Always `--clip3pPolyG yes` for NovaSeq/NextSeq data
  (auto-detected in CellRanger4 mode, but explicit is safer).
- **Chemistry**: Set `--crChemistry TRU` or `NXT` explicitly when known.
  Auto-detect can misclassify certain samples (see AALG1 autodetect bug).
- **Per-library Hamming**: Use `star_max_hamming` column when mixing short
  guides (h=1) with long barcodes (h=5). See section above.

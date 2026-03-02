# Handoff: CR vs STAR End-to-End Comparison — MSK 30polyKO (LARRY + PolyIII gRNA)

Date: 2026-02-26
Branch: `CR-Larry-perturb` (off `master` at `d230e0c`)
Prior: `docs/HANDOFF_MULTI_FEATURE_PHASE5_VALIDATION_20260226.md`

## Objective

Run Cell Ranger 9.0.1 `multi` and STAR (with multi-feature libraries) on
the same MSK 30polyKO perturb-seq dataset and produce a quantitative parity
report covering: cell barcode sets, GEX gene counts, feature (gRNA + LARRY)
assignment counts, and feature calling results.

This is the first real-data comparison of STAR's multi-feature pipeline
against Cell Ranger on a dataset with **three library types** (GEX, PolyIII
gRNA, LARRY lineage barcodes).

## Dataset: MSK 30polyKO ES/DE

| Component | Location | Details |
|-----------|----------|---------|
| FASTQ root | `/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM` | 352 FASTQ files total |
| mRNA (GEX) | `mRNA_DE_30KO_XM_IGO_16692_*` | 3 samples (S25, S36, S2), 4 lanes each = 28 R1 files |
| PolyIII (gRNA) | `PolyIII_DE_30KO_XM_IGO_16692_C_7_S55_*` | 1 sample, 4 lanes = 8 R1 files |
| LARRY (lineage) | `LARRY_DE_30KO_XM_IGO_16692_C_3_S35_*` | 1 sample, 4 lanes = 8 R1 files |
| gRNA feature ref | `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv` | 29 features, `Custom` type, pattern `(BC)GCTTTAAGGCCGGTCCTAGC` |
| LARRY feature ref | `/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv` | 245,979 features, `Custom` type, 40-nt sequences |
| Genome reference | `/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar` | GRCh38 2024-A, gencode v44, mkref 8.0.0 |
| Cell barcode whitelist | `/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt` | 3M-feb-2018 TRU |
| Cell Ranger | `/usr/local/bin/cellranger` | v9.0.1 |
| STAR | `core/legacy/source/STAR` | `multi-feature` squash merge `d230e0c` |

### Important: Feature Ref Type Compatibility

The MSK feature ref CSVs use `feature_type=Custom` for **both** gRNA and
LARRY. Cell Ranger routes `Custom` features through its generic feature
barcode pipeline (not the CRISPR Guide Capture calling pipeline).

**Decision needed before running CR**: Should the gRNA ref be rewritten to use
`feature_type=CRISPR Guide Capture` for Cell Ranger, so CR applies its
GMM-based guide calling? If so, create a separate CR-specific gRNA ref.
If not, both tools will treat gRNA as generic custom features (no guide
calling on either side), which makes the comparison simpler but less
informative for CRISPR calling parity.

**Recommendation**: Run **two CR configs** — one with `Custom` (apples-to-apples
with STAR's current ref) and one with `CRISPR Guide Capture` (to test whether
CR's guide calling adds value). STAR's pf-multi pipeline already supports both
routing modes via `feature_types` in the config.

## Scope

### Phase A: Single-Sample Downsampled Smoke (Fast, <30 min)

Pick one mRNA sample (S25), use lane L001 only, downsample to 100K reads per
library. Purpose: validate configs, catch format/path errors, produce a small
baseline for script debugging.

1. Create downsampled fixture (or reuse `create_fixture_downsampled.sh` with
   100K reads).
2. Write CR `multi` config CSV.
3. Write STAR `pfMultiConfig` CSV.
4. Run both tools.
5. Verify both complete without error.

### Phase B: Full Single-Sample Run (S25, all lanes)

Run on the full S25 sample (all 4 lanes × 3 library types). This is the
primary comparison dataset.

1. **Cell Ranger `multi`** — Run with `--localcores 16 --localmem 64`.
   Expected runtime: 30–60 min.
2. **STAR multi-feature** — Run with `--runThreadN 16`.
   Expected runtime: 15–30 min.
3. Both outputs should be in timestamped directories under
   `/storage/MSK-perturb-comparison/`.

### Phase C: Parity Comparison

Extract and compare the following from both outputs:

#### C1: Cell Barcode Parity
- Extract filtered barcode sets from both tools.
- STAR: `Solo.out/GeneFull/filtered/barcodes.tsv`
- CR: `outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz` (raw) and
  `outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz` (filtered)
- Compute Jaccard index, overlap count, STAR-only count, CR-only count.

#### C2: GEX Gene Count Parity
- For the intersection of filtered barcodes:
  - Read GEX count matrices from both tools.
  - STAR: `Solo.out/GeneFull/filtered/matrix.mtx`
  - CR: `outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix/`
    (filter to `Gene Expression` features)
  - Compute per-gene Pearson correlation (summed across cells).
  - Compute per-cell total UMI correlation.
  - Report gene-level and cell-level Pearson r.

#### C3: Feature Assignment Parity (gRNA)
- STAR: `cr_assign/CRISPR_Guide_Capture/grna_de/*/matrix.mtx` (or
  `cr_assign/Custom/grna_de/*/matrix.mtx` depending on `feature_types`)
- CR: extract the CRISPR Guide Capture (or Custom) features from CR's
  combined feature-barcode matrix.
- Compare: total UMIs per feature, total UMIs per cell, Pearson on the
  intersection barcode set.

#### C4: Feature Assignment Parity (LARRY)
- Same as C3 but for the LARRY lineage barcodes.
- Note: LARRY has 245K features; expect many zero-count features in both
  tools. Focus on features with ≥1 UMI.

#### C5: Feature Calling Parity (optional, only if CR uses CRISPR Guide Capture)
- CR's guide calling output: `outs/multi/count/crispr_analysis/`
- STAR's CRISPR calling: `cr_assign/.../crispr_analysis/` (if enabled)
- Compare guide assignments per cell.

### Phase D: Report

Produce a markdown summary with:
- Run parameters (versions, threads, wall time, peak RAM)
- Cell barcode Jaccard and overlap table
- GEX Pearson (gene-level, cell-level)
- Feature assignment Pearson per library
- Key differences and potential root causes
- Recommendations for follow-up

## Implementation Guide

### Step 1: Create CR Multi Config

Cell Ranger `multi` requires a specific CSV format. Template:

```
[gene-expression]
reference,/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar
create-bam,true

[feature]
reference,/path/to/combined_feature_ref.csv

[libraries]
fastqs,sample,library_type
/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM,mRNA_DE_30KO_XM_IGO_16692_C_9_S25,Gene Expression
/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM,PolyIII_DE_30KO_XM_IGO_16692_C_7_S55,Custom
/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM,LARRY_DE_30KO_XM_IGO_16692_C_3_S35,Custom
```

**Critical**: CR's `[feature] reference` must be a **single combined CSV**
containing both gRNA and LARRY features. Concatenate the two refs (keeping
one header). CR does not support per-library feature refs.

```bash
head -1 /mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv > /tmp/combined_feature_ref.csv
tail -n +2 /mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv >> /tmp/combined_feature_ref.csv
tail -n +2 /mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv >> /tmp/combined_feature_ref.csv
```

### Step 2: Create STAR Multi Config

```
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM,mRNA_DE_30KO_XM_IGO_16692_C_9_S25,Gene Expression,Gene Expression,TRU,,gex_s25
/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM,PolyIII_DE_30KO_XM_IGO_16692_C_7_S55,Custom,Custom,NXT,/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv,grna_s55
/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM,LARRY_DE_30KO_XM_IGO_16692_C_3_S35,Custom,Custom,TRU,/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv,larry_s35
```

STAR supports per-library `star_feature_ref`, so no combined ref is needed.

### Step 3: Run Commands

**Cell Ranger** (from a clean working directory):
```bash
OUT=/storage/MSK-perturb-comparison/cr_multi_s25_$(date +%Y%m%d_%H%M%S)
mkdir -p "$OUT" && cd "$OUT"
cellranger multi \
    --id cr_msk_s25 \
    --csv /path/to/cr_multi_config.csv \
    --localcores 16 \
    --localmem 64
```

**STAR**:
```bash
OUT=/storage/MSK-perturb-comparison/star_multi_s25_$(date +%Y%m%d_%H%M%S)
mkdir -p "$OUT"
STAR \
    --runMode alignReads \
    --runThreadN 16 \
    --genomeDir /storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star \
    --readFilesIn <GEX_R2> <GEX_R1> \
    --readFilesCommand zcat \
    --pfMultiConfig /path/to/star_multi_config.csv \
    --defaultCrCompat yes \
    --crChemistry auto \
    --outFileNamePrefix "$OUT/" \
    --outSAMtype BAM Unsorted \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist /storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures Gene GeneFull \
    --soloBarcodeReadLength 0
```

### Step 4: Comparison Script

Write `comparisons/msk_cr_vs_star_parity.py` (or `.sh`) that:
1. Reads both output directories
2. Extracts barcode sets, GEX matrices, feature matrices
3. Computes Jaccard, Pearson, per-feature UMI totals
4. Outputs a markdown parity report

Model after the existing `scripts/ucsf_parity/` infrastructure (see
`build_star_exact_vs_cr.py` for MEX reading and comparison patterns).

## Expected Challenges

1. **LARRY ref size (245K features)**: CR and STAR will both be slow on
   feature assignment for LARRY. Downsampled smoke (Phase A) should use
   the 500-feature subset from `create_fixture_downsampled.sh`. Full run
   (Phase B) must use the complete ref.

2. **Chemistry mismatch**: The PolyIII gRNA library uses NXT chemistry.
   STAR handles this via `star_chemistry=NXT`. CR infers chemistry
   automatically from the FASTQ header; verify CR auto-detects correctly.

3. **FASTQ sample naming**: CR uses the `sample` column in `[libraries]`
   to match FASTQ filenames (the `_S<N>_` token). Make sure the sample
   names in the CR config match exactly (e.g., `PolyIII_DE_30KO_XM_IGO_16692_C_7_S55`).
   If CR does not find matching FASTQs, it will fail with a confusing error.

4. **feature_type routing**: Both refs currently use `feature_type=Custom`.
   CR's `library_type` in `[libraries]` must match what the ref says
   (`Custom`). Do not mix `CRISPR Guide Capture` in the libraries CSV with
   `Custom` in the feature ref, or CR will not match features to libraries.

5. **Barcode namespace**: STAR's mRNA and LARRY are TRU; PolyIII is NXT.
   STAR normalizes at the merge boundary. CR handles this internally.
   Verify barcode sets are comparable (same namespace in output).

## Output Structure

```
/storage/MSK-perturb-comparison/
├── cr_multi_s25_YYYYMMDD_HHMMSS/
│   └── cr_msk_s25/
│       └── outs/
│           ├── multi/count/raw_feature_bc_matrix/
│           ├── per_sample_outs/*/count/sample_filtered_feature_bc_matrix/
│           └── ...
├── star_multi_s25_YYYYMMDD_HHMMSS/
│   ├── Solo.out/GeneFull/filtered/
│   ├── cr_assign/Custom/grna_s55/
│   ├── cr_assign/Custom/larry_s35/
│   └── Log.out
├── parity_report_YYYYMMDD.md
└── configs/
    ├── cr_multi_config.csv
    └── star_multi_config.csv
```

## Exit Criteria

1. Both tools complete on the full S25 sample without error.
2. Cell barcode Jaccard ≥ 0.90 (expect ≥ 0.95 based on UCSF parity).
3. GEX gene-level Pearson ≥ 0.99.
4. GEX cell-level total-UMI Pearson ≥ 0.99.
5. gRNA feature-level Pearson ≥ 0.95 (lower bar: different matching engines).
6. LARRY feature-level Pearson ≥ 0.90 (lower bar: 245K features, many sparse).
7. Parity report produced and committed to `comparisons/`.

## Existing Infrastructure to Reuse

| Tool/Script | Use |
|-------------|-----|
| `scripts/customize_cellranger_multi.sh` | Generate/customize CR multi config |
| `scripts/ucsf_parity/build_star_exact_vs_cr.py` | Pattern for MEX reading + Pearson computation |
| `scripts/ucsf_parity/inspect_barcode_feature_totals.py` | Barcode/feature total inspection |
| `tests/multi_feature/create_fixture_downsampled.sh` | Downsampled fixture for smoke testing |
| `tests/multi_feature/run_phase5_validation.sh` | Reference for STAR multi-feature invocation |

## Non-Goals

- Do not modify STAR or CR source code in this phase.
- Do not attempt multi-sample (S25+S36+S2) comparison yet; start with S25.
- Do not tune parameters for parity; use defaults on both sides.
- Do not compare BAM outputs (too large, low signal for this comparison).

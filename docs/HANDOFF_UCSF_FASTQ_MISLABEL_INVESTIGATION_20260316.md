# UCSF Perturb-seq FASTQ Mislabelling Investigation

**Date**: 2026-03-16
**Dataset**: UCSF Perturb-seq (16 samples x 2 library types, CRISPRa screen)
**Branch**: CR-Larry-perturb

## 1. Background

The UCSF Perturb-seq dataset was provided by a collaborator as 32 sample
directories (16 base samples x 2 indices: AALG1 and AALG2), organized into two
top-level directories on `/mnt/pikachu/ucsf-perturb-seq/`:

```
GEX/
  EBs1_1_AALG1/  EBs1_1_AALG2/  ...  iPSC2_3_AALG1/  iPSC2_3_AALG2/
guides/
  EBs1_1_AALG1/  EBs1_1_AALG2/  ...  iPSC2_3_AALG1/  iPSC2_3_AALG2/
```

The assumption was that `GEX/` contained gene expression data and `guides/`
contained CRISPR guide capture data, and that AALG1 and AALG2 were separate
biological samples or experimental conditions. **Both assumptions turned out to
be wrong.**

## 2. How the Problem Was Discovered

### 2.1 Initial benchmark with iPSC2_1_AALG1

A STAR benchmark was run on `iPSC2_1_AALG1`, forcing TRU chemistry
(`--crChemistry TRU`). The GEX results looked reasonable (~11,000 filtered
cells) and closely matched CellRanger. However, when we checked the **feature
(CRISPR guide capture) results**, a major discrepancy appeared:

| Metric | STAR | CellRanger |
|---|---|---|
| Cells with guide calls | 595 | 1 |
| Guide capture UMIs | Very low | Near zero |

CellRanger found effectively **zero guide reads** in the AALG1 "guide" library.
STAR found a small number, but far fewer than expected for a Perturb-seq
experiment. This prompted a deeper investigation.

### 2.2 Read content analysis

We sampled R2 reads from the AALG1 "guide" FASTQs looking for the guide RNA
scaffold pattern `GTTT[N]AGAGCTAAGC`. Result: **0% scaffold hits**. The reads
were 91bp genomic/transcriptomic sequences — gene expression data, not guide
capture.

Meanwhile, the AALG2 "GEX" directory R2 reads showed **~80% scaffold hits** —
these were actually guide capture reads mislabeled as GEX.

### 2.3 Systematic survey of all 48 directories

A comprehensive survey was run across all 24 samples x 2 library types
(48 directories total), checking scaffold content in R2 reads:

```
AALG1 "GEX" directories  → 0% scaffold  → Actual GEX data       ✓ Correct
AALG1 "guide" directories → 0% scaffold  → Actual GEX data       ✗ MISLABELED
AALG2 "GEX" directories  → ~80% scaffold → Actual guide data     ✗ MISLABELED
AALG2 "guide" directories → ~80% scaffold → Actual guide data     ✓ Correct
```

The pattern was **100% consistent** across all 16 base samples: the directories
were organized by **flowcell**, not by library content.

## 3. What Was Tried (and Why Earlier Approaches Failed)

### 3.1 Raw R1 barcode comparison (failed)

First attempt: extract the top-50 most frequent 16bp R1 barcodes from each
FASTQ and compare across AALG1/AALG2 pairs by raw sequence identity.

**Result**: Zero overlap between AALG1 and AALG2 for the same base sample.

**Why it failed**: AALG1 and AALG2 use different 10x Chromium barcode chemistries
(TRU vs NXT). The raw 16bp sequences are completely different even for the same
physical cell. Without chemistry detection and translation, raw comparison is
meaningless.

### 3.2 Processed barcodes.tsv comparison (partially worked)

Second attempt: use `barcodes.tsv` files from previous STAR batch runs (which
had already been error-corrected and matched to whitelists).

**Result**: Some overlap appeared, but inconsistent. Same-name Jaccard scores
ranged from 0.003 to 0.04 — better than raw comparison but not conclusive.

**Why it was inconsistent**: The previous batch run had applied the **wrong
whitelist** to AALG1 samples. It used the NXT whitelist for everything, but
AALG1 data is TRU chemistry. This caused poor cell calling for AALG1 samples
(~30 cells vs ~5000 expected), making the `barcodes.tsv` files unreliable for
those runs.

### 3.3 Chemistry-aware barcode matching (definitive solution)

Final approach: detect chemistry from raw R1 data, then translate to a common
namespace before matching.

**Method**:
1. Extract 500K raw R1 barcodes per FASTQ directory
2. Match against both TRU and NXT whitelists to detect chemistry
3. Translate all barcodes to TRU namespace (NXT→TRU via the
   `3M-february-2018_NXT.txt` translation file)
4. Take top-500 barcodes by count for each FASTQ set
5. Compute Jaccard similarity between all pairs

**Results**:

Chemistry detection was unambiguous:
```
All AALG1 files: TRU  (85-90% TRU whitelist hit rate, <1% NXT)
All AALG2 files: NXT  (85-90% NXT whitelist hit rate, <1% TRU)
```

After TRU normalization, within-sample 4-way Jaccard analysis:
```
Base         A1gex-A1gui  A1gex-A2gui  A1gex-A2gex  A2gui-A2gex
EBs1_1        0.0839       0.0804       0.0013       0.0693
iPSC2_1       0.0601       0.0699       0.0010       0.0686
...
```

Key observations:
- `A1gex ↔ A1gui` (Jaccard ~0.06–0.12): **Same cells** — both are GEX from
  the same TRU library
- `A1gex ↔ A2gui` (Jaccard ~0.06–0.12): **Same cells** — GEX (TRU) matches
  guide capture (NXT) after translation
- `A1gex ↔ A2gex` (Jaccard ~0.001): **Different content** — GEX vs guide data
  have different barcode frequency distributions

This confirmed that **AALG1 and AALG2 are the same cells** sequenced with
different chemistry kits for different library types.

## 4. Root Cause

The FASTQ files were organized by **sequencing flowcell**, not by library
content:

| Flowcell | Lanes | Original directory | AALG1 content | AALG2 content |
|---|---|---|---|---|
| 22MW3JLT3 | 8 (L001–L008) | `GEX/` | GEX reads (TRU) | Guide reads (NXT) |
| 22MTK2LT3 | 2 (L001–L002) | `guides/` | GEX reads (TRU) | Guide reads (NXT) |

Both flowcells carried **both sample indices** (AALG1 and AALG2). The
flowcell IDs are only visible in the FASTQ read headers
(`@A00519:111:22MW3JLT3:...` vs `@A00519:112:22MTK2LT3:...`), not in the
filenames. Whoever organized the directories assumed one flowcell = GEX and
the other = guides, which was incorrect.

The depth difference between GEX and guide files (~4x) reflects the library
pooling ratio, not flowcell allocation. GEX files are ~2.6 GB/lane while guide
files are ~640 MB/lane on the same flowcell.

## 5. Corrected Organization

Created `/mnt/pikachu/ucsf-perturb-seq-corrected/` with symlinks organizing
data by **actual library content**, preserving original filenames for
provenance:

```
ucsf-perturb-seq-corrected/
  <sample>/
    GEX/          ← ALL AALG1 files (TRU chemistry)
                     L001–L008: original GEX dir (flowcell 22MW3JLT3)
                     L009–L010: original guide dir (flowcell 22MTK2LT3, renumbered)
    guides/       ← ALL AALG2 files (NXT chemistry)
                     L001–L002: original guide dir (flowcell 22MTK2LT3)
                     L009–L016: original GEX dir (flowcell 22MW3JLT3, renumbered)
    pf_multi_config.csv
```

Lane numbers are renumbered (L009+) only where needed to avoid filename
collisions. Original sample names, S-indices, and AALG1/AALG2 labels are
preserved exactly as received from the collaborator.

Each `pf_multi_config.csv`:
```
[libraries]
fastqs,sample,library_type,feature_types
.../EBs2_2/GEX,EBs2_2_AALG1,Gene Expression,Gene Expression
.../EBs2_2/guides,EBs2_2_AALG2,CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv
```

## 6. Actual Dataset Summary

**16 biological samples** (not 32), each with paired GEX + guide capture:

| Group | Samples | Description |
|---|---|---|
| EBs1 | EBs1_1 through EBs1_5 | Embryoid bodies batch 1 (5 samples) |
| EBs2 | EBs2_1 through EBs2_5 | Embryoid bodies batch 2 (5 samples) |
| iPSC1 | iPSC1_1 through iPSC1_3 | iPSC batch 1 (3 samples) |
| iPSC2 | iPSC2_1 through iPSC2_3 | iPSC batch 2 (3 samples) |

Chemistry mapping:
- **AALG1** = TRU (v3) chemistry → Gene Expression library
- **AALG2** = NXT (v4) chemistry → CRISPR Guide Capture library (549 guides,
  pattern `(BC)GTTTNAGAGCTAAGC`)

Feature reference: `cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv`

## 7. Key Scripts

| Script | Purpose |
|---|---|
| `/tmp/detect_and_crossmatch.py` | Chemistry detection + TRU-normalized barcode cross-matching |
| `/tmp/survey_all_fastqs2.sh` | Scaffold pattern survey across all 48 directories |
| `/tmp/organize_corrected_v2.sh` | Creates the corrected symlink directory structure |

## 8. Implications for Processing

1. **Mixed chemistry**: GEX uses TRU, guides use NXT. Both STAR and CellRanger
   must handle per-library chemistry detection and barcode translation.
2. **Barcode namespace**: All output should use TRU namespace
   (`--crOutputChemistry TRU`) so GEX and guide barcodes are in the same space.
3. **Previous batch runs**: Any runs using the original directory organization
   with AALG1 guide data are invalid (they processed GEX data as guides and
   vice versa). Must be re-run with corrected organization.
4. **No data loss**: All files from both flowcells contribute additional
   sequencing depth. Each sample gets 10 lanes of GEX and 10 lanes of guides.

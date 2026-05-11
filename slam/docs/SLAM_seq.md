# SLAM-seq Quantification in STAR-SLAM

STAR-SLAM includes a fully integrated SLAM-seq quantification module ("STAR-Slam") that allows for gene-level quantification of metabolic labeling experiments (e.g., s4U labeling) directly within the aligner. This implementation provides parity with GRAND-SLAM while adding performance improvements and new features like variance-based auto-trimming and comprehensive QC reports.

## Overview

SLAM-seq (Thiol-linked Alkylation for the Metabolic sequencing of RNA) uses 4-thiouridine (4sU) incorporation into nascent RNA. Upon reverse transcription, 4sU interacts with Cytosine, leading to T-to-C transitions (in the sense strand) during sequencing. By quantifying these T->C conversions, one can distinguish newly synthesized RNA from pre-existing RNA.

STAR-Slam performs the following steps:
1.  **Alignment**: Maps reads to the reference genome.
2.  **SNP Handling**: Masks native SNPs (T->C) to avoid false positives, either using an input VCF/BED or internal detection.
3.  **Trimming**: Identifies and excludes read ends with high error rates (auto-trimming).
4.  **Quantification**: Counts T->C conversions and estimates the New-to-Total Ratio (NTR) for each gene using a binomial or EM model.

## New Features in STAR-SLAM

*   **GRAND-SLAM Compatible Output**: Generates `<prefix>SlamQuant.grandslam.tsv` with the exact column schema required by downstream tools expecting GRAND-SLAM output.
*   **Variance-based Auto-trimming**: Automatically detects optimal 5' and 3' trim points by analyzing the variance of T->C rates across read positions.
*   **Self-contained QC Reports**: Generates interactive HTML reports visualizing T->C conversion rates, error profiles, and trimming decisions.

## Usage

### Basic Command

```bash
STAR \
  --runMode alignReads \
  --genomeDir /path/to/index \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --slamQuantMode 1 \
  --slamGrandSlamOut 1 \
  --autoTrim variance \
  --slamQcReport qc_output/my_sample \
  --outFileNamePrefix output/
```

### Batch Layout + Blank-First Error Rate

For multi-fastq SLAM runs where a **blank** (e.g., no4sU) should set the background error
rate, STAR can derive a per-sample output layout automatically while keeping a common root.
This is useful for large series (0h/6h/24h) where you want a consistent directory structure.

```bash
STAR \
  --runMode alignReads \
  --genomeDir /path/to/index \
  --readFilesIn /path/to/blank_no4su_R1.fastq.gz \
  --readFilesCommand zcat \
  --slamQuantMode 1 \
  --slamGrandSlamOut 1 \
  --slamSnpBed /path/to/snps.bed \
  --trimScope first \
  --trimSource /path/to/blank_no4su_R1.fastq.gz \
  --slamErrorRateFromBlank 1 \
  --outFileNamePrefix /path/to/output_root/ \
  --outFileNamePrefixAuto 1
```

**What `--outFileNamePrefixAuto` does**

- Derives the sample name from the first read file:
  - Strips `.fastq(.gz)`/`.fq(.gz)` and trailing `_R1/_R2` tokens.
- Uses `--outFileNamePrefix` as a **root** and creates subdirectories:
  - `alignments/`, `qc/`, `counts/`, `y_removed/`
- Auto‑routes outputs:
  - SLAM QC (`slamQcJson/Html`) → `qc/`
  - SlamQuant tables → `counts/`
  - `_Y.bam` / `_noY.bam` → `y_removed/`

### SLAM Batch Mode (--slamBatchMode)

For processing multiple SLAM-seq FASTQs in a single STAR invocation with shared genome loading:

```bash
STAR \
  --runThreadN 24 \
  --genomeDir /path/to/index \
  --readFilesIn /path/to/blank.fastq.gz,/path/to/0h.fastq.gz,/path/to/6h.fastq.gz,/path/to/24h.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix /path/to/output/ \
  --outFileNamePrefixAuto 1 \
  --slamQuantMode 1 \
  --slamGrandSlamOut 1 \
  --slamSnpMaskIn /path/to/snps.bed.gz \
  --trimScope first \
  --slamErrorRateFromBlank 1 \
  --slamBatchMode 1
```

**What `--slamBatchMode 1` does**:

- Enables batch processing of multiple FASTQs with a single genome load
- When combined with `--slamErrorRateFromBlank 1`, the **first FASTQ** (blank) is used to compute the background T→C error rate, which is then applied to all subsequent samples
- Requires `--slamQuantMode 1`
- Incompatible with Solo/single-cell processing (use single-sample mode for scRNA-seq)
- Recommended with `--outFileNamePrefixAuto 1` for per-sample output routing

**V1 Limitations**:
- Per-sample outputs are separated by prefix/dirs; Log.out remains a combined run log
- For paired-end input, provide **two comma‑separated lists** (R1 and R2):
  ```
  --readFilesIn R1_1.fastq.gz,R1_2.fastq.gz R2_1.fastq.gz,R2_2.fastq.gz
  ```

### Key Parameters

#### Quantification Control
*   `--slamQuantMode 1`: Enables SLAM-seq quantification.
*   `--slamGrandSlamOut 1`: (Default: 1) Enables generation of the GRAND-SLAM compatible output file `<prefix>SlamQuant.grandslam.tsv`.
*   `--slamOutFile`: (Optional) Custom name for the standard STAR-style output (default: `<prefix>SlamQuant.out`).
*   `--slamCbOut 1`: Enables model-ready count-binomial output `<prefix>SlamQuant.cB.tsv`.
*   `--slamCbOutFile`: (Optional) Custom cB output path. A `.gz` suffix writes gzip-compressed TSV.
*   `--slamCbFormat star|ezbakr`: Selects the compact STAR-SLAM cB schema or an EZbakR-like schema.
*   `--slamErrorRateFromBlank 1`: Derive `slamErrorRate` from the detection pass (useful when the first file is a blank).
*   `--trimScope first` + `--trimSource <blank>`: Use the blank to define trimming parameters that apply to all files.

#### Auto-Trimming (`--autoTrim variance`)
Reads often exhibit higher error rates or chemical artifacts at the 5' and 3' ends, which can confound T->C conversion detection. STAR-SLAM introduces a robust **variance-based** auto-trimming method.

*   **Theory**: We calculate the standard deviation of T->C transition rates at each read position across a sample of reads. Artifact-prone ends typically show high variance, while the "clean" middle section shows low, stable variance.
*   **Method**: A segmented regression (2, 3, or 4 segments) is fitted to the smoothed standard deviation curve. The "knees" or breakpoints of this regression define the start and end of the reliable region.
*   **Parameters**:
    *   `--autoTrim variance`: Enable variance-based trimming.
    *   `--trimScope`:
        *   `first`: (Default) Compute trims from the first file and apply to all (efficient for similar libraries).
        *   `per-file`: Compute and apply trims independently for each input file (requires distinct file markers).

#### Quality Control (`--slamQcReport`)
Generates comprehensive QC reports to validate the experiment and trimming performance.

*   `--slamQcReport <prefix>`: Generates:
    *   `<prefix>.slam_qc.json`: Raw QC metrics.
    *   `<prefix>.slam_qc.html`: Interactive self-contained HTML report.
*   **Report Contents**:
    *   **T->C Rate**: T->C mismatch rate by read position.
    *   **T->C Stdev & Fit**: Visualization of the variance curve and the segmented regression lines used for trimming.
    *   **Trim Overlay**: Visual indication of the chosen 5' and 3' trim points.

#### SNP Handling
To avoid counting native SNPs as metabolic conversions, you must mask them.
*   `--slamSnpBed <file.bed>`: Provide a BED file of known SNPs (recommended if available).
*   `--slamSnpDetect 1`: Enable internal SNP detection (if no BED provided). Estimates a mismatch threshold to identify likely SNPs.

## Output Files

1.  **`<prefix>SlamQuant.grandslam.tsv`**: (New) A tab-separated file mirroring the GRAND-SLAM output format.
    *   Columns: `Gene`, `Symbol`, `Readcount`, `0.05 quantile`, `Mean`, `MAP`, `0.95 quantile`, `alpha`, `beta`, `Conversions`, `Coverage`, `Double-Hits`, `Double-Hit Coverage`, `min2`, `Length`.
    *   Note: `MAP` (Maximum A Posteriori) and `Mean` correspond to the NTR (New Transcript Ratio).
2.  **`<prefix>SlamQuant.out`**: The standard STAR-Slam output format.
3.  **`<prefix>SlamQuant.cB.tsv`**: Optional count-binomial output when `--slamCbOut 1`.
    *   STAR schema: `sample`, `feature_id`, `feature_name`, `nT`, `TC`, `n`.
    *   EZbakR-like schema: `sample`, `rname`, `sj`, `GF`, `XF`, `nT`, `TC`, `n`.
    *   `n` is the STAR-SLAM weighted read count for that `(sample, feature, nT, TC)` bin and can be fractional.
4.  **`<prefix>.slam_qc.html`**: The interactive QC report (if `--slamQcReport` is used).

## Binary Dump Format

When `--slamDumpBinary 1` / `--slamDumpWeights 1` are enabled, STAR-SLAM writes
binary dumps (`*_slam_dump.bin`) and weight sidecars (`*_slam_weights.bin`).
In batch mode with `--outFileNamePrefixAuto 1`, these live under
`alignments/<sample>/` alongside the BAMs.
The **bitwise** file format (header + records) is documented in:

```
slam/docs/SLAM_DUMP_FORMAT.md
```

## Theory: Variance-Based Trimming

Standard trimming methods often rely on base quality (Phred) or simple adapter matching. However, in SLAM-seq, chemical modifications can cause specific conversion artifacts at read ends that don't always correspond to low quality scores.

The **variance-based** approach is based on the observation that true metabolic labeling signals (T->C) should be distributed somewhat uniformly or stochastically across the read length in the aggregate, whereas artifacts often manifest as highly variable conversion rates at specific positions (usually ends).

By plotting the **Standard Deviation of T->C rates** vs **Read Position**, we typically observe:
*   High/Variable Stdev at 5' end (Artifacts).
*   Low/Flat Stdev in the middle (Clean signal).
*   High/Variable Stdev at 3' end (Artifacts).

The algorithm fits a **Piecewise Linear Regression** (segmented regression) to this curve. It uses the Bayesian Information Criterion (BIC) to select the best model (2, 3, or 4 segments) to accurately locate the "knees" where the variance stabilizes. These knees become the 5' and 3' trim coordinates.

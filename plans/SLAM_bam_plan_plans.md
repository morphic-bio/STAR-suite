# SLAM BAM Plan (Blank + BAM parity)

Goal: establish a reproducible BAM‑based SLAM pipeline that (a) uses a no4su blank to estimate baseline T>C rates and trimming, (b) applies SNP masking consistently, and (c) achieves parity vs GEDI (and STAR‑SLAM outputs) on a 1M‑read subset before scaling.

## Scope
- Inputs: existing 1M downsampled FASTQs in `/mnt/pikachu/NW-5-21/SLAM-Seq-1M`, VCF → BED mask.
- Blank: no4su (per‑group blank; global blank used only as QC).
- Parity target: Pearson r ≥ 0.9 for transition rates (GEDI vs STAR‑SLAM outputs).
- Deliverables: scripts + documented runs + artifacts locations.

## Stage 0 — Inventory + SNP mask (done/verify)
**Objective:** ensure consistent SNP masking assets and file list.
1) Confirm SNP BED mask exists and is compressed/indexed:
   - `/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz`
2) Enumerate 0h‑1 and 24h‑1 FASTQs (1M):
   - Produce a table with sample, timepoint, and R1/R2 files.
3) Document the list + counts in `plans/HANDOFF_slam_blank_stage1_2_20260201.md`.

**Exit criteria:** mask file exists; file list table written.

## Stage 1 — Blank model from no4su (done/verify)
**Objective:** estimate baseline T>C for no4su and derive trimming constants.
1) Run blank model script on no4su BAMs or STAR‑SLAM outputs (per plan).
2) Save outputs under:
   - `/mnt/pikachu/slam_blank_artifacts_20260201/blank_model_no4su/`
3) Record `p_e_blank`, `p_e_regression`, and chosen trims.

**Exit criteria:** `blank_transition_model.tsv`, `pe_estimates.tsv`, and summary JSON exist.

## Stage 2 — 1M FASTQ alignment + STAR‑SLAM transitions
**Objective:** generate STAR‑SLAM transition summaries with SNP masking and no trimming (compat mode).
1) Align 0h‑1 + 24h‑1 1M FASTQs with STAR‑SLAM:
   - `--runThreadN 24`
   - `--slamSnpMaskIn /mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz`
   - no trimming for compatibility run
2) Collect `SlamQuant.out.transitions.tsv` and per‑sample summaries.
3) Store outputs under:
   - `/mnt/pikachu/slam_blank_artifacts_20260201/fastq_1M_runs/`

**Exit criteria:** transitions TSVs exist for all 0h‑1 and 24h‑1 samples.

## Stage 3 — GEDI parity (BAM‑based)
**Objective:** compare GEDI vs STAR‑SLAM using the same SNP‑masked data.
1) Generate SNP‑filtered BAMs for GEDI:
   - `bedtools intersect -v -abam <bam> -b snps_from_vcf.bed.gz > <bam>.snpmask.bam`
2) Run GEDI with the SNP‑masked BAMs and the same reference/annotations.
3) Compare transition rates (per gene/global) and compute Pearson r.

**Exit criteria:** Pearson r ≥ 0.9 on 1M subset; discrepancy analysis logged.

## Stage 4 — Trimming using blank‑derived constants
**Objective:** run STAR‑SLAM with constant trimming derived from no4su blank.
1) Apply trim constants (from Stage 1) in STAR‑SLAM and re‑compute transitions.
2) Compare to compat mode and GEDI parity.

**Exit criteria:** trimmed results documented; parity preserved or improved.

## Stage 5 — Scale‑up and QC
**Objective:** verify stability on larger runs and QC blank vs per‑group blank.
1) Run on larger subset or full data (if needed).
2) Compare per‑group blank estimates vs global blank (QC only).
3) Document in `plans/HANDOFF_slam_blank_stage1_2_20260201.md`.

**Exit criteria:** QC checks pass and differences are explained.

## Deliverables Checklist
- Scripts:
  - FASTQ inventory + STAR‑SLAM run wrapper
  - GEDI run wrapper
  - BAM SNP‑masking helper
  - Summary table generator
- Documentation updates:
  - `plans/HANDOFF_slam_blank_stage1_2_20260201.md`
  - `tests/ARTIFACTS.md` for new outputs
- Results:
  - Parity tables (STAR‑SLAM vs GEDI)
  - Blank model outputs (no4su)

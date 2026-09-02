# MorPhiC Data Working Group 3-Dataset Slide Runbook

Date: 2026-05-15

Purpose: outline a 6-8 slide update for MorPhiC program officers, PIs,
co-Is, and analysts covering the current STAR-suite release status for:

- MSK 30KO Perturb-seq
- NW SLAM-seq paired-end production
- NW ATAC-seq Chromap-suite/RapidMACS production

The recommended deck is 8 slides. For a 6-slide version, merge slides 3/4 and
5/6.

## Slide 1 - Why This Update

Program-officer / PI message:

- Three production-facing data surfaces are now ready or close to ready for
  working-group review.
- The work is focused on making each release reproducible, inspectable, and
  usable by downstream analysts, not just producing primary files.
- Each dataset has a different maturity profile: MSK is single-cell
  h5ad/QC-focused, SLAM has completed expression DESeq2, and ATAC has a new
  faster Chromap-suite/RapidMACS peak-calling path.

Co-I / analyst detail:

- This is a release/status deck, not a biological interpretation deck.
- Emphasize exact release paths, file layouts, known caveats, and the intended
  first analysis surfaces.
- Use the release READMEs as the detailed source of truth.

## Slide 2 - Release Map Across The Three Datasets

Program-officer / PI message:

- The three releases are organized by dataset-specific user workflows:
  h5ad/QC objects for single-cell, per-sample SLAM outputs plus DESeq2, and
  per-sample ATAC tracks/peaks plus large-file archives.
- Large files are separated from analyst-facing result files where possible.

Co-I / analyst detail:

| Dataset | Current release surface | Globus / handoff shape | Main exclusions |
| --- | --- | --- | --- |
| MSK 30KO Perturb-seq | `/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/` | BAM archive root `/MSK-05-13-26-large-files/`; staged h5ad/QC tree mirrors previous MSK sample-wise layout | FASTQs, raw MEX, CellBender raw `.h5`, copied provider `.rds` metadata |
| SLAM-seq PE | `/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/` plus curated upload trees | Globus root `SLAM-seq-PE-results/prod_full_20260511T103908Z/` with `samples/`, `de/`, manifests, README | Large local BAM/FASTQ files from curated smaller handoff |
| ATAC-seq | `/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/` | Globus root `/ATAC-results-05-05-26/`; per-sample `_large_files/` bundles plus root `README.md` | FASTQs not duplicated in analyst-facing result tree |

Note for the deck: be explicit that the ATAC root README has been transferred
to `/ATAC-results-05-05-26/README.md`.

## Slide 3 - MSK 30KO: What Was Done

Program-officer / PI message:

- MSK 30KO was reprocessed from FASTQs for the nine provider reference groups:
  `ES`, `DE`, `DE_GemX`, `PP1`, `PP2`, `S5_1`, `S5_2`, `S6_1`, and `S6_2`.
- The corrected `DE_GemX` dataset is treated as the newer, correct GEM-X DE
  surface where conflicts exist.
- The production output is an AnnData-based release with GeneFull counts,
  Velocyto layers, CellBender denoising, CRISPR guide calls, LARRY barcode
  calls, and transferred cell-type labels.

Co-I / analyst detail:

- FASTQ manifest: `docs/MSK_30KO_FASTQ_MANIFEST.tsv`.
- Provider comparison reference: `/mnt/pikachu/df.meta.rds`.
- Historical h5ad files are not analysis inputs.
- No Y-chromosome removal was applied for MSK.
- BAM files were retained separately and uploaded by Globus; FASTQs were not
  uploaded.
- GEM-X chemistry uses the May 2023 GEM-X whitelist family and PolyIII
  NXT-to-TRU translation whitelist. This is the main protocol difference from
  the other eight sample groups.

## Slide 4 - MSK 30KO: Cell Typing, QC, And Release Contents

Program-officer / PI message:

- Scimilarity is used as an embedding space, not as the final cell-type caller.
- Provider labels from `/mnt/pikachu/df.meta.rds` seed a RandomForest model in
  Scimilarity space; the saved model is then applied consistently across the
  full MSK sample set.
- Provider metadata is deliberately not copied into the h5ad files to avoid
  confusing reference labels with production annotations.

Co-I / analyst detail:

- Main h5ads per sample:
  - `downstream_genefull_velocyto_cellbender/filtered_counts.h5ad`
  - `downstream_genefull_velocyto_cellbender/default_singlet_filtered_counts.h5ad`
  - feature-library h5ads under `downstream_genefull_velocyto_cellbender/feature_libraries/`
- Included fields:
  - Velocyto layers: `spliced`, `unspliced`, `ambiguous`
  - CellBender layer: `denoised`
  - CRISPR and LARRY obs feature-call columns
  - `obsm["X_scimilarity"]`
  - `rf_celltype`, `rf_subcelltype`, confidence, and margin columns
- Excluded fields:
  - copied provider `.rds` metadata
  - provisional Scimilarity kNN labels
- Model artifacts:
  - `results/provider_rf_celltyping_20260515_065032/models/rf_celltype.joblib`
  - `results/provider_rf_celltyping_20260515_065032/models/rf_subcelltype.joblib`
- Held-out RF validation:
  - celltype accuracy 95.52%, balanced accuracy 87.92%
  - subcelltype accuracy 94.44%, balanced accuracy 80.85%
- Provider comparison summary:
  - STAR non-empty recovers provider cells at 99.59% to 100.00%
  - CRISPR provider-positive confident-call recovery 92.52% to 98.06%
  - LARRY exact top-barcode concordance in filtered cells 98.63% to 99.87%
  - RF celltype/subcelltype exact match on matched provider cells 99.45% to 100.00%

## Slide 5 - SLAM-seq PE: What Was Done

Program-officer / PI message:

- The SLAM-seq panel has been rerun as a paired-end STAR-SLAM production
  workflow.
- The PE run adds model-ready count-binomial output and fixes several
  paired-end-specific methodology issues from the earlier single-end workflow.
- The release includes per-sample SLAM outputs, QC reports, TranscriptVB
  expression outputs, and a completed DESeq2 expression analysis.

Co-I / analyst detail:

- Input FASTQs: `/mnt/pikachu/SLAM-Seq`.
- Samples processed: 290 paired-end samples.
- Run root: `/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z`.
- Genome index: `/storage/autoindex_110_44/bulk_index`.
- SNP mask: `/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz`.
- Core method:
  - paired-end STAR input
  - TranscriptVB expression quantification
  - SLAM strandness `Sense`
  - minimum callable length 30
  - count-binomial output in STAR format
  - sample-wise Globus transfer and cleanup
- PE method changes:
  - mate-specific trims: mate 1 `5p=8, 3p=13`; mate 2 `5p=19, 3p=14`
  - repaired PE TranscriptVB autodetection
  - `BanSingleEnd` transcriptome output to avoid singleton quantification
  - corrected PE mismatch coordinate/overlap reporting

## Slide 6 - SLAM-seq PE: Outputs, DESeq2, And Caveats

Program-officer / PI message:

- The primary SLAM output is per-gene NTR and conversion evidence.
- The primary expression DE output is DESeq2 on STAR TranscriptVB
  tximport-compatible gene expression, not on SLAM NTR or count-binomial tables.
- DESeq2 has been run in three modes to separate primary inference,
  optimistic sensitivity, and conservative descriptive summaries.

Co-I / analyst detail:

- Key per-sample files:
  - `counts/<sample>.SlamQuant.out`
  - `counts/<sample>.SlamQuant.cB.tsv`
  - `qc/<sample>.slam_qc.html` and `.json`
  - `run/star_quant.genes.tximport.sf`
  - `run/star_Log.final.out`
- DESeq2 input surface:
  - `samples/<sample>/run/star_quant.genes.tximport.sf`
  - imported with `tximport`
- Completed DESeq2 modes:
  - `blocked_pairwise`: `design <- ~ set_id + time`
  - `unblocked_pairwise`: `design <- ~ time`
  - `collapsed_pairwise`: collapsed timepoint summaries, no DESeq2 p-values
- DESeq2 status:
  - 290 samples in DE manifest
  - 35 target labels, 34 normalized perturbed genes
  - 33 ready targets run through the full panel
  - 89 valid pairwise contrasts per mode
  - 0 failed jobs
- Important caveats:
  - six R1/R2 name mismatches were paired by shared sample index and are listed
    in `manifests/pairing_warnings.tsv`
  - `MYC` and `RAD21` were kept review-only for the full DESeq2 panel
  - SLAM cB/NTR outputs should not be treated as ordinary DESeq2 count matrices

## Slide 7 - ATAC-seq: What Was Done

Program-officer / PI message:

- The NW ATAC-seq set was processed with Chromap-suite and the integrated
  RapidMACS (`librapidmacs`) FRAG peak-calling path.
- This replaces the older external MACS2/MACS3 BAMPE callpeak workflow for this
  delivery.
- The new librapidmacs library path was approximately 8x to 11x faster in
  validation for the peak-calling stage.

Co-I / analyst detail:

- Release root: `/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/`.
- Globus root: `/ATAC-results-05-05-26/`.
- Result set:
  - 31 named targets plus `unknown`
  - 200 complete processed sample outputs
  - 1 explicit skip marker for `ATAC-NUP133-0h-1` because it had already been
    processed as the final validation sample
- Mapping/processing:
  - paired-end bulk ATAC
  - adapter trimming
  - duplicate removal
  - unique mapping
  - MAPQ 30
  - maximum insert size 2000
  - Tn5 4/4 symmetric offset convention
- Peak-calling parameters:
  - `macs3_frag_pvalue = 0.01`
  - `bdgpeakcall_cutoff_neg_log10_p = 2`
  - `macs3_frag_min_length = 200`
  - `macs3_frag_max_gap = 30`
  - `macs3_uint8_counts = 1`
  - `effective_genome_size = 2913022398`
  - `fragments_source = memory`

## Slide 8 - ATAC-seq Release Contents And Working-Group Decisions

Program-officer / PI message:

- The ATAC release gives analysts the expected tracks, fragments, peak calls,
  summits, summaries, and logs.
- Large files are uploaded separately so analysts can work from compact result
  files without pulling BAMs or split FASTQs unless needed.
- The main next decision is how much downstream biological summarization to add
  on top of the core release.

Co-I / analyst detail:

- Per-sample result files:
  - `<sample>.fragments.bed`
  - `<sample>_unstranded.bw`
  - `<sample>.peaks.narrowPeak`
  - `<sample>.peaks.summits.bed`
  - `<sample>.summary.txt`
  - `<sample>.summary.txt.macs3_frag_peaks.tsv`
  - `<sample>.y_readnames.txt`
  - Chromap stdout/stderr logs
  - Globus task metadata and monitor logs
- Large-file Globus summary:
  - 200 per-sample Globus tasks
  - 200 succeeded
  - 2,000 large files recorded
  - 1,076,632,529,803 bytes recorded
  - destination pattern: `/ATAC-results-05-05-26/<target>/<timepoint>/<sample>_large_files/`
- Older peak-calling parameter reference:
  - `macs3 callpeak -t <sample>.bam -f BAMPE -g hs -q 0.01 --keep-dup all`
  - `macs2 callpeak -t <sample>.bam -f BAMPE -g hs -q 0.01 --keep-dup all`
- Caveat:
  - `estfrip` is zero in Chromap summaries and should not be treated as the
    final peak FRiP metric for this delivery.

## Six-Slide Compression

Use this if time is tight:

1. Why this update and release map.
2. MSK 30KO status plus release contents.
3. MSK cell typing/QC details.
4. SLAM-seq PE status plus method changes.
5. SLAM-seq DESeq2 and release details.
6. ATAC-seq status, RapidMACS update, release details.

## Suggested Closing Slide Or Speaker Close

Program-officer / PI message:

- All three datasets now have reproducible production-facing outputs and
  documented release layouts.
- The immediate value is analyst access to standardized outputs with clear
  provenance and known caveats.
- The next layer is biological synthesis: target-level summaries, cross-assay
  comparisons, and decisions about which QC metrics become formal release gates.

Co-I / analyst action list:

- Confirm final public-facing Globus roots for each release.
- For MSK, verify that the intended nine-sample h5ad/QC tree, including
  `DE_GemX`, is the release surface analysts should use.
- For SLAM, review the six pairing-warning samples before interpretation.
- For ATAC, compute or agree on final FRiP/peak QC summaries rather than using
  the zero-valued `estfrip` field.
- Keep the release READMEs synchronized with any re-upload or file-layout
  changes.

## Source Documents

- `docs/RUNBOOK_MSK_30KO_GEMX_PRODUCTION_DRAFT.md`
- `/mnt/pikachu/msk-30ko-production-planning/docs/RUNBOOK_MSK_30KO_PRODUCTION_AND_CELL_TYPING.md`
- `/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/README.md`
- `/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/README.md`
- `/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z_deseq2_globus_upload_20260515T061311Z/README.md`
- `/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/README.md`

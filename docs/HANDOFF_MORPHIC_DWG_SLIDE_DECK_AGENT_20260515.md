# Handoff: MorPhiC Data Working Group Slide Deck Agent

Date: 2026-05-15

Audience: MorPhiC data working group, program officers, PIs, co-Is, and
analysts.

Task: create the actual 6-8 slide deck summarizing three STAR-suite production
data surfaces:

- MSK 30KO Perturb-seq
- NW SLAM-seq paired-end production
- NW ATAC-seq Chromap-suite/RapidMACS production

Start from the slide outline here:

```text
/mnt/pikachu/STAR-suite/docs/RUNBOOK_MORPHIC_DWG_3_DATASET_SLIDES_20260515.md
```

That runbook has the recommended 8-slide structure and a 6-slide compression.
Use this handoff for file locations, QC/analysis assets, and caveats.

## Recommended Deck Shape

Use 8 slides if there is room:

1. Why this update / release framing
2. Release map across the three datasets
3. MSK 30KO: what was done
4. MSK 30KO: cell typing, QC, release contents
5. SLAM-seq PE: what was done
6. SLAM-seq PE: outputs, DESeq2, caveats
7. ATAC-seq: what was done
8. ATAC-seq: release contents and next decisions

For 6 slides, merge slides 3/4 and 5/6.

Style target: high-level statements on the slide face; detailed paths and
method parameters in speaker notes or an appendix.

## Dataset And Release Locations

### MSK 30KO Perturb-seq

Primary local h5ad/QC staging tree:

```text
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/
```

Main documentation:

```text
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/README.md
/mnt/pikachu/msk-30ko-production-planning/docs/RUNBOOK_MSK_30KO_PRODUCTION_AND_CELL_TYPING.md
/mnt/pikachu/STAR-suite/docs/RUNBOOK_MSK_30KO_GEMX_PRODUCTION_DRAFT.md
/mnt/pikachu/STAR-suite/docs/MSK_30KO_FASTQ_MANIFEST.tsv
```

Reference provider metadata used for comparison only:

```text
/mnt/pikachu/df.meta.rds
```

Globus large-file root visible on Morphic Processing:

```text
/MSK-05-13-26-large-files/
```

Important: the visible Globus MSK root contains retained BAMs by sample, for
example `/MSK-05-13-26-large-files/ES/Aligned.out.bam`. The local h5ad/QC
release tree above is the explicit staged analyst-facing tree. Before claiming
the h5ad/QC tree is on Globus, verify that upload separately.

### SLAM-seq PE

Primary local production root:

```text
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/
```

Curated smaller handoff tree, excluding BAM/FASTQ:

```text
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z_remaining_globus_upload_20260512T184310Z/
```

DESeq2 handoff tree:

```text
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z_deseq2_globus_upload_20260515T061311Z/
```

Main documentation:

```text
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/README.md
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z_deseq2_globus_upload_20260515T061311Z/README.md
```

Globus root:

```text
SLAM-seq-PE-results/prod_full_20260511T103908Z/
```

The Globus root currently shows `README.md`, `samples/`, `de/`, and DESeq2
upload manifests.

### NW ATAC-seq

Primary local result tree:

```text
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/
```

Main documentation:

```text
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/README.md
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/README.md
```

Globus root:

```text
/ATAC-results-05-05-26/
```

The root README has already been uploaded:

```text
/ATAC-results-05-05-26/README.md
```

Per-sample large-file folders are under:

```text
/ATAC-results-05-05-26/<target>/<timepoint>/<sample>_large_files/
```

## Deck-Ready QC And Analysis Assets

### MSK 30KO Assets

Useful slide-ready visual assets are already present for all 9 samples:

- 9 `gene_quantile_histogram.png` files
- 9 `cellbender_counts.pdf` files
- 9 `cellbender_counts_report.html` files

Path pattern:

```text
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/<sample>/downstream_genefull_velocyto_cellbender/gene_quantile_histogram.png
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/<sample>/downstream_genefull_velocyto_cellbender/cellbender/cellbender_counts.pdf
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/<sample>/downstream_genefull_velocyto_cellbender/cellbender/cellbender_counts_report.html
```

Representative examples:

```text
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/DE_GemX/downstream_genefull_velocyto_cellbender/gene_quantile_histogram.png
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/DE_GemX/downstream_genefull_velocyto_cellbender/cellbender/cellbender_counts.pdf
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/DE_GemX/downstream_genefull_velocyto_cellbender/cellbender/cellbender_counts_report.html
```

Provider comparison and RF cell-typing tables:

```text
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/analysis/results/provider_comparison_rf_20260515_065652/cell_calling_concordance.tsv
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/analysis/results/provider_comparison_rf_20260515_065652/crispr_concordance.tsv
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/analysis/results/provider_comparison_rf_20260515_065652/larry_concordance.tsv
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/analysis/results/provider_comparison_rf_20260515_065652/celltype_rf_concordance.tsv
/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/analysis/results/provider_rf_celltyping_20260515_065032/rf_validation_metrics.tsv
```

Good plot ideas if the deck needs more than screenshots:

- Bar chart: provider cell recovery in STAR non-empty, CellBender, filtered
  counts, and default singlet sets by sample.
- Bar chart: CRISPR confident-call recovery by sample.
- Bar chart: LARRY exact top-barcode concordance by sample.
- Small table: RF held-out validation metrics for `rf_celltype` and
  `rf_subcelltype`.

Key numbers to use, from the production runbook:

- STAR non-empty recovers provider cells at 99.59% to 100.00%.
- CellBender recovers provider cells well; minimum observed provider recovery
  is 96.70% for `PP1`.
- CRISPR provider-positive confident-call recovery is 92.52% to 98.06%.
- LARRY exact top-barcode concordance in filtered cells is 98.63% to 99.87%.
- RF held-out validation: celltype accuracy 95.52%, balanced accuracy 87.92%;
  subcelltype accuracy 94.44%, balanced accuracy 80.85%.

Do not show provider labels as h5ad-native annotations. Provider `.rds` labels
were used for external comparison and RF training, but copied provider metadata
was intentionally not placed in the production h5ads.

### SLAM-seq PE Assets

There are 290 per-sample QC HTML reports and 290 matching JSON files.

Path pattern:

```text
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/samples/<sample>/qc/<sample>.slam_qc.html
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/samples/<sample>/qc/<sample>.slam_qc.json
```

Representative examples:

```text
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/samples/ARID1A-0h-1_S40/qc/ARID1A-0h-1_S40.slam_qc.html
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/samples/ARID1A-0h-1_S40/qc/ARID1A-0h-1_S40.slam_qc.json
```

The QC JSON includes trim and per-mate QC fields such as:

```text
version=2.1
type=comprehensive_qc
variance_histogram_mode=per_mate_separate
mate_len_mate1=50
mate_len_mate2=50
trim5p_mate1=8
trim3p_mate1=13
trim5p_mate2=19
trim3p_mate2=14
```

For a slide, screenshot one representative QC HTML report if a visual is
needed. The HTML reports summarize position-specific T-to-C behavior,
mate-specific trim behavior, mismatch/quality behavior, and related QC panels.

Core DESeq2 summary tables:

```text
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/run_inventory.tsv
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/padj_counts_by_mode_contrast.tsv
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/blocked_pairwise.contrast_summaries.tsv
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/unblocked_pairwise.contrast_summaries.tsv
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/collapsed_pairwise.contrast_summaries.collapsed.tsv
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/gene_level_consensus.tsv.gz
```

Useful DESeq2 status numbers:

- 290 input samples.
- 33 ready targets run through the full panel.
- 89 valid pairwise contrasts per mode.
- 0 failed DESeq2 jobs.
- Modes: `blocked_pairwise`, `unblocked_pairwise`, `collapsed_pairwise`.

Optional preliminary biological-analysis tables exist locally, but they were
not part of the curated DESeq2-only Globus handoff. Use only if clearly labeled
as preliminary:

```text
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/target_knockdown_self_check_summary.tsv
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/known_correlate_marker_check_summary.tsv
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/literature_signature_de_ntr_summary.tsv
/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z/de/full_modes_20260515T033955Z/combined_summaries/tier1_marker_method_concordance.tsv
```

Good plot ideas:

- Small table or bar chart from `run_inventory.tsv`.
- Bar chart from `padj_counts_by_mode_contrast.tsv`, showing total
  padj < 0.05 counts by mode/contrast.
- One representative SLAM QC HTML screenshot.

Do not use SLAM NTR or count-binomial tables as ordinary DESeq2 count matrices.
The DESeq2 surface is STAR TranscriptVB `star_quant.genes.tximport.sf`.

### ATAC-seq Assets

The ATAC release tree does not currently contain standalone `.png`, `.pdf`, or
`.html` QC plots. It does contain slide-usable track and peak files:

```text
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/<target>/<timepoint>/<sample>/<sample>_unstranded.bw
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/<target>/<timepoint>/<sample>/<sample>.peaks.narrowPeak
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/<target>/<timepoint>/<sample>/<sample>.peaks.summits.bed
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/<target>/<timepoint>/<sample>/<sample>.fragments.bed
```

Representative example:

```text
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/NUP210/6h/ATAC-NUP210-6h-1/ATAC-NUP210-6h-1_unstranded.bw
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/NUP210/6h/ATAC-NUP210-6h-1/ATAC-NUP210-6h-1.peaks.narrowPeak
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/NUP210/6h/ATAC-NUP210-6h-1/ATAC-NUP210-6h-1.summary.txt
/mnt/pikachu/NW-ATAC-SEQ-results-libmacs3-05-05-26/final/NUP210/6h/ATAC-NUP210-6h-1/ATAC-NUP210-6h-1.summary.txt.macs3_frag_peaks.tsv
```

Good plot ideas:

- Generate a simple aggregate QC chart from per-sample `*.summary.txt`:
  total reads, duplicate reads, unmapped reads, low MAPQ reads.
- Generate a peak-count histogram from all `*.peaks.narrowPeak` files.
- Use IGV or a genome-browser screenshot from representative bigWig +
  narrowPeak tracks.
- Use the README's inventory numbers as a clean status graphic:
  200 complete samples, 200 bigWigs, 200 fragments, 200 narrowPeak files,
  200 summit files, and 200 successful large-file Globus tasks.

Known ATAC caveat:

- `estfrip` is zero in Chromap summaries and should not be shown as the final
  peak FRiP metric. If a FRiP slide is desired, compute it separately from the
  final peaks and fragments/BAMs.

## Do Not Use / Avoid

- Do not use historical MSK h5ad files as analysis inputs.
- Do not frame MSK "Scimilarity labels" as final labels. The final MSK labels
  are provider-seeded RF labels in Scimilarity embedding space.
- Do not include copied provider `.rds` metadata as if it were embedded in the
  production h5ads.
- Do not use older ATAC roots such as `/NW-ATAC-Seq-results/` as the new
  RapidMACS release. Those reflect the older external MACS2/MACS3 layout.
- Do not claim SLAM exploratory literature or marker summaries are part of the
  curated DESeq2-only handoff unless that is explicitly desired and labeled.
- Do not treat ATAC `estfrip=0` as a QC failure or as a valid FRiP metric.

## Quick Narrative Points

Use these if the deck needs concise text:

- MSK: "Nine provider-defined sample groups were reprocessed from FASTQs into
  UCSF-style AnnData objects with Velocyto, CellBender, CRISPR, LARRY, and
  provider-seeded RF cell typing."
- SLAM: "The PE production run processed 290 samples with STAR-SLAM, added
  count-binomial outputs, and completed DESeq2 from TranscriptVB expression
  using blocked, unblocked, and collapsed modes."
- ATAC: "The ATAC run processed 200 complete sample outputs with Chromap-suite
  and the new integrated RapidMACS (`librapidmacs`) FRAG peak-calling path,
  approximately 8x to 11x faster for peak calling in validation."

## If More Figures Are Needed

Suggested fast follow-ups for a plotting agent:

1. MSK: produce a four-panel concordance summary from the provider comparison
   TSVs.
2. SLAM: produce a DESeq2 status bar chart from
   `padj_counts_by_mode_contrast.tsv` and a compact run-inventory table.
3. ATAC: produce an aggregate sample QC plot from all `*.summary.txt` and a
   peak-count distribution from all `*.peaks.narrowPeak`.

Keep any generated figures next to this handoff or in a named artifact folder,
for example:

```text
/mnt/pikachu/STAR-suite/docs/artifacts/morphic_dwg_20260515/
```

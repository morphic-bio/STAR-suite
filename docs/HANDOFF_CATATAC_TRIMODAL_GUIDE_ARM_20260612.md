# HANDOFF — CAT-ATAC trimodal single-pass: wire the CRISPR guide arm

**Date:** 2026-06-12 · **Audience:** STAR-suite agent (review + implement code) ·
**Status:** all data staged; needs guide-arm wiring + a downsampled smoke, then a full concordance run.

## 1. Goal and why it matters (read first)

Demonstrate a **single STAR-suite run that emits all three multiome perturb-seq
modalities — GeneFull GEX + ATAC peaks + CRISPR guide calls** — on the CAT-ATAC
dataset (GSE288996, K562 DMSO rep1).

This is strategically load-bearing: it is the **in-pocket counter to the "STAR-suite
is just a 10x re-implementation" critique** on the current Nature Methods submission.
No vendor pipeline does single-pass trimodal:

- `cellranger-arc count` = GEX + ATAC, **no feature/guide library** (its only library
  types are Gene Expression + Chromatin Accessibility).
- `cellranger multi`/`count` = GEX + CRISPR guide, **no ATAC**.

The CAT-ATAC authors themselves had to **stitch** it: `cellranger-arc count` for
GEX+ATAC **plus a separate hand-rolled `bbduk`+`awk` guide caller** (see their repo,
§5). The assay is published (Cell Rep Methods 2025); the **novelty here is the
unified single-pass engine**, which no commercial tool — and not even the assay's
authors — provides. GEX+ATAC is already benchmarked (suite ~42 min vs CR-ARC 2:28 +
Signac re-call ~3 min; see `/mnt/pikachu/catatac_gse288996/BENCHMARK.md`). The guide
arm is the missing third modality.

## 2. The ask (TL;DR)

1. Wire the **CAT-ATAC guide-capture library** into the suite's
   `process_features`/`assignBarcodes` + `call_features` path to produce a
   guide×cell matrix + per-cell guide assignments.
2. Integrate it as **library #3** in the multiome recipe so one orchestrated run
   emits GeneFull + ATAC peaks + guide calls (compose-up: add a guide profile/flag,
   defaults unchanged).
3. **Smoke** on the downsampled fixture (functional: does it run + produce a
   guide×cell matrix). Then a **full** run for concordance vs the deposited calls.

The main review question is whether `assignBarcodes` already supports this *exact*
read layout (CB/UMI/feature spread across **three different reads**); the engine and
ATAC-mode support exist, so expect wiring + minor bug-catching, not new architecture.

## 3. What is already in the suite (do NOT rebuild)

- **GEX+ATAC in one binary**: STAR built `WITH_CHROMAP=1` runs STARsolo GeneFull +
  chromap ATAC in-process. Recipe: `morphic-recipes/scripts/run_star_multiome_lane_smoke.sh`
  (+ `run_multiome_minimal.sh`), already validated on this dataset.
- **Feature/guide engine**: `core/features/process_features/` (`libprocess_features.a`,
  `pf_api.h`, `call_features.c`) + standalone `/mnt/pikachu/process_features/`
  (`assignBarcodes`, `demux_fastq`, `feature_lists/`). README states `assignBarcodes`
  does **exhaustive search for feature barcodes in ATAC-seq AND RNA-seq data** + fuzzy
  matching + UMI dedup — i.e. ATAC-side guide capture (CAT-ATAC's case) is a supported
  mode.
- **Precedent**: CRISPR guide-call concordance already done for MSK perturb-seq —
  `comparisons/msk_30polyko_full_benchmark_20260306/crispr_call_concordance.json`;
  GEX-side guide example `scripts/run_ucsf_full_compat_forward_rescue_guides.sh`
  (appends GUIDE_R1/R2 to the STAR readFilesIn + a CR-style feature ref).
- **3 simultaneous libraries** already proven by the Larry perturb example.
- **ATAC barcode handling matches the guide R2 exactly**: the chromap ATAC read
  format we used is `bc:8:23:-` (16 bp barcode at 0-based 8–23, reverse-complement) +
  `atac2gex` translation. The CAT-ATAC guide R2 barcode is the **same 10x-ARC i5
  barcode** — the paper extracts it as `revcomp(substr(seq,9,24))` (1-based 9–24 =
  0-based 8–23). **Reuse the existing chromap ATAC barcode config for the guide R2.**

## 4. CAT-ATAC guide-capture spec (the new thing to support)

Bespoke **ATAC-side** capture, three reads per fragment (`SRR32265756`, K562 DMSO rep1):

| Read | len | content |
|---|---|---|
| R1 | ~100 bp | **capture sequence** (`CS1=CAAGTTGATAACGGACTAGCC`, `CS2=CAAGTTGTAAACGGACTAGCC`) + **UMI = first 12 bp** |
| R2 | 24 bp | **cell barcode** = 10x-ARC ATAC i5 (`bc:8:23:-`; `atac2gex` translate) |
| R3 | ~100 bp | **protospacer** (the spacer; match vs the library) |

- **Dual-guide library**: each target has 2 spacers named `GENE_1`, `GENE_2`
  (treated as separate features in the deposited calls). 56 CRISPRi spacer rows
  in the deposited reference (54 unique feature names after STAR-suite
  first-definition-wins dedup in `catatac_crispri_guide_capture.csv`)
  (28 targets × 2), 19–20 bp.
- Reference (status-quo) matching logic, for parity: bbduk matches CS1/CS2 in R1
  (k=22) and the protospacer in R3 (k=20, hdist 1) vs a protospacer FASTA; UMI =
  R1[1:12]; CB = R2 (`atac2gex`); dedup UMIs per (cell, guide). Reads are pre-filtered
  to those containing CS1/CS2 in R1.

## 5. Data locations (all on this host)

- **Spacer library**: `/mnt/pikachu/catatac_gse288996/guide_ref/ps_ref_crispri.tsv`
  (cols: `target_gene  guide_name  spacer`; 56 rows) and
  `catatac_guide_library_full.tsv` (+CRISPRa). Source: paper `mmc2.xlsx`
  (`/mnt/pikachu/1-s2.0-S2667237525002589-mmc2.xlsx`), sheet `CRISPRi guide sequences_1/_2`.
- **Downsampled smoke fixture** (2M reads, in sync):
  `/mnt/pikachu/catatac_gse288996/guide_redump/fixture/guide_R{1,2,3}.fastq.gz` (+MANIFEST).
- **Full guide FASTQs** (12.6M reads): `/mnt/pikachu/catatac_gse288996/guide_redump/SRR32265756_{1,2,3}.fastq.gz`
  (`_1`=R1 capture+UMI, `_2`=R2 24bp barcode, `_3`=R3 protospacer). Re-dumped with
  `fasterq-dump --include-technical` (the original dump dropped the technical barcode read).
- **Whitelists**: ATAC `…/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt`,
  `atac2gex.tsv` (same dir); GEX `/mnt/pikachu/GEX_whitelist/737K-arc-v1.txt`.
- **Concordance target** (paper's per-cell guide calls):
  `/mnt/pikachu/catatac_gse288996/catatac_public/guide_caller_calls/DMSO1/DMSO1_calls.rds`
  + `DMSO1_thresholds.csv` (R `.rds`; read with R — the Signac container
  `catatac-signac:macs2-2.2.9.1` has R/Seurat).
- **GEX+ATAC** for the full trimodal run (2020-A refs):
  GEX `fastq/GEX/SRR32265752_{1,2}.fastq.gz`, ATAC `fastq/ATAC/SRR32265760_{1,2,3}.fastq.gz`;
  STAR GEX index `/mnt/pikachu/autoindex_98_32/pe_index`, chromap idx
  `…/catatac_gse288996/refs/GRCh38-arc.chromap.idx`, ARC genome
  `/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa`.
- **Reference (status-quo) pipeline**: `…/catatac_gse288996/catatac_public/00_catatac_pipeline/primary/`
  (`02_reads2tsv.sh`, `03_select_reads_by_guide_with_mismatch.sh`, `04_run_count.sh`).

## 6. Implementation tasks

1. **Feature reference** — convert `ps_ref_crispri.tsv` to the `process_features`
   feature-ref format (`feature_lists/*.csv`: `id,name,read,pattern,sequence,feature_type`).
   `sequence`=spacer, `read`=R3, `name`=guide_name, `feature_type`=CRISPR Guide Capture.
   Decide the `pattern` for a spacer matched within R3 (the paper matches it as a free
   kmer; `assignBarcodes` fuzzy/exhaustive search may not need a fixed anchor).
2. **assignBarcodes config for the 3-read layout** — CB from **R2** (`bc:8:23:-`,
   `atac2gex`), UMI from **R1[1:12]**, feature/spacer from **R3**. Confirm/extend
   `assignBarcodes` to take CB/UMI/feature from three *different* reads (the
   MSK/UCSF precedent was GEX-side: CB+UMI in R1, feature in R2 — **this is the new
   layout and the most likely bug surface**). Optionally pre-filter on CS1/CS2 in R1.
3. **call_features** → guide×cell matrix + per-cell assignment (dominant or GMM).
4. **Recipe integration (compose-up)** — add the guide library as an optional layer in
   `run_star_multiome_lane_smoke.sh` (e.g. `--guide-r1/--guide-r2/--guide-r3`,
   `--guide-ref`, and a `--profile` that turns it on) so one run emits GeneFull + ATAC
   peaks + guide calls; defaults stay unchanged (see `morphic-recipes/AGENTS.md`
   "Compose to the target" + the COMPOSITION block in that recipe).

## 7. Smoke + validation plan

- **Downsampled smoke** (fixture): guide arm runs end-to-end → non-empty guide×cell
  matrix containing the expected guides (`ADNP_1`, `BRD2_1`, … `HIC2_1/2`, `NTC*`).
  NB: fixture cells do **not** overlap the GEX/ATAC tiny fixture (independent
  downsamples) — that's fine for a wiring smoke; do not assert cross-modality overlap here.
- **Full concordance** (after wiring works): full guide run on real cells → compare
  per-cell guide assignments to `DMSO1_calls.rds` (agreement / Jaccard on shared
  barcodes; guide-level count correlation; assignment-rate vs the paper's ~stated rate).
  Read the `.rds` in the Signac container.

## 8. Open questions / likely bug surface

- Does `assignBarcodes` accept CB/UMI/feature across **three** reads (CB=R2, UMI=R1,
  feature=R3)? If not, this is the code to add/extend.
- Barcode orientation — should be the existing `bc:8:23:-` + `atac2gex` (verified equal
  to the paper's `revcomp(substr(R2,9,24))`); reuse, don't reinvent.
- Pre-filter reads on CS1/CS2 in R1 before counting (reduces noise; the paper does).
- Dual guides `GENE_1`/`GENE_2` are **separate features**.
- UMI dedup strategy appropriate for ATAC-side targeted capture.
- Future scaling TODO: the current mixed-length feature fallback is fine for
  this 54-guide library, but large mixed-length guide/lineage libraries should
  use length-bucketed feature prehashes. STAR-Spatial has an analogous local
  pattern in `/mnt/pikachu/STAR-Spatial/native/hd_sw_rescue_harness.cpp`
  (`group_oligos_by_length()` + length-window candidate search).

## 9. References

- Paper: Cell Reports Methods 2025, PII `S2667237525002589` (bioRxiv `2025.02.11.637716`).
  Guide library = supplementary `mmc2.xlsx`.
- `catatac_public`: github.com/ucsf-lgr/catatac_public (cloned at
  `…/catatac_gse288996/catatac_public/`) — the reference stitched pipeline.
- Suite: `core/features/process_features/README.md`,
  `docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md`, `docs/feature_barcodes.md`,
  `docs/RUNBOOK_PROCESS_FEATURES_CBQ_NATIVE.md`,
  `scripts/run_ucsf_full_compat_forward_rescue_guides.sh`.

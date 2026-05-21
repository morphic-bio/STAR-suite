# OCM Composite-Barcode scRNA-seq Runbook

Date: 2026-05-20

## Status

This runbook is the validation path for the native STAR-suite Flex-style OCM
effective barcode model. The STAR reader keeps the input FASTQ layout unchanged
(`CB16 + UMI12`) and, when `--ocmMultiBarcodeMode flex` is enabled, appends an
internal OCM TAG8 suffix to the cell barcode before inline barcode correction,
UMI collapse, EmptyDrops, Gene/GeneFull counting, and Velocyto counting.

This is Flex-style barcode semantics for a normal GEX assay, not a direct
`--flex yes` probe run. The current STAR-Flex direct-hash MEX writer is
probe-list oriented and uses a compact probe/gene key that is not the GeneFull
GEX surface. For OCM GEX production, keep `--soloInlineHashMode no` and use the
standard STAR Gene/GeneFull/Velocyto writers on the effective `CB16+TAG8`
barcode axis. The shared part with Flex is the methodological placement of the
sample tag before CB correction, UMI collapse, and cell calling.

The native C++ OCM materializer now handles the production GeneFull/Velocyto
split after STAR writes the pool MEX. `scripts/ocm_composite_adapter.py`
remains in the repo as a Cell Ranger comparator and fallback. It can still read
historical helper outputs that used a synthetic one-base `CB17` barcode, but
the harness now defaults to native STAR mode and reuses existing Cell Ranger
artifacts unless `--run-cr` is explicitly requested.

## Method

The current post-hoc OCM materializer corrects the ordinary 16 bp GEM-X cell
barcode first, then uses bases 8-9 of that corrected barcode to split samples.
That pools all OCM populations during barcode correction and cell calling. For
OCM this is the wrong prior: an abundant `CB-OB2` population can help rescue a
low-count `CB-OB1` observation.

The method treats the OCM assignment as part of the effective cell barcode
before correction and counting:

```text
effective_cb = CB16 + OCM_TAG8
UMI          = original UMI12
```

The native STAR mode uses fixed Flex-compatible eight-base tag suffixes:

| CB bases 8-9 | OCM ID | Internal TAG8 |
| --- | --- | --- |
| `GT` | `OB1` | `GTGTGTGT` |
| `CA` | `OB2` | `CACACACA` |
| `TC` | `OB3` | `TCTCTCTC` |
| `AG` | `OB4` | `AGAGAGAG` |

The effective whitelist is derived in STAR from the May-2023 GEM-X TRU
whitelist by appending the matching TAG8 to each OCM-classified barcode. This
preserves the original whitelist cardinality for classified OCM barcodes and
prevents barcode correction from crossing OCM populations.

After STAR counting, OCM outputs split by TAG8 and strip the TAG8 from
Cell Ranger-compatible per-sample barcode labels. The methodological boundary
is split before ED, not after pooled cell calling.

## Velocyto Semantics

Velocyto counting should use the same effective barcode axis as Gene/GeneFull.
That means spliced, unspliced, and ambiguous counts are accumulated after
composite barcode correction and UMI collapse, not as a late mask over a pooled
CB16 matrix. Per-sample Velocyto MEX output is then a direct split of the
already-corrected composite barcode matrix.

The smoke materializer writes `raw_velocyto_feature_bc_matrix` and
`filtered_velocyto_feature_bc_matrix` from `Solo.out/Velocyto/raw` for the
GeneFull materialization. In native mode those matrices are already on the
effective `CB16+TAG8` axis, so GeneFull totals and Velocyto totals share barcode
semantics.

When `--soloCellFilter None` is used, STAR does not create a pooled
`GeneFull/filtered` MEX. The OCM materializer treats this as the native
production path: it streams per-sample raw MEX, runs per-sample EmptyDrops, then
streams filtered GeneFull and Velocyto outputs from the same composite barcode
axis. This keeps the ED evidence separated by OCM tag instead of cell-calling a
pooled CB16 population.

## 100K Smoke

Run the full 100K gate from the repo root:

```bash
scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh \
  --read-pairs 100000 \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_composite_smoke_100k_$(date -u +%Y%m%dT%H%M%SZ) \
  --run-all
```

The harness stages the same downsample for both tools:

- original R1/R2 FASTQs for Cell Ranger 9 `multi` only when `--run-cr` is
  explicitly requested to build a new reference;
- original R1/R2 FASTQs for STAR native OCM-Flex mode;
- the original GEM-X TRU whitelist for STAR.

STAR uses:

```text
--genomeDir /storage/autoindex_110_44/bulk_index
--soloCBlen 16
--soloUMIstart 17
--soloUMIlen 12
--soloCBwhitelist /storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
--soloInlineCBCorrection yes
--soloUMIfiltering MultiGeneUMI_CR
--soloUMIdedup 1MM_CR
--soloMultiMappers Unique
--soloCellFilter None
--soloFeatures Gene GeneFull Velocyto
--soloCrMultimapRescue yes
--soloInlineHashMode no
--outSAMtype BAM Unsorted
--emitNoYBAM yes
--emitYNoYFastq yes
--emitYNoYFastqCompression gz
--outSAMattributes NH HI AS nM NM GX GN
--ocmMultiEnable auto
--ocmMultiConfig /mnt/pikachu/JAX_scRNAseq02/cellranger-logs/config.csv
--ocmMultiBarcodeMode flex
--ocmMultiBamSplit yes
```

`--ocmMultiEnable auto` is used for the comparator smoke because
`--soloCellFilter None` intentionally defers per-sample EmptyDrops to the
native OCM materializer when no pool filtered MEX is present. Production runs
should keep this split-before-ED behavior and use `--ocmMultiEnable yes` or
`auto` with a valid OCM config.

The existing Cell Ranger control uses Cell Ranger 9.0.1 with:

```text
reference,/mnt/pikachu/CR-references/refdata-gex-GRCh38-2024-A
create-bam,true
include-introns,false
```

The parity comparison intentionally uses STAR `Gene`, not `GeneFull`, because
the Cell Ranger control config has `include-introns,false`. GeneFull and
Velocyto production outputs are written by native STAR under
`star_composite/outs` and `star_composite/samples`. Python GeneFull/Velocyto
materialization is only a fallback when native outputs are absent.

Validated 100K STAR-only smoke after the CB/UB read-id guard:

```text
/tmp/ocm_flex_star_smoke_20260520T174046Z
STAR: completed
Native OCM BAM split: completed
Gene materialization: completed
GeneFull + Velocyto materialization: completed
```

The read-id guard is important for this surface: `trackReadIdsForTags` must not
force CB/UB writeback storage on a feature stream that does not carry read
indices. Otherwise a stride-2 matrix stream can be interpreted as if it had a
read-id column, creating impossible read-id to CB conflicts during BAM tag
injection.

Validated native C++ OCM materialization smoke:

```text
/mnt/pikachu/JAX_scRNAseq02_processed/ocm_native_cpp_100k_20260521T003616Z
STAR-only elapsed: 2m22s, no Cell Ranger rerun
Native OCM materialization: 00:38:16-00:38:37 UTC
CR reference reused: /mnt/pikachu/JAX_scRNAseq02_processed/ocm_materialize_opt_100k_20260520T231602Z/cellranger/25E32-L3_ocm_composite_100000
```

Validated 50M timing/parity gate:

```text
/mnt/pikachu/JAX_scRNAseq02_processed/ocm_native_cpp_50M_20260521T003926Z
STAR-only elapsed: 13m44s, no Cell Ranger rerun
Mapping: 00:40:28-00:42:43 UTC
Solo counting: 00:43:32-00:45:02 UTC
Native OCM materialization: 00:45:12-00:47:24 UTC
Tagged BAM writeback: 00:47:24-00:53:09 UTC
CR reference reused: /mnt/pikachu/JAX_scRNAseq02_processed/ocm_composite_smoke_50m_bam_20260520T164318Z/cellranger/25E32-L3_ocm_composite_50000000
```

This timing run was generated before the composite harness was corrected to
default `STAR_YREMOVE=yes`. Treat it as the native OCM materialization timing
gate only. JAX OCM production and any final smoke intended to feed delivery
artifacts must use `--emitNoYBAM yes --emitYNoYFastq yes
--emitYNoYFastqCompression gz`, which the wrapper now emits by default unless
`STAR_YREMOVE=no` or `--star-yremove no` is set explicitly.

Validated 50M native streaming Y/noY gate:

```text
/mnt/pikachu/JAX_scRNAseq02_processed/ocm_native_cpp_streambam_sample_yremove_50M_20260521T024938Z
STAR job: 02:49:38-03:01:09 UTC
Mapping: 02:50:44-02:57:21 UTC
Solo counting: finished 02:58:49 UTC
Harness done: 03:02:34 UTC
Native OCM BAM split: 4 samples, with per-sample Y/noY BAMs
CB/UB BAM tags: disabled
```

The 50M gate above was generated before the wrapper default was tightened to
omit `GX/GN`. Production OCM runs should keep BAM attributes alignment-only
(`NH HI AS nM NM`) unless a diagnostic run explicitly needs legacy gene tags.
STAR now warns when `GX/GN` are requested with `GeneFull` and/or `Velocyto`
but without `Gene`, because those tags are alignment-level annotations and are
not the final GeneFull/Velocyto UMI-collapsed counting policy.

50M Gene comparator parity:

| sample | STAR cells | CR cells | Jaccard | barcode UMI Pearson | feature UMI Pearson |
| --- | ---: | ---: | ---: | ---: | ---: |
| GCM1-Day-4 | 2751 | 2733 | 0.993457 | 0.999977 | 0.999564 |
| GRHL1-Day-4 | 2896 | 2900 | 0.998621 | 0.999958 | 0.999574 |
| OVOL1-Day-4 | 3219 | 3227 | 0.997521 | 0.999922 | 0.999614 |
| WT-PrS-20pct-Day-4 | 2214 | 2223 | 0.995951 | 0.999970 | 0.999598 |

Expected output locations:

```text
<out-root>/star_composite/run/
<out-root>/star_composite/run/Aligned.out_Y.bam
<out-root>/star_composite/run/Aligned.out_noY.bam
<out-root>/star_composite/outs/per_sample_outs/<sample_id>/count/sample_alignments.bam
<out-root>/star_composite/outs/per_sample_outs/<sample_id>/count/sample_alignments_Y.bam
<out-root>/star_composite/outs/per_sample_outs/<sample_id>/count/sample_alignments_noY.bam
<out-root>/star_composite/outs/multi/count/unassigned_alignments.bam
<out-root>/star_materialized/Gene/
<out-root>/star_materialized/GeneFull/                 # fallback only
<out-root>/parity_gene_vs_cr9.tsv
<out-root>/parity_gene_vs_cr9.json
```

Review gates:

- STAR finishes successfully and writes native OCM outputs.
- STAR writes pooled `run/Aligned.out_Y.bam` and `run/Aligned.out_noY.bam`
  when Y-removal is enabled.
- When `--ocmMultiBamSplit yes` is active, STAR writes four native per-OCM
  BAM groups under `star_composite/outs/per_sample_outs/`: unsplit,
  sample-level Y, and sample-level noY. The OCM sample and Y/noY state are
  both streaming routing determinants.
- STAR writes an unassigned OCM BAM under `star_composite/outs/multi/count/`.
- Cell Ranger reference artifacts are reused from an existing completed run.
  The harness only launches Cell Ranger with explicit `--run-cr`.
- Per-sample `sample_filtered_feature_bc_matrix` exists for STAR `Gene`.
- Per-sample GeneFull `raw_velocyto_feature_bc_matrix` and
  `filtered_velocyto_feature_bc_matrix` exist.
- `parity_gene_vs_cr9.tsv` reports per-sample cell counts, barcode overlap,
  UMI totals, barcode UMI Pearson, and feature UMI Pearson.

The native STAR OCM BAM split writes during direct unsorted BAM streaming when
final `CB/UB` tags are not requested, not by re-reading `Aligned.out.bam`:

```text
star_composite/
  run/Aligned.out.bam
  run/Aligned.out_Y.bam
  run/Aligned.out_noY.bam
  outs/per_sample_outs/<sample_id>/count/sample_alignments.bam
  outs/per_sample_outs/<sample_id>/count/sample_alignments_Y.bam
  outs/per_sample_outs/<sample_id>/count/sample_alignments_noY.bam
  outs/multi/count/unassigned_alignments.bam
```

`run/Aligned.out.bam` is the standard STAR pooled BAM. For JAX OCM production,
do not request barcode/UMI BAM tags or `GX/GN`; `STAR_BAM_CBUB_TAGS=no` and
`STAR_BAM_GXGN_TAGS=no` are the wrapper defaults and keep BAM emission
streaming with alignment-only attributes. Native per-sample BAMs are routed by
the OCM tag during that same stream, so they do not require a late effective
barcode tag or a Python post-split pass. When Y-removal is enabled, each
per-sample BAM is written as the unsplit BAM plus sample-level Y and noY BAM
outputs. If `CB/UB` tags are explicitly requested, STAR falls back to the
tagged replay path because final corrected tags are only available after
counting.

## 50M Gate

Only run the 50M gate after the 100K smoke completes and the output structure is
correct:

```bash
scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh \
  --read-pairs 50000000 \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_composite_smoke_50m_$(date -u +%Y%m%dT%H%M%SZ) \
  --run-all
```

`--run-all` reuses a completed Cell Ranger reference by default
(`CR_REUSE_RUN_DIR=auto`) and will fail rather than rerun Cell Ranger if no
matching reference exists. Keep the benchmark serialized. Do not run it in
parallel with other STAR jobs on this host.

## Interpreting Parity

This is a comparator, not a validator. STAR-suite and Cell Ranger differ in
alignment and multimapper policy, so exact equality is not expected.

Use the comparison as a method check:

- Barcode overlap should be meaningful within each OCM sample and should not
  show cross-tag rescue behavior.
- Per-feature Pearson should be interpreted on `Gene` for CR parity.
- GeneFull materialization should be inspected for production shape, not direct
  Cell Ranger equality.
- Large shifts in sample-specific cell counts should be reviewed before running
  full production.

## Native Implementation

The native STAR core implementation uses:

1. `--ocmMultiBarcodeMode flex` to select OCM-Flex barcode semantics.
2. The input `CB16` and original GEM-X whitelist; STAR derives the effective
   `CB16+OCM_TAG8` whitelist internally.
3. `--soloInlineCBCorrection yes`, required because the effective 24 bp
   barcode must use the 64-bit inline corrector rather than the legacy 32-bit
   whitelist hash.
4. Barcode correction, UMI collapse, Gene/GeneFull, and Velocyto on the
   effective barcode.
5. OCM splitting before EmptyDrops.
6. CR-compatible per-sample MEX outputs with standard CB16 barcodes.
7. `--ocmMultiBamSplit yes` to stream BAM records into per-sample unsplit,
   sample-level Y, sample-level noY, and unassigned BAMs while the standard
   STAR BAM path writes the pooled BAM/Y/noY outputs.

The adapter scripts should remain in the repo as reproducible parity fixtures
and as a fallback for investigating native implementation regressions.

## Materialization Optimization

The current smoke harness still uses the Python adapter for Gene and
GeneFull/Velocyto materialization. The adapter path now uses a streaming
MatrixMarket router for per-sample raw/filtered MEX, parallel per-sample
EmptyDrops, and one GeneFull+Velocyto production materialization command.
`Gene` remains a comparator artifact and should be skipped for production
unless explicitly requested.

The optimization plan is tracked in
[RUNBOOK_SCRNA_OCM_MATERIALIZATION_OPTIMIZATION_20260520.md](/mnt/pikachu/STAR-suite/docs/RUNBOOK_SCRNA_OCM_MATERIALIZATION_OPTIMIZATION_20260520.md).
The preferred implementation is a streaming MatrixMarket router that shares OCM
column maps, one per-sample cell-calling result, and the GeneFull barcode axis
across Velocyto raw and filtered outputs.

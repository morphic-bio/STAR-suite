# Paper Benchmark Methodology

This document records how the benchmarks in the STAR Suite manuscript and in
the top-level `README.md` were run and evaluated. **Section 1 is the current
methodology (STAR Suite 1.9.5) and matches the manuscript's Methods.** Section 2
keeps the records of earlier investigations for history; wherever Section 2
disagrees with Section 1 (thresholds, comparator pipelines, datasets, script
defaults), Section 1 applies.

## 1. STAR Suite 1.9.5 manuscript methodology

### 1.1 Benchmark set

| Assay | Dataset | Reference | STAR Suite entry point |
|---|---|---|---|
| Bulk RNA-seq | PPARG revertant WT replicate 3, KOLF2.2J iPSC (35,089,336 read pairs; GEO GSE288287), without and with Y-chromosome removal | External stepwise pipeline (1.4) | `scripts/paper/run_pe_bulk_feature_benchmark.sh --integrated-only --threads 32 --quant-vb-component-parallel 1` (`--no-yremove` for the no-Y arm) |
| scRNA-seq | 10x PBMC 10K, 3' v3 (638.9M read pairs) | Cell Ranger 9.0.1 `count` run made for this work (the outputs published with the dataset are Cell Ranger 3) | Direct STAR invocation (1.6) |
| Perturb-seq | A375 1k CRISPR 5' GEM-X (gene expression + CRISPR guides) | Cell Ranger 9.0.1 `multi` | `scripts/paper/run_a375_benchmark.sh` |
| Perturb-seq | MSK 30-KO ES sample (gene expression + 30 CRISPR guides + 245,979-barcode LARRY library, one pass) | Sum of the Cell Ranger 9.0.1 GEX+guide and GEX+LARRY runs | `scripts/paper/run_msk_30polyko_benchmark.sh` on the ES staging set |
| 10x Flex | JAX SC2300771 (4 tags, 8 lanes, 2.011B read pairs) | Cell Ranger 9.0.1 `multi` run made for this work | Direct STAR invocation (1.6) |
| 10x Flex | GSE325982 pool (4 single-tag samples, 1.121B read pairs) | Cell Ranger 9.0.1 `multi` run made for this work, default probe filtering; a second run with the submitters' `filter-probes` false reproduces their GEO deposit byte for byte | Direct STAR invocation (1.6) |
| 10x Flex | 10x 320k scFFPE GEM-X (16 tags pooled in pairs into 8 samples, 7.303B read pairs) | Cell Ranger 9.0.1 `multi` run made for this work | Direct STAR invocation (1.6) |
| SLAM-seq | GRAND-SLAM 100K-read human fixture (external BED SNP mask) | GRAND-SLAM reference NTRs (`from_nosnp` oracle) | `tests/run_slam_fixture_parity.sh` |

The paper runs used copies of the `scripts/paper/` wrappers with these
recorded corrections: an explicit 32-thread budget and TranscriptVB
component-parallel setting, and a gate that refuses any input other than the
PPARG pair (bulk); explicit selection of the native BGZF reader, a check that
the comparator's matrix header is a Cell Ranger 9 run, and `--soloStrand Reverse`
(A375; now also in the repository script); ES defaults,
ES output labels and a 28-file input-inventory check (MSK); STAR run and EM
comparison switched on (they default off), 32 threads and an explicit SNP-mask
arm (SLAM), with parity scored against the canonical `from_nosnp` oracle.
All references are GRCh38-2024-A; Flex uses probe set v1.1.0 restricted to its
included probes, Cell Ranger's default, on all three datasets and for every tool.
GSE325982's submitters deposited outputs made with filtering off; running Cell
Ranger that way reproduces their matrices byte for byte, and the two settings
call identical cells and give per-cell and per-gene correlations of 1.000000 on
the genes both report.

### 1.2 Timing protocol

- Hardware: one server, Intel i9-13900KF (24 cores / 32 threads), 126 GiB RAM,
  32 threads for every tool. The 320k Flex dataset ran on a rented AWS instance
  limited to the same size (m6id.16xlarge at 32 vCPUs and 126 GiB of memory,
  two local NVMe devices in RAID0, Ubuntu 24.04), because its Cell Ranger run
  needs about 3 TB of scratch; Cell Ranger and cyto ran on the same instance.
- Every input, output, temporary, spill and index path is resolved (following
  symlinks) and asserted to be on local NVMe before the run; a run that fails
  the check is aborted, not timed.
- One run at a time on a quiet machine (no other STAR, Cell Ranger, cyto,
  Salmon, Trim Galore, compression or staging process; load average below 3),
  with the page cache dropped immediately before each run and a fresh output
  directory.
- Wall time is the whole process tree from `/usr/bin/time -v`.
- Cell Ranger 9.0.1 runs with `--create-bam=false` and `--nosecondary`
  (`create-bam,false` / `no-secondary,true` in `multi` configs), the same thread
  count, and is timed to completion, because its output directory is delivered
  only when the pipeline finishes. Where the same work needs two Cell Ranger
  runs (MSK 30-KO), the walls are summed. If a Cell Ranger run fails, the
  failure is recorded, the same pipestance is resumed with the identical command
  and the attempts are summed. One Cell Ranger run needed attention: the
  2026-09-23 GSE325982 run crashed after writing its per-sample matrices, and a
  clean re-run with the identical command is the run reported. The time at which
  Cell Ranger writes the per-sample matrices is recorded for reference only.
- CBQ (BINSEQ) input is converted once from the FASTQ; the conversion time is
  excluded for STAR Suite and cyto alike.

### 1.3 Flex input routes

- Sequencer-delivered BGZF FASTQ (JAX, 320k): the in-process parallel reader,
  `--readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1`.
- Plain single-stream gzip (GSE325982 as re-compressed by ENA):
  `--readFilesBgzfMode off --readFilesCommand rapidgzip -d -c -P 8`
  (rapidgzip 0.16.0, an external program, not part of STAR Suite).
- CBQ: `--readFilesType Binseq PE --readFilesCbqRangeMode range`.

### 1.4 Comparators

- **Bulk external pipeline** (one process tree, each step also timed):
  1. `trim_galore --paired --quality 20 --length 20 --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --cores 8 --fastqc --fastqc_args "--threads 32"` (Trim Galore 0.6.10, cutadapt 5.1, FastQC 0.11.9)
  2. upstream STAR 2.7.11b: `--runThreadN 32 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --readFilesCommand zcat`
  3. with Y removal only: the post-hoc `awk`/`samtools`/`gzip` Y split of both BAMs and both trimmed FASTQs
  4. `salmon quant -t transcriptome.fa -l A -a Aligned.toTranscriptome.out.bam -g tx2gene.tsv --gcBias -p 32` (Salmon 1.10.3)

  Salmon's automatic library-type detection reports ISR for this stranded
  library. An earlier recording used `-l IU`, which assigns antisense reads to
  opposite-strand genes; it is superseded.
- **Cell Ranger 9.0.1** as in 1.2.
- **CellGENI STARsolo** (scRNA-seq baseline): upstream STAR 2.7.11b with the 10x
  3' options of the CellGENI wrapper (`cellgeni/STARsolo` at `fbcd9ac`)
  verbatim, including `--soloMultiMappers EM`, `--outFilterScoreMin 30`,
  `--soloFeatures Gene GeneFull Velocyto`, `--clipAdapterType CellRanger4` and
  `--soloCellFilter EmptyDrops_CR`; its `GeneFull` filtered matrix is compared.
- **cyto 0.4.7** (Flex), pinned by absolute path. Its probe-to-gene tables are
  derived from the probe set CSV of the corresponding Cell Ranger run and follow
  that run's `filter-probes` setting, which is the default (included probes) for
  every reported run.

### 1.5 Parity metrics

- **Gene level.** Spearman correlation (tied ranks averaged) and Pearson
  correlation of raw, untransformed per-gene count totals over the cells called
  by both tools, over every gene present in both annotations, zero-count genes
  included. For scRNA-seq and Perturb-seq these are `spearman_all_genes` and
  `pearson_all_genes` in the `filtered_vs_filtered` block of the per-gene
  section of `scripts/report_additional_parity_metrics.py` (run with
  `--gene-corr-min-counts 20 --gene-corr-min-cells-pct 0.01`; those thresholds
  only affect the `*_filtered_genes` fields, which are not reported). For Flex
  the same quantities are computed from STAR Suite's per-sample MEX against
  Cell Ranger's `sample_filtered_feature_bc_matrix`, with all samples pooled
  (the quoted values) and per sample, averaged over samples (Supplementary
  Table S3), by `scripts/paper/concordance_levels.py`;
  `scripts/paper/flex_gene_correlation.py` gives the per-sample gene values
  alone. A Pearson correlation, or any correlation restricted to well-expressed
  genes, is dominated by the most abundant genes and hides differences in
  low-count genes, which is where counting rules differ most; both are reported.
- **Cell calls.** Jaccard index |A ∩ B| / |A ∪ B| on the called-cell sets. For
  Flex a cell is its 16-base barcode together with its sample tag.
- **Per cell.** Two different measures, both reported for every single-cell
  benchmark and both on the cells called by both tools. *Cell Pearson*: Pearson
  correlation, across cells, of each cell's total UMI count. *Mean per-cell
  Pearson*: for each cell, Pearson correlation across genes of its
  `log(1 + x)` counts in the two outputs, over genes with at least 20 counts in
  both outputs and detected in at least 1% of shared cells, averaged over cells
  (cells whose counts do not vary are skipped).
- **Multiplexed libraries (Flex).** Every measure is computed with all samples
  pooled, which ignores the samples and reflects cell calling and counting (the
  quoted values), and per sample, averaged over samples, which also reflects any
  disagreement in assigning reads to samples (Supplementary Table S3). Dropping
  the sample tag and summing counts per barcode is not used: in the 320k
  dataset 325,410 Cell Ranger cells share 128,394 barcodes, so that comparison
  is between partitions, not cells.
- **Feature barcodes.** Per-guide Pearson correlation of UMI sums over shared
  cells; CRISPR-call agreement is the fraction of shared cells whose called guide
  set is identical.
- **Bulk RNA-seq.** Assigned read counts, Salmon `NumReads` against TranscriptVB `NumReads`, at gene
  (`quant.genes.sf`) and transcript (`quant.sf`) level, all features, on the
  no-Y-removal benchmark.
- **SLAM-seq.** NTR Pearson against GRAND-SLAM at >= 20, 50 and 100 reads.
- **Deterministic components** are compared for exact identity: trimmed FASTQ
  hashed against Cutadapt output, BAM records compared after BGZF decoding and
  coordinate sorting.
- Parity is computed on the full, non-Y-removed read sets.

### 1.6 STAR Suite options

scRNA-seq (PBMC 10K), with `STAR_SOLO_NONFLEX_HASH_BRIDGE=1` in the environment:

```
--runThreadN 32 --readFilesBgzfMode auto --bgzfReaderThreads 0
--clipAdapterType CellRanger4 --clip3pPolyG yes --alignEndsType Local --chimSegmentMin 1000000
--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
--soloBarcodeReadLength 0 --soloCBwhitelist 3M-february-2018_TRU.txt
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
--soloUMIdedup 1MM_CR --soloMultiMappers Unique --soloCellFilter EmptyDrops_CR
--soloCbUbRequireTogether no --soloStrand Forward --soloFeatures GeneFull
--soloCrGexFeature genefull --soloCrMultimapRescue yes
--dynamicThreadInterface 1 --dynamicThreadConstMapPermits 32 --dynamicThreadTelemetry 1
--outSAMtype None --soloInlineHashMode yes
```

Perturb-seq, options shared by A375 and MSK (as recorded in the 1.9.5 runs'
`Log.out`):

```
--readFilesBgzfMode auto --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 1000000
--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
--soloBarcodeReadLength 0 --soloInlineHashMode yes --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
--soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloMultiMappers Unique
--soloCellFilter EmptyDrops_CR --soloCbUbRequireTogether no --soloFeatures GeneFull
--soloCrGexFeature genefull --soloCrMultimapRescue yes --pfMultiConfig multi_config.csv
--crAssignMaxHamming 1 --crAssignFeatureOffset 0 --crAssignMinCounts 0
--crAssignMaxBarcodeMismatches 5 --crAssignFeatureN 0 --crAssignBarcodeN 1
--crAssignSearchThreads 1 --crAssignSkipQcOutputs 1
--dynamicThreadInterface 1 --dynamicThreadConstMapPermits 32 --dynamicThreadTelemetry 1
```

| Option | A375 | MSK 30-KO ES |
|---|---|---|
| `--soloStrand` | `Reverse` (5' R2-only library) | `Forward` |
| `--clip3pPolyG` | not set | `yes` |
| `--crChemistry` | `TRU` | `auto` |
| `--crMinUmi` | `10` | `2` |
| feature reader | `--crAssignBgzfMode auto` | plain gzip input |

Flex (JAX shown; GSE325982 and 320k differ only in inputs, tag lists, probe
tables, cache file and input route):

```
--runThreadN 32 --flex yes --flexPipeline yes --flexNoAlign 1
--soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17 --soloCBlen 16 --soloUMIlen 12
--soloBarcodeReadLength 0 --soloCBwhitelist 737K-fixed-rna-profiling.txt
--soloSampleWhitelist <tags> --soloFlexAllowedTags <tags> --soloSampleProbes probe-barcodes-fixed-rna-profiling-rna.txt
--soloSampleProbeOffset 68 --soloProbeList <probe genes> --soloHashScreenFile <half-probe cache>
--soloFeatures Gene --soloCellFilter None --soloMultiMappers Rescue
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
--soloStrand Unstranded --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 0
--soloKeysCompat cr --soloBucketMode ram --soloBucketCount 256
--soloRunFlexFilter yes --soloFlexCellCaller tag-aware --soloCellFilterBootstrapThreads 32
--soloFlexEdFdrThreshold 0.01 --outSAMtype None --outSJtype None --dynamicThreadInterface 1
```

Strand follows the library: 3' libraries are `Forward`; 10x 5' R2-only
libraries (Cell Ranger chemistry `SC5P-R2*`) are `Reverse`, because read 2 is
antisense to the transcript; Flex probe reads are counted `Unstranded`. The A375 arm originally ran `Unstranded`, which
drops reads where opposite-strand genes overlap and counts antisense reads
(gene Spearman 0.952 against Cell Ranger); it was re-run with `Reverse` (0.988).

The EmptyDrops Monte Carlo draws are split across bootstrap streams, so
`--soloCellFilterBootstrapThreads` is part of the reproducibility contract:
matching it reproduces a call set exactly.

## 2. Earlier investigations (history)

The records below predate release 1.9.4. They document how the canonical
settings were found and are kept for reference; their thresholds, datasets,
comparator pipelines and numbers are superseded by Section 1.

### Root Cause of the A375 Regression

The A375 paper script (`scripts/paper/run_a375_benchmark.sh`) was missing
`--soloCrMultimapRescue yes`, which the UCSF and MSK scripts both had. This
single omission caused:

- 1.7M reads shifted from unique to multi-mapped (no rescue to recover them)
- 12% fewer UMIs on common barcodes (18.7M vs 21.4M)
- Gene Pearson dropped from 0.975 to 0.943

A secondary factor was a different genome index (autoindex `bulk_index` with
`sjdbOverhang 100` and `cellrangerLegacyGtfFilter No`) vs the original CellRanger 9
pre-built reference (`sjdbOverhang 90`, legacy GTF filter). The original CR9
reference at `/storage/A375-CR-9.01/` was deleted, and the `crstar` symlink index
at `/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star/` has
15 broken symlinks pointing to that deleted path.

#### Timeline

| Date | Event |
|---|---|
| 2026-03-15 | Original A375 benchmark: STAR vs local CR9 run, Gene Pearson **0.975** |
| 2026-03-17 | Paper script created, inadvertently dropped `--soloCrMultimapRescue yes` |
| 2026-03-17 | Re-run produced Gene Pearson **0.943** (also used autoindex instead of CR9 ref) |
| 2026-03-17 | Investigation traced root cause to missing rescue flag + different index |

#### Verification

Running the old STAR build (commit c30a106, same era as the 0.975 benchmark)
against the CR9 reference with the official parity script also yields ~0.943
Gene Pearson — confirming the script/index are the variables, not the STAR build.

### Canonical STAR Parameters (All Perturb-seq Benchmarks)

Every paper benchmark script MUST include the following CR-compat parameters.
Dataset-specific values are noted.

#### Shared across all datasets

```
--clipAdapterType CellRanger4
--alignEndsType Local
--chimSegmentMin 1000000
--soloType CB_UMI_Simple
--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
--soloBarcodeReadLength 0
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
--soloUMIfiltering MultiGeneUMI_CR
--soloUMIdedup 1MM_CR
--soloMultiMappers Unique
--soloCellFilter EmptyDrops_CR
--soloCbUbRequireTogether no
--soloFeatures GeneFull
--soloCrGexFeature genefull
--soloCrMultimapRescue yes          ← CRITICAL, do not omit
--crAssignMaxHamming 1
--crAssignFeatureOffset <dataset>   ← 0 for A375, -1 for UCSF
--crAssignLimitSearch -1
--crAssignMinCounts 0
--crAssignMaxBarcodeMismatches 5
--crAssignFeatureN 0
--crAssignBarcodeN 1
--crAssignConsumerThreads -1
--crAssignSearchThreads 1
--crAssignSkipQcOutputs 1
--dynamicThreadInterface 1
--dynamicThreadConstMapPermits <threads>
--dynamicThreadTelemetry 1
```

#### Dataset-specific parameters

| Parameter | A375 | UCSF EBs2_2 | MSK 30polyKO |
|---|---|---|---|
| `--soloStrand` | `Unstranded` (superseded: `Reverse`, Section 1.6) | `Forward` | `Forward` |
| `--clip3pPolyG` | omit (not NovaSeq) | `yes` | `yes` |
| `--crChemistry` | `TRU` | `auto` | `auto` |
| `--crOutputChemistry` | n/a (TRU native) | `TRU` | n/a |
| `--crWhitelist` | single-col TRU | 2-col NXT | per-library |
| `--crMinUmi` | `10` | `3` | `2` |
| `--crAssignFeatureOffset` | `0` | `-1` | `0` |

#### BAM tags are not parity outputs

The canonical paper recipes deliberately omit `--outSAMattributes` and use
`--outSAMtype None` unless a dataset-specific diagnostic requires a BAM. An
optional BAM mode does not change the comparison surface and must not add
`GX`/`UR` unless the experiment explicitly tests raw-tag compatibility.

| Surface | Meaning | Use for count parity? |
|---|---|---|
| `GX`/`GN` | Alignment-level gene compatibility annotation | No; it is not final UMI-collapsed `GeneFull`, Velocyto, multimap-rescue, or CR-compatible policy |
| `UR` | Raw, uncorrected UMI sequence | No; it is intentionally not updated after correction |
| `UB` | Corrected UMI, when BAM tag injection is requested | No; it can support BAM diagnostics, but the matrix remains authoritative |
| `GeneFull`/CR-compatible MEX and cell calls | Final counting and calling products | Yes |

Adding `GX` or `UR` to an otherwise canonical benchmark changes BAM generation
work without improving the count comparison. It therefore invalidates a
controlled performance comparison unless BAM-tag generation is itself the
declared test variable.

#### Why `--soloCrMultimapRescue yes` matters

CellRanger 9 uses a multimap rescue algorithm that reassigns multi-mapped reads
to unique genes when one mapping is exonic and the others are intergenic or
intronic. Without this flag, STAR leaves those reads as multi-mapped and they are
excluded from Solo counting (which defaults to `--soloMultiMappers Unique`).

For A375, rescue recovered ~1.77M reads in the original benchmark:
- Gene-vs-non-gene fast path: 817K
- Exonic winner (Phase 3): 949K
- Intronic fallback (Phase 3): 3K

### Canonical Parity Script

All paper parity numbers MUST use the same script and parameters:

```bash
python3 scripts/report_additional_parity_metrics.py \
  --cr-run <CR_DIR> \
  --star-run <STAR_DIR> \
  --gene-corr-min-counts 20 \
  --gene-corr-min-cells-pct 0.01 \
  --translate <none|both>
```

- `--translate none` for TRU-native datasets (A375) where no barcode translation
  is needed.
- `--translate both` for NXT datasets (UCSF EBs2_2) that use 2-column whitelists.
- `--gene-corr-min-counts 20` and `--gene-corr-min-cells-pct 0.01` are the
  standard gene-level filtering thresholds.

#### Gene filtering logic

The script uses a strict AND condition for gene filtering
(`summarize_gene_corr`, lines 398–403):

```
left_sum >= min_counts
AND right_sum >= min_counts
AND left_cell_n >= min_cells_abs
AND right_cell_n >= min_cells_abs
```

where `min_cells_abs = ceil(len(common_barcodes) * min_cells_pct)`.

Do NOT use ad-hoc scripts with OR-based filtering — those produce inflated
Gene Pearson values by including low-expression genes that are noise-dominated.

#### CellRanger reference outputs

For parity, the CR reference must be a **CellRanger 9** run, not the
10x-published reference downloads (which may be CellRanger 8.0.0 or earlier).

For A375, the 10x-published `sample_filtered_feature_bc_matrix` was generated by
CellRanger 8.0.0 (different from CR9). Using it as the reference produces
different parity numbers (0.943 Gene Pearson vs 0.975 against actual CR9 output).
Always verify the `software_version` in the MTX header.

### Genome Index

The original A375 benchmark used the CellRanger 9 pre-built reference at
`/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A/star/` (now deleted).

Current benchmarks use `/storage/autoindex_110_44/bulk_index` (sjdbOverhang 100,
gencode v44, `cellrangerLegacyGtfFilter No`). This index produces a ~3%
lower unique mapping rate (70.3% vs 73.6%) because it includes more gene models,
shifting reads from unique to multi-mapped. With `--soloCrMultimapRescue yes`,
most of these are recovered.

**Action needed**: The CR9 reference needs to be re-downloaded or the autoindex
needs to be validated as producing equivalent parity when rescue is enabled.

### PE Bulk Benchmark

The PE bulk benchmark (`scripts/paper/run_pe_bulk_feature_benchmark.sh`) compares
integrated STAR-suite (trimming + alignment + Y-removal + internal TranscriptVB)
against an external stepwise pipeline (trimvalidate + STAR + remove_y_reads +
Salmon).

#### Quantification contract

- The production STAR-suite arm is internal TranscriptVB, enabled with
  `--quantMode TranscriptVB --quantVBgcBias 1`.
- The integrated recipe must pass the same transcriptome FASTA used by Salmon via
  `--transcriptomeFasta`, so GC-bias and effective-length calculation use the
  pinned reference sidecar/input.
- External Salmon is part of the external stepwise production baseline, because
  that arm needs a separate quantifier after STAR emits `Aligned.toTranscriptome.out.bam`.
- Integrated TranscriptomeSAM emission and integrated Salmon QC are opt-in
  parity artifacts enabled with `--parity-qc`. They are not part of the
  STAR-suite production timing.
- Internal TranscriptVB collects transcript evidence during the STAR run, then
  merges/finalizes ECs, GC/effective-length state, and VB/EM convergence after
  alignment EOF. Do not describe this as external Salmon running concurrently.
- Headline speedups must be measured with this wrapper and the same output mode
  on both arms. The integrated total is the timed STAR-suite production command;
  the external total is decompress + trimvalidate + STAR TranscriptomeSAM +
  optional remove_y_reads + Salmon. Lean sanity checks using direct Trim Galore
  or unsorted BAM are valid diagnostics, but their ratios are not
  apples-to-apples with the paper benchmark.

#### Modes

- `--yremove` (default): Both arms perform Y-chromosome removal. Integrated uses
  `--emitNoYBAM yes --emitYNoYFastq yes`; external uses `remove_y_reads` on
  trimmed FASTQs.
- `--no-yremove`: Y-chromosome removal disabled on both arms. Useful for
  datasets/analyses where Y-removal is not relevant.
- `--integrated-only`: Run only the STAR-suite production arm. Use this for
  STAR-suite timing refreshes when the external control does not need to be
  rerun. In default mode this is equivalent to `--skip-external --skip-compare`
  and does not enable parity artifacts.

#### Parity comparison

When `--parity-qc` is enabled, transcript and gene-level Pearson/Spearman
correlations are computed by `tests/transcriptvb/compare_salmon_star.py`. The
benchmark produces a `comparison_metrics.tsv` with three comparisons per stage:
1. Integrated TranscriptVB vs Integrated Salmon (internal consistency)
2. Integrated Salmon vs External Salmon (pipeline effect on same tool)
3. Integrated TranscriptVB vs External Salmon (the headline comparison)

Do not run parity-QC mode for headline wall-time claims. It deliberately emits
and rereads an integrated transcriptome BAM that normal STAR-suite production
does not require.

#### 2026-06-27 PPARG timing reruns

The corrected production-mode STAR-suite-only PPARG no-Y rerun on `/storage`
measured integrated trim+align+sorted BAM+internal TranscriptVB at `8:54.52`.
It did not emit integrated TranscriptomeSAM and did not run integrated Salmon
QC, matching the normal STAR-suite production arm. Run root:
`/storage/JAX_PE/results/pparg_prod_benchmark_no_y_20260627_172349/`.

The matched external no-Y control was then completed in the same run root:
decompress `1:18.01`, trimvalidate `6:53.71`, STAR TranscriptomeSAM `6:56.55`,
and Salmon `1:01.34`, for an external total of `16:09.61`. The production
STAR-suite arm is therefore `1.81x` faster for this PPARG no-Y wrapper run.

An earlier PPARG no-Y sanity rerun before the production-mode wrapper correction
measured STAR-suite integrated trim+align+TranscriptVB plus TranscriptomeSAM
emission at `7:18.06`; Salmon QC on the integrated transcriptome BAM added
`1:00.08`.
A lean serial comparator using Trim Galore, upstream STAR, unsorted transcriptome
BAM, and Salmon took `9:50.17`, giving `1.35x` for that diagnostic setup. This
run confirmed high current Salmon parity (NumReads Pearson `0.999979` on all
transcripts; `0.999980` at sum>=1000) but should not replace the archived paper
speedup because its serial comparator is lighter than the paper wrapper and that
earlier integrated arm was still parity-artifact enabled.

### Paper Scripts

| Dataset | Script | Status |
|---|---|---|
| A375 | `scripts/paper/run_a375_benchmark.sh` | Fixed (rescue added) |
| UCSF EBs2_2 | `scripts/paper/run_ucsf_ebs2_2_benchmark.sh` | OK |
| MSK 30polyKO | `scripts/paper/run_msk_30polyko_benchmark.sh` | OK (but crstar index broken) |
| PE Bulk | `scripts/paper/run_pe_bulk_feature_benchmark.sh` | OK (supports --yremove/--no-yremove) |

### Broken Index: crstar

`/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star/` has
15 broken symlinks pointing to the deleted `/storage/A375-CR-9.01/` directory.
The MSK paper script defaults to this path and will fail until repaired.

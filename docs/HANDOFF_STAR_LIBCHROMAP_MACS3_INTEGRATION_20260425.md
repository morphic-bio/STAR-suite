# Handoff: STAR libchromap MACS3 FRAG Peak Integration

Date: 2026-04-25

Current production note (2026-05-18): this handoff is historical. The current
multiome production recipe is in
`docs/RUNBOOK_MULTIOME_MEX_MUDATA_20260516.md`. STAR now writes the Chromap BAM
plus binary fragment sidecar locally; peaks, ATAC peak MEX, and ATAC metrics
are built after STAR with
`core/features/libchromap_contract/star_multiome_atac_peak_mex`. Do not use the
removed STAR CLI flag `chromapAtacMacs3FragPeaksSource` in new recipes.

## Goal

Complete the STAR-suite integration of the hardened `libchromap` API and the
new MACS3-compatible FRAG peak-calling path. The intended first production shape
is still file-oriented: STAR constructs Chromap parameters, Chromap writes BAM or
CRAM plus fragments, and optional MACS3 FRAG peaks are emitted as narrowPeak and
summits files.

## Repositories Involved

- `Chromap-suite`: `/mnt/pikachu/Chromap-suite`
- `STAR-suite`: `/mnt/pikachu/STAR-suite`

Both repositories currently have pre-existing dirty state. Do not assume every
dirty file was changed by this integration pass.

## What Is Already Implemented

### Chromap-suite

`libchromap` is now the source of truth for Chromap mapping execution:

- `chromap` links through `libchromap.a`.
- `libchromap.a` includes the peak-caller objects needed by the MACS3 FRAG
  workspace and pipeline.
- `libchromap::RunMapping()` owns mapping dispatch plus post-mapping MACS3 FRAG
  peak orchestration.
- `libchromap::ValidateMappingParameters()` is public and validates expected
  configuration errors before mapping starts.
- `chromap_driver.cc` and `chromap_lib_runner.cc` call the shared validator
  instead of owning separate validation paths.

Relevant dirty files in `Chromap-suite`:

- `Makefile`
- `src/chromap_driver.cc`
- `src/chromap_lib_runner.cc`
- `src/libchromap.cc`
- `src/libchromap.h`

Known untracked files that were present and should not be committed accidentally:

- `tests/__pycache__/`
- `tests/test_frag_compact_store`

### STAR-suite

The existing STAR-side contract was used rather than creating a new adapter:

- `core/features/libchromap_contract/include/star_chromap_contract.h`
- `core/features/libchromap_contract/src/star_chromap_contract.cpp`
- `core/features/libchromap_contract/tools/star_libchromap_contract_runner.cpp`
- `core/legacy/source/star_chromap_orchestration.cpp`
- `core/legacy/source/Parameters.h`
- `core/legacy/source/Parameters.cpp`
- `docs/LIBCHROMAP_CONTRACT.md`

The contract now exposes optional MACS3 FRAG peak-calling fields:

- `call_macs3_frag_peaks`
- `macs3_frag_peaks_output`
- `macs3_frag_summits_output`
- `macs3_frag_keep_intermediates_dir`
- `macs3_frag_pvalue`
- `macs3_frag_qvalue` plus q-value threshold mode
- `macs3_frag_min_length`
- `macs3_frag_max_gap`
- `macs3_frag_uint8_counts`
- `macs3_frag_peaks_source` (`file` or `memory`; contract-runner/internal
  validation only, not a current STAR CLI flag)

STAR runtime parameters were added with the `chromapAtac*` naming style:

- `chromapAtacCallMacs3FragPeaks`
- `chromapAtacMacs3FragPeaksOutput`
- `chromapAtacMacs3FragSummitsOutput`
- `chromapAtacMacs3FragKeepIntermediates`
- `chromapAtacMacs3FragPvalue`
- `chromapAtacMacs3FragQvalue` (`0` disables q-value mode)
- `chromapAtacMacs3FragMinLength`
- `chromapAtacMacs3FragMaxGap`
- `chromapAtacMacs3FragUint8Counts`

The contract runner now accepts matching flags:

- `--call-macs3-frag-peaks`
- `--macs3-frag-peaks-output`
- `--macs3-frag-summits-output`
- `--macs3-frag-pvalue`
- `--macs3-frag-qvalue`
- `--macs3-frag-min-length`
- `--macs3-frag-max-gap`
- `--macs3-frag-peaks-source file|memory` (contract-runner validation only)
- `--macs3-frag-keep-intermediates`
- `--macs3-frag-no-uint8-counts`

## Validation Already Run

From `Chromap-suite`:

```bash
make -C /mnt/pikachu/Chromap-suite chromap chromap_callpeaks chromap_lib_runner
make -C /mnt/pikachu/Chromap-suite test-peak-integration-100k test-peak-memory-source-100k
git -C /mnt/pikachu/Chromap-suite diff --check
```

Negative `chromap_lib_runner` checks were run for:

- missing reference
- mismatched read counts
- barcode whitelist without barcode reads
- `-t 0`

All returned structured validation failures.

From `STAR-suite`:

```bash
make -C /mnt/pikachu/STAR-suite star-libchromap-contract
make -C /mnt/pikachu/STAR-suite/core/legacy/source STAR WITH_CHROMAP=1
git -C /mnt/pikachu/STAR-suite diff --check
```

Negative `star_libchromap_contract_runner` checks were run for:

- `--call-macs3-frag-peaks` without peak outputs
- peak calling without secondary fragments
- invalid `--macs3-frag-peaks-source`

All returned validation errors. The `WITH_CHROMAP=1` STAR build completed.

## Important Current Limitation

The peak path is wired, but it has not yet had a positive end-to-end fixture run
through STAR with `chromapAtacCallMacs3FragPeaks 1`.

Also, the current first supported shape for peak calling is dual ATAC output:

- primary output: BAM or CRAM
- secondary fragments: `chromapAtacSecondaryFragments`
- optional peaks: narrowPeak and summits

Do not assume BED-only primary fragment output plus peak calling is validated.

## Next Step: E2E Validation

Run the existing 100K multiome/ATAC fixture through the contract runner first.
Use the paths documented in `docs/LIBCHROMAP_CONTRACT.md` as the starting point,
or the known fixture roots referenced there:

- 100K benchmark root:
  `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k`
- Prior validated libchromap/Chromap reference:
  `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/libchromap_parity_20260424_041454`

Historical validation order from the original handoff:

1. Contract runner, mapping plus secondary fragments, no peaks.
2. Contract runner, mapping plus secondary fragments, peaks enabled with
   `--macs3-frag-peaks-source file`.
3. Contract runner, same but `--macs3-frag-peaks-source memory`.
4. STAR `WITH_CHROMAP=1`, `postMapping`, peaks enabled.
5. STAR `WITH_CHROMAP=1`, `concurrent`, peaks enabled.

Current production validation order:

1. STAR `WITH_CHROMAP=1`, `concurrent`, sorted BAM plus binary sidecar.
2. Run `star_multiome_atac_peak_mex` with `--sidecar <out>/atac_fragments.bin`
   and `--call-peaks-from-sidecar` to produce peaks, ATAC peak MEX, and metrics.
3. Remote post-MEX MuData construction and layer validation.

Compare:

- fragment tuple identity against existing Chromap CLI/libchromap reference
- narrowPeak and summits against `chromap_callpeaks` or integrated Chromap output
- file-source peaks vs memory-source peaks
- STAR GEX outputs still complete successfully

## Suggested Contract Runner Command Shape

Fill in fixture-specific paths:

```bash
/mnt/pikachu/STAR-suite/core/features/libchromap_contract/star_libchromap_contract_runner \
  --ref <genome.fa> \
  --index <chromap.index> \
  --read1 '<R1.fastq.gz,...>' \
  --read2 '<R2.fastq.gz,...>' \
  --barcode '<BC.fastq.gz,...>' \
  --barcode-whitelist <whitelist.txt> \
  --output <out>/atac_possorted.bam \
  --output-format BAM \
  --sort-bam \
  --atac-fragments <out>/atac_fragments.bin \
  --summary <out>/chromap_summary.tsv \
  --threads 8 \
  --call-macs3-frag-peaks \
  --macs3-frag-peaks-output <out>/atac_peaks.narrowPeak \
  --macs3-frag-summits-output <out>/atac_summits.bed \
  --macs3-frag-peaks-source memory
```

## Suggested STAR Command Shape

Historical STAR peak-calling shape from this handoff used in-STAR peaks. New
production runs should instead use the sidecar post-materialization command
shown after this block.

Use the same fixture paths and add:

```text
--chromapAtacEnable 1
--chromapAtacStartMode postMapping
--chromapAtacReferenceFasta <genome.fa>
--chromapAtacIndex <chromap.index>
--chromapAtacRead1 <R1.csv>
--chromapAtacRead2 <R2.csv>
--chromapAtacBarcode <BC.csv>
--chromapAtacBarcodeWhitelist <whitelist.txt>
--chromapAtacOutputFormat BAM
--chromapAtacOutputFragments <out>/atac_possorted.bam
--chromapAtacSecondaryFragments <out>/atac_fragments.bin
--chromapAtacSortBam 1
--chromapAtacThreads 8
--chromapAtacHtsThreads 2
--chromapAtacCallMacs3FragPeaks 1
--chromapAtacMacs3FragPeaksOutput <out>/atac_peaks.narrowPeak
--chromapAtacMacs3FragSummitsOutput <out>/atac_summits.bed
```

After `postMapping` succeeds, repeat with:

```text
--chromapAtacStartMode concurrent
```

Current production post-materialization shape:

```bash
core/features/libchromap_contract/star_multiome_atac_peak_mex \
  --sidecar <out>/atac_fragments.bin \
  --barcode-translate <refs>/atac2gex.tsv \
  --barcode-translate-from-first \
  --call-peaks-from-sidecar \
  --peaks <out>/atac_peaks.narrowPeak \
  --summits-out <out>/atac_summits.bed \
  --out-dir <sample>/atac/peak_mex \
  --metrics-tsv <sample>/atac/atac_metrics.tsv \
  --threads 16 \
  --temp-dir <sample>/star_sample/chromap_tmp
```

## EmptyCells / Multiome Cell Calling Work

The new MACS3 peak integration only produces peaks. It does not yet connect those
peaks to the multiome empty-cell / cell-calling evidence path.

Current downstream ATAC evidence code expects ARC-style per-barcode metrics such
as:

- `atac_peak_region_cutsites`
- `atac_peak_region_fragments`
- `atac_fragments`
- `atac_peak_fraction`

Relevant existing files:

- `core/features/libscrna/tools/scrna_build_atac_evidence.cpp`
- `core/features/libscrna/tools/scrna_multiome_combine.cpp`
- `scripts/run_multiome_cell_call_from_arc.sh`
- `scripts/run_multiome_cell_call_external_gex_from_arc.sh`

What remains for EmptyCells:

1. Add a small tool that consumes Chromap fragments plus MACS3 narrowPeak and
   emits per-barcode ATAC evidence with the columns above.
2. Validate this evidence against Cell Ranger ARC-style metrics on the 100K
   fixture.
3. Wire the tool into the existing multiome cell-calling scripts.
4. Only after that, consider integrating the evidence generation directly into
   the STAR/libchromap path.

Suggested new tool name:

```text
scrna_build_atac_evidence_from_peaks
```

Historical text-fragment inputs:

- `--fragments atac_fragments.tsv.gz`
- `--peaks atac_peaks.narrowPeak`
- optional barcode whitelist / barcode suffix arguments

Current production inputs are instead the binary sidecar, peaks, and barcode
translation table consumed by `star_multiome_atac_peak_mex`.

Output:

- `atac_evidence.tsv`

Implementation detail:

- Count total fragments per barcode.
- Count fragments overlapping any MACS3 peak per barcode.
- Count cut sites within peaks per barcode if needed to match ARC semantics.
- Emit `atac_peak_fraction = peak_fragments / atac_fragments`.

## Risks / Caveats

- `libchromap` now handles expected validation failures, but deeper runtime
  failures in Chromap internals can still call `ExitWithMessage()`.
- `memory` peak source should be validated in STAR E2E specifically; it relies on
  the in-memory workspace being allocated and populated during mapping.
- Concurrent mode is already supported in STAR orchestration, but dynamic permit
  hooks are still intentionally rejected. Use fixed thread budgets.
- The STAR repo has many unrelated dirty files. Review `git status` carefully and
  stage only the libchromap contract / Chromap ATAC integration files for this
  work.

## Files To Review Before Continuing

Chromap:

- `/mnt/pikachu/Chromap-suite/src/libchromap.h`
- `/mnt/pikachu/Chromap-suite/src/libchromap.cc`
- `/mnt/pikachu/Chromap-suite/src/peak_caller/macs3_frag_peak_pipeline.cc`
- `/mnt/pikachu/Chromap-suite/src/peak_caller/macs3_frag_workspace.cc`

STAR:

- `/mnt/pikachu/STAR-suite/core/features/libchromap_contract/include/star_chromap_contract.h`
- `/mnt/pikachu/STAR-suite/core/features/libchromap_contract/src/star_chromap_contract.cpp`
- `/mnt/pikachu/STAR-suite/core/features/libchromap_contract/tools/star_libchromap_contract_runner.cpp`
- `/mnt/pikachu/STAR-suite/core/legacy/source/star_chromap_orchestration.cpp`
- `/mnt/pikachu/STAR-suite/core/legacy/source/Parameters.h`
- `/mnt/pikachu/STAR-suite/core/legacy/source/Parameters.cpp`
- `/mnt/pikachu/STAR-suite/docs/LIBCHROMAP_CONTRACT.md`

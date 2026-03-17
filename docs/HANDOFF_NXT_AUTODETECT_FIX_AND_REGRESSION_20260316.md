# Handoff: NXT Auto-detect Fix + Feature Assignment Regression (2026-03-16)

## Update After Review

The headline "feature assignment reads the full FASTQ data 10x" conclusion does
**not** hold up on re-check.

## Update After Fix

The actual regression was that raw pf-multi feature assignment had been moved
to run **inline before `soloMain.processAndOutput()`**, which removed the old
Phase A/Phase B overlap and made feature counting wait in the same serial
post-mapping section as Solo.

This is now restored:

- `STAR.cpp` launches `startPfFeatureAssignment(...)` during initialization,
  right after parameter parsing / pf preload
- raw feature assignment (`runPfMultiAssignPhase`) runs in a background thread
  again
- `processPfMultiConfig(...)` now joins the async handle only after Solo
  finishes, then runs the merge/finalize phase
- deferred `cr_assign/.../filtered` outputs still work after Solo, so the
  original dependency fix is preserved

Validation on the mixed-feature smoke:

- output dir:
  `/tmp/mixed_chem_async_regression_check_20260316_232038`
- `Log.out` timestamps:
  - `23:18:17` Phase A launched
  - `23:18:20` raw pf assignment finished
  - `23:18:37` mapping finished
  - `23:18:38` Solo finished
  - `23:18:38` pf merge/finalize started
- smoke script:
  `tests/multi_feature/run_mixed_chemistry_filtered_smoke.sh`
  passed with `11_passed / 0_failed`

One test expectation also needed updating: in the deferred post-Solo path, the
per-library filtered barcode subsets in this TRU-whitelist smoke are validated
in TRU namespace, so the smoke now checks `grna_de` filtered barcodes as `TRU`
instead of `NXT`.

What the review verified:

- the producer path in `assignBarcodes.c` still creates **one producer thread
  per lane** and each producer is initialized with `nfiles = 1`
- the benchmark log thread messages are consistent with **per-lane depth**, not
  rereading all lanes
- one representative guide lane,
  `/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/guides/EBs2_2_AALG2_S26_L001_R1_001.fastq.gz`,
  contains `26,874,904` reads by direct line count
- the benchmark log's final feature-processing total,
  `Finished processing 278480782 reads`, is therefore consistent with roughly
  ten ~27-28M-read guide lanes, not a 28M-read library reread ten times

So the apparent "thread regression" was caused by an incorrect assumption about
the total guide-library depth in the corrected `EBs2_2` fixture, not by a code
change in the producer/consumer file partitioning path.

The autodetect fixes below still stand. The remaining performance work should be
framed as optimization, not as a file-reading regression.

## Summary

This session fixed two bugs in the NXT/TRU auto-detect and barcode
normalization path, then discovered a performance regression in the feature
assignment pipeline where the production pass reads the full FASTQ data
10x (once per reader thread) instead of partitioning across threads.

## Context

Dataset: UCSF Perturb-seq `EBs2_2` (corrected organization — see
`docs/HANDOFF_UCSF_FASTQ_MISLABEL_INVESTIGATION_20260316.md`).

- GEX: TRU chemistry, 10 lanes, 445M reads, 27 GB
- Guides: NXT chemistry, 10 lanes, 28M reads, 12 GB
- 548 CRISPR features, ~13,700 filtered cells

CellRanger benchmark: 51 min wall clock.
STAR benchmark (with fixes): 66 min wall clock — **slower than CR**.

## Bug 1: Auto-detect Inversion for NXT Whitelists

**File**: `core/legacy/source/PfMultiProcess.cpp` (~line 1929)

When `--crChemistry auto` is used with a 2-column NXT whitelist
(`translation/3M-february-2018_NXT.txt`), `read_whiteList` loads only
column 1 (NXT barcodes) into the hash. For NXT guide reads, `chem_detect`
finds `RAW_MATCH` (raw barcodes match NXT hash). The old code then set
`translateNxt = false`, which is correct for a TRU whitelist but inverted
for an NXT one: NXT reads matching an NXT hash still need translation to
TRU for output.

**Fix applied**:

```cpp
bool needTranslate = (detectedMatchMode == "TRANSLATED_MATCH");
if (assignmentWhitelistNamespace == "NXT") {
    needTranslate = !needTranslate;
    P.inOut->logMain << "NOTICE: inverting translate decision for NXT-namespace whitelist"
                     << " (detect=" << detectedMatchMode
                     << ", translateNxt=" << needTranslate << ")\n";
}
runAssignOpts.translateNxt = needTranslate;
```

## Bug 2: Filtered Barcode Normalization Wrong Source Namespace

**File**: `core/legacy/source/PfMultiProcess.cpp` (~line 1810)

Solo GEX barcodes are always in TRU namespace (from `--soloCBwhitelist`
which is the 1-column TRU file). But `normalizeFilteredBarcodesForAssignNamespace`
assumed they were in `assignmentWhitelistNamespace` (NXT for 2-column files).
This caused TRU barcodes to be matched against NXT whitelist without translation,
producing massive unmatched counts (10,308/12,202 unmatched in the smoke test).

**Fix applied**:

```cpp
const string assignSourceNS =
    (sourceLabel == "solo") ? gexNormalizationChem
                            : assignmentWhitelistNamespace;
```

## Fix Validation

Smoke test (100K reads) before/after:

| Metric | Buggy | Fixed |
|--------|-------|-------|
| translateNxt | 0 | 1 |
| Filtered barcodes | 4 | 9,188 |
| Deduped counts | 19 | 57,272 |
| Normalization unmatched | 10,308 | 0 |

Full run (`EBs2_2` corrected, 2-column NXT whitelist):

| Metric | STAR | CellRanger |
|--------|------|------------|
| GEX filtered cells | 13,721 | 13,760 |
| GEX total counts (filtered) | 256M | 255M |
| Feature barcodes (filtered) | 11,345 | 13,370 |
| Feature counts (filtered) | 7.18M | 8.34M |
| Feature types | 548/548 | 548/548 |

Both fixes confirmed by `assignBarcodes.api_run.txt`:
```
translateNxt=1
barcode_normalization.translated_to_set=13721
barcode_normalization.unmatched=0
```

## Performance Regression: Feature Assignment Reads 10x Too Much Data

**This is the main issue for the next agent.**

STAR timing breakdown (EBs2_2 corrected, 32 threads, permits enabled):

| Phase | Duration |
|-------|----------|
| Genome load | 1 min |
| GEX mapping (445M reads) | 8 min |
| Solo (UMI collapse + EmptyDrops) | 8 min |
| **Feature assignment (28M reads)** | **48 min** |
| CRISPR GMM calling | 0.5 min |
| **Total** | **66 min** |

CellRanger does the same job in 51 min. STAR's GEX pipeline (17 min) is
3x faster than CR, but the feature assignment wipes out the advantage.

### Root cause

The production `runAssignBarcodes` call spawns 10 reader threads (one per
FASTQ lane pair). Each thread reads the **entire** 28M-record guide library
instead of its own lane (~2.8M records). Evidence from the log:

```
Thread 0 has read 1 million records.
...
Thread 0 has read 26 million records.
Thread 0 done reading
Thread 1 has read 1 million records.
...
Thread 1 has read 28 million records.
Thread 1 done reading
...  (same pattern for Thread 2-9)
```

20 total "done reading" events = 2 passes × 10 threads:
- Pass 1 (auto-detect probe): 10 threads, but `maxReads=10000` so exits fast
- Pass 2 (production): 10 threads × 28M each = **280M reads processed** for
  a 28M-read library

The probe pass correctly samples only 10K reads (checked: `reads_sampled=10000`
in the log) and exits quickly. The regression is in the production pass where
each reader thread reads the full dataset instead of its lane partition.

### Additional waste in probe pass

Even though the probe pass correctly limits chemistry detection to 10K reads,
it still generates **full heatmaps, histograms, and filtered outputs** for the
`.autodetect_probe` directory. These are thrown away immediately. Skipping
QC outputs for the probe pass would save time.

### Where to look

1. **Reader thread FASTQ distribution**: `core/features/process_features/src/assignBarcodes.c`
   around the `process_sample_fastqs` / producer thread setup (~line 5121-5204).
   The `sample_size` variable (number of FASTQ file pairs) is used for both
   producer thread count and data distribution. Check whether each producer
   thread is being given all files or just its assigned lane.

2. **Probe pass QC suppression**: The probe call at `PfMultiProcess.cpp:1892`
   uses `detectOpts` which should disable all output except the chemistry
   decision. Check if heatmap/histogram generation is gated on
   `autodetectChemistry` mode.

3. **Heatmap skip flag**: `mex_writer.h` has `write_mex_core_only()` that
   skips heatmaps, but there is no CLI flag to invoke it. Adding
   `--crSkipHeatmaps` or similar would help benchmarks.

### Relevant files

| File | What's there |
|------|-------------|
| `core/legacy/source/PfMultiProcess.cpp` | Two `runAssignBarcodes` calls: probe (line 1892) and production (line 2190); both bug fixes applied here |
| `core/legacy/source/PfMultiAssign.cpp` | Bridges C++ options to C `pf_config`; `maxReads` wiring (line 311) |
| `core/legacy/source/PfMultiAssign.h` | `AssignOptions` struct: `maxReads`, `autodetectChemistryReads=10000` |
| `core/features/process_features/src/assignBarcodes.c` | Reader/producer thread setup (~5121), consumer loop (~4605), `max_reads` gate (~3851), `chem_detect` sampling (~4716) |
| `core/features/process_features/include/mex_writer.h` | `write_mex_core_only()` and `mex_write_heatmaps()` |

### Expected fix

After fixing the reader thread regression, the production pass should
process 28M reads total (not 280M), bringing feature assignment down to
~5-8 minutes. Combined with GEX (17 min), total STAR time should be
~22-25 minutes — approximately 2x faster than CR's 51 minutes, consistent
with the benchmarks in the README.

## Benchmark Artifacts

| Artifact | Path |
|----------|------|
| STAR run (with fixes, permits) | `/storage/bench_EBs2_2_corrected_20260316/star/run/` |
| STAR run script | `/storage/bench_EBs2_2_corrected_20260316/run_star.sh` |
| STAR log (parallel run) | `/storage/bench_EBs2_2_corrected_20260316/star/star_run_parallel.log` |
| CellRanger run | `/storage/bench_EBs2_2_corrected_20260316/cr/EBs2_2_corrected/` |
| Parity report (sequential run) | `/storage/bench_EBs2_2_corrected_20260316/parity_fixed.log` |
| Corrected data | `/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/` |
| Smoke test (100K) | `/tmp/star_smoke_nxt/` |

## Pending After Regression Fix

1. Rerun STAR benchmark with the reader thread fix — expect ~22-25 min
2. Run parity comparison with feature correlation on common barcode set

## Update: feature QC suppression wired

The previously-unused lightweight feature-output path is now exposed in the
runtime as:

- STAR CR-compat: `--crAssignSkipQcOutputs 1`
- standalone `process_features`: `--skip_qc_outputs`

Current behavior:

- keeps core feature outputs (`barcodes.tsv`, `features.tsv`, `matrix.mtx`,
  `stats.txt`, `feature_per_cell.csv`)
- suppresses Plotly heatmaps and HTML histogram QC outputs
- still leaves the plain-text `deduped_counts_histograms.txt` sidecar, because
  that is written outside `mex_writer`

Validation artifact:

- `/tmp/msk_skip_qc_smoke_20260316_222408/`
  - `cr_assign/Custom/larry_de/assignBarcodes.api_run.txt` records
    `skipQcOutputs=1`
  - `cr_assign/Custom/larry_de/LARRY/` contains core matrix/stat outputs
  - no `Feature_*heatmap.html/json` or `feature_*histogram.html` files are
    emitted there
3. Update `README.md` benchmark table with corrected dataset numbers
   (replacing the old mislabeled-data UCSF row)
4. Update `docs/feature_barcodes.md` optimized parameters section to
   mandate the 2-column NXT whitelist for `--crChemistry auto`

## Build

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

Binary with both fixes: `core/legacy/source/STAR` (compiled 2026-03-16T17:09:25).

## Update: CPU-aware PF permit controller

The dynamic PF controller now has an optional CPU-aware secondary signal layered
on top of the existing ETA/chunked permit controller.

New STAR parameters:

- `--dynamicThreadPfControllerCpuAware 0|1`
- `--dynamicThreadPfControllerCpuSampleMs <ms>`
- `--dynamicThreadPfControllerCpuEmaAlpha <0..1>`
- `--dynamicThreadPfControllerCpuIdleThreshold <0..1>`
- `--dynamicThreadPfControllerCpuBusyThreshold <0..1>`

Implementation summary:

- `/proc/stat` deltas are sampled in `ThreadControl`
- the busy fraction is smoothed with an EMA
- `PfMultiProcess.cpp` treats this as a conservative secondary nudge:
  - sustained idle can raise the permit target by `+1`
  - sustained saturation can lower it by `-1`
  - only after two stable controller ticks, to avoid permit thrash
- the signal is only allowed for `eta` / `chunked` controller modes

Validation:

- clean rebuild succeeded:
  `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`
- mixed-feature CR-compat smoke passed with CPU-aware mode enabled:
  `/tmp/mixed_chem_cpuaware_smoke_20260316_235655`
- filtered-barcode validator passed:
  `/tmp/mixed_chem_cpuaware_smoke_20260316_235655/filtered_barcode_report.tsv`
- async overlap still holds in `Log.out`:
  - `pf-feature-assign Phase A launched (async)`
  - `finished mapping`
  - `finished Solo counting`
  - `started pf-multi merge/finalize`

Important limitation from the smoke:

- on this tiny fixture, the PF assign stages finished before enough controller
  runtime accumulated for a `500 ms` CPU sample, so
  `assignBarcodes.api_run.txt` records the CPU-aware settings but
  `cpuInitializedFinal=0`
- this is expected for short-lived PF stages and does not indicate a failure in
  the core logic

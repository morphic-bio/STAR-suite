# PE Bulk Benchmark Results (2026-03-10)

## Scope

This note captures the first complete rerun of the PE bulk benchmark after the
two PE regressions were fixed:

- integrated Y/noY FASTQ routing
- TranscriptVB PE `A` autodetect collapsing to `IU`

It also records the benchmark-side Salmon setting change required for this
dataset:

- use explicit `salmon -l IU`
- do not use `salmon -l A` on this PE benchmark sample

The benchmark compares:

- integrated STAR-suite path
  - internal adapter trimming
  - Y chromosome removal
  - TranscriptVB quantification
  - transcriptome BAM generation
- external stepwise baseline
  - `trimvalidate`
  - STAR alignment + transcriptome BAM
  - `remove_y_reads`
  - Salmon alignment-mode QC

## Benchmark Command

Script:

- `scripts/paper/run_pe_bulk_feature_benchmark.sh`

Run:

```bash
scripts/paper/run_pe_bulk_feature_benchmark.sh \
  --threads 32 \
  --out-root /tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337
```

Output root:

- `/tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337`

Run manifest:

- `/tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337/RUN_MANIFEST.txt`

Reference:

- STAR index: `/storage/autoindex_110_44/bulk_index`
- transcriptome FASTA: `/storage/autoindex_110_44/bulk_index/transcriptome.fa`
- GTF: `/storage/autoindex_110_44/bulk_index/cellranger_ref/genes.gtf`

## Important Benchmark Setting

For this PE benchmark sample, Salmon alignment-mode autodetect misclassifies
the library as stranded when invoked with `-l A`.

The paper benchmark therefore pins:

- `Salmon libtype: IU`

This is now the default in
`scripts/paper/run_pe_bulk_feature_benchmark.sh`.

## Timing Summary

All timings below are wall clock times from the benchmark summaries.

### Downsampled

- integrated STAR: `23.83s`
- integrated Salmon QC: `6.38s`
- integrated total: `30.21s`
- external decompress: `0.21s`
- external trimvalidate: `1.17s`
- external STAR: `14.05s`
- external remove_y_reads: `1.13s`
- external Salmon QC: `6.36s`
- external total: `22.92s`

### Full sample on `/storage`

- integrated STAR: `29.20s`
- integrated Salmon QC: `31.51s`
- integrated total: `60.71s`
- external decompress: `13.44s`
- external trimvalidate: `16.89s`
- external STAR: `49.94s`
- external remove_y_reads: `14.03s`
- external Salmon QC: `31.03s`
- external total: `125.33s`

### Full sample copied to `/mnt/pikachu`

- integrated STAR: `37.25s`
- integrated Salmon QC: `31.05s`
- integrated total: `68.30s`
- external decompress: `13.26s`
- external trimvalidate: `17.00s`
- external STAR: `48.61s`
- external remove_y_reads: `13.81s`
- external Salmon QC: `31.52s`
- external total: `124.20s`

Notes:

- the `pikachu` stage reused existing copied FASTQs, so copy timing is `NA`
- HDD-local input mainly slowed the integrated STAR stage
  - `/storage`: `29.20s`
  - `/mnt/pikachu`: `37.25s`

## Y/noY FASTQ Parity

After the Y-route fix, integrated and external Y/noY FASTQ counts match exactly
on all three stages.

### Downsampled

- integrated Y: `121`
- external Y: `121`
- integrated noY: `98550`
- external noY: `98550`

### Full sample on `/storage`

- integrated Y: `8265`
- external Y: `8265`
- integrated noY: `6422076`
- external noY: `6422076`

### Full sample on `/mnt/pikachu`

- integrated Y: `8265`
- external Y: `8265`
- integrated noY: `6422076`
- external noY: `6422076`

## Quant Comparison Summary

### Integrated TranscriptVB vs integrated Salmon QC

These comparisons use explicit `Salmon -l IU`.

#### Downsampled

- transcript-level Spearman all: `0.968841`
- transcript-level Pearson all: `0.990877`
- transcript-level Spearman expressed: `0.994690`
- transcript-level Pearson expressed: `0.996190`
- total reads:
  - Salmon: `84124.011`
  - TranscriptVB: `84157.001`
  - absolute difference: `32.990`

#### Full sample on `/storage`

- transcript-level Spearman all: `0.975818`
- transcript-level Pearson all: `0.994534`
- transcript-level Spearman expressed: `0.993772`
- transcript-level Pearson expressed: `0.994874`
- total reads:
  - Salmon: `5482280.993`
  - TranscriptVB: `5485631.025`
  - absolute difference: `3350.032`
- gene-level Spearman all: `0.993644`
- gene-level Pearson all: `0.996834`

#### Full sample on `/mnt/pikachu`

- transcript-level Spearman all: `0.975800`
- transcript-level Pearson all: `0.994512`
- transcript-level Spearman expressed: `0.993755`
- transcript-level Pearson expressed: `0.994881`
- total reads:
  - Salmon: `5482281.099`
  - TranscriptVB: `5485631.117`
  - absolute difference: `3350.018`
- gene-level Spearman all: `0.993550`
- gene-level Pearson all: `0.996826`

Interpretation:

- the previous large benchmark-side divergence is gone
- the remaining TranscriptVB vs Salmon difference is small and stable
  - about `3350` reads over `~5.49M`
  - about `0.06%`

### Integrated Salmon QC vs external Salmon QC

These rows isolate the STAR-suite integrated path from the external stepwise
path while holding the Salmon model fixed.

#### Downsampled

- total read difference: `0.020`
- Spearman all: `0.975058`
- Pearson all: `0.995285`

#### Full sample on `/storage`

- total read difference: `0.010`
- Spearman all: `0.996187`
- Pearson all: `0.999921`

#### Full sample on `/mnt/pikachu`

- total read difference: `0.086`
- Spearman all: `0.996489`
- Pearson all: `0.999932`

Interpretation:

- integrated and external Salmon QC are effectively identical on the corrected
  benchmark
- the residual difference is therefore between TranscriptVB and Salmon, not
  between integrated STAR-suite preprocessing and the external stepwise path

## Primary Conclusion

With the PE regressions fixed and Salmon pinned to `IU` for this benchmark:

- integrated Y/noY FASTQ output matches the external reference split exactly
- integrated Salmon QC matches external Salmon QC essentially exactly
- integrated end-to-end runtime is much faster than the external stepwise path
  on the full sample
- the remaining TranscriptVB vs Salmon residual is small enough to investigate
  separately as a model/allocation difference rather than a pipeline regression

## Source Artifacts

Stage summaries:

- `/tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337/downsampled/BENCHMARK_SUMMARY.txt`
- `/tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337/storage/BENCHMARK_SUMMARY.txt`
- `/tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337/pikachu/BENCHMARK_SUMMARY.txt`

Comparison tables:

- `/tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337/downsampled/comparison/comparison_metrics.tsv`
- `/tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337/storage/comparison/comparison_metrics.tsv`
- `/tmp/pe_bulk_feature_benchmark_full_iu32_20260310_033337/pikachu/comparison/comparison_metrics.tsv`

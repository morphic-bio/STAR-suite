# Production half-probe khash: L004 validation, 2026-09-10

The production loader now accepts a compact persisted H0 + two 25-base
Hamming-1 khash representation. A full-probe match intersects parent-ID lists;
a unique half-match carries its parent directly into the existing Hamming
scorer. No large full-variant hash or legacy seed index is built on this path.

The biological algorithm, model feature universe, gene IDs, sample/tag policy,
UMI collapse and cell caller are unchanged. This is a storage and lookup
regression benchmark. The earlier classifier microbenchmark is documented in
[FLEX_KHASH_HALF_PROBE_20260910.md](FLEX_KHASH_HALF_PROBE_20260910.md).

## Completed implementation checks

- Clean Chromap-enabled STAR build based on master `4c33c014`, with the
  uncommitted hash-storage changes captured in the artifact source manifest.
- Unit tests: legacy formats and stored general tables, half-list intersection,
  same-gene ambiguity, split probes, Hamming limit, N handling, sample-specific
  H0 lookup, compact table roundtrip, immutable output and invalid formats.
- Exhaustive conversion checks against 316,072,780 source records and every
  pair of generated half variants; no changed or newly introduced primary
  matches. Reopened compact tables match the verified arrays byte-for-byte.
- Full STAR CBQ and BGZF small fixtures: identical raw MEX values, eight cell
  calls, and aggregate read classifications against the established control.
  No reference index, BAM or decision sidecar.

## Model cache

| Property | Value |
|---|---:|
| Unique exact probes | 54,580 |
| Left half keys | 4,148,062 |
| Right half keys | 4,148,060 |
| Source records | 316,072,780 |
| General khash file | 20,615,036,928 bytes |
| Compact half file, including sample H0 mappings | 331,268,096 bytes |
| Convert, verify and write from stored general cache | 45.80 seconds |

Compact SHA256:
`801fc6143a383b6a94c55307f816bf824c7a9c685d3049405b2ffcc4d72db296`.

## L004 protocol

Serial executions on cloud instance `i-06de289faa5d78117`, 48 threads, with a
cold OS page cache before each timed run. Inputs are lane_003 CBQ or the L004
BGZF R2/R1 pair: 1,823,648,323 read pairs. The three arms are general stored
khash + CBQ control, compact half khash + CBQ, and compact half khash + BGZF.
There is one execution per arm; no repeat trials.

All use the deprecated-complete 19,068-feature model and 18,129-feature filtered
export, grouped samples, tag-aware cell calling and 100,000 simulations.
The genome directory is empty, and logs must confirm count-only no-genome.
BAM and decision sidecar outputs are disabled; caller diagnostics are enabled.

General-khash control STAR SHA256:
`dca77afeb1d6e57c36d98b24806c2d0962ffb0304a855a7ef0addd9e649e8790`.
Compact production STAR SHA256:
`6d9a3dd888248eeb7c79af021a0aad7562dd62219054831981a89e3323a4f3d1`.

## Completed L004 runtime results

| Arm | Total seconds | Peak RSS (GiB) | Called cells |
|---|---:|---:|---:|
| Stored full khash, CBQ | 524.44 | 93.36 | 249,194 |
| Stored half khashes, CBQ | 522.28 | 75.11 | 249,194 |
| Stored half khashes, BGZF | 1030.83 | 75.90 | 249,194 |

The CBQ runtime difference is 0.4%, which is not a material speed result from
one execution per arm. Peak memory fell by 18.25 GiB (19.55%). Startup was 21
vs 18 seconds and the logged read phase was 110 seconds in both CBQ runs.
The isolated matcher speedup does not imply the same end-to-end speedup.

The single-lane BGZF run used genuine BGZF input, on-demand BC/BSIZE range
readers and 48 inflater threads. Its read phase took 620 seconds. A matching
full-khash BGZF control was not run here, so its time must not be interpreted
as the effect of changing hashes.

All three completed arms have identical aggregate read classifications and
all 249,194 called cell identities. The complete raw MEX axes and coordinate
body are identical, including 499,298,264 nonzero entries and 809,101,584 UMIs
across 6,658,487 observed CB16+TAG8 combinations. All eight per-sample filtered matrices also match exactly. No full four-lane
320K execution is part of this validation.

The comparison script initially expected per-sample raw directories. This
configuration exports one global raw matrix and per-sample filtered matrices;
the script was corrected and resumed from the completed raw comparison. No
STAR run was repeated.

## Cold initialization

With the same production loader in isolated processes and a cold OS page
cache before each execution:

| Cache | Load seconds | Loader peak RSS (MiB) |
|---|---:|---:|
| full | 2.360480 | 6811.54 |
| half | 0.148786 | 297.80 |

Loading is about 15.9 times faster. These figures cover initialization, not
querying all reads. Whole STAR startup includes other initialization work,
so it decreased from 21 to 18 seconds in the CBQ runs. Total wall times also
include cleanup after STAR's final timestamp; phase timestamps have one-second
resolution.

## Concordance with the saved L004 Cell Ranger outputs

These values are identical for all three STAR arms. Full CB16+TAG8 identities
are used throughout. Expression correlations use shared called cells and the
18,129 included features; gene correlations require at least 20 UMIs and
presence in at least 1% of shared cells in both matrices.

| Sample | STAR cells | CR cells | Cell Jaccard | Cell UMI Pearson | Gene-total Pearson |
|---|---:|---:|---:|---:|---:|
| BreastCancer1_BC7-8 | 9,489 | 12,724 | 0.741513 | 0.99999464 | 0.99999940 |
| Colorectal_BC3-4 | 37,380 | 36,717 | 0.978241 | 0.99998998 | 0.99999996 |
| Endo_BC15-16 | 33,856 | 33,431 | 0.978855 | 0.99999320 | 0.99999971 |
| Glioblastoma_BC1-2 | 25,797 | 32,659 | 0.787591 | 0.99999383 | 0.99999928 |
| Kidney_BC11-12 | 27,523 | 27,440 | 0.992857 | 0.99999474 | 0.99999988 |
| LNReactive_BC9-10 | 21,872 | 18,173 | 0.821884 | 0.99999660 | 0.99999957 |
| LungCancer2_BC5-6 | 56,174 | 55,906 | 0.993313 | 0.99998838 | 0.99999987 |
| SkinMelanoma_BC13-14 | 37,103 | 36,823 | 0.981081 | 0.99999531 | 0.99999991 |
| POOLED | 249,194 | 253,873 | 0.934732 | 0.99998928 | 0.99999986 |

Pooled count-weighted Jaccard on shared cells is 0.994396. This storage change
preserves the caller's existing residual differences: BreastCancer1 and
Glioblastoma have fewer called cells than CR; LNReactive has 3,807 STAR-only
calls and 108 CR-only calls (cell Jaccard 0.821884). Those differences occur
in the control and both new runs, despite highly concordant expression counts.
No cell-caller parameters or decisions were changed by this patch.

Only CR result H5/MEX data and metadata were read. All eight filtered STAR and
CR matrices are also saved as compressed H5 CSR caches for subsequent analysis.

## Artifacts

- Local sources, build, unit tests, launch scripts and completed reports:
  `/mnt/pikachu/star_suite_paper/analysis/flex_half_production_20260910/`.
- Cloud outputs and compact model:
  `/scratch/flex_half_production_20260910_v1/`.
- S3 bucket `star-suite-320k-benchmark-alt-171440768238-us-west-2-20260904`,
  prefix `analysis-tools/flex_half_production_20260910_v1/`. The archive is
  complete: model, paired config, binary, filtered STAR/CR H5 caches, reports
  and `completed/l004_raw_mex.tar.gz` (1,562,377,804 bytes). Raw archive SHA256:
  `7b2c438781154efb460b9ee1f0999d6afbb33173894f2f71f8564f9f40453544`.
  `config/half_cache_manifest.json` records the matching gene-list hashes.
- Runtime selection: `--soloHashScreenFile model.half.khash` with the same
  `model_gene_ids.txt`. See [format and generation instructions](../FLEX_KHASH_CACHE.md).

The validated binary is installed at `core/legacy/source/STAR`. The previous
full-khash executable is preserved as `STAR.before_half` in the local artifact
directory. The code changes remain uncommitted.

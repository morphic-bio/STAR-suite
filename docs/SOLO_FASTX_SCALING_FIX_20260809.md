# STARsolo Fastx Scaling Fix (2026-08-09)

## Status

This note records the diagnosis, implementation, and local validation for the
plain-GEX STARsolo scaling regression reported against STAR Suite v1.3.0b.
The fix is implemented on `fix/solo-fastx-scaling-20260809`, based on
`60eab145837903b2979daaeed15ac142cc814025`.

The regression is not an inherent limitation of legacy STARsolo. Two STAR
Suite integrations were active in ordinary Fastx runs even though the legacy
path did not need them:

1. The generic Fastx record adapter parsed and reconstructed each record while
   STAR's global input mutex was held. That serialized the reader and starved
   mapping threads.
2. `ParametersSolo` always constructed `CbCorrector`. For a 16-base 10x 3M
   whitelist this precomputes up to 48 one-mismatch probes per barcode, about
   177 million probes, even when inline correction and inline hashing are off.
   Plain STARsolo uses `matchCBtoWL()` instead, so the structure contributed no
   matrix result but cost roughly 8.1 GiB and 47 seconds on this host.

## Implementation

### Direct Fastx reader

`Parameters::readFilesInit()` still creates and validates an
`InputSourcePlan`, including manifest normalization and Y/noY format gates.
Ordinary Fastx data then use STAR's direct chunk reader. The generic
record-at-a-time adapter is no longer installed in the mutex-protected
`readLoad()` path.

Manifest input is routed through the existing `cat` FIFO path so the normalized
manifest filenames, rather than placeholder names, reach the direct reader.
The runtime log records:

```text
Fastx input path: direct STAR chunk reader
```

The standalone Fastx contract module remains available for validation and for
a future batched integration that does parsing outside the input mutex.

### Conditional inline barcode corrector

`CbCorrector` is now constructed only when `soloInlineCBCorrection` or
`soloInlineHashMode` is enabled. Standard STARsolo retains the existing H0
whitelist hash and `matchCBtoWL()` behavior. The runtime log records the
decision:

```text
CbCorrector not required: inlineCBCorrection=0, inlineHashMode=0
```

Flex and explicit inline modes retain the corrector and their existing
behavior.

## Benchmark method

The local host is a 24-core/32-thread Intel Core i9-13900KF with 125 GiB RAM.
The scaling comparison uses 8 and 32 threads as the available analogue of the
issue's 16-to-64-thread comparison. All runs were serialized and used fresh
output and temporary directories.

- Input: 10,000,000 paired GEX reads, one million from each of ten lanes of
  UCSF `iPSC2_1/GEX`.
- Input mode: gzip FASTQ with `--readFilesCommand zcat`.
- Index: local GRCh38 Cell Ranger-compatible index at
  `/storage/autoindex_110_44/bulk_index`.
- Whitelist: 10x `3M-february-2018_TRU.txt`, 3,686,400 barcodes.
- Counting: `GeneFull`, `1MM_CR`, `MultiGeneUMI_CR`, `Rescue`, and
  `EmptyDrops_CR`.
- Parity controls: `--clip3pPolyG no`; strict upstream comparisons also use
  `--soloEmptyDropsLegacy yes`.
- Wrapper: `tests/run_ucsf_solo_gex_100k_benchmark.sh` with
  `UCSF_GEX_TOTAL_READS=10000000`.
- Artifact root: `/storage/solo_fastx_scaling_fix_20260809/`.

Comparators:

- Affected local STAR Suite v1.3.1 binary, SHA-256
  `0de93f83f3ff49b5da85c1da4de3c0824f6cef6cfeac95982861c48b33ded6b5`.
  This contains the same record-adapter and unconditional-corrector behavior
  implicated in the v1.3.0b issue.
- Official upstream STAR 2.7.11b static release binary, SHA-256
  `36e94b899a56b0ea5de5d65e722f55bc552b6713bd2d1e83a5c2abbd06c9881a`.

Times derived from STAR timestamps have one-second resolution. `Wall` and RSS
come from `/usr/bin/time -v`. `Map+Solo` is the interval that STAR uses for the
`Log.final.out` mapping-speed calculation.

## Results

| Binary and mode | Threads | Startup (s) | Mapping (s) | Solo (s) | Map+Solo (s) | Log speed (M/h) | Wall (s) | Peak RSS (GB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Affected Suite v1.3.1, default | 8 | 56 | 24 | 24.52 | 49 | 734.69 | 136.07 | 39.28 |
| Affected Suite v1.3.1, default | 32 | 57 | 20 | 24.62 | 46 | 782.61 | 133.56 | 42.52 |
| Fixed Suite, default CR9 filter | 8 | 10 | 22 | 24.80 | 48 | 750.00 | 57.79 | 30.74 |
| Fixed Suite, default CR9 filter | 32 | 9 | 12 | 24.68 | 38 | 947.37 | 47.33 | 33.97 |
| Fixed Suite, upstream-compatible filter | 8 | 10 | 22 | 7.99 | 31 | 1161.29 | 41.09 | 30.74 |
| Fixed Suite, upstream-compatible filter | 32 | 9 | 12 | 7.99 | 21 | 1714.29 | 30.73 | 33.97 |
| Official upstream STAR 2.7.11b | 8 | 8 | 21 | about 4 | 25 | 1440.00 | 33.17 | 30.42 |
| Official upstream STAR 2.7.11b | 32 | 9 | 11 | about 4 | 15 | 2400.00 | 23.97 | 33.60 |

Key observations:

- Affected mapping scaled only 1.20x from 8 to 32 threads (24 to 20 seconds).
  Fixed mapping scaled 1.83x (22 to 12 seconds), close to upstream's 1.91x
  (21 to 11 seconds).
- The unconditional corrector removal reduced peak RSS by 8.54 to 8.55
  decimal GB and eliminated 46 to 48 seconds of startup. Fixed peak RSS is
  within 0.32 to 0.37 GB of upstream.
- The fixed default `Log.final.out` rate scales only 1.26x because STAR labels
  the whole mapping-through-Solo interval as mapping time. STAR Suite's default
  CR9/libscrna EmptyDrops stage is deterministic but serial and takes about 17
  seconds here. It is not evidence that the mapping threads remain flat.
- With `--soloEmptyDropsLegacy yes`, which selects the same cell-filter family
  as upstream, the fixed end-to-end rate scales 1.48x. The residual local gap
  is mostly Suite bookkeeping and the additional Rescue matrix output; true
  mapping differs from upstream by only one second at either thread count.

## Output parity

The following comparisons passed:

- Affected Suite versus fixed Suite default: raw and filtered `barcodes.tsv`,
  `features.tsv`, and `matrix.mtx` are byte-identical.
- Fixed Suite default at 8 versus 32 threads: the same six MEX components are
  byte-identical.
- Fixed Suite with `--soloEmptyDropsLegacy yes` versus official upstream at
  each thread count: the same six raw and filtered MEX components are
  byte-identical.

STAR Suite additionally writes `raw/UniqueAndMult-Rescue.mtx`; that extension
has no upstream counterpart and was excluded from the byte comparison.

## Why the OCM run did not catch it

The archived OCM FASTQ run did exercise the affected plain Fastx path. It did
not demonstrate that the path scaled:

- it used one thread point, 16 threads;
- it used only 100,000 reads;
- it compared FASTQ with CBQ for output parity, not STAR Suite with upstream
  STAR or 16 threads with a higher thread count;
- its alignment phase was only about two seconds, while initialization,
  Velocyto, and OCM materialization dominated the run; and
- both comparison arms included shared serialized work.

It was therefore a correctness test with too little parallel work to expose a
throughput regression. Future OCM acceptance should retain the parity run but
add at least two thread counts and report mapping-phase time separately.

## Validation gates

After rebuilding the changed source, the following passed:

- `tests/run_fastx_contract_star_smoke.sh`: plain FASTQ, internal gzip, explicit
  `zcat`, comma-separated lanes, manifest input, and Y/noY format rejection.
- `tests/run_fastx_input_harness_smoke.sh`: standalone normalized input
  contract.
- `tests/run_solo_smoke.sh`: genome generation with GTF, two-thread STARsolo,
  sorted BAM, matrix output, and plain-Solo corrector gate.
- An explicit `--soloInlineCBCorrection yes` gate initialized the corrector on
  the same tiny fixture and retained raw-matrix parity with the plain run.
- `tests/run_ychrom_fastq_se_test.sh`: single-end Y/noY FASTQ emission.
- 10-million-read 8/32-thread determinism and upstream MEX parity comparisons
  described above.

## User-facing interpretation

No dynamic-thread flag is required to restore scaling. In particular,
`--dynamicThreadInterface` and `--dynamicThreadConstMapPermits` do not repair
the affected reader; the source fix does.

For strict drop-in comparison with upstream STAR 2.7.11b, also disable STAR
Suite enhancements that intentionally change upstream behavior:

```bash
--clip3pPolyG no --soloEmptyDropsLegacy yes
```

For ordinary STAR Suite scRNA-seq, keep the default CR9/libscrna cell caller.
Its fixed-cost post-mapping time should be reported separately from true
mapping throughput.

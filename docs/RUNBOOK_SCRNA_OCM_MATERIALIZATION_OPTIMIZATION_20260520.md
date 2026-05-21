# OCM GeneFull and Velocyto Materialization Optimization Runbook

Date: 2026-05-20

Status: implemented in the native STAR OCM materializer for production
GeneFull/Velocyto output, with `scripts/ocm_composite_adapter.py` retained for
Gene comparator materialization and fallback/debug runs.

## Goal

Make JAX scRNAseq02 OCM production materialization fast enough to scale without
changing the validated OCM method:

```text
effective_cb = CB16 + OCM_TAG8
UMI          = original UMI12
```

The optimized path should produce the production-shaped `GeneFull` and
Velocyto outputs together, using one per-sample cell-calling result and one set
of OCM column maps. `Gene` materialization should remain available for Cell
Ranger comparator runs, but it should not be part of the normal production path.

## Current Bottleneck

The historical smoke harness used `scripts/ocm_composite_adapter.py
materialize` twice:

1. `--feature Gene` for Cell Ranger parity.
2. `--feature GeneFull --include-velocyto` for production-shape inspection.

This is useful for validation, but inefficient for production. On the 50M OCM
smoke:

```text
STAR core run                         9m10s
old Python STAR BAM post-split        3m24s
Gene materialization                  1m19s
GeneFull + Velocyto materialization   4m09s
```

The BAM post-split is already replaced by native `--ocmMultiBamSplit yes`.
The optimized adapter now uses streaming MatrixMarket routing, parallel
per-sample EmptyDrops, and parallel gzip finalization. The remaining
materialization cost before this optimization was mostly repeated
MatrixMarket I/O:

- split `GeneFull/raw/matrix.mtx` by OCM;
- run SimpleED separately for each sample;
- scan/write each filtered GeneFull MEX;
- split Velocyto raw layers (`spliced`, `unspliced`, `ambiguous`, and sometimes
  `matrix`);
- scan/write filtered Velocyto layers again for each sample;
- gzip all MatrixMarket outputs.

Several steps scan the same matrix twice: once to count output `nnz` for the
MatrixMarket header, then again to write records. Gene and GeneFull runs also
repeat the per-sample EmptyDrops work in comparator mode.

## Production Contract

For production OCM scRNA-seq, the default should be:

```text
--soloFeatures GeneFull Velocyto
--soloCellFilter None
--ocmMultiBarcodeMode flex
--ocmMultiBamSplit yes
```

Then the native STAR OCM materializer should:

- split `Solo.out/GeneFull/raw` by OCM sample;
- call cells per OCM sample using the agreed SimpleED/EmptyDrops policy;
- write per-sample raw and filtered `GeneFull` MEX;
- write per-sample raw and filtered Velocyto MEX using the same filtered
  barcodes;
- write CR-compatible `outs/per_sample_outs/...` and
  `samples/<sample>/run/outs/...` mirrors;
- avoid producing `Gene` unless a comparator run explicitly asks for it.

Fallback helper command:

```bash
scripts/ocm_composite_adapter.py materialize-production \
  --repo-root <repo> \
  --star-run-dir <run> \
  --config <config.csv> \
  --out-dir <out> \
  --threads 4
```

The comparator contract remains separate:

```text
--materialize-gene-for-cr-parity
```

Comparator output should never be required for downstream production handoff.

## Optimization Plan

### Phase 0: Measurement Harness

Add timing and size logging around each materializer substep:

```text
load/write barcode maps
split GeneFull raw
per-sample cell calling
write GeneFull filtered
split Velocyto raw
write Velocyto filtered
link/copy downstream mirrors
```

Record:

- elapsed seconds per substep;
- input matrix dimensions and `nnz`;
- output `nnz` by sample and layer;
- peak RSS for the materializer process;
- output directory size.

Use the existing 100K and 50M downsamples as gates.

### Phase 1: Collapse Duplicate Work in Python

Implemented in the helper:

1. Add one production subcommand:

```bash
scripts/ocm_composite_adapter.py materialize-production \
  --star-run-dir <run> \
  --config <config.csv> \
  --out-dir <out> \
  --include-velocyto
```

2. Build OCM barcode maps once from `GeneFull/raw/barcodes.tsv`.
3. Run per-sample EmptyDrops once from GeneFull raw MEX.
4. Reuse those filtered barcode lists for GeneFull filtered and Velocyto
   filtered outputs.
5. Do not materialize `Gene` unless `--materialize-gene-for-cr-parity` is set.
6. Use hardlinks for repeated metadata files and downstream mirrors when
   possible.

Acceptance for Phase 1:

- output cells and counts match the current two-stage helper exactly on 100K;
- 50M production materialization time drops by removing the `Gene` comparator
  stage and duplicate ED work;
- no `star_materialized/bam` post-split is created in native mode.

### Phase 2: Stream MatrixMarket Once Per Layer

Implemented in the helper for OCM GEX and Velocyto materialization.

MatrixMarket requires `nnz` in the header, which caused the helper to do a
count pass and a write pass. Replace that with a one-parse temp-body writer:

1. Open one temporary body file per output sample/layer.
2. Stream input entries once.
3. Route each entry to the proper sample body and count `nnz`.
4. Close bodies.
5. Write final `.mtx.gz` with the header followed by the body.
6. Remove temp bodies.

For Velocyto, write raw and filtered outputs in the same stream:

```text
for each entry(row, col, value):
  sample = ocm_sample_for_col[col]
  raw_col = raw_col_map[sample][col]
  write raw body
  if col is in filtered set for sample:
      filtered_col = filtered_col_map[sample][col]
      write filtered body
```

This avoids splitting raw Velocyto and then re-reading the per-sample raw
Velocyto directories just to make filtered outputs. If STAR writes Velocyto
`spliced`, `unspliced`, and `ambiguous` without a total `matrix.mtx`, the
helper writes a synthetic total `matrix.mtx.gz` in the same streaming pass by
emitting the layer entries into the total matrix body. Readers that coalesce
MatrixMarket duplicates recover the exact summed total.

Acceptance for Phase 2:

- each input matrix layer is parsed once;
- raw and filtered Velocyto outputs are emitted from the same pass;
- output MEX is byte-order stable enough for deterministic tests, or documented
  as count-equivalent if gzip timestamps/body ordering differ.

### Phase 3: Native STAR C++ Materializer

Implemented in the native materializer:

- `core/legacy/source/OcmMultiMaterialize.cpp`
- `core/legacy/source/VelocytoMexWriter.cpp`
- shared MEX streaming helpers if needed.

Internal API shape:

```text
streamOcmMexLayers(
  inputDir,
  matrixLayers,
  sampleDefinitions,
  rawColumnMaps,
  filteredColumnMaps,
  outputLayout
)
```

The native implementation:

- share one barcode axis and feature axis across GeneFull and Velocyto;
- share one per-sample filtered-cell list;
- support union sample definitions such as `OB1|OB2`;
- write both CR-compatible `outs/...` and downstream `samples/...` mirrors;
- keep the existing `--ocmMultiEnable` output layout behavior;
- preserve `--ocmMultiBamSplit yes` as the native BAM path.

It streams MatrixMarket routing and writes temp bodies for final gzip headers;
it does not add another large in-memory sparse matrix representation. If
`Solo.out/GeneFull/filtered` is absent because `--soloCellFilter None` was used,
the materializer runs per-sample EmptyDrops from the OCM-split raw MEX and then
streams filtered GeneFull and Velocyto from the same composite barcode axis.

### Phase 4: Production Wrapper Cleanup

Production wrappers now use the native boundary:

- do not call Python `split-bam`;
- do not materialize `Gene` unless a comparator flag is requested;
- do not run SimpleED twice;
- run one native GeneFull+Velocyto materialization step;
- transfer only the production outputs needed for downstream h5ad/h5mu creation.

The preferred downstream boundary remains after:

```text
GeneFull MEX
Velocyto MEX
native OCM BAMs
QC/cell-calling metadata
```

CellBender, h5ad/h5mu construction, and heavier downstream analysis can remain
on the remote GPU server after this boundary.

## Tests

### Unit Fixtures

Add or extend tiny fixtures to cover:

- four OCM tags with one sample per tag;
- a union sample (`OB1|OB2`);
- shuffled Velocyto barcode order;
- missing/empty Velocyto layers;
- path-safe sample IDs;
- filtered barcode reuse across GeneFull and Velocyto.

Expected checks:

- raw per-sample GeneFull sums equal routed pool counts;
- filtered GeneFull is an exact subset of raw GeneFull;
- Velocyto raw and filtered barcode axes match GeneFull axes;
- `spliced`, `unspliced`, and `ambiguous` columns remain order-independent;
- no sample receives records from another OCM tag.

### Smoke Tests

Run from repo root:

```bash
scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh \
  --read-pairs 100000 \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_materialize_opt_100k_$(date -u +%Y%m%dT%H%M%SZ) \
  --run-all
```

Then run a production-shaped smoke without CR:

```bash
scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh \
  --read-pairs 100000 \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_materialize_prod_100k_$(date -u +%Y%m%dT%H%M%SZ) \
  --prepare --run-star
```

The second run is the STAR/native-materializer performance signal because it
excludes Cell Ranger and the `Gene` comparator. Add `--materialize --compare`
only when a Gene-vs-CR parity table is needed.

When a previous Cell Ranger control exists for the same deterministic
downsample, reuse it instead of rerunning CR. `CR_REUSE_RUN_DIR=auto` is the
default and searches `/mnt/pikachu/JAX_scRNAseq02_processed`; pass an explicit
path when comparing to a known reference:

```bash
CR_REUSE_RUN_DIR=/path/to/cellranger/<id> \
OCM_MATERIALIZE_THREADS=4 \
scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh \
  --read-pairs 50000000 \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_materialize_opt_50m_$(date -u +%Y%m%dT%H%M%SZ) \
  --prepare --run-star --materialize --compare
```

### 50M Gate

After 100K passes, run:

```bash
scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh \
  --read-pairs 50000000 \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_materialize_opt_50m_$(date -u +%Y%m%dT%H%M%SZ) \
  --prepare --run-star
```

Compare against the previous 50M baseline:

```text
Gene materialization                  1m19s
GeneFull + Velocyto materialization   4m09s
```

Target:

- production materialization should be closer to one GeneFull+Velocyto pass
  than the current Gene + GeneFull duplicate path;
- no separate BAM post-split time;
- no repeated ED time;
- stable or lower peak RSS than the current helper.

For the parity table, run the same output root after STAR completes:

```bash
scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh \
  --read-pairs 50000000 \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_materialize_opt_50m_<stamp> \
  --materialize --compare
```

That step materializes STAR `Gene` only for the Cell Ranger comparator and
reuses existing Cell Ranger artifacts. It skips Python GeneFull/Velocyto when
native STAR OCM outputs exist.

Validated 2026-05-21 result:

```text
out_root=/mnt/pikachu/JAX_scRNAseq02_processed/ocm_native_cpp_50M_20260521T003926Z
STAR-only elapsed=13m44s
native_Ocm_materialization=2m12s
tagged_BAM_writeback=5m45s
peak_logged_VmRSS=44.7 GB during mapping; post-materialization VmRSS=11.6 GB
```

Note: this timing run predates the composite harness fixes that default
`STAR_YREMOVE=yes` and `STAR_BAM_CBUB_TAGS=no`. It remains the native
materialization timing benchmark, but production JAX OCM runs should include
the Y/noY sidecars and should not request barcode/UMI BAM tags unless they are
needed for a specific downstream consumer.

The comparator step reused the existing CR reference and did not rerun Cell
Ranger:

```text
CR_REUSE_RUN_DIR=/mnt/pikachu/JAX_scRNAseq02_processed/ocm_composite_smoke_50m_bam_20260520T164318Z/cellranger/25E32-L3_ocm_composite_50000000
```

50M Gene parity remained in the expected range:

| sample | STAR cells | CR cells | Jaccard | barcode UMI Pearson | feature UMI Pearson |
| --- | ---: | ---: | ---: | ---: | ---: |
| GCM1-Day-4 | 2751 | 2733 | 0.993457 | 0.999977 | 0.999564 |
| GRHL1-Day-4 | 2896 | 2900 | 0.998621 | 0.999958 | 0.999574 |
| OVOL1-Day-4 | 3219 | 3227 | 0.997521 | 0.999922 | 0.999614 |
| WT-PrS-20pct-Day-4 | 2214 | 2223 | 0.995951 | 0.999970 | 0.999598 |

## Acceptance Gates

Before using the optimized path for production:

- 100K optimized output matches current helper counts and cells.
- 50M optimized output has the same parity profile as the prior run where
  relevant.
- Native per-sample BAM outputs exist under `star_composite/outs/...`.
- GeneFull h5ad/h5mu construction sees the same `obs` and velocity-layer
  barcode axes.
- Logs explicitly state whether `Gene` comparator materialization was skipped
  or requested.
- The runbook and wrapper command line document the exact STAR binary hash and
  materializer implementation version used.

## Non-Goals

- Do not change OCM methodology.
- Do not use pooled CB16 EmptyDrops for production OCM.
- Do not optimize Cell Ranger control output.
- Do not make `Gene` parity materialization part of the production handoff.
- Do not add an in-memory full sparse matrix cache to avoid writing a streaming
  implementation.

# Handoff: Solo GeneFull + OCM Materialization Memory Optimizations

**Date**: 2026-05-19
**Status**: Implemented in working tree; smoke-tested on tiny OCM fixture, 100k JAX OCM low-memory smoke, 100k JAX OCM aggressive-spill smoke (8192 buckets, binary GeneFull spool), 100k JAX OCM profile, and 2M CountingSink replay. **OCM streaming now applies CR barcode normalization** (`buildCrBarcodeLayoutForColumns`, unsuffixed-input unit test). **Not validated at full production depth (1.3B reads).**
**Branch context**: JAX scRNAseq02 OCM oracle path (`soloInlineHashMode=no`, GeneFull + Velocyto)

## Executive summary

Two production blockers were observed. The first OOM occurred **during GeneFull
post-map counting**, immediately after `nReadsInput ≈ 1.317×10⁹`. After the
CountingSink fixes, GeneFull completed, but the next production attempt reached
**~117 GB RSS during Velocyto counting** on the default stream path and was
killed before OCM materialization.

This handoff covers:

1. **P0 — GeneFull `CountingSink` peak reduction** (production blocker)
2. **Instrumentation** — `STAR_SOLO_MEMORY_PROFILE` + harness scripts
3. **P1 — OCM streaming column-subset MEX** (post-Solo materialization spike)
4. **P1 — Velocyto low-memory range spill** (bounded per-CB UMI maps)

---

## Problem diagnosis

### Production failure point

- Run died in Solo post-map, right after logging full `nReadsInput`, before downstream stages.
- OCM materializer was never reached on that run.

### Root-cause structures (GeneFull, `soloInlineHashMode=no`)

| Structure | Role | Scale concern |
|-----------|------|----------------|
| `PackedReadInfo` | 8 B / input read | ~10.5 GB at 1.3B reads (necessary for BAM CB/UB + Velocyto stream) |
| `CountingSink::perWL` | Was `vector<vector<ReadInfoRecord>>` over **full whitelist** | Tens of GB at full record count |
| `CountingSink::readToCb` | `unordered_map` conflict guard per read index | Node overhead; plausibly **largest avoidable** term |
| `rGeneUMI` | Second materialization after `perWL` | Peak while **both** buffers lived |
| Velocyto `cuTrTypes` | Per-CB `unordered_map<UMI, vector<trTypeStruct>>` | Default stream path holds all CB maps at once; production saw ~117 GB RSS |

### Second spike (OCM, after Solo)

| Old behavior | Issue |
|--------------|--------|
| `PfMultiMerge::readMex()` ×2 | Full pool raw + filtered triplets in RAM |
| `subsetMexColumns()` per sample | `oldToNew` sized to pool columns; `reserve(full pool nnz)` |
| Velocyto load + subset × layers | Six full triplet vectors + `std::map` total matrix |
| Duplicate writes | Same subset held in memory for `per_sample_outs` and `samples/.../run/outs` |

---

## Implemented optimizations

### 1. GeneFull `CountingSink` (`flex/source/SoloReadInfoSink.{h,cpp}`)

**Three changes (all required for the pathological peak; gating `readToCb` alone is insufficient):**

1. **`readToCb` gated in production**
   - Active only when `STAR_DEBUG_COUNTING_SINK_READ_TO_CB=1`.
   - Default: no hash map on the hot path.

2. **Observed-CB compact buckets**
   - Replaced `perWL[cbWLsize]` (millions of empty `std::vector` slots) with
     `unordered_map<uint32_t, vector<CountingSinkRow>>` (12 B/row: `featureId`, `umi`, `readRef`).
   - Only whitelist indices that receive records allocate storage.

3. **Per-CB bucket clear during `rGeneUMI` fill**
   - After allocating full `rGeneUMI`, fill **per-CB** and **clear each bucket** in the same loop.
   - **Caveat:** `finalize()` still allocates the full `rGeneUMI` while all compact buckets are resident; buckets shrink only during the fill loop. This removes `readToCb` and whitelist `perWL` overhead (likely enough for the observed OOM) but is **not** a strict single-buffer peak. A two-pass fill (count then allocate) would be P2.

**Related:** `friend class CountingSink` on `SoloFeature` for profile logging.

### 2. Memory profiling (`core/legacy/source/SoloMemoryProfile.{h,cpp}`)

- Enable: `export STAR_SOLO_MEMORY_PROFILE=1`
- Logs: `Solo memory: <label> | <counters> | VmRSS/VmHWM/VmPeak` to `logMain` (`Log.out` / `logs/star.log` depending on run layout).

**GeneFull / Velocyto checkpoints** — see `docs/SOLO_MEMORY_PROFILE_HARNESS.md`.

**OCM checkpoints:**

- `ocm_materialize_begin`
- `ocm_materialize_axes_loaded`
- `ocm_materialize_multi_raw_copied`
- `ocm_materialize_velocyto_axes_loaded`
- `ocm_materialize_sample_begin:<sample_id>` / `ocm_materialize_sample_done:<sample_id>`
- `ocm_materialize_all_samples_done`

### 3. CountingSink replay harness (no mapping)

- **`--runMode countingSinkStress`** (`core/legacy/source/CountingSinkStress.{h,cpp}`, `STAR.cpp`)
- **Script:** `tests/run_counting_sink_stress.sh`

Isolates post-map counting from genome index / alignment RAM.

| Env var | Default | Meaning |
|---------|---------|---------|
| `STAR_COUNTING_SINK_STRESS_NRECORDS` | 50M (or script arg) | Synthetic counted records |
| `STAR_COUNTING_SINK_STRESS_ACTIVE_CBS` | 20000 | Distinct CB indices |
| `STAR_COUNTING_SINK_STRESS_COLLAPSE` | off | Also run `collapseUMIall` |
| `STAR_DEBUG_COUNTING_SINK_READ_TO_CB` | off | Re-enable conflict map |

Requires `--genomeDir` + `--soloCBwhitelist` (transcriptome init only).

### 4. OCM streaming materialization

**GEX (`core/legacy/source/PfMultiMerge.{h,cpp}`):**

- `readMexAxes()` — features + barcodes only
- `buildColumnRemap()` — `oldToNew` only for selected columns (4 × pool_cols, not pool nnz)
- `streamMatrixColumnSubset()` — two-pass stream of `matrix.mtx` → `matrix.mtx.gz`
- `buildCrBarcodeLayoutForColumns()` — GEM `-1` suffix, lexicographic sort, NXT/TRU translation (same rules as `writeCombinedMex`)
- `writeColumnSubsetMexGz()` — per-sample streamed MEX with CR-compat barcodes
- `writeStreamedPoolMexGzCrCompat()` — pool `outs/multi/.../raw` with the full raw barcode axis + CR formatting (replaces raw `copyMexGzDir` on pool input)
- `copyMexGzDir()` — downstream mirrors only (copies already CR-formatted per-sample output); uses `createDirectory()`

**Velocyto (`core/legacy/source/VelocytoMexWriter.{h,cpp}`):**

- `readVelocytoAxes()` / `loadSoloVelocytoAxes()` — axes only
- `streamVelocytoColumnSubsetToDir()` — stream spliced/unspliced/ambiguous; build `matrix.mtx.gz` from a **sample-sized** merge map (three layer passes over disk, not six full in-memory layers)
- Filtered-without-layers: join GeneFull filtered barcodes to **raw** Velocyto matrix columns (`outputBarcodesOverride`)
- OCM runs skip pooled run-level Velocyto `outs/` materialization; per-sample
  Velocyto is streamed from `Solo.out/Velocyto`.

**Velocyto counting (`core/legacy/source/SoloFeature_countVelocyto*.cpp`):**

- `STAR_VELOCYTO_LOW_MEM=1` selects the range-spill path directly.
- Records are spilled to contiguous CB-range buckets; each bucket is loaded,
  sorted, merged, finalized into `countCellGeneUMI`, and released before the
  next bucket.
- Default low-memory bucket count is 4096 unless
  `STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS` is set; the JAX harness uses
  8192 and `STAR_VELOCYTO_UMI_RESERVE_CAP=32`.
- The JAX harness also enables `STAR_SOLO_BINARY_SPOOL=1` so the GeneFull replay
  spool is written as compact on-disk binary records instead of legacy text, and
  sets `MALLOC_ARENA_MAX=2` / `MALLOC_TRIM_THRESHOLD_=131072` to reduce allocator
  arena retention between phases.

**OCM loop (`core/legacy/source/OcmMultiMaterialize.cpp`):**

- Load axes once; tag/classify barcodes in memory
- Multi raw: `writeStreamedPoolMexGzCrCompat` (not a blind pool copy)
- Per sample: `buildCrBarcodeLayoutForColumns` + `writeColumnSubsetMexGz` → `copyMexGzDir` to downstream mirrors
- Velocyto: stream with GEX-derived `CrBarcodeLayout` (sorted/suffixed barcodes, column remap aligned with GEX)
- Velocyto: stream once per sample to `samples/<id>/run/outs/`

---

## Testing performed

### Build

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
make -C core/legacy/source -j8 ocm-multi-unit-tests   # optional
```

### Unit / fixture tests (pass)

| Test | Command | Result |
|------|---------|--------|
| OCM unit | `core/legacy/source/ocm_multi_unit_tests all` | PASS (includes `cr_barcode` unsuffixed-input suffix test) |
| OCM tiny materializer | `tests/test_ocm_mex_materializer_tiny.sh` | PASS (structure + MEX values) |
| OCM Velocyto shuffle | `tests/test_ocm_velocyto_barcode_shuffle.sh` | PASS |
| OCM streaming memory smoke | `tests/run_ocm_materializer_memory_smoke.sh` | PASS |
| JAX OCM 100k aggressive spill | `scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh --downsample-read-pairs 100000 --run-star --validate` | PASS (`spill_buckets=8192`, `STAR_SOLO_BINARY_SPOOL=1`) |

### Solo smoke with profiling

```bash
export STAR_SOLO_MEMORY_PROFILE=1
tests/run_solo_smoke.sh
```

Checkpoints appear in `tests/solo_smoke/output/Log.out` (Gene path; tiny read count).

### JAX OCM 100k profile run

```bash
STAR_SOLO_MEMORY_PROFILE=1 \
  tests/run_solo_memory_profile_harness.sh \
  --jax-ocm-smoke --downsample-read-pairs 100000 --run-star \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_memprof_100k_<stamp>
```

**Parse checkpoints from** `samples/25E32-L3/run/Log.out` (not `logs/star.log` alone — tee may omit `logMain` lines).

**100k findings (VmRSS deltas, illustrative):**

| Phase | Δ VmRSS (approx.) | Notes |
|-------|-------------------|--------|
| `sumThreads_done:GeneFull` | baseline ~10.5 GB | Genome + spool; dominates at this scale |
| `CountingSink_loader` | **+133 MB** | Largest **counting-phase** step |
| `readToCb` | 85,793 entries (= buffered records) with new binary | Would be **0** in production |
| Velocyto `cuTrTypes` | ~0 at 100k | Not yet limiting at this scale |
| OCM | N/A on this run | Materializer runs after Solo completes |

**Caution:** 100k is **~13,000×** smaller than production by read count; use for **ordering** of phases, not linear GB extrapolation.

### CountingSink replay (2M synthetic records)

```bash
STAR_SOLO_MEMORY_PROFILE=1 tests/run_counting_sink_stress.sh 2000000
```

- Counting-structure RSS growth from loader through finalize: **~80 MB** (excluding ~9 GB transcriptome/whitelist baseline in stress mode).
- `counting_sink_readToCb=0` in profile counters.

### OCM tiny fixture memory profile

`tests/run_ocm_materializer_memory_smoke.sh`:

- **Peak VmRSS ~8.6 MB** for 5 samples (4-barcode toy MEX).
- Per-sample steps: **+0–68 kB** after first sample; **flat** across remaining samples.
- Confirms streaming path does not accumulate pool triplets in RAM.

---

## Production RAM scaling (deduced)

### After CountingSink fixes (if Solo completes)

| Component | Old order of magnitude | New expectation |
|-----------|------------------------|-----------------|
| `readToCb` | Tens of GB possible | **0** (production) |
| `perWL` full whitelist + `ReadInfoRecord` | Tens of GB | **O(records × 12 B)** in buckets + streaming clear |
| `rGeneUMI` + buckets | 2× record bytes briefly | Buckets cleared during fill; brief overlap of full `rGeneUMI` + all buckets at loop start |
| `PackedReadInfo` | ~10.5 GB @ 1.3B reads | Unchanged (still required) |

**Gate before full production retry:**

```bash
STAR_SOLO_MEMORY_PROFILE=1 tests/run_counting_sink_stress.sh 50000000
# optional: STAR_COUNTING_SINK_STRESS_COLLAPSE=1
```

Require: roughly **linear** VmRSS vs `N`, flat across replay scales, `readToCb=0`.

### After OCM streaming (post-Solo)

| Component | Old | New |
|-----------|-----|-----|
| Pool GEX triplets | `2 × nnz × 12 B` in RAM | **Streamed**; axes only in RAM |
| Per-sample subset | `reserve(pool nnz)` + duplicate `MexData` | Stream + **gz copy** |
| Velocyto | 6× full layers + map on large subset | Stream layers; **map ≈ sample nnz** |
| Velocyto counting maps | Full `cuTrTypes` for all CBs | Range-spill; **one CB bucket** plus final count matrix |

Expected OCM footprint: **O(pool barcodes + genes)** fixed plus **O(sample nnz)** per sample sequential — typically **low GB**, not pool-scale tens of GB.

**Remaining Velocyto risk:** `matrix.mtx.gz` total still uses an in-memory `map` over **sample** (cell, gene) pairs after three streamed layers. Huge per-sample nnz could still be costly; it is no longer pool × layers × duplicate writes.

---

## Operational commands

### Profile full JAX OCM path

```bash
export STAR_SOLO_MEMORY_PROFILE=1
tests/run_solo_memory_profile_harness.sh \
  --jax-ocm-smoke --downsample-read-pairs 2000000 --run-star
```

Parse:

```bash
tests/run_solo_memory_profile_harness.sh --parse-log \
  <out_root>/samples/25E32-L3/run/Log.out
```

### Profile OCM only (after Solo outputs exist)

```bash
export STAR_SOLO_MEMORY_PROFILE=1
export OCM_TEST_FIXTURE_ROOT=...   # or use production run dir via unit test env
tests/run_ocm_materializer_memory_smoke.sh
```

### Debug CountingSink conflict map

```bash
export STAR_DEBUG_COUNTING_SINK_READ_TO_CB=1
```

---

## Files touched (main)

| Area | Files |
|------|--------|
| CountingSink | `flex/source/SoloReadInfoSink.{h,cpp}`, `core/legacy/source/SoloFeature.h` |
| Profiling | `core/legacy/source/SoloMemoryProfile.{h,cpp}`, hooks in `SoloFeature_countCBgeneUMI.cpp`, `SoloFeature_processRecords.cpp`, `SoloFeature_countVelocytoBridge.cpp`, `flex/source/SoloReadInfoSink.cpp` |
| Stress mode | `core/legacy/source/CountingSinkStress.{h,cpp}`, `STAR.cpp`, `parametersDefault` |
| OCM streaming | `core/legacy/source/OcmMultiMaterialize.cpp`, `PfMultiMerge.{h,cpp}`, `VelocytoMexWriter.{h,cpp}` |
| Tests / docs | `tests/run_solo_memory_profile_harness.sh`, `tests/run_counting_sink_stress.sh`, `tests/run_ocm_materializer_memory_smoke.sh`, `docs/SOLO_MEMORY_PROFILE_HARNESS.md`, `tests/ARTIFACTS.md` |

---

## Recommended next steps

1. **Monitor current first-sample production attempt**:
   `/mnt/pikachu/JAX_scRNAseq02_processed/ocm_prod_25E32-L3_aggressive_lowmem_20260519T201025Z`
   (profiling disabled; `STAR_VELOCYTO_LOW_MEM=1`,
   `STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS=8192`,
   `STAR_VELOCYTO_UMI_RESERVE_CAP=32`, `STAR_SOLO_BINARY_SPOOL=1`).
2. **If production fails**, use the CountingSink replay ladder: 50M → 200M records; confirm linear RSS and `readToCb=0`.
3. If Solo completes, confirm OCM materialization stays **flat per sample** on the multi-million-cell pool run.
4. **Optional P2** (see `docs/RUNBOOK_STARSOLO_SOLO_PHASE_OPTIMIZATION_20260319.md`): single-pass `rGeneUMI` fill, collapse in-place sort, and further Velocyto fusion with the GeneFull bridge if full-depth profiling still warrants it.
5. **OCM at scale:** consider single-pass Velocyto total matrix (avoid second read of layers for merge map) if sample nnz is huge.

---

## Related documentation

- `docs/SOLO_MEMORY_PROFILE_HARNESS.md` — checkpoint map and harness usage
- `docs/RUNBOOK_STARSOLO_SOLO_PHASE_OPTIMIZATION_20260319.md` — P0/P1/P2 Solo phase plan
- `docs/RUNBOOK_SCRNA_OCM_MULTI_MEX_MATERIALIZER_IMPLEMENTATION_20260519.md` — OCM functional spec (if present in tree)
- `tests/ARTIFACTS.md` — `ocm_memprof_*`, `ocm_oracle_smoke_*` output locations

---

## Artifact locations

- `/mnt/pikachu/JAX_scRNAseq02_processed/ocm_memprof_100k_*` — 100k JAX profile runs
- `/mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_smoke_*` — JAX OCM smoke (may predate profiling binary)
- `/tmp/counting_sink_stress_*` — CountingSink replay outputs
- `/tmp/ocm_materializer_mem.log` — tiny OCM streaming profile log from `run_ocm_materializer_memory_smoke.sh`

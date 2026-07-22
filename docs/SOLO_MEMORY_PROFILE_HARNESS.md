# Solo post-map memory profiling harness

**Date**: 2026-05-19
**Scope**: GeneFull / Velocyto counting OOM on `soloInlineHashMode=no` (CountingSink + PackedReadInfo + Velocyto `cuTrTypes`)

## Enable

```bash
export STAR_SOLO_MEMORY_PROFILE=1
# optional: STAR_SOLO_PHASE_DEBUG=1 for extra phase lines
```

Rebuild STAR after pulling profiling hooks:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

## Run

**JAX OCM smoke (matches production OCM path):**

```bash
tests/run_solo_memory_profile_harness.sh \
  --jax-ocm-smoke --downsample-read-pairs 2000000 --run-star
```

**Bisect scale:**

```bash
tests/run_solo_memory_profile_harness.sh \
  --jax-ocm-smoke --ladder 100000,500000,2000000 --run-star
```

**UCSF 100k GeneFull+Velocyto (same CountingSink path, explicit `soloInlineHashMode no`):**

```bash
tests/run_solo_memory_profile_harness.sh \
  --ucsf-velocyto-100k --threads 8 --run-star \
  --out-prefix /mnt/pikachu/JAX_scRNAseq02_processed/ucsf_velocyto_memprof_100k_TEST
```

**Parse an existing log:**

```bash
tests/run_solo_memory_profile_harness.sh --parse-log /path/to/logs/star.log
```

## Checkpoint map (read-scale structures)

| Label | Structure / phase |
|-------|-------------------|
| `sumThreads_done:GeneFull` | Per-thread `readFeat` spool merged |
| `packedReadInfo_allocated` | `PackedReadInfo.data` (~8 B/read) |
| `CountingSink_loader_begin/done` | Loader pass over spool |
| `CountingSink::finalize_before_rGeneUMI` | `perWL` buffers + `readToCb` (see counters) |
| `CountingSink::finalize_after_fill` | `rGeneUMI` / `rCBp` materialized |
| `collapseUMIall_begin/done` | Per-CB sort/copy in collapse |
| `count_phase_done:GeneFull` | End of GeneFull counting |
| `countVelocytoStreamThreads_maps_allocated` | `cuTrTypes` vector of per-CB `unordered_map` |
| `countVelocytoStreamThreads_input_done` | After stream ingest (`velocyto_umi_buckets`, `velocyto_tr_entries`) |
| `countVelocyto_finalize_done` | Velocyto collapse / matrix build |

## Interpreting results

1. Use the harness `--parse-log` output: largest **single-step VmRSS jump** pinpoints the phase.
2. Compare **counter lines** at `finalize_before_rGeneUMI` vs `after_fill`:
   - `counting_sink_buffered_records` × ~32 B ≈ `perWL` cost
   - `rGeneUMI_bytes_est` ≈ stride × records × 4 (often similar magnitude, so **double materialization** shows as two jumps)
3. Velocyto: if RSS jumps at `maps_allocated` or `input_done`, `cuTrTypes` dominates; compare `velocyto_umi_buckets` to GeneFull record counts.
4. OCM materializer runs **after** Solo; do not attribute pool OOM to `OcmMultiMaterialize` until Solo checkpoints show headroom.

## CountingSink replay (no mapping)

Synthetic post-map stress without FASTQ or genome mapping:

```bash
STAR_SOLO_MEMORY_PROFILE=1 tests/run_counting_sink_stress.sh 50000000
```

Environment:

| Variable | Default | Meaning |
|----------|---------|---------|
| `STAR_COUNTING_SINK_STRESS_NRECORDS` | arg or 50M | Synthetic counted records |
| `STAR_COUNTING_SINK_STRESS_ACTIVE_CBS` | 20000 | Distinct whitelist CB indices used |
| `STAR_COUNTING_SINK_STRESS_COLLAPSE` | off | Run `collapseUMIall` after finalize |
| `STAR_DEBUG_COUNTING_SINK_READ_TO_CB` | off | Enable read-index conflict map |

Production path uses compact observed-CB buckets and omits `readToCb` unless debug is set.

**Gate before full production:** extrapolate peak VmRSS from 50M–200M replay records; linear growth with `N` and no 3.7M-slot vector overhead is required.

## Next code fixes (when profiling confirms)

See `docs/RUNBOOK_STARSOLO_SOLO_PHASE_OPTIMIZATION_20260319.md` P2:

- Single-pass fill into `rGeneUMI` (drop `perWL` `ReadInfoRecord` buffer)
- Gate `readToCb` behind debug env
- Collapse in-place sort where safe

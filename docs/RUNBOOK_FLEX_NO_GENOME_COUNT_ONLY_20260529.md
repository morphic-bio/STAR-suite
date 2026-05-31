# Runbook: Flex No-Genome Count-Only Mode

Date drafted: 2026-05-29

## Goal

Avoid loading the STAR genome for pure STAR-Flex hash/count-only runs. The
current `--flexNoAlign 1 --outSAMtype None` path skips alignment inside the
Flex pipeline, but STAR still loads the genome and constructs `ReadAlignChunk`
objects before the Flex pipeline guard runs. For production Flex count-only
work, that costs startup wall time and peak RAM without producing useful state.

Target mode:

```text
--flex yes
--flexPipeline yes
--flexPipelineNTriage 0
--flexPipelineNSolo 0
--flexNoAlign 1
--outSAMtype None
```

This should support the existing FASTQ and CBQ fused Flex input surfaces.

## Current Behavior

Key landmarks:

- `core/legacy/source/STAR.cpp`
  - `Genome genomeMain(P, P.pGe); genomeMain.genomeLoad();` runs before mapping
    setup.
  - `samHeaders(P, *genomeMain.genomeOut.g, transcriptomeMain);` assumes a
    loaded genome.
  - `ReadAlignChunk(P, genomeMain, ...)` objects are constructed before
    `mapThreadsSpawn()`.
- `core/legacy/source/mapThreadsSpawn.cpp`
  - `flexPipelineActivationGuard(P)` runs only inside `mapThreadsSpawn()`.
  - `mapThreadsSpawnFlexPipeline(P, RAchunk)` still expects `ReadAlignChunk**`
    because it uses `RAchunk[i]->RA` for role-switch alignment and uses
    `RAchunk[0]->RA->soloRead->readFeat[...]` as the final merge target.
- `core/legacy/source/FlexPipeline.cpp`
  - `flexLaneReaderFullThread()` already skips `RA->oneReadFromPacket()` when
    `P.pSolo.flexNoAlign != 0`.
  - In fully fused no-align mode, the actual count work uses
    `SoloReadFeature`, `SoloReadBarcode`, sample detection, and the Flex hash
    cache. It does not need a loaded genome for alignment.

## Non-Goals

- Do not change the normal STAR FASTQ/CBQ mapper path.
- Do not support no-genome mode when alignment, BAM/SAM output, chimeric output,
  splice junction output, transcriptome quantification, SLAM, Chromap, or
  other genome-backed outputs are requested.
- Do not remove `--genomeDir` support globally. Existing commands can keep
  passing it; the count-only bypass should simply not load the genome.
- Do not widen this to split Flex pipeline modes yet. Start with the fully
  fused no-align mode because it is already the production benchmark surface.

## Strict Activation Predicate

Add an early predicate before `Genome genomeMain(...)` in `STAR.cpp`. It should
be stricter than the existing Flex pipeline guard:

```text
P.runMode == "alignReads"
P.pSolo.flexMode
P.pSolo.hashScreenEnabled
P.pSolo.flexPipelineStr != "no"
P.pSolo.flexNoAlign != 0
P.pSolo.flexPipelineNTriage == 0
P.pSolo.flexPipelineNSolo == 0
P.outSAMtype[0] == "None"
P.readFilesN > 0
P.runThreadN >= 1
prebuilt Flex hash cache path is provided and readable
no genome-backed outputs are enabled
```

Recommended disqualifiers:

- `P.outSAMtype[0] != "None"`
- chimeric output enabled
- `P.outSJ.yes` or non-default SJ output
- transcriptome quantification enabled
- `P.pGe.transform.outYes`
- two-pass or splice-junction insertion modes
- SLAM or other quantification modes that require genome/transcriptome state
- BAM/Y/split outputs that require SAM/BAM headers or read alignment state
- missing or unreadable Flex hash cache

If the probe list is normally auto-resolved from `genomeDir`, the bypass may
still read the small probe-list file directly, but it must not call
`genomeLoad()`. Prefer requiring an explicit `--soloProbeList` for the first
implementation if that keeps the guard simple.

The no-genome mode must require prebuilt hashes. It should not build a Flex hash
cache, infer one from genome state, or silently fall back to a slower
genome-backed setup after advertising that no-genome mode is active. If the hash
cache is missing, fail the predicate before genome load and continue through the
current normal STAR path only when the requested command is otherwise valid.

## Implementation Plan

1. Factor the Flex pipeline guard.

   Move the shared checks out of the private static function in
   `mapThreadsSpawn.cpp` into a small helper that can be called from `STAR.cpp`
   before genome loading. Keep log messages clear:

   ```text
   Flex count-only no-genome: active
   ```

   and for rejected attempts:

   ```text
   Flex count-only no-genome: not active (<reason>)
   ```

2. Add a no-genome entrypoint.

   Add a new function such as:

   ```cpp
   void runFlexNoGenomeCountOnly(Parameters &P);
   ```

   It should live close to `mapThreadsSpawnFlexPipeline()` or in a dedicated
   Flex count-only source file. The function should not require `Genome`,
   `Transcriptome`, `ReadAlign`, or `ReadAlignChunk`.

3. Reuse the fully fused Flex pipeline readers.

   `FlexPipelineState`, `flexLaneReaderFullThread()`, FASTQ lane reading, and
   CBQ lane reading are already the correct ingestion surface. For no-genome
   mode, pass `RA=nullptr` in `FlexLaneReaderArgs`; `flexLaneReaderFullThread()`
   already skips alignment drain when `noAlign` is true.

4. Provide a no-genome merge/output target.

   The current fused path creates per-thread `SoloReadFeature` objects and then
   merges them into:

   ```cpp
   RAchunk[0]->RA->soloRead->readFeat[Gene]
   ```

   The no-genome path needs an equivalent non-RA owner for the final
   `SoloReadFeature` and downstream Solo/Flex output. Use the smallest
   compatible surface:

   - create one final `SoloReadFeature(SoloFeatureTypes::Gene, P, -1)`;
   - create per-thread `SoloReadFeature(SoloFeatureTypes::Gene, P, -(200+i))`;
   - merge per-thread inline hashes into the final feature;
   - call the same existing Solo/Flex output/finalization routines used after
     the current Flex pipeline completes.

   If the existing output routines are only reachable through `Solo` or
   `ReadAlign`, first extract a small helper rather than duplicating output
   logic.

5. Preserve existing accounting.

   Carry over the existing fused pipeline accounting:

   - `g_statsAll.resetN()`
   - `timeStartMap`
   - `g_statsAll.addStats(fusedStats[i])`
   - `g_statsAll.readN = state.counters.readsTotal`
   - lane read counts
   - Flex `PIPELINE_STATS`

6. Keep fallback behavior unchanged.

   If the strict predicate fails, execute the current STAR path unchanged:
   load genome, construct `ReadAlignChunk`, and call `mapThreadsSpawn()`.

## Validation Plan

Build hygiene:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
```

Smoke tests:

```bash
tests/run_cbq_flex_tiny_public_smoke.sh
```

Then run the SC2300771 100K Flex fixture in four count-only modes if inputs are
available:

- FASTQ gzip, current path
- FASTQ gzip, no-genome path
- CBQ compressed, current path
- CBQ compressed or level-0, no-genome path

Required checks:

- `Log.out` for no-genome mode must not contain the normal genome loading
  messages.
- `Log.out` must contain `Flex count-only no-genome: active`.
- `Log.out` must record the prebuilt Flex hash cache path used by the run.
- `Solo.out/Gene` must match the current count-only path byte-for-byte.
- `per_sample` Flex outputs must match byte-for-byte.
- `Barcodes.stats` must match byte-for-byte.
- `Log.final.out` read counts and hash-screen counters must match.

Production benchmark:

- Use the existing SC2300771 production no-align command from
  `docs/RUNBOOK_FLEX_CBQ_INPUT.md`.
- Use fresh output directories for every timed run.
- Do not run other benchmarks concurrently on the same host.
- Compare wall time, `Log.final.out` mapping time, and peak RSS with
  `/usr/bin/time -v`.

Expected benefit:

- Startup savings should be roughly the genome load cost, plausibly 30-60
  seconds depending on disk/cache state.
- Peak RSS should drop by the loaded genome/index footprint.
- Count output runtime should otherwise remain at parity with the current fused
  Flex pipeline.

## Acceptance Criteria

- No-genome mode is only active for strict Flex hash/count-only runs.
- No-genome mode requires a prebuilt, readable Flex hash cache and never builds
  hashes itself.
- Current Flex, STAR, CBQ, SLAM, process_features, and multiome paths are
  unchanged when the predicate is false.
- FASTQ and CBQ Flex count-only outputs match the current count-only path on
  smoke and 100K fixtures.
- Production SC2300771 no-align benchmark shows lower startup wall time and
  lower peak RSS, with identical counters.
- Documentation records benchmark paths, commit hash, command lines, wall time,
  and peak RSS.

## Completed Production Results

The implemented no-genome path was validated on 2026-05-29 with the full
SC2300771 FLEX count-only production surface:

```text
--flex yes
--flexPipeline yes
--flexPipelineNTriage 0
--flexPipelineNSolo 0
--flexNoAlign 1
--outSAMtype None
--soloFeatures Gene
```

FASTQ.gz recipes for this surface must use STAR's default internal gzip path:
omit `--readFilesCommand zcat`. For `.fastq.gz` inputs, the default
`--readFilesCommand -` logs:

```text
NOTE: .gz Fastx input detected with --readFilesCommand - ; using internal zlib streaming path.
```

Use `--readFilesCommand zcat` only for explicit legacy compatibility
experiments; it is not the default benchmark or production recipe.

Topline no-genome results on the SSD:

| Input | Wall time | Max RSS | Mapping speed | Output parity |
| --- | ---: | ---: | ---: | --- |
| FASTQ.gz, internal gzip | `10:17.97` | `44,071,412 KB` | `12,107.14M reads/hour` | PASS vs CBQ |
| Level-0 CBQ, whole-lane readers | `8:38.52` | `43,378,292 KB` | `14,538.29M reads/hour` | PASS vs FASTQ |
| Level-0 CBQ, indexed range readers | `7:22.46` | `47,982,080 KB` | `17,156.56M reads/hour` | counters match; parity assumed |

All no-genome runs processed `2,011,130,186` reads with identical Flex hash-screen
counters:

```text
triageKeep=1681459858
triageDeny=16111757
triageMiss=313558571
```

The FASTQ no-genome run was `1:39.45` slower than the original whole-lane
level-0 CBQ result, about `1.19x` the CBQ wall time. The indexed range reader
result is `1:16.06` faster than the original level-0 CBQ result (`1.17x`,
`14.7%` wall reduction) and `2:55.51` faster than the FASTQ.gz internal-gzip
result (`1.40x`, `28.4%` wall reduction). The FASTQ and original CBQ
no-genome outputs were byte-identical for `Solo.out/Gene`,
`Solo.out/Barcodes.stats`, and `per_sample_filtered`; the indexed range run was
accepted on matching counters for this benchmark checkpoint.

Artifact roots:

- FASTQ no-genome: `/tmp/star_flex_fastq_full_ssd_internalgzip_20260531T180321Z`
- CBQ no-genome/current comparison:
  `/tmp/star_flex_no_genome_full_ssd_20260529T213833Z`
- CBQ indexed range:
  `/tmp/star_flex_cbq_range_full_ssd_20260531T111048Z`

For comparison, the same SSD CBQ command forced through the current genome-load
path with `STAR_DISABLE_FLEX_NO_GENOME=1` took `9:19.20` and peaked at
`84,395,816 KB`. The no-genome CBQ path therefore saved `40.68` seconds and
about `39.1 GiB` peak RSS while preserving byte-identical count outputs.

## Risk Notes

- The main risk is Solo/Flex output ownership. The current pipeline relies on
  `ReadAlignChunk`/`ReadAlign` as the owner of the final Solo read feature. The
  no-genome path must replace that ownership cleanly, not duplicate output
  logic.
- Do not loosen the guard to cover any path that can emit alignments or BAM/SAM
  headers. Those paths need genome-derived header state.
- If startup code has hidden assumptions that `samHeaders()` has already run,
  keep no-genome output helpers narrow and add explicit log checks in smoke
  tests.

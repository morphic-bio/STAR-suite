# Handoff: Flex No-Genome Count-Only Optimization

Date: 2026-05-29

## Assignment

Implement a strict STAR-Flex count-only path that avoids loading the STAR genome
when the run is hash/count-only:

```text
--flex yes
--flexPipeline yes
--flexPipelineNTriage 0
--flexPipelineNSolo 0
--flexNoAlign 1
--outSAMtype None
```

Use `docs/RUNBOOK_FLEX_NO_GENOME_COUNT_ONLY_20260529.md` as the implementation
runbook.

## Starting Point

Current behavior:

- `core/legacy/source/STAR.cpp` loads the genome before Flex pipeline dispatch.
- `core/legacy/source/mapThreadsSpawn.cpp` decides whether the Flex pipeline is
  active only after `ReadAlignChunk` objects have been built.
- `core/legacy/source/FlexPipeline.cpp` already skips alignment in no-align
  mode, but the startup path has already paid the genome load and
  `ReadAlignChunk` cost.

Important code landmarks:

- `STAR.cpp`: genome load around the `Genome genomeMain(P, P.pGe);` block.
- `mapThreadsSpawn.cpp`: `flexPipelineActivationGuard()` and
  `mapThreadsSpawnFlexPipeline()`.
- `FlexPipeline.cpp`: `flexLaneReaderFullThread()`, `processOneLane()`,
  `processOneCbqLane()`.
- `SoloReadFeature.cpp`: constructor and inline hash ownership.

## Recommended Branch

Use a clean worktree, not the dirty local checkout, if branch operations or
large benchmarks are needed.

Suggested branch:

```bash
git checkout master
git pull --ff-only
git checkout -b feature/flex-no-genome-count-only
```

Do a clean rebuild before debugging any crash or regression:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
```

## Implementation Shape

Add an early predicate before genome loading in `STAR.cpp`. If it passes, call a
new no-genome Flex count-only entrypoint and exit cleanly after outputs are
written. If it fails, continue through the current STAR path unchanged.

The first implementation should support only fully fused no-align mode. In that
mode:

- require a prebuilt, readable Flex hash cache before activating no-genome mode;
- create per-thread `SoloReadFeature` objects directly;
- pass `RA=nullptr` to `flexLaneReaderFullThread()`;
- reuse existing FASTQ and CBQ lane readers;
- merge per-thread inline hashes into a final no-genome-owned
  `SoloReadFeature`;
- call existing Solo/Flex output/finalization helpers rather than duplicating
  output logic.

Do not try to support alignment misses, BAM/SAM, or split pipeline workers in
this pass. Do not build hash caches in no-genome mode; this optimization is for
the already-built-hashes production count surface.

## Validation Checklist

Use fresh output directories for each timed run.

1. Clean build succeeds.
2. Existing Flex CBQ tiny smoke still passes:

   ```bash
   tests/run_cbq_flex_tiny_public_smoke.sh
   ```

3. 100K SC2300771 Flex count-only FASTQ current path vs no-genome path:

   - `Solo.out/Gene` byte-identical.
   - `per_sample` byte-identical.
   - `Barcodes.stats` byte-identical.
   - `Log.final.out` counters identical.
   - no-genome `Log.out` does not contain normal genome loading messages.
   - no-genome `Log.out` records the prebuilt Flex hash cache path.

4. 100K SC2300771 Flex count-only CBQ current path vs no-genome path with the
   same parity checks.

5. Production SC2300771 CBQ level-0 no-align benchmark:

   - record `/usr/bin/time -v` wall time and max RSS;
   - compare with the current documented benchmark in
     `docs/RUNBOOK_FLEX_CBQ_INPUT.md`;
   - confirm counters match exactly.

## Expected Result

The count-only run should save the fixed genome-load startup cost and reduce
peak RSS by the loaded genome/index footprint. The Flex read/hash/count phase
should remain at parity with the current fused pipeline.

## Completion Notes

Implementation completed on 2026-05-29.

The no-genome count-only path activates before `Genome genomeMain(...)`, logs
`Flex count-only no-genome: active`, runs the fully fused Flex lane readers with
no `ReadAlign`/`ReadAlignChunk` objects, and reuses the existing Solo/Flex
output path through prefilled per-thread `SoloReadFeature` objects. Commands
that do not satisfy the strict count-only predicate fall through to the
existing STAR genome-loading path unchanged.

Full SC2300771 SSD production results:

| Input | Wall time | Max RSS | Mapping speed |
| --- | ---: | ---: | ---: |
| FASTQ.gz no-genome | `12:01.95` | `43,281,404 KB` | `10,538.67M reads/hour` |
| Level-0 CBQ no-genome | `8:38.52` | `43,378,292 KB` | `14,538.29M reads/hour` |
| Level-0 CBQ current path | `9:19.20` | `84,395,816 KB` | `13,869.86M reads/hour` |

FASTQ no-genome and CBQ no-genome produced byte-identical `Solo.out/Gene`,
`Solo.out/Barcodes.stats`, and `per_sample_filtered` outputs. Both no-genome
runs reported:

```text
Number of input reads = 2011130186
Flex pipeline complete: total=2011130186, triageKeep=1681459858, triageDeny=16111757, triageMiss=313558571
```

Completed artifact roots:

- `/tmp/star_flex_fastq_full_ssd_20260529T222432Z`
- `/tmp/star_flex_no_genome_full_ssd_20260529T213833Z`

## Stop Conditions

Stop and ask for review if:

- the no-genome path requires duplicating large Solo output routines;
- any alignment or BAM/SAM output case starts using the no-genome path;
- the implementation tries to generate Flex hashes inside no-genome mode;
- no-genome outputs differ from the current count-only path on the 100K fixture;
- the implementation needs broad changes outside `STAR.cpp`,
  `mapThreadsSpawn.cpp`, `FlexPipeline.*`, or small Solo output helper
  extraction.

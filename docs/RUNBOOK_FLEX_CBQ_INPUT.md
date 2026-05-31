# FLEX CBQ Input Runbook

## Goal

Support CBQ/BINSEQ input for STAR-Flex without changing the optimized FASTQ
path. The first production-capable path should be small and conservative:
decode CBQ records through the shared C++ reader, fill the same STAR read
buffers that FASTQ fills, and let the existing FLEX/Solo/alignment consumers run
unchanged.

## Integration Points

`STAR --readFilesType Binseq PE --flex yes ...` enters the normal
`ReadAlignChunk::processChunks()` path. For CBQ input, that path calls
`ReadAlignChunk::mapCbqChunk()`, which uses `CbqStarAdapter` to fill:

- `Read0`: per-mate ASCII sequence.
- `Read1`: per-mate STAR numeric bases.
- `Qual0`: per-mate PHRED quality.
- `readNameMates`, `readNameExtra`, `readLength`, `readLengthOriginal`,
  `readFilter`, and `readFilesIndex`.

After those fields are populated, `ReadAlign::oneReadLoaded()` runs the same
FLEX hash-screen, barcode extraction, and alignment code used by FASTQ.

The production no-align benchmark path uses the fully fused `FlexPipeline`
instead. In that mode each lane is read through the shared `CbqInputModule`,
materialized into the same bounded per-read scratch buffers used by the FASTQ
lane reader, and then passed directly to the existing FLEX hash-screen and Solo
barcode/counting routines. This avoids temporary FASTQ materialization and
keeps the optimized FASTQ reader unchanged.

## FLEX-Specific Consumers

FLEX cannot be treated as a pure two-bit alignment consumer yet. The existing
consumer surface still uses text sequence and qualities in several places:

- `SoloReadBarcode::getCBandUMI()` reads CB/UMI sequence and quality from
  `Read0` / `Qual0`.
- `ReadAlign::detectSampleFromRawR2()` consumes raw R2 text sequence for sample
  detection.
- `FlexHashScreenCache::classifyRead()` consumes raw R2 text sequence for the
  H0/H1 hash-screen path.
- Y/noY and unmapped FASTQ emission paths use `Read0` if those outputs are
  enabled.

For this phase, the CBQ adapter must therefore fill both the STAR numeric bases
and the bounded per-read ASCII buffers. This is still cheaper and safer than
materializing temporary FASTQ files.

## FlexPipeline Gate

The optimized `FlexPipeline` path closes STAR's normal read handles and reopens
lane files directly. FASTQ lanes use `gzopen()`. CBQ lanes must use
`CbqInputModule`.

CBQ is enabled only for the fully fused pipeline shape:

```text
--flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0
```

If CBQ is requested with split triage/Solo pipeline workers, the guard disables
the pipeline and logs:

```text
Flex pipeline: not active (CBQ/Binseq input currently requires fully-fused mode: --flexPipelineNTriage 0 --flexPipelineNSolo 0)
```

Non-pipeline FLEX CBQ runs continue to use the standard STAR CBQ adapter path.

## Smoke Plan

1. Build STAR plus CBQ utilities:

   ```bash
   make -C core/legacy/source -j8 STAR cbq-ordered-encoder cbq-star-adapter-harness
   ```

2. Re-run the generic CBQ STAR mapper smoke:

   ```bash
   tests/run_cbq_star_input_smoke.sh
   ```

3. Run the public tiny FLEX CBQ smoke:

   ```bash
   tests/run_cbq_flex_tiny_public_smoke.sh
   ```

   This builds the generated public tiny FLEX fixture, encodes cDNA R2 plus
   barcode R1 into an ordered CBQ, runs FASTQ and CBQ STAR-Flex, compares
   `Barcodes.stats`, and compares BAM bodies when `samtools` is available.
   This tiny fixture uses the standard non-hash-screen path, so the key log
   assertion is `readFilesType Binseq PE` plus `--flex yes`, not a hash-screen
   pipeline activation message.

4. For production-scale validation, encode a FLEX FASTQ fixture with
   `cbq_ordered_encoder`, using the same
   mate order as the working FLEX command:

   ```bash
   core/legacy/source/cbq_ordered_encoder \
     --readFilesIn flex_R2.fastq.gz flex_R1.fastq.gz \
     --outFile flex_R2_R1.cbq
   ```

5. Run the FLEX command twice, once from FASTQ and once from CBQ:

   ```bash
   STAR ... --flex yes --readFilesIn flex_R2.fastq.gz flex_R1.fastq.gz ...
   STAR ... --flex yes --readFilesType Binseq PE --readFilesIn flex_R2_R1.cbq ...
   ```

6. Compare stable outputs by read/count identity. Use source-order-preserving
   CBQ when testing chunk-boundary parity.

## Host-Local Validation

On 2026-05-29, the host-local SC2300771 100K FLEX fixture passed FASTQ-vs-CBQ
parity through the standard STAR CBQ adapter path:

- Config: `/storage/SC2300771_filtered_100K/cellranger/outs/config.csv`
- FASTQs: `/storage/downsampled_100K/SC2300771`
- Output root: `/tmp/star_suite_cbq_flex_100k_20260529T123024Z`
- CBQ encoding: eight ordered lane CBQs, cDNA R2 then barcode R1, 100,000
  records per lane.
- Threads: 4
- FASTQ wall time: `elapsed=1:00.79 user=141.71 sys=22.70 maxrss_kb=48047308`
- CBQ wall time: `elapsed=0:57.66 user=137.87 sys=21.12 maxrss_kb=48460408`
- `Solo.out/Gene` directory: byte-identical.
- `per_sample` FLEX outputs: byte-identical.
- `Solo.out/Barcodes.stats`: byte-identical.
- Unsorted BAM body: record order differs, but sorted SAM body is
  byte-identical.

The standard-adapter log confirms CBQ uses the non-pipeline path. Depending on
the surrounding Flex options, the reason can be the ordinary hash-screen gate or
the CBQ-specific fully-fused requirement:

```text
Flex pipeline: not active (CBQ/Binseq input uses the standard STAR CBQ adapter path)
```

The fully fused CBQ pipeline was also smoke-tested on eight 100K lane CBQs:

- Output root: `/tmp/star_suite_cbq_flex_pipeline_smoke_20260529T125208Z`
- Command surface: `--readFilesType Binseq PE`, `--flexPipeline yes`,
  `--flexPipelineNTriage 0`, `--flexPipelineNSolo 0`, `--flexNoAlign 1`,
  `--outSAMtype None`.
- Threads: 32
- Input reads: `800000`
- Hash screen evaluated: `655916`
- KEEP: `649698` (`99.05%`)
- DENY: `6218` (`0.95%`)
- Wall time: `elapsed=0:48.84 user=90.31 sys=15.75 maxrss_kb=38657748`
- The log confirms the fully fused CBQ pipeline:

  ```text
  Flex pipeline: runThreadN=32, nLanes=8, triage=0 (fully fused + lane steal + role switch), soloConsumers=0, fusedThreads=32, dedicatedWorkers=0, noAlign=ON (alignment skipped)
  ```

### Indexed CBQ range smoke

On 2026-05-31, the no-genome count-only CBQ path was extended to partition
indexed CBQ input by logical record range. This is intentionally narrower than
the general STAR CBQ adapter: it activates only for strict no-genome,
no-alignment, no-BAM FLEX count-only runs. If indexed range setup fails, the
pipeline falls back to the existing whole-lane CBQ readers.

Fresh same-binary 800K parity pair:

- Artifact root: `/tmp/star_suite_cbq_flex_range_100k_20260531T105047Z`
- FASTQ.gz run:
  `/tmp/star_suite_cbq_flex_range_100k_20260531T105047Z/run_fastq_no_genome_parity_20260531T105739Z`
- CBQ range run:
  `/tmp/star_suite_cbq_flex_range_100k_20260531T105047Z/run_cbq_range_no_genome_clean_20260531T105506Z`
- Threads: 32
- Input reads: `800000`
- FASTQ.gz wall/RSS: `0:33.19`, `3,412,488 KB`
- CBQ range wall/RSS: `0:33.13`, `3,619,508 KB`
- `Solo.out/Gene`: byte-identical.
- `Solo.out/Barcodes.stats`: byte-identical.
- `per_sample_filtered`: byte-identical.
- Hash-screen counters matched:
  `total=800000, triageKeep=649698, triageDeny=6218, triageMiss=144084`.

The CBQ log confirms both the no-genome guard and indexed range planner:

```text
Flex count-only no-genome: active
Flex CBQ range: active (32 ranges across 8 lanes and 800000 records)
```

This 800K run is a correctness gate, not a production performance claim; the
input is small and cache effects dominate. Use the full SC2300771 production
surface below for topline timing.

## Full Production Benchmark

### Topline no-genome FLEX results

On 2026-05-31, SC2300771 was benchmarked on the full production FLEX
count-only surface with the no-genome path active and FASTQ.gz using STAR's
internal gzip path. These are the current topline FLEX results for production
count-only runs:

For FASTQ.gz inputs, production recipes should omit `--readFilesCommand zcat`.
The default `--readFilesCommand -` uses STAR's internal gzip streaming path;
zcat is only a legacy compatibility override.

| Input | Wall time | Max RSS | Mapping speed | Input reads |
| --- | ---: | ---: | ---: | ---: |
| FASTQ.gz, internal gzip | `10:17.97` | `44,071,412 KB` | `12,107.14M reads/hour` | `2,011,130,186` |
| Level-0 CBQ | `8:38.52` | `43,378,292 KB` | `14,538.29M reads/hour` | `2,011,130,186` |

The level-0 CBQ run was `1:39.45` faster than FASTQ.gz on the same SSD-hosted
production surface, about `1.19x` faster by wall time, while using effectively
the same memory. This is the clearest production use case for CBQ in FLEX:
count-only no-genome runs avoid STAR genome load, and CBQ removes the remaining
FASTQ gzip decode bottleneck.

FASTQ no-genome and CBQ no-genome produced byte-identical outputs for:

- `Solo.out/Gene`
- `Solo.out/Barcodes.stats`
- `per_sample_filtered`

Both runs produced identical hash-screen counters:

```text
Flex pipeline complete: total=2011130186, triageKeep=1681459858, triageDeny=16111757, triageMiss=313558571
```

Artifact roots:

- FASTQ no-genome:
  `/tmp/star_flex_fastq_full_ssd_internalgzip_20260531T180321Z`
- CBQ no-genome and current-path comparison:
  `/tmp/star_flex_no_genome_full_ssd_20260529T213833Z`

For an apples-to-apples no-genome implementation check, the same level-0 CBQ
command was also forced through the current genome-load path with
`STAR_DISABLE_FLEX_NO_GENOME=1`. That run took `9:19.20` and peaked at
`84,395,816 KB`; the no-genome CBQ run took `8:38.52` and peaked at
`43,378,292 KB`.

### Full CBQ indexed-range benchmark

On 2026-05-31, the full SC2300771 level-0 CBQ no-genome count-only benchmark
was rerun from the existing SSD-staged CBQs with indexed CBQ range readers
enabled.

- Output root: `/tmp/star_flex_cbq_range_full_ssd_20260531T111048Z`
- Input list:
  `/tmp/star_suite_cbq_flex_full_level0_p8_20260529T125516Z/cbq_list.txt`
- Threads: 32
- Range planner: `39 ranges across 8 lanes and 2011130186 records`
- Wall time: `7:22.46`
- Max RSS: `47,982,080 KB`
- STAR mapping speed: `17,156.56M reads/hour`
- Input reads: `2,011,130,186`
- Counter sanity check:
  `triageKeep=1681459858, triageDeny=16111757, triageMiss=313558571`

Compared with the previous full level-0 CBQ no-genome result (`8:38.52`,
`43,378,292 KB`, `14,538.29M reads/hour`), indexed range reading improved wall
time by `1:16.06` (`1.17x`, `14.7%` reduction) and mapping speed by `18.0%`,
with about `4.39 GiB` higher max RSS. FASTQ baseline was not rerun for this
checkpoint; the corrected internal-gzip FASTQ rerun above is now the current
FASTQ topline.

### Earlier full production CBQ benchmark

On 2026-05-29, SC2300771 was benchmarked on the full production FLEX no-align
surface from uncompressed level-0 CBQ on the SSD.

Run artifacts:

- CBQ root: `/tmp/star_suite_cbq_flex_full_level0_p8_20260529T125516Z`
- STAR run:
  `/tmp/star_suite_cbq_flex_full_level0_p8_20260529T125516Z/runs/SC2300771_cbq_level0_no_bam_noalign_20260529T130646Z`
- FASTQ baseline:
  `/mnt/pikachu/JAX_scRNAseq01_2024_star_suite_20260524_080700/SC2300771_hash_noalign`
- STAR binary: `core/legacy/source/STAR`
- Genome index: `/storage/flex_filtered_reference_2024/star_index`
- Threads: 32
- Command surface: `--readFilesType Binseq PE`, `--flex yes`,
  `--flexPipeline yes`, `--flexPipelineNTriage 0`, `--flexPipelineNSolo 0`,
  `--flexNoAlign 1`, `--outSAMtype None`, `--outSAMattributes None`.

Preprocessing note: the eight source FASTQ lane pairs were converted to ordered
level-0 CBQ files in parallel with `cbq_ordered_encoder`; per-lane encode wall
times were about 10:31-10:58. Encoding time is not included in the STAR
benchmark wall time.

Input size comparison:

| Input representation | Files | Bytes | GiB | GB |
| --- | ---: | ---: | ---: | ---: |
| Source FASTQ.gz | 16 | 150,990,427,659 | 140.621 | 150.990 |
| Ordered level-0 CBQ | 8 | 106,937,671,967 | 99.593 | 106.938 |

The uncompressed level-0 CBQ set is 70.82% of the gzipped FASTQ size, saving
44,052,755,692 bytes (41.027 GiB, 44.053 GB), or 29.18%.

Timing comparison:

| Run | Total wall | STAR log job wall | STAR mapping phase | Mapping speed |
| --- | ---: | ---: | ---: | ---: |
| FASTQ.gz baseline | not captured by `/usr/bin/time` | 11:46 | 11:18 | 10,678.57M reads/hour |
| Level-0 CBQ | 9:04.62 | 8:58 | 8:31 | 14,168.43M reads/hour |

The CBQ run reduced the STAR log job wall time by 2:48 and the mapping phase by
2:47. The last `PIPELINE_STATS` line before lane completion was at 370 seconds,
so the fused CBQ read/hash phase itself finished in about 6:10; the remaining
wall time is genome load and final Solo/FLEX output bookkeeping.

Production counter parity against the FASTQ baseline:

| Counter | FASTQ.gz baseline | Level-0 CBQ |
| --- | ---: | ---: |
| Number of input reads | 2,011,130,186 | 2,011,130,186 |
| Hash screen reads evaluated | 1,697,571,615 | 1,697,571,615 |
| Hash screen KEEP | 1,681,459,858 | 1,681,459,858 |
| Hash screen KEEP % | 99.05% | 99.05% |
| Hash screen DENY | 16,111,757 | 16,111,757 |
| Hash screen DENY % | 0.95% | 0.95% |
| Hash screen PASS | 0 | 0 |

## Follow-Up

- Reduce the ASCII surface later by adding direct packed-sequence consumers for
  sample detection and hash screening. Keep this separate from initial
  correctness work.
- Add a repeatable production benchmark wrapper if this becomes a release-gate
  benchmark rather than a one-off host-local validation.

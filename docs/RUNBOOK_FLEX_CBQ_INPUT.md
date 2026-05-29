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

The standard-adapter log confirms CBQ uses the non-pipeline path:

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

## Full Production Benchmark

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

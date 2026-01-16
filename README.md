# STAR-suite

STAR-suite reorganizes STAR into module-focused directories while keeping a
single source of truth for shared code. Build outputs remain compatible with
existing STAR workflows, and the new top-level `Makefile` exposes module targets.

## Modules

- STAR-core (`core/`): Legacy STAR (indexing, bulk, Solo) plus shared utilities.
  Build: `make core` (binary at `core/legacy/source/STAR`).
- STAR-Flex (`flex/`): FlexFilter pipeline and Flex-specific integrations.
  Build tools: `make flex` or `make flex-tools`.
- STAR-SLAM (`slam/`): SLAM-seq quantification, SNP masking, trimming/QC.
  Build tools: `make slam` or `make slam-tools`.

## Flags (high-level)

- Core (legacy): Standard STAR flags. See `core/legacy/README.md`.
- Flex: `--flex yes` to enable Flex mode; `--soloRunFlexFilter yes` to run the
  filter-only pipeline on a MEX; `--soloType`, `--soloCB*`, `--soloUMI*` as in
  STARsolo. See `flex/README_flex.md`.
- SLAM: `--slamQuantMode`, `--slamSnpBed`, `--slamSnpDetect`,
  `--slamCompatMode gedi`, `--slamCompatOverlapWeight`, and trimming options
  like `--clip3pAdapterSeq`/`--clip3pAdapterMMp`. See `docs/SLAM_COMPATIBILITY_MODE.md`.

## Sample Command

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix out/ \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM MD \
  --slamQuantMode 1 \
  --slamSnpBed /path/to/snps.bed
```

Remove the `--slam*` flags to run a core-only alignment.

## More Detail

- Core usage: `core/legacy/README.md`
- Flex pipeline: `flex/README_flex.md`
- SLAM compatibility: `docs/SLAM_COMPATIBILITY_MODE.md`
- SLAM tools: `slam/tools/slam_requant/README.md`, `slam/tools/pileup_snp/README.md`

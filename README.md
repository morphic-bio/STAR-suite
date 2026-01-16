# STAR-suite

STAR-suite reorganizes STAR into module-focused directories while keeping a
single source of truth for shared code. Build outputs remain compatible with
existing STAR workflows, and the new top-level `Makefile` exposes module targets.

## Folder Structure

```
core/
  legacy/                # Upstream STAR layout (single source of truth)
  features/              # Shared feature overlays (vbem, yremove, bamsort)
flex/                    # Flex-specific code + tools
slam/                    # SLAM-seq code + tools
build/                   # Modular make fragments
docs/                    # Suite-level docs
tests/                   # Suite-level tests
tools/                   # Suite-level scripts/utilities
```

## Modules

- STAR-core (`core/`): Legacy STAR (indexing, bulk, Solo) plus shared utilities.
  Build: `make core` (binary at `core/legacy/source/STAR`).
- STAR-Flex (`flex/`): FlexFilter pipeline and Flex-specific integrations.
  Build tools: `make flex` or `make flex-tools`.
- STAR-SLAM (`slam/`): SLAM-seq quantification, SNP masking, trimming/QC.
  Build tools: `make slam` or `make slam-tools`.

## Flags (high-level)

- Core (legacy): Standard STAR flags (`--runMode`, `--genomeDir`,
  `--readFilesIn`, `--outSAMtype`, `--outSAMattributes`, etc.).
  See `core/legacy/README.md`.
- Flex: `--flex yes` to enable Flex mode; `--soloRunFlexFilter yes` to run the
  filter-only pipeline on a MEX; `--soloType`, `--soloCB*`, `--soloUMI*` as in
  STARsolo. See `flex/README_flex.md`.
- SLAM: `--slamQuantMode`, `--slamSnpBed`, `--slamSnpDetect`,
  `--slamCompatMode gedi`, `--slamCompatOverlapWeight`, and trimming options
  like `--clip3pAdapterSeq`/`--clip3pAdapterMMp`. See `docs/SLAM_COMPATIBILITY_MODE.md`.

## Sample Commands

Core alignment:

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix out/ \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM MD
```

Flex mode (inline filtering):

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix out/ \
  --outSAMtype BAM SortedByCoordinate \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloCBwhitelist /path/to/whitelist.txt \
  --flex yes
```

SLAM mode:

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

## More Detail

- Core usage: `core/legacy/README.md`
- Flex pipeline: `flex/README_flex.md`
- SLAM compatibility: `docs/SLAM_COMPATIBILITY_MODE.md`
- SLAM tools: `slam/tools/slam_requant/README.md`, `slam/tools/pileup_snp/README.md`

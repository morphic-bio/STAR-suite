# STAR Suite Compile Instructions

This document provides explicit source compile commands for each module, the
full suite, and clean rebuild workflows.

All commands are run from repo root.

## Ubuntu 24.04 Build Prerequisites

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential \
  gcc \
  g++ \
  make \
  xxd \
  cmake \
  pkg-config \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-gnutls-dev \
  libssl-dev \
  libglib2.0-dev \
  libhts-dev \
  git \
  ca-certificates
```

Notes:
- Python is not required for compiling core suite binaries.
- `libhts-dev` is needed by feature barcode tooling (`demux_bam` includes htslib headers).

## Parallel Jobs

Use 8 to 16 threads for practical compile time:

```bash
export MAKE_JOBS=8
```

## Compile Single Module/Target

```bash
# Core STAR
make -j"${MAKE_JOBS}" core

# Flex (builds core + flex-tools)
make -j"${MAKE_JOBS}" flex

# SLAM (builds core + slam-tools)
make -j"${MAKE_JOBS}" slam

# Feature barcode tools
make -j"${MAKE_JOBS}" feature-barcodes-tools

# Common default set
make -j"${MAKE_JOBS}" default

# Full suite
make -j"${MAKE_JOBS}" all
```

## Clean Rebuild Workflows

`make clean` removes core and tool build outputs (`core-clean` + `tools-clean`).

### Clean rebuild of a specific module

```bash
make clean
make -j"${MAKE_JOBS}" core
```

```bash
make clean
make -j"${MAKE_JOBS}" flex
```

```bash
make clean
make -j"${MAKE_JOBS}" slam
```

### Clean rebuild of entire suite

```bash
make clean
make -j"${MAKE_JOBS}" all
```

## Output Locations

- `core/legacy/source/STAR`
- `core/legacy/source/star_feature_call`
- `flex/tools/flexfilter/run_flexfilter_mex`
- `slam/tools/slam_requant/slam_requant`
- `slam/tools/pileup_snp/pileup_snp`
- `core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads`

## Clean Ubuntu 24.04 Validation

The repository Docker builder stage validates compilation in a clean Ubuntu 24.04 environment:

```bash
docker build --no-cache --target builder -f docker/Dockerfile --build-arg MAKE_JOBS=8 .
```

The builder stage runs:
- `make core`
- `make flex`
- `make slam`
- `make feature-barcodes-tools`
- `make default`
- `make all`

Optional explicit clean-rebuild validation inside the clean builder image:

```bash
docker build --target builder -f docker/Dockerfile --build-arg MAKE_JOBS=8 -t star-suite-builder-compilecheck .
docker run --rm star-suite-builder-compilecheck bash -lc 'cd /build && make clean && make -j8 core && make clean && make -j8 all'
```

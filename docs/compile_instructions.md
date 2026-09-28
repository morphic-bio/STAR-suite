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
  libdeflate-dev \
  git \
  ca-certificates
```

Notes:
- Python is not required for compiling core suite binaries.
- `libhts-dev` is needed only by standalone feature barcode tooling
  (`demux_bam` includes HTSlib headers) and by `make core HTSLIB=external`.
  `make core` uses STAR's bundled HTSlib.

## HTSlib Selection

`make core` builds STAR with its bundled HTSlib and needs no other checkout;
STAR Suite links no other suite. `make core-portable` is kept as an alias.

```bash
make core-clean
make -j8 core
core/legacy/source/STAR --version
```

`HTSLIB=external` compiles and links against an installed HTSlib found with
`pkg-config --cflags/--libs htslib`. Programs that embed STAR together with
other libraries built against an installed HTSlib use this setting so that one
HTSlib ABI is used throughout the executable. For a custom
installation, set `PKG_CONFIG_PATH=/path/to/htslib/lib/pkgconfig`. Without a
`.pc` file, supply both sides explicitly:

```bash
make -j8 core HTSLIB=external HTSLIB_CFLAGS=-I/path/to/htslib/include \
  HTSLIB_LIBS=/path/to/htslib/lib/libhts.so
```

`CPPFLAGS` and `CXXFLAGSextra` are honored during dependency scanning as well as
compilation. Set the runtime library search path separately when installing a
shared HTSlib outside the system loader's paths. Never mix the older bundled
STAR HTSlib headers with an external library.

With `HTSLIB=external` the build checks the installed HTSlib before scanning
source dependencies and reports missing prerequisites explicitly. A failed scan
leaves any existing `Depend.list` intact instead of retaining partially
generated dependencies.

Use a clean build when changing HTSlib installations or the `HTSLIB` setting.
Release installer archives and `.deb` packages contain prebuilt binaries and
do not require compilation.

Before 1.10.0, `make core` linked Chromap Suite for multiome runs. That
integration, including the `--chromapAtac*` and `--multiomeAtac*` parameters,
moved to Multiomics Suite, which builds the multiome binary.

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

The builder stage validates:
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

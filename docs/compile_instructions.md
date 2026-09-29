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
- `libhts-dev` is needed by the default Chromap-enabled core build and standalone
  feature barcode tooling (`demux_bam` includes HTSlib headers). The explicit
  `make core-portable` target uses STAR's bundled HTSlib instead.

## Chromap and HTSlib Discovery

`make core` remains Chromap-enabled; it does not silently disable ATAC when a
dependency is missing. It expects Chromap Suite (with initialized submodules) at
`../Chromap-suite`, or at the path supplied with `CHROMAP_SUITE_DIR`:

```bash
git clone --recursive https://github.com/morphic-bio/Chromap-suite.git ../Chromap-suite
make -j8 core CHROMAP_SUITE_DIR="$(realpath ../Chromap-suite)"
```

The core build uses `pkg-config --cflags/--libs htslib` for external HTSlib.
For a custom installation, set `PKG_CONFIG_PATH=/path/to/htslib/lib/pkgconfig`.
Without a `.pc` file, supply both sides explicitly:

```bash
make -j8 core CHROMAP_HTSLIB_CFLAGS=-I/path/to/htslib/include \
  CHROMAP_SYS_HTS=/path/to/htslib/lib/libhts.so
```

`CPPFLAGS` and `CXXFLAGSextra` are honored during dependency scanning as well as
compilation. Set the runtime library search path separately when installing a
shared HTSlib outside the system loader's paths. Never mix the older bundled
STAR HTSlib headers with the external library used by Chromap.

The build checks external HTSlib before scanning source dependencies and reports
missing prerequisites explicitly. A failed scan leaves any existing `Depend.list`
intact instead of retaining partially generated dependencies.

For RNA/Flex/SLAM without Chromap, including STARsolo poly-G trimming:

```bash
make core-clean
make -j8 core-portable
core/legacy/source/STAR --version
```

Use a clean build when changing HTSlib installations or `WITH_CHROMAP` mode.
Release installer archives and `.deb` packages contain prebuilt binaries and
do not require compilation.

For the native ATAC peak-matrix helper, runtime capability checks, and the
Launchpad UI, follow [Multiome source build and launch](LAUNCHPAD_MULTIOME.md).

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

The builder stage runs with `STAR_WITH_CHROMAP=0` by default because the
single-repo Docker context does not include the sibling Chromap-suite checkout.
Local builds use the Chromap-enabled multiome binary by default with
`make core`; use `make core-portable` or `make core WITH_CHROMAP=0` only for an
explicit no-Chromap compatibility build.

The builder stage validates:
- `make core WITH_CHROMAP=${STAR_WITH_CHROMAP}`
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

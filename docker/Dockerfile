# STAR-suite Docker image (Ubuntu 24.04)
# Multi-stage:
#   - builder   : compile and validate Makefile targets
#   - suite-base: minimal executable suite runtime (no Python)
#   - test-tier-a / test-tier-b: helper images for smoke tests

# -----------------------------------------------------------------------------
# Stage 1: Builder
# -----------------------------------------------------------------------------
FROM ubuntu:24.04 AS builder

ARG MAKE_JOBS=16
ARG STAR_WITH_CHROMAP=0
ENV DEBIAN_FRONTEND=noninteractive
ENV MAKE_JOBS=$MAKE_JOBS
ENV STAR_WITH_CHROMAP=$STAR_WITH_CHROMAP

# Build dependencies
# - core/legacy uses vendored htslib
# - feature_barcodes/demux_bam includes <htslib/sam.h>, so libhts-dev is required
RUN apt-get update && apt-get install -y --no-install-recommends \
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
    ca-certificates \
    samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy source (layer ordering: deps before source for cache)
COPY . .

# Remove host-generated artifacts that may contain host paths (Depend.list, etc.)
RUN rm -f core/legacy/source/Depend.list core/legacy/source/*.o core/legacy/source/htslib/*.a 2>/dev/null || true

# Validate Makefile targets in order. This single-repo Docker context does not
# include the sibling Chromap-suite checkout, so it uses the explicit portable
# compatibility build unless STAR_WITH_CHROMAP=1 is provided with Chromap-suite
# available in the image.
RUN make -j"${MAKE_JOBS}" core WITH_CHROMAP="${STAR_WITH_CHROMAP}" \
    && make -j"${MAKE_JOBS}" flex WITH_CHROMAP="${STAR_WITH_CHROMAP}" \
    && make -j"${MAKE_JOBS}" slam WITH_CHROMAP="${STAR_WITH_CHROMAP}" \
    && make -j"${MAKE_JOBS}" feature-barcodes-tools \
    && make default WITH_CHROMAP="${STAR_WITH_CHROMAP}" \
    && make all WITH_CHROMAP="${STAR_WITH_CHROMAP}"

# Pre-build SLAM SNP EM test helper so Tier A can run without a compiler.
RUN g++ -std=c++11 -O2 \
    -I/build/slam/source -I/build/slam/source/libem \
    /build/tests/slam/test_slam_snp_em.cpp \
    /build/slam/source/libem/slam_snp_em.cpp \
    -o /build/tests/slam/test_slam_snp_em \
    -lm

# -----------------------------------------------------------------------------
# Stage 2: Suite Base Runtime (no test helpers, no Python)
# -----------------------------------------------------------------------------
FROM ubuntu:24.04 AS suite-base

ENV DEBIAN_FRONTEND=noninteractive

# Minimal runtime deps for core suite executables on Ubuntu 24.04.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 \
    libgomp1 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4t64 \
    libssl3t64 \
    libglib2.0-0t64 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy suite executables
COPY --from=builder /build/core/legacy/source/STAR /usr/local/bin/STAR
COPY --from=builder /build/core/legacy/source/star_feature_call /usr/local/bin/star_feature_call
COPY --from=builder /build/slam/tools/slam_requant/slam_requant /usr/local/bin/slam_requant
COPY --from=builder /build/slam/tools/pileup_snp/pileup_snp /usr/local/bin/pileup_snp
COPY --from=builder /build/flex/tools/flexfilter/run_flexfilter_mex /usr/local/bin/run_flexfilter_mex
COPY --from=builder /build/flex/tools/flexfilter/run_flexfilter_mex /usr/local/bin/flexfilter
COPY --from=builder /build/core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads /usr/local/bin/remove_y_reads

WORKDIR /workspace

# Default behavior: executable suite image.
CMD ["/usr/local/bin/STAR", "--version"]

# -----------------------------------------------------------------------------
# Stage 3: Tier A Test Image (self-contained smoke tests)
# -----------------------------------------------------------------------------
FROM suite-base AS test-tier-a

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    git-lfs \
    python3 \
    python3-pandas \
    samtools \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/tests /build/tests
COPY --from=builder /build/scripts /build/scripts
COPY --from=builder /build/flex/scripts /build/flex/scripts
WORKDIR /build

# -----------------------------------------------------------------------------
# Stage 4: Tier B Test Image (fixture-backed tests with extra helpers)
# -----------------------------------------------------------------------------
FROM suite-base AS test-tier-b

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    python3 \
    python3-pandas \
    bc \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/tests /build/tests
COPY --from=builder /build/scripts /build/scripts
COPY --from=builder /build/flex/scripts /build/flex/scripts
WORKDIR /build

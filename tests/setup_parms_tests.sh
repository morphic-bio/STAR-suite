#!/usr/bin/env bash
# Preflight setup for Parameters regression test suite.
# This script sets up environment variables for test execution.
# It should be sourced (not executed) so variables persist to the calling shell.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CORE_DIR="${ROOT_DIR}/core/legacy/source"

# Load external fixture paths (if present) so defaults are centralized.
EXTERNAL_FIXTURES_ENV="${ROOT_DIR}/tests/external_fixtures_env.sh"
if [[ -f "$EXTERNAL_FIXTURES_ENV" ]]; then
    # shellcheck source=/dev/null
    source "$EXTERNAL_FIXTURES_ENV"
fi

# Tool paths (can be overridden via environment)
export STAR_BIN="${STAR_BIN:-$CORE_DIR/STAR}"
export REMOVE_Y_READS_BIN="${REMOVE_Y_READS_BIN:-$ROOT_DIR/core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads}"
export TXIMPORT_BIN="${TXIMPORT_BIN:-$ROOT_DIR/core/features/vbem/tools/tximport_compat/tximport_compat}"
export TXIMPORT_DIR="${TXIMPORT_DIR:-$ROOT_DIR/core/features/vbem/tools/tximport_compat}"

# Fixture paths with fallbacks
SLAM_FIXTURE_ROOT="${SLAM_FIXTURE_ROOT:-$ROOT_DIR/test/fixtures/slam}"
FLEX_FIXTURE_ROOT="${FLEX_FIXTURE_ROOT:-$ROOT_DIR/reference/tests/100K}"
ALT_SLAM_FIXTURE="/mnt/pikachu/STAR-Flex/test/fixtures/slam"
ALT_FLEX_FIXTURE="/mnt/pikachu/STAR-Flex/reference/tests/100K"

warn() { echo "WARN: $*" >&2; }
info() { echo "INFO: $*"; }

# Check and set SLAM fixture root with fallback
if [[ ! -d "$SLAM_FIXTURE_ROOT" ]]; then
    if [[ -d "$ALT_SLAM_FIXTURE" ]]; then
        info "SLAM fixtures not found at $SLAM_FIXTURE_ROOT; using fallback: $ALT_SLAM_FIXTURE"
        export SLAM_FIXTURE_ROOT="$ALT_SLAM_FIXTURE"
    else
        warn "SLAM fixtures not found at $SLAM_FIXTURE_ROOT or $ALT_SLAM_FIXTURE"
        export SLAM_FIXTURE_ROOT="$SLAM_FIXTURE_ROOT"  # Export anyway for test scripts to check
    fi
else
    export SLAM_FIXTURE_ROOT="$SLAM_FIXTURE_ROOT"
fi

# Check and set FLEX fixture root with fallback
if [[ ! -d "$FLEX_FIXTURE_ROOT" ]]; then
    if [[ -d "$ALT_FLEX_FIXTURE" ]]; then
        info "Flex fixtures not found at $FLEX_FIXTURE_ROOT; using fallback: $ALT_FLEX_FIXTURE"
        export FLEX_FIXTURE_ROOT="$ALT_FLEX_FIXTURE"
    else
        warn "Flex fixtures not found at $FLEX_FIXTURE_ROOT or $ALT_FLEX_FIXTURE"
        export FLEX_FIXTURE_ROOT="$FLEX_FIXTURE_ROOT"  # Export anyway for test scripts to check
    fi
else
    export FLEX_FIXTURE_ROOT="$FLEX_FIXTURE_ROOT"
fi

# Check for required tools and warn if missing
MISSING_TOOLS=()

if ! command -v python3 >/dev/null 2>&1; then
    warn "python3 not found; some tests may fail"
    MISSING_TOOLS+=("python3")
fi

if ! command -v gzip >/dev/null 2>&1 && ! command -v zcat >/dev/null 2>&1; then
    warn "gzip/zcat not found; FASTQ tests may fail"
    MISSING_TOOLS+=("gzip/zcat")
fi

if ! command -v bedtools >/dev/null 2>&1; then
    warn "bedtools not found; some SLAM tests may fail"
    MISSING_TOOLS+=("bedtools")
fi

if ! command -v samtools >/dev/null 2>&1; then
    warn "samtools not found; Y-chrom tests will fail"
    MISSING_TOOLS+=("samtools")
fi

if ! command -v Rscript >/dev/null 2>&1; then
    warn "Rscript not found; tximport parity test will skip"
    MISSING_TOOLS+=("Rscript")
else
    if ! Rscript -e "if (!requireNamespace('tximport', quietly=TRUE)) quit(status=1)" >/dev/null 2>&1; then
        warn "R package 'tximport' not available; tximport parity test will skip"
    fi
fi

# Verify tool binaries exist (warn only, don't fail)
if [[ ! -x "$STAR_BIN" ]]; then
    warn "STAR binary not found at $STAR_BIN (will be built if needed)"
fi

if [[ ! -x "$REMOVE_Y_READS_BIN" ]]; then
    warn "remove_y_reads binary not found at $REMOVE_Y_READS_BIN (will be built if needed)"
fi

if [[ ! -x "$TXIMPORT_BIN" ]]; then
    warn "tximport_compat binary not found at $TXIMPORT_BIN (will be built if needed)"
fi

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    info "Missing optional tools: ${MISSING_TOOLS[*]} (tests may skip or fail)"
fi

info "Setup complete - environment variables exported"

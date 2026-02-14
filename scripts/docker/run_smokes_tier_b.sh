#!/bin/bash
# Run Tier B smoke tests (fixture-backed; requires /storage mount)
# Usage from host: ./scripts/docker/run_smokes_tier_b.sh
#   Mount fixtures: -v /path/to/storage:/storage
#   Or set STORAGE=/path before running (script will use it for docker -v)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BASE_IMAGE_TAG="${IMAGE_TAG:-biodepot/star-suite:latest}"
TEST_IMAGE_TAG="${TEST_IMAGE_TAG:-biodepot/star-suite:test-tier-b}"
MAKE_JOBS="${MAKE_JOBS:-16}"
# Default fixture root expected by most STAR-suite integration runs.
# Can be overridden with STORAGE=/path.
STORAGE="${STORAGE:-/storage}"

cd "$REPO_ROOT"

# Ensure Tier B test image exists (built from suite-base).
if ! docker image inspect "$TEST_IMAGE_TAG" &>/dev/null; then
    echo "Tier B test image $TEST_IMAGE_TAG not found. Building it..."
    docker build \
        --target test-tier-b \
        -f docker/Dockerfile \
        -t "$TEST_IMAGE_TAG" \
        --build-arg MAKE_JOBS="$MAKE_JOBS" \
        .
fi

# Check for fixtures
if [ ! -d "$STORAGE" ]; then
    echo "SKIP: Tier B fixtures not found at $STORAGE"
    echo "  Mount fixtures: docker run -v /path/to/data:/storage ..."
    echo "  Or create $STORAGE with expected layout (see plans/docker_plan.md)"
    exit 0
fi

echo "Running Tier B smoke tests..."
echo "  Base image tag: $BASE_IMAGE_TAG"
echo "  Test image tag: $TEST_IMAGE_TAG"
echo "  Storage: $STORAGE"
echo ""

failed=0
passed=0
skipped=0

# run_cr_compat_integration_smoke requires ASAN build + run_solo_smoke index - skip in standard image
echo "=== run_cr_compat_integration_smoke ==="
echo "SKIP: run_cr_compat_integration_smoke requires ASAN build (Tier B diagnostic path)"
skipped=$((skipped + 1))
echo ""
echo "=== run_cr_multi_smoke ==="
echo "SKIP: run_cr_multi_smoke depends on external downsample scripts not packaged in Docker test image"
skipped=$((skipped + 1))
echo ""

run_one() {
    local name="$1"
    local script="$2"
    local log_file
    log_file="$(mktemp /tmp/tierb_${name}.XXXXXX.log)"
    echo "=== $name ==="
    set +e
    docker run --rm \
        -v "${STORAGE}:/storage" \
        -e STAR_BIN=/usr/local/bin/STAR \
        -e CR_COMPAT_TEST_OUTPREFIX="/tmp/${name}_$(date +%Y%m%d_%H%M%S)/" \
        -e CBUB_SAMPLE_PROBES="${CBUB_SAMPLE_PROBES:-/storage/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}" \
        -e CBUB_TMP_BASE="${CBUB_TMP_BASE:-/tmp/cbub_regress_tmp}" \
        -e FLEX_SAMPLE_PROBES="${FLEX_SAMPLE_PROBES:-/storage/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}" \
        -w /build \
        "$TEST_IMAGE_TAG" \
        bash -c "$script" 2>&1 | tee "$log_file"
    local rc=${PIPESTATUS[0]}
    set -e

    if [ "$rc" -eq 0 ]; then
        if grep -q "SKIP:" "$log_file"; then
            echo "SKIP: $name"
            skipped=$((skipped + 1))
        else
            echo "PASS: $name"
            passed=$((passed + 1))
        fi
    else
        echo "FAIL: $name"
        failed=1
    fi
    rm -f "$log_file"
    echo ""
}

run_one "test_cr_compat_crispr_calling" "tests/test_cr_compat_crispr_calling.sh"
run_one "run_unsorted_cbub_smoke_test" "tests/run_unsorted_cbub_smoke_test.sh"
run_one "run_flex_smoke" "tests/flex_smoke/run_flex_smoke.sh"
run_one "run_cbub_regression_test" "tests/run_cbub_regression_test.sh"

if [ $failed -eq 0 ]; then
    echo "Tier B summary: pass=$passed skip=$skipped fail=0."
else
    echo "Tier B summary: pass=$passed skip=$skipped fail>=1."
    exit 1
fi

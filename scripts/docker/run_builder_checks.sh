#!/bin/bash
# Run Makefile validation inside the builder stage (no image save)
# Usage: ./scripts/docker/run_builder_checks.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
IMAGE_TAG="${IMAGE_TAG:-biodepot/star-suite:latest}"
MAKE_JOBS="${MAKE_JOBS:-16}"
STAR_SUITE_COMMIT_SHA="${STAR_SUITE_COMMIT_SHA:-$(git -C "${REPO_ROOT}" rev-parse HEAD)}"

cd "$REPO_ROOT"

echo "Running builder validation (make targets) in Docker..."
echo "  MAKE_JOBS: $MAKE_JOBS"
echo ""

# Build only the builder stage and run validation
docker build \
    --target builder \
    -f docker/Dockerfile \
    -t "${IMAGE_TAG}-builder" \
    --build-arg MAKE_JOBS="$MAKE_JOBS" \
    --build-arg STAR_SUITE_COMMIT_SHA="$STAR_SUITE_COMMIT_SHA" \
    .

echo ""
echo "Builder validation passed."

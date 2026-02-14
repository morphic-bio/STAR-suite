#!/bin/bash
# Build STAR-suite base Docker image (suite runtime only; no test helpers)
# Usage: ./scripts/docker/build_image.sh
# Override: IMAGE_TAG=biodepot/star-suite:v1 ./scripts/docker/build_image.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
IMAGE_TAG="${IMAGE_TAG:-biodepot/star-suite:latest}"
MAKE_JOBS="${MAKE_JOBS:-16}"

cd "$REPO_ROOT"

echo "Building STAR-suite base image: $IMAGE_TAG"
echo "  Repo root: $REPO_ROOT"
echo "  MAKE_JOBS: $MAKE_JOBS"
echo ""

docker build \
    --target suite-base \
    -f docker/Dockerfile \
    -t "$IMAGE_TAG" \
    --build-arg MAKE_JOBS="$MAKE_JOBS" \
    .

echo ""
echo "Built: $IMAGE_TAG"
echo "  Run: docker run --rm $IMAGE_TAG"
echo "  Or:  docker run --rm -it $IMAGE_TAG bash"

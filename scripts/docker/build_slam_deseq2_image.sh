#!/usr/bin/env bash
# Build and verify the pinned STAR-suite SLAM DESeq2 analysis container.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

IMAGE_TAG="${IMAGE_TAG:-star-suite/slam-deseq2:bioc3.23-deseq2-1.52.0-tximport-1.40.0}"
BIOCONDUCTOR_IMAGE="${BIOCONDUCTOR_IMAGE:-bioconductor/bioconductor_docker:RELEASE_3_23-r-4.6.0}"
BIOCONDUCTOR_VERSION="${BIOCONDUCTOR_VERSION:-3.23}"
DESEQ2_VERSION="${DESEQ2_VERSION:-1.52.0}"
TXIMPORT_VERSION="${TXIMPORT_VERSION:-1.40.0}"

cd "$REPO_ROOT"

echo "Building SLAM DESeq2 image: $IMAGE_TAG"
echo "  Base image: $BIOCONDUCTOR_IMAGE"
echo "  Bioconductor: $BIOCONDUCTOR_VERSION"
echo "  DESeq2: $DESEQ2_VERSION"
echo "  tximport: $TXIMPORT_VERSION"

docker build \
  -f docker/slam-deseq2/Dockerfile \
  -t "$IMAGE_TAG" \
  --build-arg BIOCONDUCTOR_IMAGE="$BIOCONDUCTOR_IMAGE" \
  --build-arg BIOCONDUCTOR_VERSION="$BIOCONDUCTOR_VERSION" \
  --build-arg DESEQ2_VERSION="$DESEQ2_VERSION" \
  --build-arg TXIMPORT_VERSION="$TXIMPORT_VERSION" \
  .

echo
echo "Verifying pinned R/Bioconductor/DESeq2 versions"
docker run --rm "$IMAGE_TAG"

echo
echo "Image built and verified: $IMAGE_TAG"
echo "Record image digest with:"
echo "  docker image inspect '$IMAGE_TAG' --format '{{json .RepoDigests}}'"

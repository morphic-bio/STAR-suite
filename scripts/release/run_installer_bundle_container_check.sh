#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DOCKER_IMAGE=""
BASE_IMAGE="ubuntu:24.04"
BUNDLE=""
EXPECTED_LABEL=""
EXPECTED_VERSION="2.7.11b"
MANIFEST_OUT=""
BUILD_RUNTIME_IMAGE=1

usage() {
  cat <<USAGE
Usage: $0 --bundle <path> --expected-label <label> [--manifest-out <path>] [--expected-version <version>]
          [--docker-image <tag>] [--base-image <image>] [--skip-image-build]
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bundle)
      BUNDLE="$2"
      shift 2
      ;;
    --expected-label)
      EXPECTED_LABEL="$2"
      shift 2
      ;;
    --manifest-out)
      MANIFEST_OUT="$2"
      shift 2
      ;;
    --expected-version)
      EXPECTED_VERSION="$2"
      shift 2
      ;;
    --docker-image)
      DOCKER_IMAGE="$2"
      shift 2
      ;;
    --base-image)
      BASE_IMAGE="$2"
      shift 2
      ;;
    --skip-image-build)
      BUILD_RUNTIME_IMAGE=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$BUNDLE" || -z "$EXPECTED_LABEL" ]]; then
  echo "ERROR: --bundle and --expected-label are required" >&2
  exit 1
fi

if [[ -z "$DOCKER_IMAGE" ]]; then
  sanitized="$(printf '%s' "$BASE_IMAGE" | tr '/:@' '---')"
  DOCKER_IMAGE="star-suite-runtime-check:${sanitized}"
fi

if [[ "$BUILD_RUNTIME_IMAGE" -eq 1 ]]; then
  docker build \
    --build-arg BASE_IMAGE="$BASE_IMAGE" \
    -f "${REPO_ROOT}/scripts/release/docker/Dockerfile.runtime-check" \
    -t "$DOCKER_IMAGE" \
    "$REPO_ROOT"
fi

artifacts_dir="$(cd "$(dirname "$BUNDLE")" && pwd)"
output_dir="/tmp"
manifest_arg=()
if [[ -n "$MANIFEST_OUT" ]]; then
  mkdir -p "$(dirname "$MANIFEST_OUT")"
  output_dir="$(cd "$(dirname "$MANIFEST_OUT")" && pwd)"
  manifest_arg=(--manifest-out "/out/$(basename "$MANIFEST_OUT")")
fi

docker run --rm \
  -v "${artifacts_dir}:/artifacts:ro" \
  -v "${output_dir}:/out" \
  "$DOCKER_IMAGE" \
  /usr/local/bin/container_check_installer_bundle.sh \
    --bundle "/artifacts/$(basename "$BUNDLE")" \
    --expected-label "$EXPECTED_LABEL" \
    --expected-version "$EXPECTED_VERSION" \
    "${manifest_arg[@]}"

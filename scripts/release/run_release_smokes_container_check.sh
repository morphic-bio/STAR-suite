#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DOCKER_IMAGE=""
BASE_IMAGE="ubuntu:24.04"
MODE=""
TARBALL=""
BUNDLE=""
EXPECTED_LABEL=""
EXPECTED_VERSION="1.3.0"
BUILD_RUNTIME_IMAGE=1

usage() {
  cat <<USAGE
Usage:
  $0 --mode tarball --tarball <path> [--docker-image <tag>] [--base-image <image>] [--expected-version <version>] [--skip-image-build]
  $0 --mode bundle --bundle <path> --expected-label <label> [--docker-image <tag>] [--base-image <image>] [--expected-version <version>] [--skip-image-build]
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --tarball) TARBALL="$2"; shift 2 ;;
    --bundle) BUNDLE="$2"; shift 2 ;;
    --expected-label) EXPECTED_LABEL="$2"; shift 2 ;;
    --expected-version) EXPECTED_VERSION="$2"; shift 2 ;;
    --docker-image) DOCKER_IMAGE="$2"; shift 2 ;;
    --base-image) BASE_IMAGE="$2"; shift 2 ;;
    --skip-image-build) BUILD_RUNTIME_IMAGE=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "${MODE}" ]] || { echo "ERROR: --mode is required" >&2; exit 1; }

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

artifacts_dir=""
mode_args=(--mode "${MODE}" --expected-version "${EXPECTED_VERSION}" --repo-root /repo)
case "${MODE}" in
  tarball)
    [[ -f "${TARBALL}" ]] || { echo "ERROR: tarball not found: ${TARBALL}" >&2; exit 1; }
    artifacts_dir="$(cd "$(dirname "${TARBALL}")" && pwd)"
    mode_args+=(--tarball "/artifacts/$(basename "${TARBALL}")")
    ;;
  bundle)
    [[ -f "${BUNDLE}" ]] || { echo "ERROR: bundle not found: ${BUNDLE}" >&2; exit 1; }
    [[ -n "${EXPECTED_LABEL}" ]] || { echo "ERROR: --expected-label required for bundle mode" >&2; exit 1; }
    artifacts_dir="$(cd "$(dirname "${BUNDLE}")" && pwd)"
    mode_args+=(--bundle "/artifacts/$(basename "${BUNDLE}")" --expected-label "${EXPECTED_LABEL}")
    ;;
  *)
    echo "ERROR: --mode must be tarball or bundle" >&2
    exit 1
    ;;
esac

docker run --rm \
  -v "${artifacts_dir}:/artifacts:ro" \
  -v "${REPO_ROOT}:/repo:ro" \
  "$DOCKER_IMAGE" \
  /usr/local/bin/container_check_release_smokes.sh \
    "${mode_args[@]}"

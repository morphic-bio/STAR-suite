#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SNAPSHOT_DIR=""
STAGE_DIR=""

VERSION=""
OUT_DIR="${REPO_ROOT}/dist/release"
MAKE_JOBS="${MAKE_JOBS:-8}"
DOCKER_IMAGE=""
NATIVE_BUILD=0

cleanup() {
  if [[ -n "${STAGE_DIR}" && -d "${STAGE_DIR}" ]]; then
    rm -rf "${STAGE_DIR}"
  fi
  if [[ -n "${SNAPSHOT_DIR}" && -d "${SNAPSHOT_DIR}" ]]; then
    rm -rf "${SNAPSHOT_DIR}"
  fi
}

trap cleanup EXIT

container_usage_note() {
  cat <<EOF

Optional:
  --docker-image <image>  build from a clean committed HEAD snapshot inside a container
                          (useful for controlling the glibc baseline of the release tarball)
  --native-build          internal flag used by container mode recursion
EOF
}

usage() {
  cat <<EOF
Usage: $0 --version <version> [--out-dir <dir>] [--make-jobs <n>] [--docker-image <image>]

Builds STAR static binary and packages it as:
  STAR-static-<version>-linux-<amd64|arm64>.tar.gz
EOF
  container_usage_note
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --version)
      VERSION="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --make-jobs)
      MAKE_JOBS="$2"
      shift 2
      ;;
    --docker-image)
      DOCKER_IMAGE="$2"
      shift 2
      ;;
    --native-build)
      NATIVE_BUILD=1
      shift 1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${VERSION}" ]]; then
  echo "ERROR: --version is required" >&2
  usage
  exit 1
fi

resolve_commit_sha() {
  if git -C "${REPO_ROOT}" rev-parse --short HEAD >/dev/null 2>&1; then
    git -C "${REPO_ROOT}" rev-parse --short HEAD
  elif [[ -n "${STAR_SUITE_COMMIT_SHA:-}" ]]; then
    printf '%s\n' "${STAR_SUITE_COMMIT_SHA}"
  else
    printf 'unknown\n'
  fi
}

run_container_build() {
  if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: docker is required for --docker-image builds" >&2
    exit 1
  fi

  if ! git -C "${REPO_ROOT}" rev-parse --verify HEAD >/dev/null 2>&1; then
    echo "ERROR: --docker-image mode requires a git checkout with a committed HEAD" >&2
    exit 1
  fi

  local out_dir_abs commit_sha image_name
  SNAPSHOT_DIR="$(mktemp -d)"

  mkdir -p "${OUT_DIR}"
  out_dir_abs="$(cd "${OUT_DIR}" && pwd)"

  commit_sha="$(resolve_commit_sha)"
  image_name="${DOCKER_IMAGE}"

  git -C "${REPO_ROOT}" archive --format=tar HEAD | tar -xf - -C "${SNAPSHOT_DIR}"

  docker run --rm \
    -v "${SNAPSHOT_DIR}:/src-ro:ro" \
    -v "${out_dir_abs}:/out" \
    "${image_name}" \
    bash -lc "
      set -euo pipefail
      export DEBIAN_FRONTEND=noninteractive
      apt-get update
      apt-get install -y --no-install-recommends \
        build-essential \
        xxd \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libglib2.0-dev \
        ca-certificates \
        git
      mkdir -p /src
      cp -a /src-ro/. /src/
      cd /src
      export STAR_SUITE_COMMIT_SHA='${commit_sha}'
      export STAR_SUITE_BUILD_IMAGE='${image_name}'
      scripts/release/build_static_tarball.sh \
        --native-build \
        --version '${VERSION}' \
        --out-dir /out \
        --make-jobs '${MAKE_JOBS}'
    "
}

if [[ -n "${DOCKER_IMAGE}" && "${NATIVE_BUILD}" -eq 0 ]]; then
  mkdir -p "${OUT_DIR}"
  run_container_build
  exit 0
fi

arch_raw="$(uname -m)"
case "${arch_raw}" in
  x86_64) arch="amd64" ;;
  aarch64|arm64) arch="arm64" ;;
  *)
    echo "ERROR: unsupported architecture '${arch_raw}'" >&2
    exit 1
    ;;
esac

cd "${REPO_ROOT}"
mkdir -p "${OUT_DIR}"
OUT_DIR="$(cd "${OUT_DIR}" && pwd)"

echo "Building static STAR (jobs=${MAKE_JOBS})..."
make -j"${MAKE_JOBS}" core-static

if [[ ! -x core/legacy/source/STAR ]]; then
  echo "ERROR: expected binary missing: core/legacy/source/STAR" >&2
  exit 1
fi

STAGE_DIR="$(mktemp -d)"

mkdir -p "${STAGE_DIR}/bin"
cp core/legacy/source/STAR "${STAGE_DIR}/bin/STAR"

cat > "${STAGE_DIR}/README.txt" <<EOF
STAR-suite static release artifact
Version: ${VERSION}
Architecture: ${arch}
Commit: $(resolve_commit_sha)
Built at: $(date -u +"%Y-%m-%dT%H:%M:%SZ")
Build environment: ${STAR_SUITE_BUILD_IMAGE:-native-host}
EOF

tarball="${OUT_DIR}/STAR-static-${VERSION}-linux-${arch}.tar.gz"
tar -C "${STAGE_DIR}" -czf "${tarball}" .

echo "Created ${tarball}"

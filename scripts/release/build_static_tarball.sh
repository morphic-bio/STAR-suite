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
ASSET_PREFIX="STAR-suite"
COMPAT_LABEL=""
GLIBC_BASELINE=""

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
  cat <<USAGE

Optional:
  --docker-image <image>   build from a clean committed HEAD snapshot inside a container
                           (useful for controlling the glibc baseline of the release tarball)
  --compat-label <label>   compatibility label appended to the asset name (for example glibc234)
  --glibc-baseline <ver>   documented minimum glibc level for this asset (for example 2.34)
  --asset-prefix <prefix>  asset filename prefix (default: STAR-suite)
  --native-build           internal flag used by container mode recursion
USAGE
}

usage() {
  cat <<USAGE
Usage: $0 --version <version> [--out-dir <dir>] [--make-jobs <n>] [--docker-image <image>]

Builds a STAR-suite binary tarball release artifact.
USAGE
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
    --compat-label)
      COMPAT_LABEL="$2"
      shift 2
      ;;
    --glibc-baseline)
      GLIBC_BASELINE="$2"
      shift 2
      ;;
    --asset-prefix)
      ASSET_PREFIX="$2"
      shift 2
      ;;
    --native-build)
      NATIVE_BUILD=1
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

if [[ -z "${VERSION}" ]]; then
  echo "ERROR: --version is required" >&2
  usage >&2
  exit 1
fi

resolve_commit_sha() {
  if git -C "${REPO_ROOT}" rev-parse HEAD >/dev/null 2>&1; then
    git -C "${REPO_ROOT}" rev-parse HEAD
  elif [[ -n "${STAR_SUITE_COMMIT_SHA:-}" ]]; then
    printf '%s\n' "${STAR_SUITE_COMMIT_SHA}"
  else
    return 1
  fi
}

require_commit_sha() {
  local commit_sha
  if ! commit_sha="$(resolve_commit_sha)" \
      || [[ ! "${commit_sha}" =~ ^[0-9a-fA-F]{40}$ ]]; then
    echo "ERROR: release builds require an exact 40-character source commit SHA" >&2
    exit 1
  fi
  printf '%s\n' "${commit_sha,,}"
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

  local out_dir_abs commit_sha image_name compat_arg glibc_arg prefix_arg
  SNAPSHOT_DIR="$(mktemp -d)"

  mkdir -p "${OUT_DIR}"
  out_dir_abs="$(cd "${OUT_DIR}" && pwd)"

  commit_sha="$(require_commit_sha)"
  image_name="${DOCKER_IMAGE}"
  compat_arg=""
  glibc_arg=""
  prefix_arg=""

  if [[ -n "${COMPAT_LABEL}" ]]; then
    compat_arg="--compat-label '${COMPAT_LABEL}'"
  fi
  if [[ -n "${GLIBC_BASELINE}" ]]; then
    glibc_arg="--glibc-baseline '${GLIBC_BASELINE}'"
  fi
  if [[ -n "${ASSET_PREFIX}" ]]; then
    prefix_arg="--asset-prefix '${ASSET_PREFIX}'"
  fi

  git -C "${REPO_ROOT}" archive --format=tar HEAD | tar -xf - -C "${SNAPSHOT_DIR}"

  docker run --rm \
    -v "${SNAPSHOT_DIR}:/src-ro:ro" \
    -v "${out_dir_abs}:/out" \
    "${image_name}" \
    bash -lc "
      set -euo pipefail
      export DEBIAN_FRONTEND=noninteractive
      if [[ \"\${STAR_SUITE_BUILD_CONTAINER_READY:-0}\" != \"1\" ]]; then
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
      fi
      mkdir -p /src
      cp -a /src-ro/. /src/
      cd /src
      export STAR_SUITE_COMMIT_SHA='${commit_sha}'
      export STAR_SUITE_BUILD_IMAGE='${image_name}'
      scripts/release/build_static_tarball.sh \
        --native-build \
        --version '${VERSION}' \
        --out-dir /out \
        --make-jobs '${MAKE_JOBS}' \
        ${compat_arg} \
        ${glibc_arg} \
        ${prefix_arg}
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
commit_sha="$(require_commit_sha)"
export STAR_SUITE_COMMIT_SHA="${commit_sha}"

if [[ ! -x scripts/release/install_binary_tarball.sh ]]; then
  echo "ERROR: missing installer script source: scripts/release/install_binary_tarball.sh" >&2
  exit 1
fi

echo "Building STAR-suite binary tarball (jobs=${MAKE_JOBS})..."
make -j"${MAKE_JOBS}" core-static
make -j"${MAKE_JOBS}" molecule-first-resolver
make -j"${MAKE_JOBS}" vbem-transcriptvb-finalize trim-qc-tools

if [[ ! -x core/legacy/source/STAR ]]; then
  echo "ERROR: expected binary missing: core/legacy/source/STAR" >&2
  exit 1
fi
for tool in molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize; do
  if [[ ! -x "flex/tools/molecule_first_resolver/${tool}" ]]; then
    echo "ERROR: expected binary missing: flex/tools/molecule_first_resolver/${tool}" >&2
    exit 1
  fi
done
for tool in transcriptvb_finalize trim_qc_fastq trim_qc_merge; do
  if [[ ! -x "core/legacy/source/${tool}" ]]; then
    echo "ERROR: expected binary missing: core/legacy/source/${tool}" >&2
    exit 1
  fi
done

STAGE_DIR="$(mktemp -d)"
mkdir -p "${STAGE_DIR}/bin"
cp core/legacy/source/STAR "${STAGE_DIR}/bin/STAR"
for tool in molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize; do
  cp "flex/tools/molecule_first_resolver/${tool}" "${STAGE_DIR}/bin/${tool}"
done
for tool in transcriptvb_finalize trim_qc_fastq trim_qc_merge; do
  cp "core/legacy/source/${tool}" "${STAGE_DIR}/bin/${tool}"
done
cp scripts/release/install_binary_tarball.sh "${STAGE_DIR}/install.sh"
chmod 0755 "${STAGE_DIR}/bin/"* "${STAGE_DIR}/install.sh"

asset_name="${ASSET_PREFIX}-${VERSION}-linux-${arch}"
if [[ -n "${COMPAT_LABEL}" ]]; then
  asset_name+="-${COMPAT_LABEL}"
fi
asset_name+=".tar.gz"

cat > "${STAGE_DIR}/release-metadata.env" <<METADATA
VERSION=${VERSION}
ARCH=${arch}
COMPAT_LABEL=${COMPAT_LABEL}
GLIBC_BASELINE=${GLIBC_BASELINE}
COMMIT_SHA=${commit_sha}
BUILD_ENVIRONMENT=${STAR_SUITE_BUILD_IMAGE:-native-host}
ASSET_NAME=${asset_name}
METADATA

compat_note=""
if [[ -n "${COMPAT_LABEL}" ]]; then
  compat_note+="Compatibility label: ${COMPAT_LABEL}"$'\n'
fi
if [[ -n "${GLIBC_BASELINE}" ]]; then
  compat_note+="Documented glibc baseline: ${GLIBC_BASELINE}+"$'\n'
fi

cat > "${STAGE_DIR}/README.txt" <<README
STAR-suite binary release artifact
Version: ${VERSION}
Architecture: ${arch}
Commit: ${commit_sha}
Built at: $(date -u +"%Y-%m-%dT%H:%M:%SZ")
Build environment: ${STAR_SUITE_BUILD_IMAGE:-native-host}
${compat_note}This tarball includes:
  - bin/STAR
  - bin/molecule_first_resolver
  - bin/molecule_first_bam_ledger
  - bin/molecule_first_materialize
  - bin/transcriptvb_finalize
  - bin/trim_qc_fastq
  - bin/trim_qc_merge
  - install.sh for optional local installation
  - release-metadata.env for compatibility metadata

Note:
  Linux decides whether a binary may run on a given system.
  If an operating system rejects a binary built for a newer runtime environment,
  use a lower-compatibility tarball or the installer bundle.
README

tarball="${OUT_DIR}/${asset_name}"
tar -C "${STAGE_DIR}" -czf "${tarball}" .

echo "Created ${tarball}"

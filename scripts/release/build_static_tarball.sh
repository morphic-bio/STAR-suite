#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

VERSION=""
OUT_DIR="${REPO_ROOT}/dist/release"
MAKE_JOBS="${MAKE_JOBS:-8}"

usage() {
  cat <<EOF
Usage: $0 --version <version> [--out-dir <dir>] [--make-jobs <n>]

Builds STAR static binary and packages it as:
  STAR-static-<version>-linux-<amd64|arm64>.tar.gz
EOF
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

arch_raw="$(uname -m)"
case "${arch_raw}" in
  x86_64) arch="amd64" ;;
  aarch64|arm64) arch="arm64" ;;
  *)
    echo "ERROR: unsupported architecture '${arch_raw}'" >&2
    exit 1
    ;;
esac

mkdir -p "${OUT_DIR}"
cd "${REPO_ROOT}"

echo "Building static STAR (jobs=${MAKE_JOBS})..."
make -j"${MAKE_JOBS}" core-static

if [[ ! -x core/legacy/source/STAR ]]; then
  echo "ERROR: expected binary missing: core/legacy/source/STAR" >&2
  exit 1
fi

stage_dir="$(mktemp -d)"
trap 'rm -rf "${stage_dir}"' EXIT

mkdir -p "${stage_dir}/bin"
cp core/legacy/source/STAR "${stage_dir}/bin/STAR"

cat > "${stage_dir}/README.txt" <<EOF
STAR-suite static release artifact
Version: ${VERSION}
Architecture: ${arch}
Commit: $(git rev-parse --short HEAD)
Built at: $(date -u +"%Y-%m-%dT%H:%M:%SZ")
EOF

tarball="${OUT_DIR}/STAR-static-${VERSION}-linux-${arch}.tar.gz"
tar -C "${stage_dir}" -czf "${tarball}" .

echo "Created ${tarball}"

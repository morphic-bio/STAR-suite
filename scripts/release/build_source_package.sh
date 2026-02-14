#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/dist/release/source"
ALLOW_MISSING_DEPS=0

usage() {
  cat <<EOF
Usage: $0 [--out-dir <dir>] [--allow-missing-build-deps]

Build a Debian source package from a clean HEAD snapshot and copy source artifacts to out-dir.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --allow-missing-build-deps)
      ALLOW_MISSING_DEPS=1
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

cd "${REPO_ROOT}"

if [[ ! -f debian/control ]]; then
  echo "ERROR: debian/control not found. Debian packaging skeleton is required." >&2
  exit 1
fi

if ! command -v dh >/dev/null 2>&1; then
  cat >&2 <<EOF
ERROR: debhelper (dh) is required but not installed.
Install with: sudo apt-get install -y xxd debhelper devscripts dpkg-dev fakeroot
EOF
  exit 1
fi

build_args=(-S -sa -us -uc)
if [[ "${ALLOW_MISSING_DEPS}" -eq 1 ]]; then
  build_args+=(-d)
fi

source_name="$(dpkg-parsechangelog -S Source)"
version="$(dpkg-parsechangelog -S Version)"
upstream_version="$(dpkg-parsechangelog -S Version | sed 's/-[^-]*$//')"
base="${source_name}_${version}"

tmp_root="$(mktemp -d)"
trap 'rm -rf "${tmp_root}"' EXIT
work_src="${tmp_root}/${source_name}-${version}"
mkdir -p "${work_src}"

echo "Preparing clean source snapshot from HEAD..."
git archive --format=tar HEAD | tar -xf - -C "${work_src}"

# Source format "3.0 (quilt)" requires an orig tarball in the parent directory.
orig_tar="${tmp_root}/${source_name}_${upstream_version}.orig.tar.gz"
git archive --format=tar --prefix="${source_name}-${upstream_version}/" HEAD | gzip -n > "${orig_tar}"

echo "Building Debian source package..."
cd "${work_src}"
if ! dpkg-buildpackage "${build_args[@]}"; then
  if [[ "${ALLOW_MISSING_DEPS}" -eq 0 ]]; then
    cat >&2 <<EOF
ERROR: dpkg-buildpackage failed.
Hint: install build deps or retry with --allow-missing-build-deps for local smoke runs.
EOF
  fi
  exit 1
fi

mkdir -p "${OUT_DIR}"

copy_if_exists() {
  local f="$1"
  if [[ -f "${f}" ]]; then
    cp -f "${f}" "${OUT_DIR}/"
  fi
}

copy_if_exists "${tmp_root}/${base}.dsc"
copy_if_exists "${tmp_root}/${base}.debian.tar.xz"
copy_if_exists "${tmp_root}/${base}.debian.tar.gz"
copy_if_exists "${tmp_root}/${base}_source.changes"
copy_if_exists "${tmp_root}/${base}_source.buildinfo"
copy_if_exists "${tmp_root}/${source_name}_${upstream_version}.orig.tar.gz"
copy_if_exists "${tmp_root}/${source_name}_${upstream_version}.orig.tar.xz"

echo "Source artifacts copied to ${OUT_DIR}"

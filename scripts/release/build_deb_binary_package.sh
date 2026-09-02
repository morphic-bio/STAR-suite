#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/dist/release/deb"
ALLOW_MISSING_DEPS=0

usage() {
  cat <<EOF
Usage: $0 [--out-dir <dir>] [--allow-missing-build-deps]

Build Debian binary package(s) from a clean HEAD snapshot and copy artifacts to out-dir.
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

OUT_DIR="$(cd "${REPO_ROOT}" && mkdir -p "${OUT_DIR}" && cd "${OUT_DIR}" && pwd)"

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

build_args=(-b -us -uc)
if [[ "${ALLOW_MISSING_DEPS}" -eq 1 ]]; then
  build_args+=(-d)
fi

source_name="$(dpkg-parsechangelog -S Source)"
version="$(dpkg-parsechangelog -S Version)"
arch="$(dpkg --print-architecture)"
base="${source_name}_${version}_${arch}"

tmp_root="$(mktemp -d)"
trap 'rm -rf "${tmp_root}"' EXIT
work_src="${tmp_root}/${source_name}-${version}"
mkdir -p "${work_src}"

echo "Preparing clean source snapshot from HEAD..."
git archive --format=tar HEAD | tar -xf - -C "${work_src}"
commit_sha="$(git rev-parse HEAD)"
printf '%s\n' "${commit_sha}" > "${work_src}/debian/source-revision"

echo "Building Debian binary package..."
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
for ext in deb changes buildinfo; do
  for f in "${tmp_root}/${base}.${ext}"; do
    [[ -f "${f}" ]] && cp -f "${f}" "${OUT_DIR}/"
  done
done

echo "Binary package artifacts copied to ${OUT_DIR}"

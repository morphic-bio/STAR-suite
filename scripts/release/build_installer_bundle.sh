#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
STAGE_DIR=""
EXTRACT_DIR=""

REQUESTED_VERSION=""
OUT_DIR="${REPO_ROOT}/dist/release"
ASSET_PREFIX="STAR-suite"
declare -a TARBALLS=()

cleanup() {
  if [[ -n "${STAGE_DIR}" && -d "${STAGE_DIR}" ]]; then
    rm -rf "${STAGE_DIR}"
  fi
  if [[ -n "${EXTRACT_DIR}" && -d "${EXTRACT_DIR}" ]]; then
    rm -rf "${EXTRACT_DIR}"
  fi
}
trap cleanup EXIT

usage() {
  cat <<USAGE
Usage: $0 --version <version> --tarball <path> [--tarball <path> ...] [--out-dir <dir>]

Build a compatibility installer bundle that auto-selects the best bundled STAR-suite binary.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --version)
      REQUESTED_VERSION="$2"
      shift 2
      ;;
    --tarball)
      TARBALLS+=("$2")
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --asset-prefix)
      ASSET_PREFIX="$2"
      shift 2
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

if [[ -z "${REQUESTED_VERSION}" ]]; then
  echo "ERROR: --version is required" >&2
  exit 1
fi

if [[ "${#TARBALLS[@]}" -lt 1 ]]; then
  echo "ERROR: at least one --tarball is required" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"
OUT_DIR="$(cd "${OUT_DIR}" && pwd)"
STAGE_DIR="$(mktemp -d)"
EXTRACT_DIR="$(mktemp -d)"

mkdir -p "${STAGE_DIR}/variants"
cp "${SCRIPT_DIR}/install_compat_bundle.sh" "${STAGE_DIR}/install.sh"
chmod 0755 "${STAGE_DIR}/install.sh"
MANIFEST_PATH="${STAGE_DIR}/compat-manifest.tsv"
: > "${MANIFEST_PATH}"

bundle_version=""
bundle_arch=""
bundle_commit=""

for tarball in "${TARBALLS[@]}"; do
  tarball="$(realpath "${tarball}")"
  if [[ ! -f "${tarball}" ]]; then
    echo "ERROR: tarball not found: ${tarball}" >&2
    exit 1
  fi

  extract_subdir="${EXTRACT_DIR}/$(basename "${tarball}" .tar.gz)"
  mkdir -p "${extract_subdir}"
  tar -xzf "${tarball}" -C "${extract_subdir}"

  metadata_file="${extract_subdir}/release-metadata.env"
  if [[ ! -f "${metadata_file}" ]]; then
    echo "ERROR: tarball missing release metadata: ${tarball}" >&2
    exit 1
  fi

  unset VERSION ARCH COMPAT_LABEL GLIBC_BASELINE BUILD_ENVIRONMENT ASSET_NAME COMMIT_SHA
  # shellcheck disable=SC1090
  source "${metadata_file}"

  if [[ -z "${VERSION:-}" || -z "${ARCH:-}" || -z "${COMPAT_LABEL:-}" || -z "${GLIBC_BASELINE:-}" ]]; then
    echo "ERROR: tarball missing compatibility metadata: ${tarball}" >&2
    exit 1
  fi
  if [[ ! "${COMMIT_SHA:-}" =~ ^[0-9a-fA-F]{40}$ ]]; then
    echo "ERROR: tarball has invalid commit metadata: ${tarball}" >&2
    exit 1
  fi

  if [[ "${VERSION}" != "${REQUESTED_VERSION}" ]]; then
    echo "ERROR: tarball version ${VERSION} does not match requested bundle version ${REQUESTED_VERSION}" >&2
    exit 1
  fi

  if [[ -z "${bundle_version}" ]]; then
    bundle_version="${VERSION}"
    bundle_arch="${ARCH}"
    bundle_commit="${COMMIT_SHA,,}"
  else
    if [[ "${VERSION}" != "${bundle_version}" ]]; then
      echo "ERROR: mixed versions in installer bundle: ${bundle_version} vs ${VERSION}" >&2
      exit 1
    fi
    if [[ "${ARCH}" != "${bundle_arch}" ]]; then
      echo "ERROR: mixed architectures in installer bundle: ${bundle_arch} vs ${ARCH}" >&2
      exit 1
    fi
    if [[ "${COMMIT_SHA,,}" != "${bundle_commit}" ]]; then
      echo "ERROR: mixed commits in installer bundle: ${bundle_commit} vs ${COMMIT_SHA,,}" >&2
      exit 1
    fi
  fi

  variant_dir="${STAGE_DIR}/variants/${COMPAT_LABEL}"
  mkdir -p "${variant_dir}/bin"
  cp "${extract_subdir}/bin/STAR" "${variant_dir}/bin/STAR"
  for tool in molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize transcriptvb_finalize trim_qc_fastq trim_qc_merge; do
    cp "${extract_subdir}/bin/${tool}" "${variant_dir}/bin/${tool}"
  done
  cp "${extract_subdir}/README.txt" "${variant_dir}/README.txt"
  cp "${metadata_file}" "${variant_dir}/release-metadata.env"
  chmod 0755 "${variant_dir}/bin/"*

  printf '%s\t%s\t%s\t%s\n' \
    "${COMPAT_LABEL}" \
    "${GLIBC_BASELINE}" \
    "variants/${COMPAT_LABEL}/bin/STAR" \
    "${BUILD_ENVIRONMENT:-compatibility build}" \
    >> "${MANIFEST_PATH}"
done

if [[ -z "${bundle_arch}" ]]; then
  echo "ERROR: installer bundle assembly produced no variants" >&2
  exit 1
fi

sort -t $'\t' -k2,2V "${MANIFEST_PATH}" -o "${MANIFEST_PATH}"

compat_lines=""
while IFS=$'\t' read -r label baseline _path description; do
  [[ -z "${label}" ]] && continue
  compat_lines+="  - ${label} (glibc ${baseline}+, ${description})"$'\n'
done < "${MANIFEST_PATH}"

cat > "${STAGE_DIR}/README.txt" <<README
STAR-suite compatibility installer bundle
Version: ${bundle_version}
Architecture: ${bundle_arch}
Commit: ${bundle_commit}

This bundle contains multiple STAR-suite binaries built for different Linux compatibility levels.
Use ./install.sh to auto-select the best bundled binary for this machine.

Recommended install order:
  1. Use the .deb package on Ubuntu/Debian when possible.
  2. Otherwise use this installer bundle.
  3. Use direct compatibility tarballs only when you need manual control.

Bundled compatibility levels:
${compat_lines}
README

bundle_tarball="${OUT_DIR}/${ASSET_PREFIX}-${bundle_version}-linux-${bundle_arch}-installer.tar.gz"
tar -C "${STAGE_DIR}" -czf "${bundle_tarball}" .

echo "Created ${bundle_tarball}"

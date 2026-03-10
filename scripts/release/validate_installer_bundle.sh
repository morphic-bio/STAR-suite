#!/usr/bin/env bash

set -euo pipefail

BUNDLE_PATH=""
EXPECTED_LABEL=""
GLIBC_VERSION=""
EXTRACT_DIR=""
PREFIX_DIR=""

usage() {
  cat <<USAGE
Usage: $0 --bundle <path> [--glibc-version <version>] [--expected-label <label>]

Validate that the installer bundle selects a compatible binary and installs it successfully.
USAGE
}

cleanup() {
  if [[ -n "${EXTRACT_DIR}" && -d "${EXTRACT_DIR}" ]]; then
    rm -rf "${EXTRACT_DIR}"
  fi
}
trap cleanup EXIT

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bundle)
      BUNDLE_PATH="$2"
      shift 2
      ;;
    --glibc-version)
      GLIBC_VERSION="$2"
      shift 2
      ;;
    --expected-label)
      EXPECTED_LABEL="$2"
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

if [[ -z "${BUNDLE_PATH}" ]]; then
  echo "ERROR: --bundle is required" >&2
  exit 1
fi

BUNDLE_PATH="$(realpath "${BUNDLE_PATH}")"
if [[ ! -f "${BUNDLE_PATH}" ]]; then
  echo "ERROR: bundle not found: ${BUNDLE_PATH}" >&2
  exit 1
fi

EXTRACT_DIR="$(mktemp -d)"
PREFIX_DIR="${EXTRACT_DIR}/prefix"
mkdir -p "${EXTRACT_DIR}/unpack"
tar -xzf "${BUNDLE_PATH}" -C "${EXTRACT_DIR}/unpack"

selection_cmd=("${EXTRACT_DIR}/unpack/install.sh" --print-selection)
install_cmd=("${EXTRACT_DIR}/unpack/install.sh" --prefix "${PREFIX_DIR}" --force)
if [[ -n "${GLIBC_VERSION}" ]]; then
  selection_cmd+=(--glibc-version "${GLIBC_VERSION}")
  install_cmd+=(--glibc-version "${GLIBC_VERSION}")
fi

selection="$(${selection_cmd[@]})"
IFS=$'\t' read -r selected_label selected_baseline _selected_relpath _selected_description <<< "${selection}"

if [[ -n "${EXPECTED_LABEL}" && "${selected_label}" != "${EXPECTED_LABEL}" ]]; then
  echo "ERROR: expected label ${EXPECTED_LABEL}, got ${selected_label}" >&2
  exit 1
fi

"${install_cmd[@]}"

installed_bin="${PREFIX_DIR}/bin/STAR"
if [[ ! -x "${installed_bin}" ]]; then
  echo "ERROR: installed binary missing: ${installed_bin}" >&2
  exit 1
fi

version_output="$(${installed_bin} --version)"
echo "Selected label: ${selected_label}"
echo "Selected glibc baseline: ${selected_baseline}"
echo "Installed version: ${version_output}"

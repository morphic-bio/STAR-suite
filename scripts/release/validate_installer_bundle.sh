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

version_ge() {
  local lhs="$1"
  local rhs="$2"
  printf '%s\n%s\n' "$rhs" "$lhs" | sort -V -C
}

detect_actual_glibc() {
  if command -v getconf >/dev/null 2>&1; then
    local out
    out="$(getconf GNU_LIBC_VERSION 2>/dev/null || true)"
    if [[ "${out}" =~ ([0-9]+\.[0-9]+) ]]; then
      printf '%s\n' "${BASH_REMATCH[1]}"
      return 0
    fi
  fi

  if command -v ldd >/dev/null 2>&1; then
    local first_line
    first_line="$(ldd --version 2>&1 | head -n1)"
    if [[ "${first_line}" =~ ([0-9]+\.[0-9]+) ]]; then
      printf '%s\n' "${BASH_REMATCH[1]}"
      return 0
    fi
  fi

  printf 'unknown\n'
}

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
installed_resolver="${PREFIX_DIR}/bin/molecule_first_resolver"
if [[ ! -x "${installed_bin}" ]]; then
  echo "ERROR: installed binary missing: ${installed_bin}" >&2
  exit 1
fi
if [[ ! -x "${installed_resolver}" ]]; then
  echo "ERROR: installed molecule-first resolver missing: ${installed_resolver}" >&2
  exit 1
fi
for tool in molecule_first_bam_ledger molecule_first_materialize; do
  if [[ ! -x "${PREFIX_DIR}/bin/${tool}" ]]; then
    echo "ERROR: installed molecule-first companion missing: ${PREFIX_DIR}/bin/${tool}" >&2
    exit 1
  fi
done
for tool in transcriptvb_finalize trim_qc_fastq trim_qc_merge; do
  if [[ ! -x "${PREFIX_DIR}/bin/${tool}" ]]; then
    echo "ERROR: installed release companion missing: ${PREFIX_DIR}/bin/${tool}" >&2
    exit 1
  fi
done

actual_host_glibc="$(detect_actual_glibc)"
echo "Selected label: ${selected_label}"
echo "Selected glibc baseline: ${selected_baseline}"

if [[ "${actual_host_glibc}" == "unknown" ]]; then
  echo "NOTE: unable to determine actual host glibc; skipping installed binary runtime check."
  exit 0
fi

echo "Actual host glibc: ${actual_host_glibc}"

if ! version_ge "${actual_host_glibc}" "${selected_baseline}"; then
  echo "NOTE: skipping installed binary runtime check because actual host glibc ${actual_host_glibc} is below selected baseline ${selected_baseline}."
  echo "NOTE: selection validation passed; runtime validation requires a host or container compatible with the selected variant."
  exit 0
fi

version_output="$(${installed_bin} --version)"
echo "Installed version: ${version_output}"

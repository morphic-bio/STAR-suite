#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_BIN="${SCRIPT_DIR}/bin/STAR"
COMPANION_TOOLS=(
  molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize
  transcriptvb_finalize trim_qc_fastq trim_qc_merge
)
METADATA_FILE="${SCRIPT_DIR}/release-metadata.env"
SRC_SHARE="${SCRIPT_DIR}/share/star-suite"

PREFIX=""
BINDIR=""
TARGET_NAME="STAR"
FORCE=0
PRINT_TARGET=0

usage() {
  cat <<USAGE
Usage: $0 [--prefix <dir>] [--bindir <dir>] [--name <name>] [--force] [--print-target]

Install the bundled STAR-suite binary from this tarball.

Defaults:
  non-root install prefix: \$HOME/.local
  root install prefix:     /usr/local
USAGE
}

version_ge() {
  local lhs="$1"
  local rhs="$2"
  printf '%s\n%s\n' "$rhs" "$lhs" | sort -V -C
}

warn_if_path_missing() {
  local dir="$1"
  case ":${PATH}:" in
    *":${dir}:"*) ;;
    *) echo "NOTE: ${dir} is not on PATH." ;;
  esac
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      PREFIX="$2"
      shift 2
      ;;
    --bindir)
      BINDIR="$2"
      shift 2
      ;;
    --name)
      TARGET_NAME="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --print-target)
      PRINT_TARGET=1
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

if [[ ! -x "${SRC_BIN}" ]]; then
  echo "ERROR: bundled binary not found: ${SRC_BIN}" >&2
  exit 1
fi
for tool in "${COMPANION_TOOLS[@]}"; do
  if [[ ! -x "${SCRIPT_DIR}/bin/${tool}" ]]; then
    echo "ERROR: bundled binary not found: ${SCRIPT_DIR}/bin/${tool}" >&2
    exit 1
  fi
done
if [[ ! -f "${SRC_SHARE}/SNAPSHOTS.json" \
    || ! -f "${SRC_SHARE}/catalogs/official/catalog.yaml" \
    || ! -f "${SRC_SHARE}/evidence/official/schema/record-v1.schema.json" ]]; then
  echo "ERROR: bundled official recipe/provenance snapshots are incomplete" >&2
  exit 1
fi

if [[ -f "${METADATA_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${METADATA_FILE}"
fi

if [[ -z "${PREFIX}" ]]; then
  if [[ "$(id -u)" -eq 0 ]]; then
    PREFIX="/usr/local"
  else
    PREFIX="${HOME}/.local"
  fi
fi

if [[ -z "${BINDIR}" ]]; then
  BINDIR="${PREFIX}/bin"
fi

TARGET_PATH="${BINDIR}/${TARGET_NAME}"

if [[ "${PRINT_TARGET}" -eq 1 ]]; then
  printf '%s\n' "${TARGET_PATH}"
  exit 0
fi

mkdir -p "${BINDIR}"

if [[ -e "${TARGET_PATH}" && "${FORCE}" -ne 1 ]]; then
  if cmp -s "${SRC_BIN}" "${TARGET_PATH}"; then
    echo "STAR binary already installed at ${TARGET_PATH}; verifying companion tools"
  else
    echo "ERROR: target already exists: ${TARGET_PATH} (use --force to overwrite)" >&2
    exit 1
  fi
fi
for tool in "${COMPANION_TOOLS[@]}"; do
  target="${BINDIR}/${tool}"
  if [[ -e "${target}" && "${FORCE}" -ne 1 ]] && ! cmp -s "${SCRIPT_DIR}/bin/${tool}" "${target}"; then
    echo "ERROR: target already exists: ${target} (use --force to overwrite)" >&2
    exit 1
  fi
done

SHAREDIR="${PREFIX}/share/star-suite"
if [[ -d "${SHAREDIR}" && "${FORCE}" -ne 1 ]]; then
  if ! diff -qr "${SRC_SHARE}" "${SHAREDIR}" >/dev/null; then
    echo "ERROR: installed STAR Suite data differs at ${SHAREDIR} (use --force to update)" >&2
    exit 1
  fi
else
  mkdir -p "${SHAREDIR}"
  cp -a "${SRC_SHARE}/." "${SHAREDIR}/"
fi

tmp_target="${TARGET_PATH}.tmp.$$"
cp "${SRC_BIN}" "${tmp_target}"
chmod 0755 "${tmp_target}"
mv -f "${tmp_target}" "${TARGET_PATH}"

for tool in "${COMPANION_TOOLS[@]}"; do
  target="${BINDIR}/${tool}"
  temporary="${target}.tmp.$$"
  cp "${SCRIPT_DIR}/bin/${tool}" "${temporary}"
  chmod 0755 "${temporary}"
  mv -f "${temporary}" "${target}"
done

if [[ -n "${COMPAT_LABEL:-}" ]]; then
  echo "Installed STAR-suite (${COMPAT_LABEL}) to ${TARGET_PATH}"
else
  echo "Installed STAR-suite to ${TARGET_PATH}"
fi
if [[ -n "${GLIBC_BASELINE:-}" ]]; then
  echo "Binary compatibility baseline: glibc ${GLIBC_BASELINE}+"
fi
warn_if_path_missing "${BINDIR}"
echo "Official recipes: ${SHAREDIR}/catalogs/official/catalog.yaml"
echo "Official evidence: ${SHAREDIR}/evidence/official"

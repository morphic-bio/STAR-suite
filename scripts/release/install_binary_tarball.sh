#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_BIN="${SCRIPT_DIR}/bin/STAR"
MOLECULE_FIRST_SRC_BIN="${SCRIPT_DIR}/bin/molecule_first_resolver"
METADATA_FILE="${SCRIPT_DIR}/release-metadata.env"

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
if [[ ! -x "${MOLECULE_FIRST_SRC_BIN}" ]]; then
  echo "ERROR: bundled binary not found: ${MOLECULE_FIRST_SRC_BIN}" >&2
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
MOLECULE_FIRST_TARGET="${BINDIR}/molecule_first_resolver"

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
if [[ -e "${MOLECULE_FIRST_TARGET}" && "${FORCE}" -ne 1 ]] && ! cmp -s "${MOLECULE_FIRST_SRC_BIN}" "${MOLECULE_FIRST_TARGET}"; then
  echo "ERROR: target already exists: ${MOLECULE_FIRST_TARGET} (use --force to overwrite)" >&2
  exit 1
fi

tmp_target="${TARGET_PATH}.tmp.$$"
cp "${SRC_BIN}" "${tmp_target}"
chmod 0755 "${tmp_target}"
mv -f "${tmp_target}" "${TARGET_PATH}"

tmp_molecule_first="${MOLECULE_FIRST_TARGET}.tmp.$$"
cp "${MOLECULE_FIRST_SRC_BIN}" "${tmp_molecule_first}"
chmod 0755 "${tmp_molecule_first}"
mv -f "${tmp_molecule_first}" "${MOLECULE_FIRST_TARGET}"

if [[ -n "${COMPAT_LABEL:-}" ]]; then
  echo "Installed STAR-suite (${COMPAT_LABEL}) to ${TARGET_PATH}"
else
  echo "Installed STAR-suite to ${TARGET_PATH}"
fi
if [[ -n "${GLIBC_BASELINE:-}" ]]; then
  echo "Binary compatibility baseline: glibc ${GLIBC_BASELINE}+"
fi
warn_if_path_missing "${BINDIR}"

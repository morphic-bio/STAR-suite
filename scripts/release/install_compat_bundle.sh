#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST_FILE="${SCRIPT_DIR}/compat-manifest.tsv"
MOLECULE_FIRST_TOOLS=(molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize)

PREFIX=""
BINDIR=""
TARGET_NAME="STAR"
FORCE=0
PRINT_SELECTION=0
LIST_ONLY=0
REQUESTED_LABEL=""
OVERRIDE_GLIBC="${STAR_SUITE_INSTALLER_GLIBC_VERSION:-}"

usage() {
  cat <<USAGE
Usage: $0 [--prefix <dir>] [--bindir <dir>] [--name <name>] [--force]
          [--compat-label <label>] [--glibc-version <version>]
          [--print-selection] [--list]

Auto-select and install the most compatible STAR-suite binary bundled in this tarball.

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

list_variants() {
  while IFS=$'\t' read -r label baseline relpath description; do
    [[ -z "${label}" || "${label}" == \#* ]] && continue
    printf '%s\tglibc %s+\t%s\n' "${label}" "${baseline}" "${description}"
  done < "${MANIFEST_FILE}"
}

detect_glibc_version() {
  if [[ -n "${OVERRIDE_GLIBC}" ]]; then
    printf '%s\n' "${OVERRIDE_GLIBC}"
    return 0
  fi

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

  echo "ERROR: unable to determine glibc version on this system." >&2
  exit 1
}

select_variant() {
  local host_glibc="$1"
  local best_label=""
  local best_baseline=""
  local best_path=""
  local best_description=""
  local found_requested=0

  while IFS=$'\t' read -r label baseline relpath description; do
    [[ -z "${label}" || "${label}" == \#* ]] && continue

    if [[ -n "${REQUESTED_LABEL}" && "${label}" != "${REQUESTED_LABEL}" ]]; then
      continue
    fi

    found_requested=1

    if ! version_ge "${host_glibc}" "${baseline}"; then
      continue
    fi

    if [[ -z "${best_baseline}" ]] || version_ge "${baseline}" "${best_baseline}"; then
      best_label="${label}"
      best_baseline="${baseline}"
      best_path="${relpath}"
      best_description="${description}"
    fi
  done < "${MANIFEST_FILE}"

  if [[ -n "${REQUESTED_LABEL}" && "${found_requested}" -eq 0 ]]; then
    echo "ERROR: requested compatibility label not found: ${REQUESTED_LABEL}" >&2
    exit 1
  fi

  if [[ -z "${best_label}" ]]; then
    echo "ERROR: no compatible STAR-suite binary found for glibc ${host_glibc}." >&2
    echo "Available variants:" >&2
    list_variants >&2
    exit 1
  fi

  printf '%s\t%s\t%s\t%s\n' "${best_label}" "${best_baseline}" "${best_path}" "${best_description}"
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
    --compat-label)
      REQUESTED_LABEL="$2"
      shift 2
      ;;
    --glibc-version)
      OVERRIDE_GLIBC="$2"
      shift 2
      ;;
    --print-selection)
      PRINT_SELECTION=1
      shift
      ;;
    --list)
      LIST_ONLY=1
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

if [[ ! -f "${MANIFEST_FILE}" ]]; then
  echo "ERROR: compatibility manifest not found: ${MANIFEST_FILE}" >&2
  exit 1
fi

if [[ "${LIST_ONLY}" -eq 1 ]]; then
  list_variants
  exit 0
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

HOST_GLIBC="$(detect_glibc_version)"
selection="$(select_variant "${HOST_GLIBC}")"
IFS=$'\t' read -r selected_label selected_baseline selected_relpath selected_description <<< "${selection}"
selected_bin="${SCRIPT_DIR}/${selected_relpath}"

if [[ ! -x "${selected_bin}" ]]; then
  echo "ERROR: selected bundled binary not found: ${selected_bin}" >&2
  exit 1
fi
for tool in "${MOLECULE_FIRST_TOOLS[@]}"; do
  if [[ ! -x "$(dirname "${selected_bin}")/${tool}" ]]; then
    echo "ERROR: selected bundled binary not found: $(dirname "${selected_bin}")/${tool}" >&2
    exit 1
  fi
done

if [[ "${PRINT_SELECTION}" -eq 1 ]]; then
  printf '%s\t%s\t%s\t%s\n' "${selected_label}" "${selected_baseline}" "${selected_relpath}" "${selected_description}"
  exit 0
fi

mkdir -p "${BINDIR}"
TARGET_PATH="${BINDIR}/${TARGET_NAME}"

if [[ -e "${TARGET_PATH}" && "${FORCE}" -ne 1 ]]; then
  if cmp -s "${selected_bin}" "${TARGET_PATH}"; then
    echo "STAR binary already installed at ${TARGET_PATH}; verifying companion tools"
  else
    echo "ERROR: target already exists: ${TARGET_PATH} (use --force to overwrite)" >&2
    exit 1
  fi
fi
for tool in "${MOLECULE_FIRST_TOOLS[@]}"; do
  target="${BINDIR}/${tool}"
  if [[ -e "${target}" && "${FORCE}" -ne 1 ]] && ! cmp -s "$(dirname "${selected_bin}")/${tool}" "${target}"; then
    echo "ERROR: target already exists: ${target} (use --force to overwrite)" >&2
    exit 1
  fi
done

tmp_target="${TARGET_PATH}.tmp.$$"
cp "${selected_bin}" "${tmp_target}"
chmod 0755 "${tmp_target}"
mv -f "${tmp_target}" "${TARGET_PATH}"

for tool in "${MOLECULE_FIRST_TOOLS[@]}"; do
  target="${BINDIR}/${tool}"
  temporary="${target}.tmp.$$"
  cp "$(dirname "${selected_bin}")/${tool}" "${temporary}"
  chmod 0755 "${temporary}"
  mv -f "${temporary}" "${target}"
done

echo "Installed STAR-suite to ${TARGET_PATH}"
echo "Selected compatibility level: ${selected_label} (glibc ${selected_baseline}+)"
echo "Host glibc detected: ${HOST_GLIBC}"
warn_if_path_missing "${BINDIR}"

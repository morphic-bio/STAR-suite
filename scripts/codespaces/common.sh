#!/usr/bin/env bash
set -euo pipefail

CODESPACES_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${CODESPACES_SCRIPT_DIR}/../.." && pwd)"
DEMO_ROOT="${DEMO_ROOT:-${REPO_ROOT}/.codespaces-demo}"
CACHE_ROOT="${CACHE_ROOT:-${DEMO_ROOT}/cache}"
DATA_ROOT="${DATA_ROOT:-${DEMO_ROOT}/data}"
INDEX_ROOT="${INDEX_ROOT:-${DEMO_ROOT}/indices}"
RUN_ROOT="${RUN_ROOT:-${DEMO_ROOT}/runs}"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"

mkdir -p "${DEMO_ROOT}" "${CACHE_ROOT}" "${DATA_ROOT}" "${INDEX_ROOT}" "${RUN_ROOT}"

log() {
  printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

ensure_star_built() {
  if [[ ! -x "${STAR_BIN}" ]]; then
    log "STAR binary not found at ${STAR_BIN}; building core/legacy/source/STAR"
    make -C "${REPO_ROOT}/core/legacy/source" -j"$(nproc)" STAR
  fi
}

write_command_script() {
  local path="$1"
  shift
  {
    echo '#!/usr/bin/env bash'
    echo 'set -euo pipefail'
    printf '%q ' "$@"
    printf '\n'
  } > "${path}"
  chmod +x "${path}"
}

resolve_tiny_gtf() {
  local ref_root="$1"
  if [[ -f "${ref_root}/genes/genes.gtf" ]]; then
    echo "${ref_root}/genes/genes.gtf"
  elif [[ -f "${ref_root}/genes/genes.gtf.gz" ]]; then
    echo "${ref_root}/genes/genes.gtf.gz"
  else
    die "Could not find genes.gtf or genes.gtf.gz under ${ref_root}"
  fi
}

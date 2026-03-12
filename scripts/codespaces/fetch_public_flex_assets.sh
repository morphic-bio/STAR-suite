#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
# shellcheck source=public_demo_sources.sh
source "${SCRIPT_DIR}/public_demo_sources.sh"

OUTDIR="${OUTDIR:-${DATA_ROOT}/public_flex_k562}"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [--outdir DIR]

Fetch the small official metadata assets for the public 10x Flex dataset that
can replace the current synthetic probe overlay in Codespaces demos.

Outputs:
  <outdir>/config.csv
  <outdir>/probe_set.csv
  <outdir>/MANIFEST.txt
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

need_cmd curl
mkdir -p "${OUTDIR}"

curl -L --fail --silent --show-error \
  "${CODESPACES_PUBLIC_FLEX_CONFIG_URL}" \
  -o "${OUTDIR}/config.csv"
curl -L --fail --silent --show-error \
  "${CODESPACES_PUBLIC_FLEX_PROBE_SET_URL}" \
  -o "${OUTDIR}/probe_set.csv"

{
  printf 'dataset_url=%s\n' "${CODESPACES_PUBLIC_FLEX_DATASET_URL}"
  printf 'fastqs_tar_url=%s\n' "${CODESPACES_PUBLIC_FLEX_FASTQS_TAR_URL}"
  printf 'fastqs_tar_bytes=%s\n' "${CODESPACES_PUBLIC_FLEX_FASTQS_TAR_BYTES}"
  printf 'config_url=%s\n' "${CODESPACES_PUBLIC_FLEX_CONFIG_URL}"
  printf 'probe_set_url=%s\n' "${CODESPACES_PUBLIC_FLEX_PROBE_SET_URL}"
  printf 'fetched_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "${OUTDIR}/MANIFEST.txt"

log "Prepared ${OUTDIR}"

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
# shellcheck source=public_demo_sources.sh
source "${SCRIPT_DIR}/public_demo_sources.sh"

OUTDIR="${OUTDIR:-${DATA_ROOT}/public_perturb_k562}"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [--outdir DIR]

Fetch the small official metadata assets for the public 10x perturb dataset that
can replace the current synthetic guide companion in Codespaces demos.

Outputs:
  <outdir>/config.csv
  <outdir>/feature_reference.csv
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
  "${CODESPACES_PUBLIC_PERTURB_CONFIG_URL}" \
  -o "${OUTDIR}/config.csv"
curl -L --fail --silent --show-error \
  "${CODESPACES_PUBLIC_PERTURB_FEATURE_REF_URL}" \
  -o "${OUTDIR}/feature_reference.csv"

{
  printf 'dataset_url=%s\n' "${CODESPACES_PUBLIC_PERTURB_DATASET_URL}"
  printf 'fastqs_tar_url=%s\n' "${CODESPACES_PUBLIC_PERTURB_FASTQS_TAR_URL}"
  printf 'fastqs_tar_bytes=%s\n' "${CODESPACES_PUBLIC_PERTURB_FASTQS_TAR_BYTES}"
  printf 'config_url=%s\n' "${CODESPACES_PUBLIC_PERTURB_CONFIG_URL}"
  printf 'feature_ref_url=%s\n' "${CODESPACES_PUBLIC_PERTURB_FEATURE_REF_URL}"
  printf 'fetched_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "${OUTDIR}/MANIFEST.txt"

log "Prepared ${OUTDIR}"

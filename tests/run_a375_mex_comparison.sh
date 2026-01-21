#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1091
  source "${EXTERNAL_ENV}"
fi

CR_MEX_DIR="${CR_MEX_DIR:-/storage/A375/outputs/unpacked/sample_filtered_feature_bc_matrix}"
STAR_GENE_DIR="${STAR_GENE_DIR:-/storage/A375/star_multi_smoke_cpp/tier_1000000/outs/raw_feature_bc_matrix}"
STAR_GENEFULL_DIR="${STAR_GENEFULL_DIR:-/storage/A375/star_multi_smoke_cpp/tier_1000000/Solo.out/GeneFull/raw}"

usage() {
  cat <<'EOF'
Usage: tests/run_a375_mex_comparison.sh [--gene-only] [--genefull-only]

Environment overrides:
  CR_MEX_DIR         Cell Ranger MEX directory (default: /storage/A375/outputs/unpacked/sample_filtered_feature_bc_matrix)
  STAR_GENE_DIR      STAR Gene MEX directory
  STAR_GENEFULL_DIR  STAR GeneFull MEX directory
EOF
}

run_gene=1
run_genefull=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gene-only)
      run_genefull=0
      shift
      ;;
    --genefull-only)
      run_gene=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ ! -d "${CR_MEX_DIR}" ]]; then
  echo "Missing CR MEX dir: ${CR_MEX_DIR}" >&2
  exit 1
fi

if (( run_gene )); then
  if [[ ! -d "${STAR_GENE_DIR}" ]]; then
    echo "Missing STAR Gene MEX dir: ${STAR_GENE_DIR}" >&2
    exit 1
  fi
  echo "== Comparing Gene =="
  python3 "${REPO_ROOT}/tests/compare_a375_star_mex.py" "${CR_MEX_DIR}" "${STAR_GENE_DIR}"
fi

if (( run_genefull )); then
  if [[ ! -d "${STAR_GENEFULL_DIR}" ]]; then
    echo "Missing STAR GeneFull MEX dir: ${STAR_GENEFULL_DIR}" >&2
    exit 1
  fi
  echo "== Comparing GeneFull =="
  python3 "${REPO_ROOT}/tests/compare_a375_star_mex.py" "${CR_MEX_DIR}" "${STAR_GENEFULL_DIR}"
fi

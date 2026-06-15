#!/usr/bin/env bash
set -euo pipefail
# Tiered CITE-seq / protein Feature Barcode validation suite.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RUN_SYNTHETIC=1
RUN_DOGMA=1
RUN_PBMC1K=0
PBMC1K_DOWNLOAD=0

usage() {
  sed -n '2,/^$/s/^# \?//p' "$0"
  cat <<USAGE

Usage:
  tests/run_citeseq_validation_suite.sh [options]

Default runs the no-download quick tier:
  - process_features synthetic ADT MEX
  - process_features synthetic hash demux
  - pf-multi synthetic ADT/protein arm
  - DOGMA-HIV HTO smoke, skipped by its script if local data are unavailable

Options:
  --quick          Run only the default quick tier
  --pbmc1k        Also run the public 10x PBMC 1k TotalSeq-B benchmark
  --download      Allow the PBMC 1k benchmark to download/stage public assets
  --synthetic     Run only synthetic tests
  --dogma         Run only DOGMA-HIV smoke
  --all           Run quick tier plus PBMC 1k benchmark
  -h, --help      Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --quick)     RUN_SYNTHETIC=1; RUN_DOGMA=1; RUN_PBMC1K=0; shift ;;
    --pbmc1k)    RUN_PBMC1K=1; shift ;;
    --download)  PBMC1K_DOWNLOAD=1; shift ;;
    --synthetic) RUN_SYNTHETIC=1; RUN_DOGMA=0; RUN_PBMC1K=0; shift ;;
    --dogma)     RUN_SYNTHETIC=0; RUN_DOGMA=1; RUN_PBMC1K=0; shift ;;
    --all)       RUN_SYNTHETIC=1; RUN_DOGMA=1; RUN_PBMC1K=1; shift ;;
    -h|--help)   usage; exit 0 ;;
    *)           echo "FAIL: unknown argument: $1" >&2; exit 1 ;;
  esac
done

run_step() {
  local label="$1"
  shift
  echo ""
  echo "=== ${label} ==="
  "$@"
}

if [[ "${RUN_SYNTHETIC}" -eq 1 ]]; then
  run_step "process_features ADT MEX synthetic" \
    bash "${REPO_ROOT}/core/features/process_features/tests/test_adt_mex.sh"
  run_step "process_features HTO/CMO hash demux synthetic" \
    bash "${REPO_ROOT}/core/features/process_features/tests/test_hash_demux_mex.sh"
  run_step "pf-multi ADT/protein feature arm synthetic" \
    bash "${REPO_ROOT}/tests/multi_feature/test_adt_protein_multifeature_arm.sh"
  run_step "pf-multi HTO/CMO hash demux arm synthetic" \
    bash "${REPO_ROOT}/tests/multi_feature/test_hto_cmo_hash_demux_arm.sh"
fi

if [[ "${RUN_DOGMA}" -eq 1 ]]; then
  run_step "DOGMA-HIV native HTO demux smoke" \
    bash "${REPO_ROOT}/tests/multi_feature/test_hiv_dogma_hto_demux_smoke.sh"
fi

if [[ "${RUN_PBMC1K}" -eq 1 ]]; then
  args=()
  [[ "${PBMC1K_DOWNLOAD}" -eq 1 ]] && args+=(--download)
  run_step "10x PBMC 1k TotalSeq-B process_features benchmark" \
    bash "${REPO_ROOT}/tests/run_citeseq_pbmc1k_pf_benchmark.sh" "${args[@]}"
fi

echo ""
echo "PASS: CITE-seq validation suite completed"

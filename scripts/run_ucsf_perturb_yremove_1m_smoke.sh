#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
OUT_ROOT="${OUT_ROOT:-/mnt/pikachu/ucsf-perturb-yremove-smoke_1m_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-32}"
SAMPLES="${SAMPLES:-iPSC2_1_AALG2,iPSC2_2_AALG2}"

extra_args=("$@")
has_out_root=0
for ((i=0; i<${#extra_args[@]}; i++)); do
  if [[ "${extra_args[$i]}" == "--out-root" ]]; then
    has_out_root=1
    break
  fi
done

cmd=(
  "${SCRIPT_DIR}/run_ucsf_perturb_yremove_batch.sh"
  --star-bin "${STAR_BIN}"
  --samples "${SAMPLES}"
  --downsample-reads 1000000
  --threads "${THREADS}"
)

if [[ "${has_out_root}" == "0" ]]; then
  cmd+=(--out-root "${OUT_ROOT}")
fi

exec "${cmd[@]}" "${extra_args[@]}"

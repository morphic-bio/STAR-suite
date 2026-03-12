#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_flex_cr_config.sh"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
CR_CONFIG="${CR_CONFIG:-/storage/SC2300771_filtered_100K/cellranger/outs/config.csv}"
OUT_BASE="${OUT_BASE:-/tmp/flex_cr_config_smoke_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-4}"

if [[ ! -f "${CR_CONFIG}" ]]; then
  echo "SKIP: missing Flex Cell Ranger config: ${CR_CONFIG}"
  exit 0
fi

"${RUN_SCRIPT}" \
  --cr-config "${CR_CONFIG}" \
  --out-base "${OUT_BASE}" \
  --run-id smoke \
  --threads "${THREADS}"

run_root="${OUT_BASE}/smoke"
[[ -f "${run_root}/RUN_MANIFEST.txt" ]]
[[ -f "${run_root}/sample_whitelist.from_cr.tsv" ]]
[[ -f "${run_root}/sample_probes.from_cr.tsv" ]]
[[ -f "${run_root}/probe_list.from_cr.txt" ]]
[[ -f "${run_root}/Solo.out/Gene/raw/matrix.mtx" ]]
[[ -f "${run_root}/per_sample/flexfilter_summary.tsv" ]]
[[ -f "${run_root}/per_sample/BC004/Gene/filtered/matrix.mtx" ]]
[[ -f "${run_root}/per_sample/BC008/Gene/filtered/matrix.mtx" ]]
grep -F "cr_config=${CR_CONFIG}" "${run_root}/RUN_MANIFEST.txt" >/dev/null

echo "PASS: ${run_root}"

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_flex_cr_config.sh"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
BASE_CONFIG="${BASE_CONFIG:-/storage/SC2300771_filtered_100K/cellranger/outs/config.csv}"
STAR_PROBE_CATALOG="${STAR_PROBE_CATALOG:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
OUT_BASE="${OUT_BASE:-/tmp/flex_cr_config_star_section_smoke_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-4}"

if [[ ! -f "${BASE_CONFIG}" ]]; then
  echo "SKIP: missing Flex Cell Ranger config: ${BASE_CONFIG}"
  exit 0
fi
if [[ ! -f "${STAR_PROBE_CATALOG}" ]]; then
  echo "SKIP: missing Flex sample probe catalog: ${STAR_PROBE_CATALOG}"
  exit 0
fi

tmpdir="$(mktemp -d /tmp/flex_cr_config_star_section_cfg_XXXXXX)"
trap 'rm -rf "${tmpdir}"' EXIT
config="${tmpdir}/config.csv"
cp "${BASE_CONFIG}" "${config}"
cat >> "${config}" <<EOF

[star]
sample-probe-catalog,${STAR_PROBE_CATALOG}
sample-probe-offset,68
EOF

"${RUN_SCRIPT}" \
  --cr-config "${config}" \
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
grep -F "cr_config=${config}" "${run_root}/RUN_MANIFEST.txt" >/dev/null
grep -F "sample_probe_catalog=${STAR_PROBE_CATALOG}" "${run_root}/RUN_MANIFEST.txt" >/dev/null
grep -F "sample_probe_offset=68" "${run_root}/RUN_MANIFEST.txt" >/dev/null

echo "PASS: ${run_root}"

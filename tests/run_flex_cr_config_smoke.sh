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

# HASH_CACHE: optional half-probe cache path, reused if present and built there if not.
CACHE_ARGS=()
[[ -n "${HASH_CACHE:-}" ]] && CACHE_ARGS=(--hash-cache "${HASH_CACHE}")

"${RUN_SCRIPT}" \
  --cr-config "${CR_CONFIG}" \
  "${CACHE_ARGS[@]}" \
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
# Every configured sample gets an output directory; on a small downsample only some
# samples call cells, so require at least one filtered matrix rather than naming samples.
sample_ids="$(sed -n 's/^sample_ids=//p' "${run_root}/RUN_MANIFEST.txt")"
[[ -n "${sample_ids}" ]]
IFS=',' read -r -a samples <<< "${sample_ids}"
for sample in "${samples[@]}"; do
  [[ -d "${run_root}/per_sample/${sample}/Gene/filtered" ]]
done
compgen -G "${run_root}/per_sample/*/Gene/filtered/matrix.mtx" >/dev/null
grep -E "^hash_cache_source=(given|generated)$" "${run_root}/RUN_MANIFEST.txt" >/dev/null
grep -F "Flex probe route: half-probe H1X2 (1.9.4 default)" "${run_root}/Log.out" >/dev/null
# The input helper records the config's resolved path.
grep -F "cr_config=$(realpath -- "${CR_CONFIG}")" "${run_root}/RUN_MANIFEST.txt" >/dev/null

echo "PASS: ${run_root}"

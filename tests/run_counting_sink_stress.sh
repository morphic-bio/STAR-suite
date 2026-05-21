#!/usr/bin/env bash
# Replay synthetic counted records into CountingSink without mapping.
# Gate before full production: peak RSS should scale ~linearly with N records.
#
# Examples:
#   STAR_SOLO_MEMORY_PROFILE=1 tests/run_counting_sink_stress.sh 50000000
#   STAR_COUNTING_SINK_STRESS_COLLAPSE=1 tests/run_counting_sink_stress.sh 10000000
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
# shellcheck disable=SC1090
source "${REPO_ROOT}/tests/external_fixtures_env.sh" 2>/dev/null || true

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
NRECORDS="${1:-50000000}"
ACTIVE_CBS="${STAR_COUNTING_SINK_STRESS_ACTIVE_CBS:-20000}"
GENOME_DIR="${STAR_COUNTING_SINK_STRESS_GENOME_DIR:-${UCSF_100K_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}}"
WHITELIST="${STAR_COUNTING_SINK_STRESS_WHITELIST:-${UCSF_100K_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}}"
OUT_PREFIX="${STAR_COUNTING_SINK_STRESS_OUT:-/tmp/counting_sink_stress_$(date -u +%Y%m%dT%H%M%SZ)/}"

export STAR_SOLO_MEMORY_PROFILE="${STAR_SOLO_MEMORY_PROFILE:-1}"
export STAR_COUNTING_SINK_STRESS_NRECORDS="${NRECORDS}"
export STAR_COUNTING_SINK_STRESS_ACTIVE_CBS="${ACTIVE_CBS}"

mkdir -p "${OUT_PREFIX}"

[[ -x "${STAR_BIN}" ]] || { echo "STAR missing: ${STAR_BIN}" >&2; exit 1; }
[[ -d "${GENOME_DIR}" ]] || { echo "genomeDir missing: ${GENOME_DIR}" >&2; exit 1; }
[[ -f "${WHITELIST}" ]] || { echo "whitelist missing: ${WHITELIST}" >&2; exit 1; }

echo "counting_sink_stress: nRecords=${NRECORDS} activeCBs=${ACTIVE_CBS} out=${OUT_PREFIX}"

"${STAR_BIN}" \
  --runMode countingSinkStress \
  --runThreadN 1 \
  --genomeDir "${GENOME_DIR}" \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloCBwhitelist "${WHITELIST}" \
  --soloFeatures GeneFull \
  --outFileNamePrefix "${OUT_PREFIX}"

LOG="${OUT_PREFIX}Log.out"
if [[ -f "${LOG}" ]]; then
  "${REPO_ROOT}/tests/run_solo_memory_profile_harness.sh" --parse-log "${LOG}"
else
  echo "WARN: missing ${LOG}" >&2
  exit 1
fi

#!/usr/bin/env bash
# Build a small GEX-only debug fixture from corrected UCSF EBs2_2.
# The fixture preserves per-lane FASTQ filenames but downscales each FASTQ to
# approximately TOTAL_READS / lane_count reads.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DOWNSAMPLE="${SCRIPT_DIR}/downsample_fastq_gz.sh"

SOURCE_GEX_DIR="${UCSF_VELOCYTO_GEX_SOURCE_DIR:-/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/GEX}"
FIXTURE_ROOT="${UCSF_VELOCYTO_GEX_100K_FIXTURE_ROOT:-/storage/ucsf-velocyto-validation/fixtures/ebs2_2_gexonly_100k}"
TOTAL_READS="${UCSF_VELOCYTO_GEX_100K_TOTAL_READS:-100000}"
MODE="${UCSF_VELOCYTO_GEX_100K_MODE:-head}"
SEED="${UCSF_VELOCYTO_GEX_100K_SEED:-1}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-gex-dir) SOURCE_GEX_DIR="${2:?}"; shift 2 ;;
    --fixture-root) FIXTURE_ROOT="${2:?}"; shift 2 ;;
    --total-reads) TOTAL_READS="${2:?}"; shift 2 ;;
    --mode) MODE="${2:?}"; shift 2 ;;
    --seed) SEED="${2:?}"; shift 2 ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

die() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ -x "${DOWNSAMPLE}" ]] || die "Missing helper: ${DOWNSAMPLE}"
[[ -d "${SOURCE_GEX_DIR}" ]] || die "Missing source GEX dir: ${SOURCE_GEX_DIR}"

FIXTURE_GEX_DIR="${FIXTURE_ROOT}/GEX"
mkdir -p "${FIXTURE_GEX_DIR}"

lane_count=$(find -L "${SOURCE_GEX_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | wc -l | tr -d ' ')
[[ "${lane_count}" -gt 0 ]] || die "No R1 FASTQs under ${SOURCE_GEX_DIR}"
reads_per_fastq=$(( (TOTAL_READS + lane_count - 1) / lane_count ))

while IFS= read -r fastq; do
  base="$(basename "${fastq}")"
  dest="${FIXTURE_GEX_DIR}/${base}"
  if [[ -s "${dest}" ]]; then
    continue
  fi
  bash "${DOWNSAMPLE}" \
    --in "${fastq}" \
    --out "${dest}" \
    --reads "${reads_per_fastq}" \
    --seed "${SEED}" \
    --mode "${MODE}"
done < <(find -L "${SOURCE_GEX_DIR}" -maxdepth 1 -type f -name '*.fastq.gz' | sort)

MANIFEST="${FIXTURE_ROOT}/FIXTURE_MANIFEST.txt"
{
  echo "date=$(date -Iseconds)"
  echo "source_gex_dir=${SOURCE_GEX_DIR}"
  echo "fixture_gex_dir=${FIXTURE_GEX_DIR}"
  echo "total_reads_target=${TOTAL_READS}"
  echo "lane_count=${lane_count}"
  echo "reads_per_fastq=${reads_per_fastq}"
  echo "mode=${MODE}"
  echo "seed=${SEED}"
} > "${MANIFEST}"

echo "Fixture ready: ${FIXTURE_GEX_DIR}"
echo "Manifest: ${MANIFEST}"

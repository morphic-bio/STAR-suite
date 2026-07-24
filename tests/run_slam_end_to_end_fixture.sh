#!/usr/bin/env bash
set -euo pipefail

# Deprecated compatibility entrypoint.
#
# Fixture parity has one implementation: run_slam_parity_smoke.sh. This shim
# translates the historical fixture-E2E interface so external callers do not
# silently retain the old, incorrect alignment contract.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

WORK_DIR="${WORK_DIR:-}"
OVERWRITE=0
NO_CLEANUP=0
FORCE_NO_GEDI=0
PASSTHROUGH=()

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Deprecated compatibility wrapper for run_slam_parity_smoke.sh.

Legacy options:
  --work-dir DIR        Alias for --outdir DIR
  --overwrite           Accepted; canonical parity always runs in its outdir
  --no-cleanup          Accepted; canonical parity keeps its artifacts
  --no-gedi             Accepted; canonical parity uses the pinned reference

Canonical options such as --threads, --star-bin, --genome-dir, --fastq,
--snp-bed, --reference, and --fixture-manifest are passed through.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;
    --no-cleanup) NO_CLEANUP=1; shift ;;
    --no-gedi) FORCE_NO_GEDI=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) PASSTHROUGH+=("$1"); shift ;;
  esac
done

if [[ "${TRIM5:-0}" != "0" || "${TRIM3:-0}" != "0" ]]; then
  echo "ERROR: TRIM5/TRIM3 are not part of the pinned SLAM fixture contract." >&2
  exit 2
fi

# Historical variable aliases are translated only when their canonical names
# are unset. The fixture uses one FASTQ, not synthetic 0h/6h duplicate runs.
export FASTQ="${FASTQ:-${H6_FASTQ:-${WT_FASTQ:-}}}"
export SNPS_BED="${SNPS_BED:-${MASK_BED:-}}"
export REF_TSV="${REF_TSV:-${FIXTURE_REF:-}}"

if [[ -n "${WORK_DIR}" ]]; then
  PASSTHROUGH=(--outdir "${WORK_DIR}" "${PASSTHROUGH[@]}")
fi

echo "NOTICE: run_slam_end_to_end_fixture.sh is deprecated; delegating to run_slam_parity_smoke.sh." >&2
if [[ "${OVERWRITE}" == "1" || "${NO_CLEANUP}" == "1" || "${FORCE_NO_GEDI}" == "1" ]]; then
  echo "NOTICE: legacy execution flags are already the canonical runner behavior." >&2
fi

exec bash "${SCRIPT_DIR}/run_slam_parity_smoke.sh" "${PASSTHROUGH[@]}"

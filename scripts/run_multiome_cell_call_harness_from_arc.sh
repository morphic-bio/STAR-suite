#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_multiome_cell_call_harness_from_arc.sh [options]

Required:
  --arc-run-dir DIR                Completed Cell Ranger ARC run dir
  --out-root DIR                   Root directory for harness outputs

Optional:
  --profile NAME                   Profile to run; may be provided multiple times
                                   Defaults:
                                     default_full
                                     bootstrap_full
                                     grayzone_bootstrap_rawp
  --scripts-dir DIR                Override STAR-suite scripts dir
  -h, --help                       Show this help

Outputs under --out-root:
  profiles.tsv                     Profile definitions used by this harness
  harness_summary.tsv              Collected comparison metrics per profile
  <profile>/...                    Full output tree from run_multiome_cell_call_external_gex_from_arc.sh
EOF
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

ARC_RUN_DIR=""
OUT_ROOT=""
SCRIPTS_DIR="${SCRIPT_DIR}"
declare -a PROFILES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --arc-run-dir) ARC_RUN_DIR="$2"; shift 2;;
    --out-root) OUT_ROOT="$2"; shift 2;;
    --profile) PROFILES+=("$2"); shift 2;;
    --scripts-dir) SCRIPTS_DIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 2;;
  esac
done

if [[ -z "${ARC_RUN_DIR}" || -z "${OUT_ROOT}" ]]; then
  usage
  exit 2
fi

if [[ ! -d "${ARC_RUN_DIR}" ]]; then
  echo "Missing ARC run dir: ${ARC_RUN_DIR}" >&2
  exit 2
fi

RUNNER="${SCRIPTS_DIR}/run_multiome_cell_call_external_gex_from_arc.sh"
if [[ ! -x "${RUNNER}" ]]; then
  echo "Missing runner: ${RUNNER}" >&2
  exit 2
fi

if [[ ${#PROFILES[@]} -eq 0 ]]; then
  PROFILES=(
    default_full
    bootstrap_full
    grayzone_bootstrap_rawp
  )
fi

profile_args() {
  case "$1" in
    default_full)
      ;;
    bootstrap_full)
      printf '%s\n' --use-bootstrap
      ;;
    grayzone_bootstrap_rawp)
      printf '%s\n' \
        --use-bootstrap \
        --lower-testing-bound 1 \
        --ambient-umi-max 50 \
        --use-raw-pvalue \
        --raw-pvalue 0.05
      ;;
    *)
      echo "Unknown profile: $1" >&2
      return 1
      ;;
  esac
}

metric_value() {
  local summary_tsv="$1"
  local metric_name="$2"
  awk -F'\t' -v metric_name="${metric_name}" '$1 == metric_name { print $2; found=1 } END { if (!found) print "" }' "${summary_tsv}"
}

mkdir -p "${OUT_ROOT}"

PROFILES_TSV="${OUT_ROOT}/profiles.tsv"
{
  printf 'profile\targs\n'
  for profile in "${PROFILES[@]}"; do
    args_str="$(profile_args "${profile}" | paste -sd' ' -)"
    printf '%s\t%s\n' "${profile}" "${args_str}"
  done
} > "${PROFILES_TSV}"

SUMMARY_TSV="${OUT_ROOT}/harness_summary.tsv"
printf 'profile\trun_dir\tarc_positive\texternal_positive\ttrue_positive\ttrue_negative\tfalse_positive\tfalse_negative\taccuracy\tsensitivity\tspecificity\tprecision\n' > "${SUMMARY_TSV}"

for profile in "${PROFILES[@]}"; do
  profile_dir="${OUT_ROOT}/${profile}"
  mkdir -p "${profile_dir}"

  echo "[harness] Running profile=${profile} out_dir=${profile_dir}" >&2
  mapfile -t extra_args < <(profile_args "${profile}")

  cmd=(
    "${RUNNER}"
    --out-dir "${profile_dir}"
    --arc-run-dir "${ARC_RUN_DIR}"
  )
  if [[ ${#extra_args[@]} -gt 0 ]]; then
    cmd+=("${extra_args[@]}")
  fi

  "${cmd[@]}"

  comparison_summary="${profile_dir}/comparison_summary.tsv"
  if [[ ! -f "${comparison_summary}" ]]; then
    echo "Missing comparison summary for profile=${profile}: ${comparison_summary}" >&2
    exit 1
  fi

  arc_positive="$(metric_value "${comparison_summary}" arc_positive)"
  external_positive="$(metric_value "${comparison_summary}" external_positive_matched)"
  true_positive="$(metric_value "${comparison_summary}" true_positive)"
  true_negative="$(metric_value "${comparison_summary}" true_negative)"
  false_positive="$(metric_value "${comparison_summary}" false_positive)"
  false_negative="$(metric_value "${comparison_summary}" false_negative)"
  accuracy="$(metric_value "${comparison_summary}" accuracy)"
  sensitivity="$(metric_value "${comparison_summary}" sensitivity)"
  specificity="$(metric_value "${comparison_summary}" specificity)"
  precision="$(metric_value "${comparison_summary}" precision)"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${profile}" \
    "${profile_dir}" \
    "${arc_positive}" \
    "${external_positive}" \
    "${true_positive}" \
    "${true_negative}" \
    "${false_positive}" \
    "${false_negative}" \
    "${accuracy}" \
    "${sensitivity}" \
    "${specificity}" \
    "${precision}" >> "${SUMMARY_TSV}"
done

echo "Harness complete." >&2
echo "profiles=${PROFILES_TSV}" >&2
echo "summary=${SUMMARY_TSV}" >&2

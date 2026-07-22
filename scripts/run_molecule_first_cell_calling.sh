#!/usr/bin/env bash
set -euo pipefail

star=""
policy_root=""
out_root=""
policy="all"

usage() {
  cat <<'USAGE'
Usage: run_molecule_first_cell_calling.sh --star STAR --policy-root DIR --out-root DIR
       [--policy strict|hard|gated_hard|all]

Run EmptyDrops_CR after molecule-first UMI collapse on integer policy MEX
products. The real-valued soft_expected product is intentionally rejected.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --star) star="$2"; shift 2 ;;
    --policy-root) policy_root="$2"; shift 2 ;;
    --out-root) out_root="$2"; shift 2 ;;
    --policy) policy="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "${policy}" in
  all) policies=(strict hard gated_hard) ;;
  strict|hard|gated_hard) policies=("${policy}") ;;
  soft_expected)
    echo "ERROR: soft_expected is real-valued and must not be rounded for EmptyDrops_CR" >&2
    exit 2
    ;;
  *) echo "ERROR: unsupported policy: ${policy}" >&2; exit 2 ;;
esac

[[ -x "${star}" ]] || { echo "ERROR: STAR is not executable: ${star}" >&2; exit 2; }
[[ -d "${policy_root}" ]] || { echo "ERROR: policy root is missing: ${policy_root}" >&2; exit 2; }
[[ -n "${out_root}" ]] || { echo "ERROR: --out-root is required" >&2; exit 2; }

mkdir -p "${out_root}"
for current in "${policies[@]}"; do
  mex="${policy_root}/${current}/raw"
  matrix="${mex}/matrix.mtx"
  [[ -s "${matrix}" ]] || { echo "ERROR: missing policy matrix: ${matrix}" >&2; exit 2; }
  if ! head -n 1 "${matrix}" | grep -qx '%%MatrixMarket matrix coordinate integer general'; then
    echo "ERROR: refusing non-integer cell-calling input: ${matrix}" >&2
    exit 2
  fi
  output="${out_root}/${current}"
  [[ ! -e "${output}" ]] || { echo "ERROR: refusing to reuse output: ${output}" >&2; exit 2; }
  mkdir -p "${output}"
  "${star}" --runMode soloCellFiltering "${mex}" "${output}/" \
    --soloCellFilter EmptyDrops_CR --soloEmptyDropsLegacy no \
    --outFileNamePrefix "${output}/" \
    >"${output}/stdout.log" 2>"${output}/stderr.log"
  [[ -f "${output}/barcodes.tsv" ]] || { echo "ERROR: EmptyDrops output missing for ${current}" >&2; exit 2; }
done

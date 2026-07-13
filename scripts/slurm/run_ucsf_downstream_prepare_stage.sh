#!/usr/bin/env bash
set -euo pipefail

RUN_DIR=""
OUTPUT_DIR=""
TOOLS_DIR="/opt/star-suite/scripts"
MIN_GENES="${MIN_GENES:-200}"
MAX_GENES="${MAX_GENES:-2500}"
MT_PCT_CUTOFF="${MT_PCT_CUTOFF:-5}"
N_MAD="${N_MAD:-3}"

usage() {
  cat <<'EOF'
Usage: run_ucsf_downstream_prepare_stage.sh --run-dir DIR --output-dir DIR [options]

Build the raw-backed AnnData and QC views consumed by a separate CellBender job.

Options:
  --tools-dir DIR       STAR Suite helper directory
  --min-genes N         Minimum detected genes (default: 200)
  --max-genes N         Non-adaptive maximum genes fallback (default: 2500)
  --mt-pct-cutoff N     Mitochondrial percentage floor (default: 5)
  --n-mad N             Adaptive threshold MAD multiplier (default: 3)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --tools-dir) TOOLS_DIR="$2"; shift 2 ;;
    --min-genes) MIN_GENES="$2"; shift 2 ;;
    --max-genes) MAX_GENES="$2"; shift 2 ;;
    --mt-pct-cutoff) MT_PCT_CUTOFF="$2"; shift 2 ;;
    --n-mad) N_MAD="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "${RUN_DIR}" ]] || { echo "ERROR: --run-dir is required" >&2; exit 2; }
[[ -n "${OUTPUT_DIR}" ]] || { echo "ERROR: --output-dir is required" >&2; exit 2; }

RUN_DIR="$(realpath "${RUN_DIR}")"
OUTPUT_DIR="$(realpath -m "${OUTPUT_DIR}")"
TOOLS_DIR="$(realpath "${TOOLS_DIR}")"

BUILD_COUNTS="${TOOLS_DIR}/build_gene_full_velocyto_h5ad.py"
DOUBLET_SCRIPT="${TOOLS_DIR}/run_star_cell_doublets.R"
COMPUTE_ADAPTIVE_QC="${TOOLS_DIR}/compute_adaptive_qc_threshold.py"
APPLY_ADAPTIVE_MT="${TOOLS_DIR}/apply_adaptive_mt_filter.py"
GENERATE_QC_HISTOGRAM="${TOOLS_DIR}/generate_qc_histogram_mt_adaptive.py"
POSTPROCESS_FILTERS="${TOOLS_DIR}/postprocess_downstream_filters.py"
COMBINE_FILTERS="${TOOLS_DIR}/combine_ucsf_filters.py"

for helper in \
  "${BUILD_COUNTS}" \
  "${DOUBLET_SCRIPT}" \
  "${COMPUTE_ADAPTIVE_QC}" \
  "${APPLY_ADAPTIVE_MT}" \
  "${GENERATE_QC_HISTOGRAM}" \
  "${POSTPROCESS_FILTERS}" \
  "${COMBINE_FILTERS}"; do
  [[ -f "${helper}" ]] || { echo "ERROR: missing helper ${helper}" >&2; exit 1; }
done

[[ -d "${RUN_DIR}/outs/filtered_feature_bc_matrix" ]] || {
  echo "ERROR: missing filtered GeneFull MEX under ${RUN_DIR}" >&2
  exit 1
}
[[ -d "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix" ]] || {
  echo "ERROR: missing raw velocyto MEX under ${RUN_DIR}" >&2
  exit 1
}

mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/.numba" "${OUTPUT_DIR}/.matplotlib"
export NUMBA_CACHE_DIR="${OUTPUT_DIR}/.numba"
export NUMBA_DISABLE_JIT=1
export MPLCONFIGDIR="${OUTPUT_DIR}/.matplotlib"
export min_genes="${MIN_GENES}"
export max_genes="${MAX_GENES}"
export mt_pct_cutoff="${MT_PCT_CUTOFF}"

COUNTS_H5AD="${OUTPUT_DIR}/counts.h5ad"
UNFILTERED_H5AD="${OUTPUT_DIR}/unfiltered_counts.h5ad"
FILTERED_H5AD="${OUTPUT_DIR}/filtered_counts.h5ad"
SINGLET_H5AD="${OUTPUT_DIR}/default_singlet_filtered_counts.h5ad"
ADAPTIVE_QC_JSON="${OUTPUT_DIR}/adaptive_qc_threshold.json"

python3 "${BUILD_COUNTS}" --run-dir "${RUN_DIR}" --output-h5ad "${COUNTS_H5AD}"
Rscript "${DOUBLET_SCRIPT}" "${COUNTS_H5AD}"
python3 "${COMBINE_FILTERS}" \
  --input-file "${COUNTS_H5AD}" \
  --non-empty-barcodes "${OUTPUT_DIR}/non_empty_barcodes.txt" \
  --doublet-barcodes "${OUTPUT_DIR}/doublet_barcodes.txt" \
  --output-dir "${OUTPUT_DIR}" \
  --min-genes "${MIN_GENES}" \
  --max-genes "${MAX_GENES}" \
  --mt-pct-cutoff "${MT_PCT_CUTOFF}"

python3 "${COMPUTE_ADAPTIVE_QC}" \
  --counts-h5ad "${COUNTS_H5AD}" \
  --non-empty-barcodes "${OUTPUT_DIR}/non_empty_barcodes.txt" \
  --doublet-barcodes "${OUTPUT_DIR}/doublet_barcodes.txt" \
  --min-genes "${MIN_GENES}" \
  --n-mad "${N_MAD}" \
  --output-json "${ADAPTIVE_QC_JSON}" >/dev/null

python3 "${APPLY_ADAPTIVE_MT}" \
  --input-h5ad "${UNFILTERED_H5AD}" \
  --threshold-json "${ADAPTIVE_QC_JSON}" \
  --mt-floor "${MT_PCT_CUTOFF}" \
  --n-mad "${N_MAD}"

python3 "${GENERATE_QC_HISTOGRAM}" \
  --input-h5ad "${UNFILTERED_H5AD}" \
  --output-dir "${OUTPUT_DIR}" \
  --threshold-json "${ADAPTIVE_QC_JSON}"

python3 "${POSTPROCESS_FILTERS}" \
  --unfiltered-h5ad "${UNFILTERED_H5AD}" \
  --qc-output-h5ad "${FILTERED_H5AD}" \
  --default-singlet-output-h5ad "${SINGLET_H5AD}"

python3 - "${OUTPUT_DIR}" "${RUN_DIR}" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad

output_dir = Path(sys.argv[1])
run_dir = Path(sys.argv[2])
artifacts = {}
for name in (
    "counts.h5ad",
    "unfiltered_counts.h5ad",
    "filtered_counts.h5ad",
    "default_singlet_filtered_counts.h5ad",
    "adaptive_qc_threshold.json",
):
    path = output_dir / name
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"missing prepare artifact: {path}")
    artifacts[name] = {"path": str(path), "bytes": path.stat().st_size}

adata = ad.read_h5ad(output_dir / "unfiltered_counts.h5ad", backed="r")
shape = [int(adata.n_obs), int(adata.n_vars)]
adata.file.close()
payload = {
    "schema": "star_suite.ucsf_downstream_prepare/v1",
    "status": "complete",
    "completed_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    "run_dir": str(run_dir),
    "shape": shape,
    "artifacts": artifacts,
}
(output_dir / "PREP_COMPLETE.json").write_text(json.dumps(payload, indent=2) + "\n")
PY

echo "PASS: UCSF downstream prepare stage: ${OUTPUT_DIR}"

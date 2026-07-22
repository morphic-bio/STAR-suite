#!/usr/bin/env bash
set -euo pipefail

RUN_DIR=""
PREP_DIR=""
CELLBENDER_H5=""
OUTPUT_DIR=""
TOOLS_DIR="/opt/star-suite/scripts"
LAYER_NAME="denoised"

usage() {
  cat <<'EOF'
Usage: run_ucsf_downstream_finalize_stage.sh --run-dir DIR --prep-dir DIR \
  --cellbender-h5 FILE --output-dir DIR [options]

Options:
  --tools-dir DIR
  --layer-name NAME
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="$2"; shift 2 ;;
    --prep-dir) PREP_DIR="$2"; shift 2 ;;
    --cellbender-h5) CELLBENDER_H5="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --tools-dir) TOOLS_DIR="$2"; shift 2 ;;
    --layer-name) LAYER_NAME="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in RUN_DIR PREP_DIR CELLBENDER_H5 OUTPUT_DIR; do
  [[ -n "${!value}" ]] || { echo "ERROR: ${value} is required" >&2; exit 2; }
done

RUN_DIR="$(realpath "${RUN_DIR}")"
PREP_DIR="$(realpath "${PREP_DIR}")"
CELLBENDER_H5="$(realpath "${CELLBENDER_H5}")"
OUTPUT_DIR="$(realpath -m "${OUTPUT_DIR}")"
TOOLS_DIR="$(realpath "${TOOLS_DIR}")"

ADD_LAYER="${TOOLS_DIR}/add_cellbender_layer_from_h5.py"
PROPAGATE_LAYER="${TOOLS_DIR}/propagate_anndata_layer.py"
INTEGRATE_FEATURES="${TOOLS_DIR}/integrate_feature_library.py"
for helper in "${ADD_LAYER}" "${PROPAGATE_LAYER}" "${INTEGRATE_FEATURES}"; do
  [[ -f "${helper}" ]] || { echo "ERROR: missing helper ${helper}" >&2; exit 1; }
done

mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/feature_libraries"
for name in counts.h5ad unfiltered_counts.h5ad filtered_counts.h5ad default_singlet_filtered_counts.h5ad; do
  [[ -s "${PREP_DIR}/${name}" ]] || { echo "ERROR: missing prepare artifact ${PREP_DIR}/${name}" >&2; exit 1; }
  cp -f "${PREP_DIR}/${name}" "${OUTPUT_DIR}/${name}"
done
for name in adaptive_qc_threshold.json non_empty_barcodes.txt filtered_barcodes_with_scores.txt; do
  [[ -s "${PREP_DIR}/${name}" ]] || { echo "ERROR: missing prepare artifact ${PREP_DIR}/${name}" >&2; exit 1; }
  cp -f "${PREP_DIR}/${name}" "${OUTPUT_DIR}/${name}"
done
[[ -f "${PREP_DIR}/doublet_barcodes.txt" ]] || { echo "ERROR: missing prepare artifact ${PREP_DIR}/doublet_barcodes.txt" >&2; exit 1; }
cp -f "${PREP_DIR}/doublet_barcodes.txt" "${OUTPUT_DIR}/doublet_barcodes.txt"

for name in counts.h5ad unfiltered_counts.h5ad; do
  python3 "${ADD_LAYER}" \
    --cellbender-h5 "${CELLBENDER_H5}" \
    --input-h5ad "${OUTPUT_DIR}/${name}" \
    --output-h5ad "${OUTPUT_DIR}/${name}" \
    --layer-name "${LAYER_NAME}"
done

for name in filtered_counts.h5ad default_singlet_filtered_counts.h5ad; do
  python3 "${PROPAGATE_LAYER}" \
    --source-h5ad "${OUTPUT_DIR}/unfiltered_counts.h5ad" \
    --target-h5ad "${OUTPUT_DIR}/${name}" \
    --output-h5ad "${OUTPUT_DIR}/${name}" \
    --layer-name "${LAYER_NAME}"
done
cp -f "${OUTPUT_DIR}/unfiltered_counts.h5ad" "${OUTPUT_DIR}/final_counts.h5ad"

mapfile -t FEATURE_LIBRARY_DIRS < <(
  find "${RUN_DIR}/cr_assign" -type f -name pf_library_provenance.tsv -print 2>/dev/null \
    | sed 's#/pf_library_provenance.tsv$##' | sort
)
CRISPR_CALLS="${RUN_DIR}/outs/crispr_analysis/protospacer_calls_per_cell.csv"
for feature_library in "${FEATURE_LIBRARY_DIRS[@]}"; do
  args=(
    python3 "${INTEGRATE_FEATURES}"
    --library-dir "${feature_library}"
    --feature-output-root "${OUTPUT_DIR}/feature_libraries"
  )
  for name in counts.h5ad unfiltered_counts.h5ad filtered_counts.h5ad default_singlet_filtered_counts.h5ad final_counts.h5ad; do
    args+=(--counts-h5ad "${OUTPUT_DIR}/${name}")
  done
  if [[ "${feature_library}" == *CRISPR_Guide_Capture* && -f "${CRISPR_CALLS}" ]]; then
    args+=(--calls-csv "${CRISPR_CALLS}" --set-generic-aliases)
  fi
  "${args[@]}"
done

python3 - "${OUTPUT_DIR}" "${RUN_DIR}" "${CELLBENDER_H5}" "${LAYER_NAME}" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad

output_dir, run_dir, cellbender_h5 = map(Path, sys.argv[1:4])
layer_name = sys.argv[4]
artifacts = {}
allow_empty = {"doublet_barcodes.txt"}
for name in (
    "counts.h5ad",
    "unfiltered_counts.h5ad",
    "filtered_counts.h5ad",
    "default_singlet_filtered_counts.h5ad",
    "final_counts.h5ad",
    "adaptive_qc_threshold.json",
    "non_empty_barcodes.txt",
    "doublet_barcodes.txt",
    "filtered_barcodes_with_scores.txt",
):
    path = output_dir / name
    if not path.is_file() or (path.stat().st_size == 0 and name not in allow_empty):
        raise SystemExit(f"missing finalized artifact: {path}")
    artifacts[name] = {"path": str(path), "bytes": path.stat().st_size}

adata = ad.read_h5ad(output_dir / "final_counts.h5ad", backed="r")
if layer_name not in adata.layers:
    raise SystemExit(f"missing finalized layer: {layer_name}")
shape = [int(adata.n_obs), int(adata.n_vars)]
obs_columns = [str(name) for name in adata.obs.columns]
adata.file.close()
payload = {
    "schema": "star_suite.ucsf_downstream_finalize/v1",
    "status": "complete",
    "completed_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    "run_dir": str(run_dir),
    "cellbender_h5": str(cellbender_h5),
    "layer": layer_name,
    "shape": shape,
    "obs_columns": obs_columns,
    "artifacts": artifacts,
}
(output_dir / "FINALIZE_COMPLETE.json").write_text(json.dumps(payload, indent=2) + "\n")
PY

echo "PASS: UCSF downstream finalize stage: ${OUTPUT_DIR}"

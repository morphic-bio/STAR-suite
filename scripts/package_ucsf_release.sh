#!/usr/bin/env bash
# Package UCSF production run into a clean distribution layout for Globus upload.
#
# Creates a flat per-sample directory tree with hardlinks (no data duplication)
# from the authoritative production run root.  Generates samples.tsv manifest
# and copies the README into the output root.
#
# Usage:
#   scripts/package_ucsf_release.sh [options]
#
# Options:
#   --run-root DIR    Production run root (default: authoritative UCSF root)
#   --out-root DIR    Output directory for the release tree (required)
#   --copy            Use real copies instead of hardlinks
#   --dry-run         Print actions without creating files
#   -h, --help        Show help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RUN_ROOT="${UCSF_RUN_ROOT:-/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009}"
OUT_ROOT=""
DOWNSTREAM_SUFFIX="downstream_genefull_velocyto_cellbender"
LINK_MODE="hardlink"
DRY_RUN=0

SAMPLES=(
  EBs1_1 EBs1_2 EBs1_3 EBs1_4 EBs1_5
  EBs2_1 EBs2_2 EBs2_3 EBs2_4 EBs2_5
  iPSC1_1 iPSC1_2 iPSC1_3
  iPSC2_1 iPSC2_2 iPSC2_3
)

usage() {
  sed -n '2,/^$/s/^# \?//p' "$0"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-root)  RUN_ROOT="$2";  shift 2 ;;
    --out-root)  OUT_ROOT="$2";  shift 2 ;;
    --downstream-suffix) DOWNSTREAM_SUFFIX="$2"; shift 2 ;;
    --copy)      LINK_MODE="copy"; shift ;;
    --dry-run)   DRY_RUN=1;     shift ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -z "${OUT_ROOT}" ]]; then
  echo "ERROR: --out-root is required" >&2
  usage >&2
  exit 1
fi

[[ -d "${RUN_ROOT}/samples" ]] || { echo "ERROR: ${RUN_ROOT}/samples not found" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

link_file() {
  local src="$1" dst="$2"
  if [[ ! -f "${src}" ]]; then
    echo "  WARNING: missing ${src}" >&2
    return 1
  fi
  if (( DRY_RUN )); then
    echo "  ${LINK_MODE}: ${src} -> ${dst}"
    return 0
  fi
  if [[ "${LINK_MODE}" == "hardlink" ]]; then
    ln -f "${src}" "${dst}"
  else
    cp -f "${src}" "${dst}"
  fi
}

link_dir() {
  local src="$1" dst="$2"
  if [[ ! -d "${src}" ]]; then
    echo "  WARNING: missing ${src}" >&2
    return 1
  fi
  if (( DRY_RUN )); then
    echo "  ${LINK_MODE}-dir: ${src} -> ${dst}"
    return 0
  fi
  if [[ "${LINK_MODE}" == "hardlink" ]]; then
    cp -rl "${src}" "${dst}"
  else
    cp -r "${src}" "${dst}"
  fi
}

# ---------------------------------------------------------------------------
# Build the release tree
# ---------------------------------------------------------------------------

echo "=== UCSF Release Packaging ==="
echo "Run root:     ${RUN_ROOT}"
echo "Out root:     ${OUT_ROOT}"
echo "Downstream:   ${DOWNSTREAM_SUFFIX}"
echo "Link mode:    ${LINK_MODE}"
echo "Samples:      ${#SAMPLES[@]}"
echo ""

(( DRY_RUN )) || mkdir -p "${OUT_ROOT}"

for sample in "${SAMPLES[@]}"; do
  echo "--- ${sample} ---"
  src_ds="${RUN_ROOT}/samples/${sample}/${DOWNSTREAM_SUFFIX}"
  src_run="${RUN_ROOT}/samples/${sample}/run"
  dst="${OUT_ROOT}/${sample}"

  (( DRY_RUN )) || mkdir -p \
    "${dst}/qc" \
    "${dst}/features/crispr" \
    "${dst}/cellbender" \
    "${dst}/cellranger_compat"

  # -- Primary h5ads --
  link_file "${src_ds}/filtered_counts.h5ad"   "${dst}/filtered_counts.h5ad"
  link_file "${src_ds}/unfiltered_counts.h5ad" "${dst}/unfiltered_counts.h5ad"

  # -- QC --
  link_file "${src_ds}/summary.txt"                  "${dst}/qc/summary.txt"
  link_file "${src_ds}/adaptive_qc_threshold.json"   "${dst}/qc/adaptive_qc_threshold.json"
  link_file "${src_ds}/gene_quantile_histogram.png"  "${dst}/qc/gene_quantile_histogram.png"
  link_file "${src_ds}/gene_quantile_histogram.html" "${dst}/qc/gene_quantile_histogram.html" || true

  # -- Features / CRISPR --
  for f in protospacer_calls_per_cell.csv \
           protospacer_calls_summary.csv \
           protospacer_umi_thresholds.csv \
           protospacer_umi_thresholds.json; do
    link_file "${src_run}/outs/crispr_analysis/${f}" "${dst}/features/crispr/${f}"
  done

  # Feature library h5ads (one CRISPR library per sample)
  fl_dir="$(ls -d "${src_ds}/feature_libraries/${sample}_"* 2>/dev/null | head -1)"
  if [[ -n "${fl_dir}" ]]; then
    link_file "${fl_dir}/filtered_feature_library.h5ad" "${dst}/features/crispr/filtered_feature_library.h5ad"
    link_file "${fl_dir}/raw_feature_library.h5ad"      "${dst}/features/crispr/raw_feature_library.h5ad"
  else
    echo "  WARNING: no feature_libraries dir for ${sample}" >&2
  fi

  # -- CellBender --
  link_file "${src_ds}/cellbender/cellbender_counts.h5"          "${dst}/cellbender/cellbender_counts.h5"
  link_file "${src_ds}/cellbender/cellbender_counts.pdf"         "${dst}/cellbender/cellbender_counts.pdf"
  link_file "${src_ds}/cellbender/cellbender_counts_report.html" "${dst}/cellbender/cellbender_counts_report.html"
  link_file "${src_ds}/cellbender/cellbender_counts_metrics.csv" "${dst}/cellbender/cellbender_counts_metrics.csv"

  # -- CellRanger-compatible matrices --
  link_dir "${src_run}/outs/filtered_feature_bc_matrix"          "${dst}/cellranger_compat/filtered_feature_bc_matrix"
  link_dir "${src_run}/outs/raw_feature_bc_matrix"               "${dst}/cellranger_compat/raw_feature_bc_matrix"
  link_dir "${src_run}/outs/filtered_velocyto_feature_bc_matrix" "${dst}/cellranger_compat/filtered_velocyto_feature_bc_matrix"
  link_dir "${src_run}/outs/raw_velocyto_feature_bc_matrix"      "${dst}/cellranger_compat/raw_velocyto_feature_bc_matrix"
done

# ---------------------------------------------------------------------------
# Generate samples.tsv
# ---------------------------------------------------------------------------

echo ""
echo "=== Generating samples.tsv ==="

if (( DRY_RUN )); then
  echo "  (skipped in dry-run mode)"
else
  python3 - "${OUT_ROOT}" "${RUN_ROOT}" "${DOWNSTREAM_SUFFIX}" <<'PYEOF'
import json, sys
from pathlib import Path

import h5py

out_root = Path(sys.argv[1])
run_root = Path(sys.argv[2])
ds_suffix = sys.argv[3]

samples = [
    "EBs1_1","EBs1_2","EBs1_3","EBs1_4","EBs1_5",
    "EBs2_1","EBs2_2","EBs2_3","EBs2_4","EBs2_5",
    "iPSC1_1","iPSC1_2","iPSC1_3",
    "iPSC2_1","iPSC2_2","iPSC2_3",
]

header = [
    "sample_id",
    "filtered_cells",
    "singlet_cells",
    "star_called_cells",
    "effective_max_genes",
    "median_n_genes",
    "mad_n_genes",
    "min_genes",
    "n_mad",
]

rows = []
for sample in samples:
    ds = run_root / "samples" / sample / ds_suffix

    # Filtered cell count from h5ad shape
    with h5py.File(ds / "filtered_counts.h5ad", "r") as f:
        filtered_cells = f["obs"]["_index"].shape[0]

    # Adaptive QC thresholds
    with open(ds / "adaptive_qc_threshold.json") as fh:
        qc = json.load(fh)

    rows.append([
        sample,
        str(filtered_cells),
        str(qc["singlet_cells"]),
        str(qc["star_called_cells"]),
        str(qc["effective_max_genes"]),
        str(qc["median_n_genes"]),
        str(qc["mad_n_genes"]),
        str(qc["min_genes"]),
        str(qc["n_mad"]),
    ])

tsv = out_root / "samples.tsv"
with open(tsv, "w") as f:
    f.write("\t".join(header) + "\n")
    for row in rows:
        f.write("\t".join(row) + "\n")

print(f"  wrote {tsv}")
for row in rows:
    print(f"  {row[0]:10s}  filtered={row[1]:>6s}  singlet={row[2]:>6s}  star={row[3]:>6s}  max_genes={row[4]:>6s}")
PYEOF
fi

# ---------------------------------------------------------------------------
# Copy README
# ---------------------------------------------------------------------------

echo ""
echo "=== Copying README ==="
README_SRC="${REPO_ROOT}/scripts/ucsf_release_README.md"
if [[ -f "${README_SRC}" ]]; then
  if (( DRY_RUN )); then
    echo "  copy: ${README_SRC} -> ${OUT_ROOT}/README.md"
  else
    cp "${README_SRC}" "${OUT_ROOT}/README.md"
    echo "  wrote ${OUT_ROOT}/README.md"
  fi
else
  echo "  WARNING: ${README_SRC} not found; skipping README" >&2
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

echo ""
if (( DRY_RUN )); then
  echo "=== DRY RUN complete (no files created) ==="
else
  echo "=== Packaging complete ==="
  echo "Output: ${OUT_ROOT}"
  du -sh "${OUT_ROOT}" 2>/dev/null || true
fi

#!/usr/bin/env bash
# Paper parity: STAR-Flex vs CellRanger (CR7 / CR9) on SC2300771 full dataset.
#
# Dataset:  SC2300771 Flex (4 probe-barcode tags: BC004/BC006/BC007/BC008)
# Reads:    2.011B paired-end across 8 lanes
# Samples:  WT-Day-7 (BC004), PAX6-PTC-D9-Day7 (BC006),
#           WT-Day-8 (BC007), PAX6-PTC-D9-Day8 (BC008)
#
# Runs parity analysis (barcode Jaccard, per-barcode UMI Pearson/Spearman,
# per-gene UMI Pearson/Spearman) for each tag, then writes a combined
# summary TSV suitable for inclusion in the paper.
#
# Usage:
#   bash scripts/paper/run_flex_parity.sh                      # defaults
#   STAR_MEX_ROOT=/path/to/star/-  bash scripts/paper/run_flex_parity.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

PARITY_SCRIPT="${SCRIPT_DIR}/compute_parity_metrics.py"
[[ -f "${PARITY_SCRIPT}" ]] || { echo "ERROR: ${PARITY_SCRIPT} not found" >&2; exit 1; }

# ── Configurable paths ──────────────────────────────────────────────
# STAR per-sample filtered MEX (default: written by flexfilter to <cwd>/-/)
STAR_MEX_ROOT="${STAR_MEX_ROOT:-${REPO_ROOT}/-}"

# CellRanger 9.0.1 full-scale reference (March 19 2026 run, 32 cores, ~59 min)
CR9_ROOT="${CR9_ROOT:-/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs}"

# CellRanger 7.0.0 full-scale reference (Dec 2023 production run, 32 cores)
CR7_ROOT="${CR7_ROOT:-/home/lhhung/cellranger-multi/per_sample_outs}"

OUTDIR="${FLEX_PARITY_OUTDIR:-${REPO_ROOT}/comparisons/flex_parity_$(date +%Y%m%d_%H%M%S)}"

# ── Tag → sample-name mapping (from cellranger multi config) ────────
declare -A TAG_TO_CR_SAMPLE=(
    [BC004]="WT-Day-7"
    [BC006]="PAX6-PTC-D9-Day7"
    [BC007]="WT-Day-8"
    [BC008]="PAX6-PTC-D9-Day8"
)

TAGS="BC004 BC006 BC007 BC008"

# ── Validate inputs ─────────────────────────────────────────────────
for tag in ${TAGS}; do
    star_dir="${STAR_MEX_ROOT}/${tag}/Gene/filtered"
    [[ -d "${star_dir}" ]] || { echo "ERROR: STAR MEX not found: ${star_dir}" >&2; exit 1; }
done

mkdir -p "${OUTDIR}"

echo "================================================================="
echo "  Flex Parity: STAR-Flex vs CellRanger (SC2300771, 2.011B reads)"
echo "================================================================="
echo "STAR MEX root:   ${STAR_MEX_ROOT}"
echo "CR9 root:        ${CR9_ROOT}"
echo "CR7 root:        ${CR7_ROOT}"
echo "Output dir:      ${OUTDIR}"
echo ""

COMBINED_TSV="${OUTDIR}/flex_parity_combined.tsv"
HEADER_WRITTEN=false

# ── Per-tag parity: STAR vs CR9 ────────────────────────────────────
for tag in ${TAGS}; do
    sample="${TAG_TO_CR_SAMPLE[${tag}]}"
    star_dir="${STAR_MEX_ROOT}/${tag}/Gene/filtered"

    # CR9
    cr9_dir="${CR9_ROOT}/${sample}/count/sample_filtered_feature_bc_matrix"
    if [[ -d "${cr9_dir}" ]]; then
        echo "── ${tag} (${sample}) vs CR9 ──"
        tag_tsv="${OUTDIR}/${tag}_vs_cr9.tsv"
        python3 "${PARITY_SCRIPT}" \
            --cr-mex "${cr9_dir}" \
            --star-mex "${star_dir}" \
            --label "${tag}_vs_CR9" \
            --feature-types "Gene Expression" \
            --barcode-length 16 \
            --tsv-out "${tag_tsv}" \
            2>&1 | tee "${OUTDIR}/${tag}_vs_cr9.log"

        if [[ "${HEADER_WRITTEN}" == "false" ]]; then
            head -1 "${tag_tsv}" >> "${COMBINED_TSV}"
            HEADER_WRITTEN=true
        fi
        tail -n +2 "${tag_tsv}" >> "${COMBINED_TSV}"
        echo ""
    else
        echo "WARN: CR9 MEX not found for ${sample}: ${cr9_dir}"
    fi

    # CR7 — try BC tag name first, then biological sample name
    cr7_mex=""
    for cr7_name in "${tag}" "${sample}"; do
        for subpath in "count/sample_filtered_feature_bc_matrix" "filtered_feature_bc_matrix"; do
            candidate="${CR7_ROOT}/${cr7_name}/${subpath}"
            if [[ -d "${candidate}" ]]; then
                cr7_mex="${candidate}"
                break 2
            fi
        done
    done

    if [[ -n "${cr7_mex}" && -d "${cr7_mex}" ]]; then
        echo "── ${tag} (${sample}) vs CR7 ──"
        tag_tsv="${OUTDIR}/${tag}_vs_cr7.tsv"
        python3 "${PARITY_SCRIPT}" \
            --cr-mex "${cr7_mex}" \
            --star-mex "${star_dir}" \
            --label "${tag}_vs_CR7" \
            --feature-types "Gene Expression" \
            --barcode-length 16 \
            --tsv-out "${tag_tsv}" \
            2>&1 | tee "${OUTDIR}/${tag}_vs_cr7.log"

        tail -n +2 "${tag_tsv}" >> "${COMBINED_TSV}"
        echo ""
    else
        echo "INFO: CR7 per-tag MEX not found for ${tag}/${sample}, skipping CR7 comparison"
    fi
done

# ── Summary ─────────────────────────────────────────────────────────
echo "================================================================="
echo "  Combined parity summary: ${COMBINED_TSV}"
echo "================================================================="
if [[ -f "${COMBINED_TSV}" ]]; then
    column -t -s$'\t' "${COMBINED_TSV}"
fi
echo ""
echo "Per-tag logs and TSVs saved to: ${OUTDIR}/"

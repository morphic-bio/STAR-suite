#!/usr/bin/env bash
set -euo pipefail

# UCSF call-only parity recovery wrapper.
# This script intentionally gates only the call-only parity path first.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

DRIVER_SCRIPT="${DRIVER_SCRIPT:-${REPO_ROOT}/comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/run_callonly_gmm_parity.sh}"
STAR_FEATURE_CALL_BIN="${STAR_FEATURE_CALL_BIN:-${REPO_ROOT}/core/legacy/source/star_feature_call}"
CR_MEX_DIR="${CR_MEX_DIR:-/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/filtered_feature_bc_matrix}"
CR_CALLS_DIR="${CR_CALLS_DIR:-/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/crispr_analysis}"
OUT_BASE="${OUT_BASE:-/tmp/ucsf_call_parity_recovery_$(date +%Y%m%d_%H%M%S)}"
MIN_UMI_LIST="${MIN_UMI_LIST:-10 3}"
RUN_BUILD_HYGIENE="${RUN_BUILD_HYGIENE:-0}"
GATE_MIN_UMI="${GATE_MIN_UMI:-3}"

REAL_MISMATCH_MAX="${REAL_MISMATCH_MAX:-1}"
SET_EQUIV_PCT_MIN="${SET_EQUIV_PCT_MIN:-99.98}"
DELTA_UMI_COMMON_MAX="${DELTA_UMI_COMMON_MAX:-10}"
PEARSON_ALL_MIN="${PEARSON_ALL_MIN:-0.99999}"
SPEARMAN_ALL_MIN="${SPEARMAN_ALL_MIN:-0.99999}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --out-base DIR                  Output directory (default: ${OUT_BASE})
  --min-umi-list "10 3"           Space-separated min-UMI values (default: ${MIN_UMI_LIST})
  --star-feature-call-bin FILE    star_feature_call binary path
  --cr-mex-dir DIR                CR filtered_feature_bc_matrix directory
  --cr-calls-dir DIR              CR crispr_analysis directory
  --run-build-hygiene 0|1         If 1: make clean + rebuild STAR/star_feature_call
  -h, --help                      Show this help

Environment overrides are also supported for all options above.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-base)
      OUT_BASE="$2"
      shift 2
      ;;
    --min-umi-list)
      MIN_UMI_LIST="$2"
      shift 2
      ;;
    --star-feature-call-bin)
      STAR_FEATURE_CALL_BIN="$2"
      shift 2
      ;;
    --cr-mex-dir)
      CR_MEX_DIR="$2"
      shift 2
      ;;
    --cr-calls-dir)
      CR_CALLS_DIR="$2"
      shift 2
      ;;
    --run-build-hygiene)
      RUN_BUILD_HYGIENE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ ! -f "${DRIVER_SCRIPT}" ]]; then
  echo "ERROR: missing driver script: ${DRIVER_SCRIPT}" >&2
  exit 2
fi
if [[ ! -x "${STAR_FEATURE_CALL_BIN}" ]]; then
  echo "ERROR: star_feature_call binary not executable: ${STAR_FEATURE_CALL_BIN}" >&2
  exit 2
fi
if [[ ! -d "${CR_MEX_DIR}" ]]; then
  echo "ERROR: missing CR MEX directory: ${CR_MEX_DIR}" >&2
  exit 2
fi
if [[ ! -d "${CR_CALLS_DIR}" ]]; then
  echo "ERROR: missing CR calls directory: ${CR_CALLS_DIR}" >&2
  exit 2
fi
for req in \
  "${CR_MEX_DIR}/features.tsv.gz" \
  "${CR_MEX_DIR}/barcodes.tsv.gz" \
  "${CR_MEX_DIR}/matrix.mtx.gz" \
  "${CR_CALLS_DIR}/protospacer_calls_per_cell.csv" \
  "${CR_CALLS_DIR}/protospacer_umi_thresholds.csv"; do
  if [[ ! -f "${req}" ]]; then
    echo "ERROR: missing required input: ${req}" >&2
    exit 2
  fi
done

mkdir -p "${OUT_BASE}"
RUN_LOG="${OUT_BASE}/run_recovery.log"

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "repo_root=${REPO_ROOT}"
  echo "driver_script=${DRIVER_SCRIPT}"
  echo "star_feature_call_bin=${STAR_FEATURE_CALL_BIN}"
  echo "cr_mex_dir=${CR_MEX_DIR}"
  echo "cr_calls_dir=${CR_CALLS_DIR}"
  echo "out_base=${OUT_BASE}"
  echo "min_umi_list=${MIN_UMI_LIST}"
  echo "run_build_hygiene=${RUN_BUILD_HYGIENE}"
  echo "threshold_real_mismatch_max=${REAL_MISMATCH_MAX}"
  echo "threshold_set_equiv_pct_min=${SET_EQUIV_PCT_MIN}"
  echo "threshold_delta_umi_common_max=${DELTA_UMI_COMMON_MAX}"
  echo "threshold_pearson_all_min=${PEARSON_ALL_MIN}"
  echo "threshold_spearman_all_min=${SPEARMAN_ALL_MIN}"
  echo "gate_min_umi=${GATE_MIN_UMI}"
} > "${RUN_LOG}"

if [[ "${RUN_BUILD_HYGIENE}" == "1" ]]; then
  echo "[INFO] Running build hygiene (clean + rebuild)." | tee -a "${RUN_LOG}"
  make -C "${REPO_ROOT}/core/legacy/source" clean | tee -a "${RUN_LOG}"
  make -C "${REPO_ROOT}/core/legacy/source" -j8 STAR star_feature_call | tee -a "${RUN_LOG}"
fi

echo "[INFO] Running canonical call-only parity driver." | tee -a "${RUN_LOG}"
{
  STAR_FEATURE_CALL_BIN="${STAR_FEATURE_CALL_BIN}" \
  CR_MEX_DIR="${CR_MEX_DIR}" \
  CR_CALLS_DIR="${CR_CALLS_DIR}" \
  OUT_BASE="${OUT_BASE}" \
  MIN_UMI_LIST="${MIN_UMI_LIST}" \
  bash "${DRIVER_SCRIPT}"
} | tee -a "${RUN_LOG}"

for min_umi in ${MIN_UMI_LIST}; do
  for req in \
    "${OUT_BASE}/parity_minumi${min_umi}.tsv" \
    "${OUT_BASE}/star_feature_call_minumi${min_umi}.log" \
    "${OUT_BASE}/call_only_minumi${min_umi}/crispr_analysis/protospacer_calls_per_cell.csv"; do
    if [[ ! -f "${req}" ]]; then
      echo "ERROR: missing expected output after driver run: ${req}" >&2
      exit 2
    fi
  done
done

python3 - << 'PY' \
  "${OUT_BASE}" \
  "${CR_CALLS_DIR}" \
  "${MIN_UMI_LIST}" \
  "${GATE_MIN_UMI}" \
  "${REAL_MISMATCH_MAX}" \
  "${SET_EQUIV_PCT_MIN}" \
  "${DELTA_UMI_COMMON_MAX}" \
  "${PEARSON_ALL_MIN}" \
  "${SPEARMAN_ALL_MIN}"
import csv
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


def safe_int(value: str) -> int:
    v = (value or "").strip()
    if not v:
        return 0
    try:
        return int(v)
    except ValueError:
        try:
            return int(float(v))
        except ValueError:
            return 0


def parse_num_umis(value: str) -> int:
    v = (value or "").strip()
    if not v:
        return 0
    return sum(safe_int(tok) for tok in v.split("|"))


def norm_call(value: str) -> str:
    v = (value or "").strip()
    if not v or v.lower() == "none":
        return "None"
    parts = sorted(p.strip() for p in v.split("|") if p.strip())
    return "|".join(parts) if parts else "None"


def read_calls(path: Path) -> dict:
    out = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            return out
        key_col = "cell_barcode" if "cell_barcode" in reader.fieldnames else reader.fieldnames[0]
        for row in reader:
            bc = (row.get(key_col, "") or "").strip()
            if not bc:
                continue
            out[bc] = {
                "call": norm_call(row.get("feature_call", "")),
                "num_features": safe_int(row.get("num_features", "")),
                "umis": parse_num_umis(row.get("num_umis", "")),
            }
    return out


def read_parity_metrics(path: Path) -> dict:
    out = {}
    with open(path) as fh:
        header = fh.readline()
        if not header.startswith("metric\t"):
            raise RuntimeError(f"Unexpected parity header in {path}: {header.strip()}")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            out[parts[0]] = parts[1]
    return out


def rankdata(values):
    pairs = sorted((v, i) for i, v in enumerate(values))
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        current = pairs[i][0]
        while j < len(values) and pairs[j][0] == current:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[pairs[k][1]] = avg_rank
        i = j
    return ranks


def pearson(xs, ys):
    n = len(xs)
    if n == 0:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0.0 or den_y == 0.0:
        return float("nan")
    return num / (den_x * den_y)


def spearman(xs, ys):
    if not xs:
        return float("nan")
    return pearson(rankdata(xs), rankdata(ys))


def parse_runtime_real_seconds(path: Path):
    text = path.read_text()
    matches = re.findall(r"^real\s+([0-9.]+)$", text, flags=re.MULTILINE)
    if not matches:
        return float("nan")
    return float(matches[-1])


def fmt_float(v):
    if math.isnan(v):
        return "NA"
    return f"{v:.10f}"


out_base = Path(sys.argv[1])
cr_calls_dir = Path(sys.argv[2])
min_umi_list = [tok for tok in sys.argv[3].split() if tok.strip()]
gate_min_umi = str(sys.argv[4]).strip()
thr_real_mm = int(sys.argv[5])
thr_seteq_pct = float(sys.argv[6])
thr_delta = int(sys.argv[7])
thr_pearson = float(sys.argv[8])
thr_spearman = float(sys.argv[9])

if not min_umi_list:
    raise RuntimeError("MIN_UMI_LIST resolved to empty list")

cr_calls = read_calls(cr_calls_dir / "protospacer_calls_per_cell.csv")

summary_rows = []
for min_umi in min_umi_list:
    parity_path = out_base / f"parity_minumi{min_umi}.tsv"
    star_calls_path = out_base / f"call_only_minumi{min_umi}" / "crispr_analysis" / "protospacer_calls_per_cell.csv"
    runtime_log = out_base / f"star_feature_call_minumi{min_umi}.log"

    parity = read_parity_metrics(parity_path)
    star_calls = read_calls(star_calls_path)

    common = sorted(set(cr_calls) & set(star_calls))
    star_called = [bc for bc in common if star_calls[bc]["call"] != "None"]
    set_equiv_bcs = [bc for bc in common if cr_calls[bc]["call"] == star_calls[bc]["call"]]

    cr_umis_common = [cr_calls[bc]["umis"] for bc in common]
    star_umis_common = [star_calls[bc]["umis"] for bc in common]
    cr_umis_star_called = [cr_calls[bc]["umis"] for bc in star_called]
    star_umis_star_called = [star_calls[bc]["umis"] for bc in star_called]
    cr_umis_seteq = [cr_calls[bc]["umis"] for bc in set_equiv_bcs]
    star_umis_seteq = [star_calls[bc]["umis"] for bc in set_equiv_bcs]

    cr_total_common = sum(cr_umis_common)
    star_total_common = sum(star_umis_common)
    delta_common = star_total_common - cr_total_common
    delta_common_pct = (100.0 * delta_common / cr_total_common) if cr_total_common else float("nan")

    row = {
        "min_umi": str(min_umi),
        "rows_cr": int(parity.get("rows_cr", len(cr_calls))),
        "rows_star": int(parity.get("rows_star", len(star_calls))),
        "common": int(parity.get("common", len(common))),
        "only_star": int(parity.get("only_star", len(set(star_calls) - set(cr_calls)))),
        "set_equiv": int(parity.get("set_equivalent", len(set_equiv_bcs))),
        "set_equiv_pct": float(parity.get("set_equivalent_pct", "nan")),
        "real_mismatch_count": int(parity.get("real_mismatch_count", "0")),
        "cr_total_umis_common": cr_total_common,
        "star_total_umis_common": star_total_common,
        "delta_common": delta_common,
        "delta_common_pct": delta_common_pct,
        "star_total_umis_allrows": sum(v["umis"] for v in star_calls.values()),
        "pearson_all": pearson(cr_umis_common, star_umis_common),
        "spearman_all": spearman(cr_umis_common, star_umis_common),
        "n_star_called": len(star_called),
        "pearson_star_called": pearson(cr_umis_star_called, star_umis_star_called),
        "spearman_star_called": spearman(cr_umis_star_called, star_umis_star_called),
        "n_set_equiv": len(set_equiv_bcs),
        "pearson_set_equiv": pearson(cr_umis_seteq, star_umis_seteq),
        "spearman_set_equiv": spearman(cr_umis_seteq, star_umis_seteq),
        "runtime_real_s": parse_runtime_real_seconds(runtime_log),
        "real_mismatch_barcodes": parity.get("real_mismatch_barcodes", ""),
        "only_star_non_none_detail": parity.get("only_star_non_none_detail", ""),
        "threshold_mismatch_detail": parity.get("threshold_mismatch_detail", ""),
    }
    summary_rows.append(row)

summary_cols = [
    "min_umi",
    "rows_cr",
    "rows_star",
    "common",
    "only_star",
    "set_equiv",
    "set_equiv_pct",
    "cr_total_umis_common",
    "star_total_umis_common",
    "delta_common",
    "delta_common_pct",
    "star_total_umis_allrows",
    "pearson_all",
    "spearman_all",
    "n_star_called",
    "pearson_star_called",
    "spearman_star_called",
    "n_set_equiv",
    "pearson_set_equiv",
    "spearman_set_equiv",
    "runtime_real_s",
]

umi_summary_path = out_base / "umi_correlation_summary.tsv"
with open(umi_summary_path, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(summary_cols)
    for r in summary_rows:
        writer.writerow([r[c] for c in summary_cols])

target = None
for r in summary_rows:
    if r["min_umi"] == gate_min_umi:
        target = r
        break
if target is None:
    target = summary_rows[-1]

status = "PASS"
reasons = []
if target["real_mismatch_count"] > thr_real_mm:
    status = "FAIL"
    reasons.append(f"real_mismatch_count {target['real_mismatch_count']} > {thr_real_mm}")
if target["set_equiv_pct"] < thr_seteq_pct:
    status = "FAIL"
    reasons.append(f"set_equiv_pct {target['set_equiv_pct']:.6f} < {thr_seteq_pct}")
if abs(target["delta_common"]) > thr_delta:
    status = "FAIL"
    reasons.append(f"abs(delta_common) {abs(target['delta_common'])} > {thr_delta}")
if (math.isnan(target["pearson_all"])) or target["pearson_all"] < thr_pearson:
    status = "FAIL"
    reasons.append(f"pearson_all {target['pearson_all']} < {thr_pearson}")
if (math.isnan(target["spearman_all"])) or target["spearman_all"] < thr_spearman:
    status = "FAIL"
    reasons.append(f"spearman_all {target['spearman_all']} < {thr_spearman}")

summary_tsv = out_base / "PARITY_RECOVERY_SUMMARY.tsv"
with open(summary_tsv, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["key", "value"])
    writer.writerow(["status", status])
    writer.writerow(["requested_gate_min_umi", gate_min_umi])
    writer.writerow(["evaluated_min_umi", target["min_umi"]])
    writer.writerow(["out_base", str(out_base)])
    writer.writerow(["cr_calls_dir", str(cr_calls_dir)])
    writer.writerow(["threshold_real_mismatch_max", thr_real_mm])
    writer.writerow(["threshold_set_equiv_pct_min", thr_seteq_pct])
    writer.writerow(["threshold_delta_umi_common_max", thr_delta])
    writer.writerow(["threshold_pearson_all_min", thr_pearson])
    writer.writerow(["threshold_spearman_all_min", thr_spearman])
    writer.writerow(["target_real_mismatch_count", target["real_mismatch_count"]])
    writer.writerow(["target_set_equiv_pct", target["set_equiv_pct"]])
    writer.writerow(["target_delta_common", target["delta_common"]])
    writer.writerow(["target_pearson_all", target["pearson_all"]])
    writer.writerow(["target_spearman_all", target["spearman_all"]])
    writer.writerow(["target_real_mismatch_barcodes", target["real_mismatch_barcodes"]])
    writer.writerow(["target_only_star_non_none_detail", target["only_star_non_none_detail"]])
    writer.writerow(["target_threshold_mismatch_detail", target["threshold_mismatch_detail"]])
    writer.writerow(["failure_reasons", " | ".join(reasons) if reasons else ""])

summary_md = out_base / "PARITY_RECOVERY_SUMMARY.md"
with open(summary_md, "w") as fh:
    fh.write("# UCSF Call-Only Parity Recovery Summary\n\n")
    fh.write(f"- date_utc: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
    fh.write(f"- status: **{status}**\n")
    fh.write(f"- requested_gate_min_umi: `{gate_min_umi}`\n")
    fh.write(f"- evaluated_min_umi: `{target['min_umi']}`\n")
    fh.write(f"- artifacts_dir: `{out_base}`\n")
    fh.write(f"- umi_correlation_summary: `{umi_summary_path}`\n\n")

    fh.write(f"## Acceptance Gate (`min-umi {target['min_umi']}`)\n\n")
    fh.write(f"- real_mismatch_count: `{target['real_mismatch_count']}` (<= `{thr_real_mm}`)\n")
    fh.write(f"- set_equivalent_pct: `{target['set_equiv_pct']:.6f}` (>= `{thr_seteq_pct}`)\n")
    fh.write(f"- delta_umi_common_rows: `{target['delta_common']}` (abs <= `{thr_delta}`)\n")
    fh.write(f"- pearson_all_cr_rows: `{fmt_float(target['pearson_all'])}` (>= `{thr_pearson}`)\n")
    fh.write(f"- spearman_all_cr_rows: `{fmt_float(target['spearman_all'])}` (>= `{thr_spearman}`)\n\n")

    fh.write("## Per-Run Metrics\n\n")
    fh.write("| min_umi | rows_cr | rows_star | common | set_equiv | set_equiv_pct | delta_common | pearson_all | spearman_all | runtime_real_s |\n")
    fh.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in summary_rows:
        fh.write(
            f"| {r['min_umi']} | {r['rows_cr']} | {r['rows_star']} | {r['common']} | "
            f"{r['set_equiv']} | {r['set_equiv_pct']:.6f} | {r['delta_common']} | "
            f"{fmt_float(r['pearson_all'])} | {fmt_float(r['spearman_all'])} | "
            f"{r['runtime_real_s'] if not math.isnan(r['runtime_real_s']) else 'NA'} |\n"
        )
    fh.write("\n")

    if status == "FAIL":
        fh.write("## Failure Reasons\n\n")
        for reason in reasons:
            fh.write(f"- {reason}\n")
        fh.write("\n")

    fh.write(f"## Actionable Diff Details (`min-umi {target['min_umi']}`)\n\n")
    fh.write(f"- real_mismatch_barcodes: `{target['real_mismatch_barcodes'] or 'NA'}`\n")
    fh.write(f"- only_star_non_none_detail: `{target['only_star_non_none_detail'] or 'NA'}`\n")
    fh.write(f"- threshold_mismatch_detail: `{target['threshold_mismatch_detail'] or 'NA'}`\n")

print(f"PASS_FAIL={status}")
print(f"SUMMARY_MD={summary_md}")
print(f"SUMMARY_TSV={summary_tsv}")
print(f"UMI_SUMMARY_TSV={umi_summary_path}")
if reasons:
    print("FAILURE_REASONS=" + " | ".join(reasons))
PY

echo "[INFO] Done. Artifacts in: ${OUT_BASE}"
echo "[INFO] Summary: ${OUT_BASE}/PARITY_RECOVERY_SUMMARY.md"
echo "[INFO] Summary TSV: ${OUT_BASE}/PARITY_RECOVERY_SUMMARY.tsv"
echo "[INFO] Correlation TSV: ${OUT_BASE}/umi_correlation_summary.tsv"

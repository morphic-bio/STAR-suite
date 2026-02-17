#!/usr/bin/env bash
set -euo pipefail

# Reproduce UCSF iPSC2_1_AALG2 call-only GMM parity checks against Cell Ranger
# using a CRISPR-only MEX extracted from sample_filtered_feature_bc_matrix.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

STAR_FEATURE_CALL_BIN="${STAR_FEATURE_CALL_BIN:-${REPO_ROOT}/core/legacy/source/star_feature_call}"
CR_MEX_DIR="${CR_MEX_DIR:-/mnt/pikachu/ucsf-perturb-seq/iPSC2_1_AALG2_autoindex11044_run3_20260216/outs/per_sample_outs/iPSC2_1_AALG2_autoindex11044_run3_20260216/count/sample_filtered_feature_bc_matrix}"
CR_CALLS_DIR="${CR_CALLS_DIR:-/mnt/pikachu/ucsf-perturb-seq/iPSC2_1_AALG2_autoindex11044_run3_20260216/outs/per_sample_outs/iPSC2_1_AALG2_autoindex11044_run3_20260216/count/crispr_analysis}"
OUT_BASE="${OUT_BASE:-/tmp/ipsc2_callonly_mex_parity_$(date +%Y%m%d_%H%M%S)}"
MIN_UMI_LIST="${MIN_UMI_LIST:-10 3}"

if [[ ! -x "${STAR_FEATURE_CALL_BIN}" ]]; then
  echo "ERROR: star_feature_call binary not executable: ${STAR_FEATURE_CALL_BIN}" >&2
  exit 1
fi
for d in "${CR_MEX_DIR}" "${CR_CALLS_DIR}"; do
  if [[ ! -d "${d}" ]]; then
    echo "ERROR: missing directory: ${d}" >&2
    exit 1
  fi
done

mkdir -p "${OUT_BASE}/crispr_mex"

# Step 1: Build CRISPR-only MEX from sample_filtered_feature_bc_matrix
python3 - << 'PY' "${CR_MEX_DIR}" "${OUT_BASE}/crispr_mex"
import gzip
import tempfile
from pathlib import Path
import sys

mex_in = Path(sys.argv[1])
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)

f_features = mex_in / "features.tsv.gz"
f_barcodes = mex_in / "barcodes.tsv.gz"
f_matrix = mex_in / "matrix.mtx.gz"

row_map = {}
new_features = []
with gzip.open(f_features, "rt") as fi:
    old_row = 0
    for line in fi:
        old_row += 1
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3 and parts[2] == "CRISPR Guide Capture":
            row_map[old_row] = len(new_features) + 1
            new_features.append(line.rstrip("\n"))

with open(out / "features.tsv", "w") as fo:
    for line in new_features:
        fo.write(line + "\n")

with gzip.open(f_barcodes, "rt") as fi, open(out / "barcodes.tsv", "w") as fo:
    for line in fi:
        fo.write(line)

with gzip.open(f_matrix, "rt") as fi:
    comments = []
    header = None
    for line in fi:
        if line.startswith("%"):
            comments.append(line)
            continue
        header = line.strip()
        break

    if header is None:
        raise RuntimeError("matrix.mtx.gz missing Matrix Market header")

    nrows, ncols, _ = map(int, header.split())
    body_tmp = tempfile.NamedTemporaryFile("w", delete=False, dir=str(out), prefix="matrix_body_", suffix=".tmp")
    kept = 0
    for line in fi:
        line = line.strip()
        if not line:
            continue
        r, c, v = line.split()
        r = int(r)
        nr = row_map.get(r)
        if nr is None:
            continue
        body_tmp.write(f"{nr} {c} {v}\n")
        kept += 1
    body_tmp_path = Path(body_tmp.name)
    body_tmp.close()

with open(out / "matrix.mtx", "w") as fo:
    for c in comments:
        fo.write(c)
    fo.write(f"{len(new_features)} {ncols} {kept}\n")
    with open(body_tmp_path, "r") as bi:
        for line in bi:
            fo.write(line)

body_tmp_path.unlink(missing_ok=True)
print(f"crispr_features={len(new_features)}")
print(f"matrix_nnz={kept}")
PY

# Step 2: Run call-only GMM and compare to CR for each min-umi
for min_umi in ${MIN_UMI_LIST}; do
  run_dir="${OUT_BASE}/call_only_minumi${min_umi}"
  log_file="${OUT_BASE}/star_feature_call_minumi${min_umi}.log"
  mkdir -p "${run_dir}"

  /usr/bin/time -p "${STAR_FEATURE_CALL_BIN}" \
    --call-only \
    --compat-perturb \
    --min-umi "${min_umi}" \
    --mex-dir "${OUT_BASE}/crispr_mex" \
    --output-dir "${run_dir}" \
    > "${log_file}" 2>&1

  python3 - << 'PY' "${CR_CALLS_DIR}" "${run_dir}/crispr_analysis" "${OUT_BASE}/parity_minumi${min_umi}.tsv"
import csv
import sys
from pathlib import Path

cr = Path(sys.argv[1])
st = Path(sys.argv[2])
out_tsv = Path(sys.argv[3])

def read_csv(path):
    with open(path, newline="") as fh:
        r = csv.DictReader(fh)
        rows = list(r)
        return r.fieldnames, rows

def norm_call(v):
    v = (v or "").strip()
    if v == "" or v.lower() == "none":
        return "None"
    return "|".join(sorted([x.strip() for x in v.split("|") if x.strip()]))

h_cr, r_cr = read_csv(cr / "protospacer_calls_per_cell.csv")
h_st, r_st = read_csv(st / "protospacer_calls_per_cell.csv")
key_cr = h_cr[0]
key_st = h_st[0]
call_col = "feature_call"

cr_map = {r[key_cr]: r for r in r_cr}
st_map = {r[key_st]: r for r in r_st}
common = sorted(set(cr_map) & set(st_map))
only_cr = sorted(set(cr_map) - set(st_map))
only_st = sorted(set(st_map) - set(cr_map))

row_exact = sum(1 for bc in common if cr_map[bc] == st_map[bc])
call_exact = sum(1 for bc in common if cr_map[bc][call_col] == st_map[bc][call_col])
set_equiv = sum(1 for bc in common if norm_call(cr_map[bc][call_col]) == norm_call(st_map[bc][call_col]))
real_mm = [bc for bc in common if norm_call(cr_map[bc][call_col]) != norm_call(st_map[bc][call_col])]

only_st_none = [bc for bc in only_st if norm_call(st_map[bc][call_col]) == "None"]
only_st_non_none = [bc for bc in only_st if bc not in only_st_none]

th_cr_h, th_cr_r = read_csv(cr / "protospacer_umi_thresholds.csv")
th_st_h, th_st_r = read_csv(st / "protospacer_umi_thresholds.csv")
th_key_cr = th_cr_h[0]
th_key_st = th_st_h[0]
cr_th = {r[th_key_cr]: r for r in th_cr_r}
st_th = {r[th_key_st]: r for r in th_st_r}
th_common = sorted(set(cr_th) & set(st_th))
th_only_cr = sorted(set(cr_th) - set(st_th))
th_only_st = sorted(set(st_th) - set(cr_th))
th_mm = [k for k in th_common if cr_th[k] != st_th[k]]

rows = [
    ("rows_cr", len(cr_map)),
    ("rows_star", len(st_map)),
    ("common", len(common)),
    ("only_cr", len(only_cr)),
    ("only_star", len(only_st)),
    ("row_exact", row_exact),
    ("call_exact", call_exact),
    ("set_equivalent", set_equiv),
    ("set_equivalent_pct", f"{(100.0 * set_equiv / len(common)):.4f}" if common else "NA"),
    ("real_mismatch_count", len(real_mm)),
    ("only_star_none", len(only_st_none)),
    ("only_star_non_none", len(only_st_non_none)),
    ("threshold_rows_cr", len(cr_th)),
    ("threshold_rows_star", len(st_th)),
    ("threshold_common", len(th_common)),
    ("threshold_only_cr", len(th_only_cr)),
    ("threshold_only_star", len(th_only_st)),
    ("threshold_mismatch_rows", len(th_mm)),
]

with open(out_tsv, "w") as fo:
    fo.write("metric\tvalue\n")
    for k, v in rows:
        fo.write(f"{k}\t{v}\n")
    if real_mm:
        fo.write("real_mismatch_barcodes\t" + ",".join(real_mm) + "\n")
    if only_st_non_none:
        detail = ";".join([f"{bc}:{st_map[bc][call_col]}" for bc in only_st_non_none])
        fo.write("only_star_non_none_detail\t" + detail + "\n")
    if th_mm:
        detail = ";".join([f"{k}:CR={cr_th[k]['UMI threshold']},STAR={st_th[k]['UMI threshold']}" for k in th_mm])
        fo.write("threshold_mismatch_detail\t" + detail + "\n")
    if th_only_cr:
        fo.write("threshold_only_cr_detail\t" + ",".join(th_only_cr) + "\n")
PY

done

echo "OUT_BASE=${OUT_BASE}"
echo "Artifacts:"
ls -1 "${OUT_BASE}" | sed 's/^/  - /'

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_SUITE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

ASSIGN_BIN="${ASSIGN_BIN:-${STAR_SUITE_DIR}/core/features/feature_barcodes/assignBarcodes}"
MEX_STUB="${MEX_STUB:-${STAR_SUITE_DIR}/tools/feature_barcodes/assignbarcodes_mex_stub.py}"

FASTQ_DIR="${FASTQ_DIR:-/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr/downsampled_100000_v2}"
WHITELIST="${WHITELIST:-/storage/A375/3M-5pgex-jan-2023.txt}"
FEATURE_REF="${FEATURE_REF:-/storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv}"
OUT_ROOT="${OUT_ROOT:-/tmp/a375_feature_sweep}"
FILTERED_BARCODES="${FILTERED_BARCODES:-}"

BARCODE_LEN="${BARCODE_LEN:-16}"
UMI_LEN="${UMI_LEN:-12}"
STRINGENCY="${STRINGENCY:-1}"
MIN_COUNTS="${MIN_COUNTS:-1}"

THREADS="${THREADS:-1}"
SEARCH_THREADS="${SEARCH_THREADS:-4}"
CONSUMER_THREADS="${CONSUMER_THREADS:-1}"
MAX_BARCODE_MISMATCHES="${MAX_BARCODE_MISMATCHES:-}"
MAX_BARCODE_MISMATCHES_LIST="${MAX_BARCODE_MISMATCHES_LIST:-}"
MIN_POSTERIOR="${MIN_POSTERIOR:-}"
CR_MATRIX_DIR="${CR_MATRIX_DIR:-}"
CR_FEATURE_TYPE="${CR_FEATURE_TYPE:-CRISPR Guide Capture}"

# Sweep parameters (space-separated lists). Defaults start at strict settings.
# NOTE: maxHammingDistance is inclusive, so 1 permits one mismatch.
MAX_HAMMING_LIST="${MAX_HAMMING_LIST:-0 1}"
FEATURE_N_LIST="${FEATURE_N_LIST:-0 1}"
BARCODE_N_LIST="${BARCODE_N_LIST:-1}"
FEATURE_OFFSET_LIST="${FEATURE_OFFSET_LIST:-0}"
LIMIT_SEARCH_LIST="${LIMIT_SEARCH_LIST:-"-1"}"

if [[ -z "${MAX_BARCODE_MISMATCHES_LIST}" ]]; then
    if [[ -n "${MAX_BARCODE_MISMATCHES}" ]]; then
        MAX_BARCODE_MISMATCHES_LIST="${MAX_BARCODE_MISMATCHES}"
    else
        MAX_BARCODE_MISMATCHES_LIST="0 1 2 3"
    fi
fi

if [[ ! -x "${ASSIGN_BIN}" ]]; then
    echo "ERROR: assignBarcodes not found or not executable: ${ASSIGN_BIN}" >&2
    echo "Run: make feature-barcodes-tools" >&2
    exit 1
fi
if [[ ! -f "${MEX_STUB}" ]]; then
    echo "ERROR: MEX stub not found: ${MEX_STUB}" >&2
    exit 1
fi
if [[ ! -d "${FASTQ_DIR}" ]]; then
    echo "ERROR: FASTQ dir not found: ${FASTQ_DIR}" >&2
    exit 1
fi
if [[ ! -f "${WHITELIST}" ]]; then
    echo "ERROR: whitelist not found: ${WHITELIST}" >&2
    exit 1
fi
if [[ ! -f "${FEATURE_REF}" ]]; then
    echo "ERROR: feature reference not found: ${FEATURE_REF}" >&2
    exit 1
fi

mkdir -p "${OUT_ROOT}"
SUMMARY="${OUT_ROOT}/summary.tsv"
printf "run_id\tmax_hamming\tfeature_n\tbarcode_n\tmax_barcode_mismatches\tfeature_offset\tlimit_search\tbarcodes\tfeatures\ttotal_counts\twhitelisted_barcodes\tunmatched_reads\tpercent_assigned\tcr_barcodes\tcr_features\tcr_total_counts\tcounts_delta\tcounts_ratio\tassign_out\n" > "${SUMMARY}"

FILTERED_BARCODES_CLEAN=""
if [[ -n "${FILTERED_BARCODES}" ]]; then
    if [[ ! -f "${FILTERED_BARCODES}" ]]; then
        echo "ERROR: filtered barcodes file not found: ${FILTERED_BARCODES}" >&2
        exit 1
    fi
    FILTERED_BARCODES_CLEAN="${OUT_ROOT}/filtered_barcodes.clean.txt"
    if [[ "${FILTERED_BARCODES}" == *.gz ]]; then
        gzip -dc "${FILTERED_BARCODES}" | awk -F',' '{print $NF}' | tr -d '\r' | sed 's/-1$//' | awk 'length($0) > 0' > "${FILTERED_BARCODES_CLEAN}"
    else
        awk -F',' '{print $NF}' "${FILTERED_BARCODES}" | tr -d '\r' | sed 's/-1$//' | awk 'length($0) > 0' > "${FILTERED_BARCODES_CLEAN}"
    fi
fi

shopt -s nullglob
barcode_files=("${FASTQ_DIR}"/*_R1_*.fastq.gz)
forward_files=("${FASTQ_DIR}"/*_R2_*.fastq.gz)
shopt -u nullglob

if [[ ${#barcode_files[@]} -eq 0 ]]; then
    echo "ERROR: no R1 FASTQs found in ${FASTQ_DIR}" >&2
    exit 1
fi
if [[ ${#forward_files[@]} -eq 0 ]]; then
    echo "ERROR: no R2 FASTQs found in ${FASTQ_DIR}" >&2
    exit 1
fi

tagify_limit_search() {
    local value="$1"
    if [[ "${value}" == -* ]]; then
        echo "n${value#-}"
    else
        echo "${value}"
    fi
}

find_assign_out() {
    local root="$1"
    local matrix_path
    matrix_path=$(find "${root}" -maxdepth 4 -type f -name matrix.mtx | head -n 1 || true)
    if [[ -z "${matrix_path}" ]]; then
        echo ""
        return
    fi
    local dir
    dir="$(dirname "${matrix_path}")"
    if [[ "$(basename "${dir}")" == "filtered" ]]; then
        dir="$(dirname "${dir}")"
    fi
    echo "${dir}"
}

count_matrix_sum() {
    local matrix_path="$1"
    awk 'BEGIN {sum=0; header=0} /^%/ {next} {if (header==0) {header=1; next} sum+=$3} END {printf "%.0f", sum}' "${matrix_path}"
}

compute_cr_stats() {
    local matrix_path="$1"
    local features_path="$2"
    local barcodes_path="$3"
    local feature_type="$4"
    python3 - "$matrix_path" "$features_path" "$barcodes_path" "$feature_type" <<'PY'
import gzip
import sys

matrix_path, features_path, barcodes_path, feature_type = sys.argv[1:5]

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

with open_maybe_gz(barcodes_path) as f:
    n_barcodes = sum(1 for _ in f if _.strip())

feature_rows = []
with open_maybe_gz(features_path) as f:
    for idx, line in enumerate(f, start=1):
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            parts = line.split(",")
        if len(parts) >= 3 and parts[2] == feature_type:
            feature_rows.append(idx)

feature_rows_set = set(feature_rows)
total_counts = 0

with open_maybe_gz(matrix_path) as f:
    for line in f:
        if line.startswith("%"):
            continue
        parts = line.strip().split()
        if len(parts) == 3:
            row = int(parts[0])
            if row in feature_rows_set:
                total_counts += int(parts[2])

print(f"{n_barcodes}\t{len(feature_rows)}\t{total_counts}")
PY
}

CR_BARCODES="NA"
CR_FEATURES="NA"
CR_TOTAL_COUNTS="NA"
if [[ -n "${CR_MATRIX_DIR}" ]]; then
    cr_matrix=""
    cr_features=""
    cr_barcodes=""
    if [[ -f "${CR_MATRIX_DIR}/matrix.mtx.gz" ]]; then
        cr_matrix="${CR_MATRIX_DIR}/matrix.mtx.gz"
    elif [[ -f "${CR_MATRIX_DIR}/matrix.mtx" ]]; then
        cr_matrix="${CR_MATRIX_DIR}/matrix.mtx"
    fi
    if [[ -f "${CR_MATRIX_DIR}/features.tsv.gz" ]]; then
        cr_features="${CR_MATRIX_DIR}/features.tsv.gz"
    elif [[ -f "${CR_MATRIX_DIR}/features.tsv" ]]; then
        cr_features="${CR_MATRIX_DIR}/features.tsv"
    fi
    if [[ -f "${CR_MATRIX_DIR}/barcodes.tsv.gz" ]]; then
        cr_barcodes="${CR_MATRIX_DIR}/barcodes.tsv.gz"
    elif [[ -f "${CR_MATRIX_DIR}/barcodes.tsv" ]]; then
        cr_barcodes="${CR_MATRIX_DIR}/barcodes.tsv"
    fi
    if [[ -n "${cr_matrix}" && -n "${cr_features}" && -n "${cr_barcodes}" ]]; then
        read -r CR_BARCODES CR_FEATURES CR_TOTAL_COUNTS < <(compute_cr_stats "${cr_matrix}" "${cr_features}" "${cr_barcodes}" "${CR_FEATURE_TYPE}")
    else
        echo "WARNING: CR matrix files not found under ${CR_MATRIX_DIR}; skipping CR comparison." >&2
    fi
fi

for max_hamming in ${MAX_HAMMING_LIST}; do
    for feature_n in ${FEATURE_N_LIST}; do
        for barcode_n in ${BARCODE_N_LIST}; do
            for max_barcode_mismatches in ${MAX_BARCODE_MISMATCHES_LIST}; do
                for feature_offset in ${FEATURE_OFFSET_LIST}; do
                    for limit_search in ${LIMIT_SEARCH_LIST}; do
                        limit_tag="$(tagify_limit_search "${limit_search}")"
                        run_id="m${max_hamming}_fn${feature_n}_bn${barcode_n}_mbm${max_barcode_mismatches}_off${feature_offset}_ls${limit_tag}"
                        run_out="${OUT_ROOT}/${run_id}"
                        rm -rf "${run_out}"
                        mkdir -p "${run_out}"

                        log_file="${run_out}/assignBarcodes.log"
                        cmd_file="${run_out}/command.txt"

                        {
                            echo "${ASSIGN_BIN} \\"
                            echo "  --whitelist ${WHITELIST} \\"
                            echo "  --featurelist ${FEATURE_REF} \\"
                            echo "  --directory ${run_out} \\"
                            echo "  --barcode_fastq_pattern _R1_ \\"
                            echo "  --forward_fastq_pattern _R2_ \\"
                            echo "  --barcode_length ${BARCODE_LEN} \\"
                            echo "  --umi_length ${UMI_LEN} \\"
                            echo "  --maxHammingDistance ${max_hamming} \\"
                            echo "  --feature_n ${feature_n} \\"
                            echo "  --barcode_n ${barcode_n} \\"
                            echo "  --max_barcode_mismatches ${max_barcode_mismatches} \\"
                            echo "  --feature_constant_offset ${feature_offset} \\"
                            echo "  --limit_search ${limit_search} \\"
                            echo "  --stringency ${STRINGENCY} \\"
                            echo "  --min_counts ${MIN_COUNTS} \\"
                            if [[ -n "${MIN_POSTERIOR}" ]]; then
                                echo "  --min_posterior ${MIN_POSTERIOR} \\"
                            fi
                            if [[ -n "${FILTERED_BARCODES_CLEAN}" ]]; then
                                echo "  --filtered_barcodes ${FILTERED_BARCODES_CLEAN} \\"
                            fi
                            echo "  --threads ${THREADS} \\"
                            echo "  --search_threads ${SEARCH_THREADS} \\"
                            echo "  --consumer_threads_per_set ${CONSUMER_THREADS} \\"
                            echo "  ${FASTQ_DIR}"
                        } > "${cmd_file}"

                        echo "Running ${run_id}..."
                        "${ASSIGN_BIN}" \
                            --whitelist "${WHITELIST}" \
                            --featurelist "${FEATURE_REF}" \
                            --directory "${run_out}" \
                            --barcode_fastq_pattern "_R1_" \
                            --forward_fastq_pattern "_R2_" \
                            --barcode_length "${BARCODE_LEN}" \
                            --umi_length "${UMI_LEN}" \
                            --maxHammingDistance "${max_hamming}" \
                            --feature_n "${feature_n}" \
                            --barcode_n "${barcode_n}" \
                            --max_barcode_mismatches "${max_barcode_mismatches}" \
                            --feature_constant_offset "${feature_offset}" \
                            --limit_search "${limit_search}" \
                            --stringency "${STRINGENCY}" \
                            --min_counts "${MIN_COUNTS}" \
                            ${MIN_POSTERIOR:+--min_posterior "${MIN_POSTERIOR}"} \
                            ${FILTERED_BARCODES_CLEAN:+--filtered_barcodes "${FILTERED_BARCODES_CLEAN}"} \
                            --threads "${THREADS}" \
                            --search_threads "${SEARCH_THREADS}" \
                            --consumer_threads_per_set "${CONSUMER_THREADS}" \
                            "${FASTQ_DIR}" \
                            > "${log_file}" 2>&1

                        assign_out="$(find_assign_out "${run_out}")"
                        if [[ -z "${assign_out}" ]]; then
                            echo "WARNING: no matrix.mtx found for ${run_id}" >&2
                            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\t%s\t%s\t%s\tNA\tNA\tNA\n" \
                                "${run_id}" "${max_hamming}" "${feature_n}" "${barcode_n}" "${max_barcode_mismatches}" \
                                "${feature_offset}" "${limit_search}" "${CR_BARCODES}" "${CR_FEATURES}" "${CR_TOTAL_COUNTS}" >> "${SUMMARY}"
                            continue
                        fi

                        python3 "${MEX_STUB}" \
                            --assign-out "${assign_out}" \
                            --feature-csv "${FEATURE_REF}" \
                            --force \
                            >> "${log_file}" 2>&1 || true

                        matrix_path="$(find "${assign_out}" -maxdepth 2 -type f -name matrix.mtx | head -n 1 || true)"
                        if [[ -z "${matrix_path}" ]]; then
                            matrix_path="$(find "${assign_out}/filtered" -maxdepth 2 -type f -name matrix.mtx | head -n 1 || true)"
                        fi
                        if [[ -z "${matrix_path}" ]]; then
                            echo "WARNING: no matrix.mtx found for ${run_id} after MEX stub" >&2
                            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\t%s\t%s\t%s\tNA\tNA\t%s\n" \
                                "${run_id}" "${max_hamming}" "${feature_n}" "${barcode_n}" "${max_barcode_mismatches}" \
                                "${feature_offset}" "${limit_search}" "${CR_BARCODES}" "${CR_FEATURES}" "${CR_TOTAL_COUNTS}" "${assign_out}" >> "${SUMMARY}"
                            continue
                        fi

                        matrix_dir="$(dirname "${matrix_path}")"
                        barcodes_path="${matrix_dir}/barcodes.tsv"
                        features_path="${matrix_dir}/features.tsv"
                        if [[ ! -f "${barcodes_path}" ]]; then
                            barcodes_path="${matrix_dir}/barcodes.txt"
                        fi
                        if [[ ! -f "${features_path}" ]]; then
                            features_path="${matrix_dir}/features.txt"
                        fi

                        barcodes_count="NA"
                        features_count="NA"
                        if [[ -f "${barcodes_path}" ]]; then
                            barcodes_count="$(wc -l < "${barcodes_path}" | tr -d ' ')"
                        fi
                        if [[ -f "${features_path}" ]]; then
                            features_count="$(wc -l < "${features_path}" | tr -d ' ')"
                        fi

                        stats_path="$(find "${assign_out}" -maxdepth 2 -type f -name stats.txt | head -n 1 || true)"
                        whitelisted_barcodes="NA"
                        unmatched_reads="NA"
                        percent_assigned="NA"
                        if [[ -f "${stats_path}" ]]; then
                            whitelisted_barcodes="$(awk '/^Total whitelisted barcodes/ {print $4}' "${stats_path}")"
                            unmatched_reads="$(awk '/^Total_unmatched_reads/ {print $2}' "${stats_path}")"
                            percent_assigned="$(awk '/^Percentage reads assigned to barcode/ {print $6}' "${stats_path}")"
                        fi

                        total_counts="$(count_matrix_sum "${matrix_path}")"
                        counts_delta="NA"
                        counts_ratio="NA"
                        if [[ "${CR_TOTAL_COUNTS}" != "NA" ]]; then
                            counts_delta="$((total_counts-CR_TOTAL_COUNTS))"
                            if [[ "${CR_TOTAL_COUNTS}" -gt 0 ]]; then
                                counts_ratio="$(python3 - <<PY
cr=${CR_TOTAL_COUNTS}
star=${total_counts}
print(f"{star/cr:.6f}")
PY
)"
                            fi
                        fi

                        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                            "${run_id}" "${max_hamming}" "${feature_n}" "${barcode_n}" "${max_barcode_mismatches}" \
                            "${feature_offset}" "${limit_search}" "${barcodes_count}" "${features_count}" "${total_counts}" \
                            "${whitelisted_barcodes}" "${unmatched_reads}" "${percent_assigned}" \
                            "${CR_BARCODES}" "${CR_FEATURES}" "${CR_TOTAL_COUNTS}" "${counts_delta}" "${counts_ratio}" "${assign_out}" >> "${SUMMARY}"
                    done
                done
            done
        done
    done
done

echo "Sweep complete."
echo "Summary: ${SUMMARY}"

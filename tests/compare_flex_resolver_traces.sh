#!/usr/bin/env bash

set -euo pipefail

if [[ "$#" -ne 2 ]]; then
    echo "usage: $0 SPATIAL_TRACE ORDINARY_TRACE" >&2
    exit 2
fi

spatial_trace="$1"
ordinary_trace="$2"

for trace in "${spatial_trace}" "${ordinary_trace}"; do
    if [[ ! -s "${trace}" ]]; then
        echo "FAIL: missing or empty resolver trace: ${trace}" >&2
        exit 2
    fi
done

work_dir="$(mktemp -d /tmp/star-flex-resolver-parity.XXXXXX)"
trap 'rm -rf "${work_dir}"' EXIT

normalize_terminal_rows() {
    local trace="$1"
    local output="$2"
    awk -F '\t' '
        BEGIN { OFS="\t" }
        NR == 1 { next }
        $9 == "KEEP_HASH" || $9 == "RESOLVER_DROP" {
            if (++seen[$2] > 1) {
                printf "duplicate terminal resolver row for iRead=%s in %s\n",
                       $2, FILENAME > "/dev/stderr"
                duplicate=1
            }
            print $2, $7, $8, $9, $10
        }
        END { if (duplicate) exit 3 }
    ' "${trace}" | LC_ALL=C sort -t $'\t' -k1,1 > "${output}"
}

normalize_terminal_rows "${spatial_trace}" "${work_dir}/spatial.tsv"
normalize_terminal_rows "${ordinary_trace}" "${work_dir}/ordinary.tsv"

spatial_count="$(wc -l < "${work_dir}/spatial.tsv")"
ordinary_count="$(wc -l < "${work_dir}/ordinary.tsv")"

if [[ "${spatial_count}" -eq 0 ]]; then
    echo "FAIL: spatial trace has no terminal resolver rows" >&2
    exit 1
fi
if [[ "${ordinary_count}" -eq 0 ]]; then
    echo "FAIL: ordinary trace has no terminal resolver rows" >&2
    exit 1
fi

cut -f1 "${work_dir}/spatial.tsv" > "${work_dir}/spatial.ids"
cut -f1 "${work_dir}/ordinary.tsv" > "${work_dir}/ordinary.ids"
comm -23 "${work_dir}/ordinary.ids" "${work_dir}/spatial.ids" \
    > "${work_dir}/ordinary_only.ids"
if [[ -s "${work_dir}/ordinary_only.ids" ]]; then
    echo "FAIL: ordinary resolver invoked for reads absent from the spatial resolver trace" >&2
    head -20 "${work_dir}/ordinary_only.ids" >&2
    exit 1
fi

join -t $'\t' -j 1 "${work_dir}/spatial.tsv" "${work_dir}/ordinary.tsv" \
    > "${work_dir}/joined.tsv"
awk -F '\t' '
    BEGIN { OFS="\t" }
    $2 != $6 || $3 != $7 || $4 != $8 || $5 != $9 { print }
' "${work_dir}/joined.tsv" > "${work_dir}/mismatches.tsv"

mismatch_count="$(wc -l < "${work_dir}/mismatches.tsv")"
if [[ "${mismatch_count}" -ne 0 ]]; then
    echo "FAIL: ${mismatch_count} per-read resolver decisions differ" >&2
    echo "iRead spatial_isProbe spatial_gene spatial_reason spatial_extra ordinary_isProbe ordinary_gene ordinary_reason ordinary_extra" >&2
    head -20 "${work_dir}/mismatches.tsv" >&2
    exit 1
fi

echo "spatial_terminal_reads=${spatial_count}"
echo "ordinary_terminal_reads=${ordinary_count}"
echo "ordinary_resolver_coverage_fraction=$(awk -v n="${ordinary_count}" -v d="${spatial_count}" 'BEGIN { printf "%.9f", n/d }')"
echo "resolver_decision_mismatches=0"
echo "PASS: ordinary inline Flex and native spatial Flex use identical terminal resolver decisions on every shared hash-miss read"

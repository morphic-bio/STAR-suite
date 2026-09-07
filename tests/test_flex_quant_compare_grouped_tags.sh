#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp_root=$(mktemp -d "${TMPDIR:-/tmp}/flex_quant_compare_tags.XXXXXX")
trap 'rm -rf -- "$tmp_root"' EXIT

tool="$tmp_root/flex_quant_compare"
g++ -O3 -DNDEBUG -std=c++17 -fopenmp \
    "$repo_root/docs/benchmarks/jax_matrix_20260904/flex_quant_compare.cpp" \
    -lz -o "$tool"

write_single_cell_mex() {
    local directory=$1
    local barcode=$2
    local count=$3
    mkdir -p "$directory"
    printf 'ENSG1\tGENE1\tGene Expression\n' >"$directory/features.tsv"
    printf '%s\n' "$barcode" >"$directory/barcodes.tsv"
    printf '%s\n' \
        '%%MatrixMarket matrix coordinate integer general' \
        '%' \
        '1 1 1' \
        "1 1 $count" >"$directory/matrix.mtx"
}

cell_barcode=AAAAAAAAAAAAAAAA
write_single_cell_mex \
    "$tmp_root/query/BC001/Gene/filtered" "$cell_barcode" 30
write_single_cell_mex \
    "$tmp_root/query/BC002/Gene/filtered" "$cell_barcode" 40

cr_mex="$tmp_root/cr/SAMPLE/count/sample_filtered_feature_bc_matrix"
mkdir -p "$cr_mex"
printf 'ENSG1\tGENE1\tGene Expression\n' >"$cr_mex/features.tsv"
printf '%s\n' \
    "${cell_barcode}ACTTTAGG-1" \
    "${cell_barcode}AACGGGAA-1" >"$cr_mex/barcodes.tsv"
printf '%s\n' \
    '%%MatrixMarket matrix coordinate integer general' \
    '%' \
    '1 2 2' \
    '1 1 30' \
    '1 2 40' >"$cr_mex/matrix.mtx"

printf '%s\n' \
    '[samples]' \
    'sample_id,probe_barcode_ids' \
    'SAMPLE,BC001|BC002' >"$tmp_root/cellranger.config.csv"
printf 'BC001\tACTTTAGG\nBC002\tAACGGGAA\n' >"$tmp_root/tag_map.tsv"

mkdir -p "$tmp_root/star_run/per_sample" "$tmp_root/star_run/Solo.out/Gene/raw"
ln -s "$tmp_root/query/BC001" "$tmp_root/star_run/per_sample/BC001"
ln -s "$tmp_root/query/BC002" "$tmp_root/star_run/per_sample/BC002"
printf '%s\n' \
    "${cell_barcode}ACTTTAGG" \
    "${cell_barcode}AACGGGAA" \
    >"$tmp_root/star_run/Solo.out/Gene/raw/barcodes.tsv"
printf '%s\n' \
    '%%MatrixMarket matrix coordinate integer general' \
    '%' \
    '1 2 2' \
    '1 1 30' \
    '1 2 40' >"$tmp_root/star_run/Solo.out/Gene/raw/matrix.mtx"

out_prefix="$tmp_root/report"
"$tool" \
    --query-root "$tmp_root/query" \
    --input-kind star \
    --label grouped-tag-test \
    --cr-root "$tmp_root/cr" \
    --cr-config "$tmp_root/cellranger.config.csv" \
    --tag-map "$tmp_root/tag_map.tsv" \
    --star-run "$tmp_root/star_run" \
    --sample-whitelist "$tmp_root/tag_map.tsv" \
    --out-prefix "$out_prefix" \
    --threads 1

raw_row=$(awk '$1 == "SAMPLE" {print $1, $2, $3, $4}' \
    "$out_prefix.raw-concordance.txt")
if [[ "$raw_row" != 'SAMPLE BC001+BC002 2/2 1.0000' ]]; then
    echo "FAIL: grouped cells were not retained independently: $raw_row" >&2
    exit 1
fi

paper_row=$(awk '$1 == "SAMPLE" {print $1, $2, $3, $5}' \
    "$out_prefix.paper-protocol.txt")
if [[ "$paper_row" != 'SAMPLE BC001+BC002 2 1.0000' ]]; then
    echo "FAIL: common-cell count or Jaccard is wrong: $paper_row" >&2
    exit 1
fi

calling_row=$(awk '$3 == "STAR-Flex" {print $1, $2, $3, $4, $5, $6}' \
    "$out_prefix.cell-calling-pr.txt")
if [[ "$calling_row" != 'SAMPLE BC001+BC002 STAR-Flex 2 2 2' ]]; then
    echo "FAIL: grouped cell-calling identities were merged: $calling_row" >&2
    exit 1
fi

if "$tool" \
    --query-root "$tmp_root/query" \
    --input-kind star \
    --label missing-map \
    --cr-root "$tmp_root/cr" \
    --cr-config "$tmp_root/cellranger.config.csv" \
    --out-prefix "$tmp_root/missing-map" \
    --threads 1 >"$tmp_root/missing-map.stdout" 2>"$tmp_root/missing-map.stderr"; then
    echo 'FAIL: grouped comparison unexpectedly succeeded without --tag-map' >&2
    exit 1
fi
grep -q 'requires --tag-map' "$tmp_root/missing-map.stderr" || {
    echo 'FAIL: missing grouped-tag map did not produce a useful diagnostic' >&2
    cat "$tmp_root/missing-map.stderr" >&2
    exit 1
}

echo 'PASS: native comparator keys grouped cells by CB16|tag and requires a tag map'

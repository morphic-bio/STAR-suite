#!/usr/bin/env bash
set -euo pipefail

tool=${1:-./flex_mex_cell_stats}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

mkdir -p "$work/mex" "$work/mex_plain" "$work/out"

cat >"$work/mex/features.tsv" <<'EOF'
ENSG000001	GENE1	Gene Expression
ENSG000002	GENE2	Gene Expression
ENSG000003	GENE3	Gene Expression
EOF

# The first two cells deliberately have the same CB16 and different TAG8s.
cat >"$work/barcodes.tsv" <<'EOF'
AAAAAAAAAAAAAAAAACGTACGT-1
AAAAAAAAAAAAAAAATGCATGCA-1
CCCCCCCCCCCCCCCCACGTACGT-1
EOF
gzip -c "$work/barcodes.tsv" >"$work/mex/barcodes.tsv.gz"

cat >"$work/matrix.mtx" <<'EOF'
%%MatrixMarket matrix coordinate integer general
% deterministic tag-aware fixture
3 3 4
1 1 5
2 1 2
3 2 11
1 3 3
EOF
gzip -c "$work/matrix.mtx" >"$work/mex/matrix.mtx.gz"
cp "$work/mex/features.tsv" "$work/mex_plain/features.tsv"
cp "$work/barcodes.tsv" "$work/mex_plain/barcodes.tsv"
cp "$work/matrix.mtx" "$work/mex_plain/matrix.mtx"

cat >"$work/tag_a.txt" <<'EOF'
AAAAAAAAAAAAAAAAACGTACGT-1
CCCCCCCCCCCCCCCCACGTACGT-1
GGGGGGGGGGGGGGGGACGTACGT-1
EOF
cat >"$work/tag_b.txt" <<'EOF'
AAAAAAAAAAAAAAAATGCATGCA-1
EOF
cat >"$work/both.txt" <<'EOF'
AAAAAAAAAAAAAAAAACGTACGT-1
AAAAAAAAAAAAAAAATGCATGCA-1
EOF

"$tool" \
    --mex-dir "$work/mex_plain" \
    --cohort "tag_a=$work/tag_a.txt" \
    --cohort "tag_b=$work/tag_b.txt" \
    --cohort "both=$work/both.txt" \
    --input-label fixture \
    --out-prefix "$work/out/first"

"$tool" \
    --mex-dir "$work/mex" \
    --cohort "tag_a=$work/tag_a.txt" \
    --cohort "tag_b=$work/tag_b.txt" \
    --cohort "both=$work/both.txt" \
    --input-label fixture \
    --out-prefix "$work/out/second"

cmp "$work/out/first.cells.tsv" "$work/out/second.cells.tsv"
cmp "$work/out/first.summary.tsv" "$work/out/second.summary.tsv"

cat >"$work/expected.cells.tsv" <<'EOF'
input	cohort	barcode	status	mex_column	total_counts	detected_genes	nnz
fixture	tag_a	AAAAAAAAAAAAAAAAACGTACGT	found	1	7	2	2
fixture	tag_a	CCCCCCCCCCCCCCCCACGTACGT	found	3	3	1	1
fixture	tag_a	GGGGGGGGGGGGGGGGACGTACGT	missing	NA	NA	NA	NA
fixture	tag_b	AAAAAAAAAAAAAAAATGCATGCA	found	2	11	1	1
fixture	both	AAAAAAAAAAAAAAAAACGTACGT	found	1	7	2	2
fixture	both	AAAAAAAAAAAAAAAATGCATGCA	found	2	11	1	1
EOF
cmp "$work/expected.cells.tsv" "$work/out/first.cells.tsv"

cat >"$work/expected.summary.tsv" <<'EOF'
input	cohort	expected_cells	found_cells	missing_cells	cells_with_counts	total_counts	total_detected_genes	total_nnz	counts_min	counts_q25	counts_median	counts_q75	counts_p90	counts_p95	counts_p99	counts_max	genes_min	genes_q25	genes_median	genes_q75	genes_p90	genes_p95	genes_p99	genes_max
fixture	tag_a	3	2	1	2	10	3	3	3	3	3	7	7	7	7	7	1	1	1	2	2	2	2	2
fixture	tag_b	1	1	0	1	11	1	1	11	11	11	11	11	11	11	11	1	1	1	1	1	1	1	1
fixture	both	2	2	0	2	18	3	3	7	7	7	11	11	11	11	11	1	1	1	2	2	2	2	2
EOF
cmp "$work/expected.summary.tsv" "$work/out/first.summary.tsv"

printf 'flex_mex_cell_stats fixture: PASS\n'

#!/usr/bin/env bash
set -euo pipefail

resolver="${1:?resolver path required}"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

cat >"$tmp/input.tsv" <<'EOF'
read_id	feature_id	raw_umi	candidate	log_sequence_likelihood	exact_read_count
r1	gene	AAAA	A	0	8
r1	gene	AAAA	B	0	0
r2	gene	AAAA	A	0	8
EOF

"$resolver" --input "$tmp/input.tsv" --out-dir "$tmp/out"

{
  head -n 1 "$tmp/input.tsv"
  tail -n +2 "$tmp/input.tsv" | LC_ALL=C sort -t $'\t' -k2,2 -k1,1 -k4,4
} >"$tmp/input.feature-sorted.tsv"
"$resolver" --input "$tmp/input.feature-sorted.tsv" --input-feature-sorted \
  --out-dir "$tmp/out-streaming"
for table in read_cliques.tsv hard_call_audit.tsv strict_molecules.tsv \
  hard_molecules.tsv gated_hard_molecules.tsv soft_expected_molecules.tsv; do
  diff -u \
    <(tail -n +2 "$tmp/out/$table" | LC_ALL=C sort) \
    <(tail -n +2 "$tmp/out-streaming/$table" | LC_ALL=C sort)
done
grep -q $'^execution_mode\tfeature_sorted_streaming$' \
  "$tmp/out-streaming/resolved_config.tsv"

cat >"$tmp/parallel.tsv" <<'EOF'
read_id	feature_id	raw_umi	candidate	log_sequence_likelihood	exact_read_count
r1	geneA	AAAA	A	0	8
r2	geneB	AAAT	B	0	7
r3	geneC	AATT	C	0	6
EOF
"$resolver" --input "$tmp/parallel.tsv" --input-feature-sorted \
  --out-dir "$tmp/out-serial-features"
"$resolver" --input "$tmp/parallel.tsv" --input-feature-sorted --threads 3 \
  --out-dir "$tmp/out-parallel-features"
for table in read_cliques.tsv hard_call_audit.tsv strict_molecules.tsv \
  hard_molecules.tsv gated_hard_molecules.tsv soft_expected_molecules.tsv; do
  cmp "$tmp/out-serial-features/$table" "$tmp/out-parallel-features/$table"
done
grep -q $'^feature_threads\t3$' "$tmp/out-parallel-features/resolved_config.tsv"

if "$resolver" --input "$tmp/parallel.tsv" --threads 2 \
  --out-dir "$tmp/out-invalid-threads" 2>/dev/null; then
  echo "parallel global input mode unexpectedly accepted" >&2
  exit 1
fi

expected='gated_hard_molecules.tsv
hard_call_audit.tsv
hard_molecules.tsv
read_cliques.tsv
resolved_config.tsv
soft_expected_molecules.tsv
strict_molecules.tsv
summary.tsv'
actual="$(find "$tmp/out" -maxdepth 1 -type f -printf '%f\n' | LC_ALL=C sort)"
[[ "$actual" == "$expected" ]]
grep -q $'^1mm_cr.hard_count\t1$' "$tmp/out/summary.tsv"
grep -q $'^prior_application\tonce_per_read_clique$' "$tmp/out/resolved_config.tsv"

mkdir "$tmp/nonempty"
touch "$tmp/nonempty/existing"
if "$resolver" --input "$tmp/input.tsv" --out-dir "$tmp/nonempty" 2>/dev/null; then
  echo "nonempty output directory unexpectedly accepted" >&2
  exit 1
fi

cat >"$tmp/inconsistent.tsv" <<'EOF'
read_id	feature_id	raw_umi	candidate	log_sequence_likelihood	exact_read_count
r1	gene	AAAA	A	0	8
r2	gene	AAAT	A	0	9
EOF
if "$resolver" --input "$tmp/inconsistent.tsv" --out-dir "$tmp/inconsistent-out" 2>/dev/null; then
  echo "inconsistent exact-read prior unexpectedly accepted" >&2
  exit 1
fi

cat >"$tmp/nonfinite.tsv" <<'EOF'
read_id	feature_id	raw_umi	candidate	log_sequence_likelihood	exact_read_count
r1	gene	AAAA	A	nan	8
EOF
if "$resolver" --input "$tmp/nonfinite.tsv" --out-dir "$tmp/nonfinite-out" 2>/dev/null; then
  echo "non-finite likelihood unexpectedly accepted" >&2
  exit 1
fi

cat >"$tmp/out-of-order.tsv" <<'EOF'
read_id	feature_id	raw_umi	candidate	log_sequence_likelihood	exact_read_count
r1	geneB	AAAA	A	0	8
r2	geneA	AAAT	A	0	8
EOF
if "$resolver" --input "$tmp/out-of-order.tsv" --input-feature-sorted \
  --out-dir "$tmp/out-of-order-out" 2>/dev/null; then
  echo "out-of-order feature streaming input unexpectedly accepted" >&2
  exit 1
fi

echo "molecule-first CLI tests passed"

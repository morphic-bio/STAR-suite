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

echo "molecule-first CLI tests passed"

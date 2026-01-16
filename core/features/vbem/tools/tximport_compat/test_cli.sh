#!/bin/bash
# test_cli.sh - Smoke test for tximport_compat CLI

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLI="$SCRIPT_DIR/tximport_compat"
TEST_DIR="/tmp/tximport_cli_test_$$"

cleanup() {
    rm -rf "$TEST_DIR"
}
trap cleanup EXIT

mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

echo "=== Creating test fixtures ==="

# Create test quant.sf
cat > quant.sf <<'EOF'
Name	Length	EffectiveLength	TPM	NumReads
tx1	1000	800	100.000000	80.000
tx2	500	400	50.000000	20.000
tx3	2000	1800	200.000000	360.000
tx4	1500	1300	100.000000	130.000
EOF

# Create test tx2gene.tsv
cat > tx2gene.tsv <<'EOF'
tx1	gene1
tx2	gene1
tx3	gene2
tx4	gene2
EOF

echo "=== Test 1: Basic invocation ==="
$CLI --quant quant.sf --tx2gene tx2gene.tsv --output genes.sf --stats
cat genes.sf

echo ""
echo "=== Test 2: Sorted output ==="
$CLI --quant quant.sf --tx2gene tx2gene.tsv --output genes_sorted.sf --sort
cat genes_sorted.sf

echo ""
echo "=== Test 3: Verify gene count ==="
GENE_COUNT=$(tail -n +2 genes.sf | wc -l)
if [ "$GENE_COUNT" -eq 2 ]; then
    echo "PASS: Got expected 2 genes"
else
    echo "FAIL: Expected 2 genes, got $GENE_COUNT"
    exit 1
fi

echo ""
echo "=== Test 4: Verify total counts preserved ==="
# Total raw counts = 80 + 20 + 360 + 130 = 590
# Output counts should sum to ~590
TOTAL=$(awk 'NR>1 {sum += $5} END {printf "%.1f", sum}' genes.sf)
if [ "$TOTAL" = "590.0" ]; then
    echo "PASS: Total counts preserved ($TOTAL)"
else
    echo "FAIL: Expected total 590.0, got $TOTAL"
    exit 1
fi

echo ""
echo "=== Test 5: Help message ==="
$CLI --help | head -3

echo ""
echo "=== All CLI tests PASSED ==="



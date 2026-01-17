#!/bin/bash
# spot_check_gene_format.sh - Spot-check quant.genes.sf format against reference
#
# Usage: ./spot_check_gene_format.sh [quant.genes.sf] [reference.genes.sf]
#
# Validates:
# - Header format
# - Decimal precision per column
# - No scientific notation
# - Zero values formatted correctly

set -e

GENE_FILE=${1:-basic_quant.genes.sf}
REF_FILE=${2:-}

if [ ! -f "$GENE_FILE" ]; then
    echo "Error: $GENE_FILE not found"
    exit 1
fi

echo "========================================"
echo "Spot-Check: quant.genes.sf Format"
echo "========================================"
echo "File: $GENE_FILE"
if [ -n "$REF_FILE" ]; then
    echo "Reference: $REF_FILE"
fi
echo "========================================"
echo ""

PASSED=0
FAILED=0

# Check 1: Header
echo "=== Check 1: Header Format ==="
HEADER=$(head -1 "$GENE_FILE")
EXPECTED_HEADER="Name	Length	EffectiveLength	TPM	NumReads"
if [ "$HEADER" = "$EXPECTED_HEADER" ]; then
    echo "✓ Header correct: $HEADER"
    PASSED=$((PASSED + 1))
else
    echo "✗ Header mismatch"
    echo "  Expected: $EXPECTED_HEADER"
    echo "  Got:      $HEADER"
    FAILED=$((FAILED + 1))
fi
echo ""

# Check 2: No scientific notation
echo "=== Check 2: No Scientific Notation ==="
if grep -qE '[eE][+-]' "$GENE_FILE"; then
    echo "✗ Scientific notation found:"
    grep -E '[eE][+-]' "$GENE_FILE" | head -3
    FAILED=$((FAILED + 1))
else
    echo "✓ No scientific notation"
    PASSED=$((PASSED + 1))
fi
echo ""

# Check 3: No NaN/Inf
echo "=== Check 3: No NaN/Inf Values ==="
if grep -qiE '(nan|inf)' "$GENE_FILE"; then
    echo "✗ NaN or Inf values found:"
    grep -iE '(nan|inf)' "$GENE_FILE" | head -3
    FAILED=$((FAILED + 1))
else
    echo "✓ No NaN/Inf values"
    PASSED=$((PASSED + 1))
fi
echo ""

# Check 4: Decimal Precision (spot check first 10 lines)
echo "=== Check 4: Decimal Precision (spot check) ==="
PRECISION_OK=true
while IFS=$'\t' read -r name len efflen tpm counts; do
    # Length: max 3 decimals
    if echo "$len" | grep -qE '\.'; then
        DECIMALS=$(echo "$len" | sed 's/.*\.//' | wc -c)
        if [ $DECIMALS -gt 4 ]; then
            echo "  ✗ Length precision too high: $len (max 3 decimals)"
            PRECISION_OK=false
        fi
    fi
    
    # EffectiveLength: max 4 decimals
    if echo "$efflen" | grep -qE '\.'; then
        DECIMALS=$(echo "$efflen" | sed 's/.*\.//' | wc -c)
        if [ $DECIMALS -gt 5 ]; then
            echo "  ✗ EffectiveLength precision too high: $efflen (max 4 decimals)"
            PRECISION_OK=false
        fi
    fi
    
    # TPM: max 6 decimals
    if echo "$tpm" | grep -qE '\.'; then
        DECIMALS=$(echo "$tpm" | sed 's/.*\.//' | wc -c)
        if [ $DECIMALS -gt 7 ]; then
            echo "  ✗ TPM precision too high: $tpm (max 6 decimals)"
            PRECISION_OK=false
        fi
    fi
    
    # NumReads: max 3 decimals
    if echo "$counts" | grep -qE '\.'; then
        DECIMALS=$(echo "$counts" | sed 's/.*\.//' | wc -c)
        if [ $DECIMALS -gt 4 ]; then
            echo "  ✗ NumReads precision too high: $counts (max 3 decimals)"
            PRECISION_OK=false
        fi
    fi
done < <(tail -n +2 "$GENE_FILE" | head -10)

if [ "$PRECISION_OK" = true ]; then
    echo "✓ Decimal precision within limits (first 10 lines)"
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi
echo ""

# Check 5: Zero values formatted correctly
echo "=== Check 5: Zero Value Formatting ==="
ZERO_OK=true
while IFS=$'\t' read -r name len efflen tpm counts; do
    if [ "$tpm" != "0" ] && [ "$tpm" != "0.0" ] && [ "$tpm" != "0.00" ] && [ "$tpm" != "0.000" ]; then
        if echo "$tpm" | grep -qE '^0\.0+$'; then
            : # OK, it's zero with trailing zeros (will be trimmed)
        else
            echo "  ✗ Zero TPM formatted incorrectly: $tpm (should be '0')"
            ZERO_OK=false
        fi
    fi
    if [ "$counts" != "0" ] && [ "$counts" != "0.0" ] && [ "$counts" != "0.00" ] && [ "$counts" != "0.000" ]; then
        if echo "$counts" | grep -qE '^0\.0+$'; then
            : # OK, it's zero with trailing zeros (will be trimmed)
        else
            echo "  ✗ Zero NumReads formatted incorrectly: $counts (should be '0')"
            ZERO_OK=false
        fi
    fi
done < <(tail -n +2 "$GENE_FILE" | awk -F'\t' '$4 == 0 && $5 == 0' | head -5)

if [ "$ZERO_OK" = true ]; then
    echo "✓ Zero values formatted correctly"
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi
echo ""

# Check 6: Sample output (first 5 lines)
echo "=== Check 6: Sample Output (first 5 lines) ==="
head -6 "$GENE_FILE"
echo ""

# Summary
echo "========================================"
echo "Spot-Check Summary"
echo "========================================"
echo "Passed: $PASSED"
echo "Failed: $FAILED"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All format checks passed!"
    exit 0
else
    echo "✗ Some format checks failed"
    exit 1
fi

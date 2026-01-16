#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURES="$SCRIPT_DIR/fixtures/cellranger_format"
PERL_SCRIPT="/mnt/pikachu/Bulk-RNA-seq/star_salmon_aws/widgets/star_salmon_aws/10x_format_fa_gtf/Dockerfiles/format-fa-gtf.pl"
CPP_TOOL="$SCRIPT_DIR/test_cellranger_format"

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

if [ ! -f "$PERL_SCRIPT" ]; then
    echo "ERROR: Perl script not found: $PERL_SCRIPT"
    exit 1
fi

if [ ! -x "$CPP_TOOL" ]; then
    echo "ERROR: C++ tool not found or not executable: $CPP_TOOL"
    echo "Build with: cd source && make test_cellranger_format"
    exit 1
fi

echo "=============================================="
echo "CellRangerFormatter Parity Test"
echo "=============================================="
echo "Perl script: $PERL_SCRIPT"
echo "C++ tool: $CPP_TOOL"
echo "Fixtures: $FIXTURES"
echo ""

# 1. Run Perl script
echo ">>> Running Perl script..."
INPUTFA="$FIXTURES/test_genome.fa" OUTPUTFA="$TMPDIR/perl_genome.fa" \
INPUTGTF="$FIXTURES/test_genes.gtf" OUTPUTGTF="$TMPDIR/perl_genes.gtf" \
overwrite=1 \
perl "$PERL_SCRIPT" 2>&1 | grep -v "^working on" || true

# 2. Run C++ formatter
echo ">>> Running C++ formatter..."
"$CPP_TOOL" "$FIXTURES/test_genome.fa" "$FIXTURES/test_genes.gtf" \
    "$TMPDIR/cpp_genome.fa" "$TMPDIR/cpp_genes.gtf"

# 3. Compare outputs
echo ""
echo ">>> Comparing outputs..."

PASS_COUNT=0
FAIL_COUNT=0

# Compare FASTA
if diff -q "$TMPDIR/perl_genome.fa" "$TMPDIR/cpp_genome.fa" > /dev/null 2>&1; then
    echo "    ✓ FASTA: IDENTICAL"
    PASS_COUNT=$((PASS_COUNT + 1))
else
    echo "    ✗ FASTA: DIFFERS"
    echo "      First diff:"
    diff "$TMPDIR/perl_genome.fa" "$TMPDIR/cpp_genome.fa" | head -10 || true
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# Compare GTF
if diff -q "$TMPDIR/perl_genes.gtf" "$TMPDIR/cpp_genes.gtf" > /dev/null 2>&1; then
    echo "    ✓ GTF: IDENTICAL"
    PASS_COUNT=$((PASS_COUNT + 1))
else
    echo "    ✗ GTF: DIFFERS"
    echo "      First diff:"
    diff "$TMPDIR/perl_genes.gtf" "$TMPDIR/cpp_genes.gtf" | head -10 || true
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# Summary
echo ""
echo "=============================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=============================================="

if [ "$FAIL_COUNT" -eq 0 ]; then
    echo "✓ PARITY TEST PASSED"
    exit 0
else
    echo "✗ PARITY TEST FAILED"
    exit 1
fi


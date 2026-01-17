#!/bin/bash
# tximport_parity_test.sh - Verify parity between tximport_compat and R tximport
#
# Usage: ./tximport_parity_test.sh [--skip-synthetic] [--skip-transcriptvb]
#
# Runs:
# 1. Synthetic fixture parity (small, deterministic)
# 2. TranscriptVB golden dataset parity

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_FLEX="$(cd "$SCRIPT_DIR/../.." && pwd)"
CLI="$STAR_FLEX/tools/tximport_compat/tximport_compat"
RENV_DIR="$STAR_FLEX/tools/tximport_compat"
REPORT_DIR="${REPORT_DIR:-/storage/production/tximport_parity_$(date +%Y%m%d_%H%M%S)}"

SKIP_SYNTHETIC=false
SKIP_TRANSCRIPTVB=false

for arg in "$@"; do
    case $arg in
        --skip-synthetic) SKIP_SYNTHETIC=true ;;
        --skip-transcriptvb) SKIP_TRANSCRIPTVB=true ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# Create report directory
mkdir -p "$REPORT_DIR"
cd "$REPORT_DIR"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
error_exit() { echo "[ERROR] $*" >&2; exit 1; }

# Helper: Compare two quant files numerically (values + ordering)
compare_quant_files() {
    local expected="$1"
    local actual="$2"
    local tolerance="${3:-1e-6}"
    
    python3 << PYEOF
import sys

def parse_sf_ordered(path):
    """Parse quant.sf preserving row order"""
    order = []
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            name = fields[0]
            order.append(name)
            data[name] = {h: float(fields[i]) if i > 0 else fields[i] 
                         for i, h in enumerate(header)}
    return order, data

exp_order, expected = parse_sf_ordered("$expected")
act_order, actual = parse_sf_ordered("$actual")
tolerance = $tolerance

# Check same genes
if set(expected.keys()) != set(actual.keys()):
    missing = set(expected.keys()) - set(actual.keys())
    extra = set(actual.keys()) - set(expected.keys())
    if missing:
        print(f"Missing genes: {list(missing)[:5]}...")
    if extra:
        print(f"Extra genes: {list(extra)[:5]}...")
    sys.exit(1)

# Check gene ordering matches
if exp_order != act_order:
    # Find first mismatch
    for i, (e, a) in enumerate(zip(exp_order, act_order)):
        if e != a:
            print(f"FAIL: Gene order mismatch at position {i}")
            print(f"  Expected: {e}")
            print(f"  Actual: {a}")
            sys.exit(1)
    # Length mismatch
    print(f"FAIL: Gene count mismatch: expected {len(exp_order)}, actual {len(act_order)}")
    sys.exit(1)

# Compare values including Length column
max_rel_diff = 0
max_diff_gene = None
max_diff_col = None

for gene in expected:
    for col in ['Length', 'EffectiveLength', 'TPM', 'NumReads']:
        exp_val = expected[gene][col]
        act_val = actual[gene][col]
        
        if abs(exp_val) < 1e-10 and abs(act_val) < 1e-10:
            continue  # both essentially zero
        
        if abs(exp_val) < 1e-10:
            rel_diff = abs(act_val)
        else:
            rel_diff = abs(exp_val - act_val) / max(abs(exp_val), 1e-10)
        
        if rel_diff > max_rel_diff:
            max_rel_diff = rel_diff
            max_diff_gene = gene
            max_diff_col = col

if max_rel_diff > tolerance:
    print(f"FAIL: Max relative difference {max_rel_diff:.2e} exceeds tolerance {tolerance}")
    print(f"  Gene: {max_diff_gene}, Column: {max_diff_col}")
    print(f"  Expected: {expected[max_diff_gene][max_diff_col]}")
    print(f"  Actual: {actual[max_diff_gene][max_diff_col]}")
    sys.exit(1)
else:
    print(f"PASS: Max relative difference {max_rel_diff:.2e} within tolerance {tolerance}")
    print(f"  Genes compared: {len(expected)}, ordering verified")
    sys.exit(0)
PYEOF
}

# ============================================================================
# Test 1: Synthetic fixture parity
# ============================================================================
if [ "$SKIP_SYNTHETIC" = false ]; then
    log "=== Test 1: Synthetic Fixture Parity ==="
    
    SYNTH_DIR="$REPORT_DIR/synthetic"
    mkdir -p "$SYNTH_DIR"
    cd "$SYNTH_DIR"
    
    # Copy fixtures
    cp "$SCRIPT_DIR/fixtures/synthetic_quant.sf" quant.sf
    cp "$SCRIPT_DIR/fixtures/synthetic_tx2gene.tsv" tx2gene.tsv
    
    log "Running R tximport..."
    cd "$RENV_DIR"
    Rscript - "$SYNTH_DIR" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
synth_dir <- args[1]

library(tximport)

# Read inputs
quant_path <- file.path(synth_dir, "quant.sf")
tx2gene <- read.delim(file.path(synth_dir, "tx2gene.tsv"), header=FALSE, 
                      col.names=c("tx","gene"), stringsAsFactors=FALSE)

# Run tximport
txi <- tximport(quant_path, type="salmon", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM", dropInfReps=TRUE)

# Write output
output <- data.frame(
    Name = rownames(txi$counts),
    Length = round(txi$length[,1], 3),
    EffectiveLength = round(txi$length[,1], 3),
    TPM = txi$abundance[,1],
    NumReads = round(txi$counts[,1], 3)
)
write.table(output, file.path(synth_dir, "expected.sf"), 
            sep="\t", quote=FALSE, row.names=FALSE)

cat("R tximport completed\n")
cat("Genes:", nrow(output), "\n")
RSCRIPT
    
    cd "$SYNTH_DIR"
    log "Running tximport_compat CLI..."
    $CLI --quant quant.sf --tx2gene tx2gene.tsv --output actual.sf --stats
    
    log "Comparing outputs..."
    compare_quant_files expected.sf actual.sf 1e-6
    
    log "Synthetic test PASSED"
fi

# ============================================================================
# Test 2: TranscriptVB golden dataset parity
# ============================================================================
if [ "$SKIP_TRANSCRIPTVB" = false ]; then
    log "=== Test 2: TranscriptVB Golden Dataset Parity ==="
    
    TXVB_DIR="$REPORT_DIR/transcriptvb"
    mkdir -p "$TXVB_DIR"
    cd "$TXVB_DIR"
    
    # Copy inputs
    cp "$STAR_FLEX/tests/transcriptvb/golden/vb_quant.sf" quant.sf
    cp "$STAR_FLEX/tests/transcriptvb/tx2gene.tsv" tx2gene.tsv
    
    log "Running R tximport..."
    cd "$RENV_DIR"
    Rscript - "$TXVB_DIR" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
txvb_dir <- args[1]

library(tximport)

# Read inputs
quant_path <- file.path(txvb_dir, "quant.sf")
tx2gene <- read.delim(file.path(txvb_dir, "tx2gene.tsv"), header=FALSE, 
                      col.names=c("tx","gene"), stringsAsFactors=FALSE)

# Run tximport
txi <- tximport(quant_path, type="salmon", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM", dropInfReps=TRUE)

# Write output
output <- data.frame(
    Name = rownames(txi$counts),
    Length = round(txi$length[,1], 3),
    EffectiveLength = round(txi$length[,1], 3),
    TPM = txi$abundance[,1],
    NumReads = round(txi$counts[,1], 3)
)
write.table(output, file.path(txvb_dir, "expected.sf"), 
            sep="\t", quote=FALSE, row.names=FALSE)

cat("R tximport completed\n")
cat("Genes:", nrow(output), "\n")
cat("Total counts:", sum(txi$counts[,1]), "\n")
RSCRIPT
    
    cd "$TXVB_DIR"
    log "Running tximport_compat CLI..."
    $CLI --quant quant.sf --tx2gene tx2gene.tsv --output actual.sf --stats
    
    log "Comparing outputs..."
    compare_quant_files expected.sf actual.sf 1e-6
    
    log "TranscriptVB test PASSED"
fi

# ============================================================================
# Summary
# ============================================================================
log "=== All Parity Tests PASSED ==="
log "Report directory: $REPORT_DIR"

# Create summary markdown
cat > "$REPORT_DIR/SUMMARY.md" << EOF
# tximport Parity Test Results

**Date**: $(date)
**tximport_compat**: $CLI
**R tximport version**: 1.38.2

## Tests Run

1. **Synthetic Fixture Parity**: $( [ "$SKIP_SYNTHETIC" = true ] && echo "SKIPPED" || echo "PASSED" )
2. **TranscriptVB Golden Dataset**: $( [ "$SKIP_TRANSCRIPTVB" = true ] && echo "SKIPPED" || echo "PASSED" )

## Verification

All tests compare numerical values with relative tolerance of 1e-6.

Output files preserved in subdirectories:
- \`synthetic/\`: Synthetic test fixtures and outputs
- \`transcriptvb/\`: TranscriptVB golden dataset test

EOF

log "Summary written to $REPORT_DIR/SUMMARY.md"


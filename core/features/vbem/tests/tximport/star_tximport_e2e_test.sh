#!/bin/bash
# star_tximport_e2e_test.sh - End-to-end test for STAR's tximport gene quantification
#
# Runs STAR with --quantVBgenesMode Tximport and compares output to R tximport

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_FLEX="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR_BIN="${STAR_BIN:-$STAR_FLEX/source/STAR}"
CLI="$STAR_FLEX/tools/tximport_compat/tximport_compat"
RENV_DIR="$STAR_FLEX/tools/tximport_compat"

# Use existing test genome and data from transcriptvb
GENOME_DIR="${GENOME_DIR:-/tmp/star_vb_test/star_new_index}"
READS_1="${READS_1:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz}"
READS_2="${READS_2:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz}"
GTF="${GTF:-/mnt/pikachu/test-datasets-rnaseq/reference/genes.gtf}"

OUTDIR="${OUTDIR:-/storage/production/tximport_e2e_$(date +%Y%m%d_%H%M%S)}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
error_exit() { echo "[ERROR] $*" >&2; exit 1; }

# Validate prerequisites
[ -x "$STAR_BIN" ] || error_exit "STAR binary not found: $STAR_BIN"
[ -d "$GENOME_DIR" ] || error_exit "Genome directory not found: $GENOME_DIR"
[ -f "$READS_1" ] || error_exit "Reads file not found: $READS_1"
[ -f "$GTF" ] || error_exit "GTF file not found: $GTF"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

log "=== STAR tximport End-to-End Test ==="
log "STAR: $STAR_BIN"
log "Genome: $GENOME_DIR"
log "Output: $OUTDIR"

# Step 1: Run STAR with TranscriptVB and tximport mode
log "Step 1: Running STAR TranscriptVB with tximport mode..."

$STAR_BIN --runMode alignReads \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READS_1" "$READS_2" \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --quantVBgenesMode Tximport \
    --outSAMtype None \
    --runThreadN 4 \
    --outFileNamePrefix star_

if [ ! -f star_quant.sf ]; then
    error_exit "STAR did not produce quant.sf"
fi

if [ ! -f star_quant.genes.tximport.sf ]; then
    error_exit "STAR did not produce quant.genes.tximport.sf"
fi

log "STAR outputs generated:"
ls -la star_quant*.sf

# Step 2: Generate tx2gene from GTF
log "Step 2: Generating tx2gene from GTF..."
"$STAR_FLEX/tests/transcriptvb/make_gene_map_from_gtf.sh" \
    --gtf "$GTF" --out tx2gene.tsv

# Step 3: Run CLI tool on STAR output
log "Step 3: Running tximport_compat CLI on STAR quant.sf..."
$CLI --quant star_quant.sf --tx2gene tx2gene.tsv --output cli_genes.sf --stats

# Step 4: Run R tximport
log "Step 4: Running R tximport for reference..."
cd "$RENV_DIR"
Rscript - "$OUTDIR" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

library(tximport)

quant_path <- file.path(outdir, "star_quant.sf")
tx2gene <- read.delim(file.path(outdir, "tx2gene.tsv"), header=FALSE, 
                      col.names=c("tx","gene"), stringsAsFactors=FALSE)

txi <- tximport(quant_path, type="salmon", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM", dropInfReps=TRUE)

output <- data.frame(
    Name = rownames(txi$counts),
    Length = round(txi$length[,1], 3),
    EffectiveLength = round(txi$length[,1], 3),
    TPM = txi$abundance[,1],
    NumReads = round(txi$counts[,1], 3)
)
write.table(output, file.path(outdir, "r_tximport_genes.sf"), 
            sep="\t", quote=FALSE, row.names=FALSE)

cat("R tximport completed\n")
cat("Genes:", nrow(output), "\n")
cat("Total counts:", sum(txi$counts[,1]), "\n")
RSCRIPT

cd "$OUTDIR"

# Step 5: Compare STAR tximport output vs R tximport
log "Step 5: Comparing STAR tximport output vs R tximport..."

python3 << 'PYEOF'
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

tolerance = 1e-6

star_order, star_tximport = parse_sf_ordered("star_quant.genes.tximport.sf")
r_order, r_tximport = parse_sf_ordered("r_tximport_genes.sf")

# Check same genes
if set(star_tximport.keys()) != set(r_tximport.keys()):
    missing = set(r_tximport.keys()) - set(star_tximport.keys())
    extra = set(star_tximport.keys()) - set(r_tximport.keys())
    if missing:
        print(f"Missing genes in STAR: {list(missing)[:5]}...")
    if extra:
        print(f"Extra genes in STAR: {list(extra)[:5]}...")
    print("FAIL: Gene sets don't match")
    sys.exit(1)

# Check gene ordering matches
if star_order != r_order:
    for i, (s, r) in enumerate(zip(star_order, r_order)):
        if s != r:
            print(f"FAIL: Gene order mismatch at position {i}")
            print(f"  STAR: {s}")
            print(f"  R: {r}")
            sys.exit(1)
    print(f"FAIL: Gene count mismatch: STAR {len(star_order)}, R {len(r_order)}")
    sys.exit(1)

# Compare values including Length column
max_rel_diff = 0
max_diff_gene = None
max_diff_col = None

for gene in star_tximport:
    for col in ['Length', 'EffectiveLength', 'TPM', 'NumReads']:
        star_val = star_tximport[gene][col]
        r_val = r_tximport[gene][col]
        
        if abs(star_val) < 1e-10 and abs(r_val) < 1e-10:
            continue
        
        if abs(r_val) < 1e-10:
            rel_diff = abs(star_val)
        else:
            rel_diff = abs(star_val - r_val) / max(abs(r_val), 1e-10)
        
        if rel_diff > max_rel_diff:
            max_rel_diff = rel_diff
            max_diff_gene = gene
            max_diff_col = col

if max_rel_diff > tolerance:
    print(f"FAIL: Max relative difference {max_rel_diff:.2e} exceeds tolerance {tolerance}")
    print(f"  Gene: {max_diff_gene}, Column: {max_diff_col}")
    print(f"  STAR: {star_tximport[max_diff_gene][max_diff_col]}")
    print(f"  R: {r_tximport[max_diff_gene][max_diff_col]}")
    sys.exit(1)
else:
    print(f"PASS: Max relative difference {max_rel_diff:.2e} within tolerance {tolerance}")
    print(f"  Genes compared: {len(star_tximport)}, ordering verified, Length column verified")
    sys.exit(0)
PYEOF

PARITY_RESULT=$?

# Step 6: Also compare CLI output vs R tximport
log "Step 6: Comparing CLI tool output vs R tximport..."

python3 << 'PYEOF'
import sys

def parse_sf(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            name = fields[0]
            data[name] = {h: float(fields[i]) if i > 0 else fields[i] 
                         for i, h in enumerate(header)}
    return data

tolerance = 1e-6

cli = parse_sf("cli_genes.sf")
r_tximport = parse_sf("r_tximport_genes.sf")

max_rel_diff = 0
for gene in cli:
    if gene not in r_tximport:
        continue
    for col in ['TPM', 'NumReads']:
        cli_val = cli[gene][col]
        r_val = r_tximport[gene][col]
        if abs(r_val) > 1e-10:
            rel_diff = abs(cli_val - r_val) / abs(r_val)
            max_rel_diff = max(max_rel_diff, rel_diff)

if max_rel_diff > tolerance:
    print(f"FAIL: CLI vs R max diff {max_rel_diff:.2e}")
    sys.exit(1)
else:
    print(f"PASS: CLI vs R max diff {max_rel_diff:.2e}")
    sys.exit(0)
PYEOF

CLI_RESULT=$?

# Summary
log "=== Summary ==="
if [ $PARITY_RESULT -eq 0 ] && [ $CLI_RESULT -eq 0 ]; then
    log "ALL TESTS PASSED"
    log "Output files in: $OUTDIR"
    
    # Create summary
    cat > "$OUTDIR/SUMMARY.md" << EOF
# STAR tximport End-to-End Test Results

**Date**: $(date)
**STAR**: $STAR_BIN
**tximport_compat**: $CLI

## Results

- STAR tximport vs R tximport: **PASSED**
- CLI tool vs R tximport: **PASSED**

## Files

- \`star_quant.sf\`: STAR transcript-level quantification
- \`star_quant.genes.sf\`: STAR gene-level (Legacy mode)
- \`star_quant.genes.tximport.sf\`: STAR gene-level (tximport mode)
- \`cli_genes.sf\`: CLI tool output
- \`r_tximport_genes.sf\`: R tximport reference

EOF
    exit 0
else
    log "TESTS FAILED"
    exit 1
fi



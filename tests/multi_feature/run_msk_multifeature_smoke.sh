#!/usr/bin/env bash
# E2E smoke test: MSK 30polyKO multi-feature (mRNA + PolyIII + LARRY)
# Validates that STAR pf-multi can process multiple feature libraries with
# per-library star_feature_ref, star_chemistry, and data-driven routing.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR="$REPO_ROOT/core/legacy/source/STAR"
GENOME="${MSK_MULTI_GENOME:-/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star}"
FIXTURE=${1:-/tmp/msk_multi_fixture}
OUT_DIR=${2:-/tmp/msk_multifeature_smoke_$(date +%Y%m%d_%H%M%S)}

if [ ! -x "$STAR" ]; then
    echo "FATAL: STAR binary not found at $STAR"
    exit 1
fi
if [ ! -d "$FIXTURE/mRNA" ] || [ ! -d "$FIXTURE/PolyIII" ] || [ ! -d "$FIXTURE/LARRY" ]; then
    echo "FATAL: fixture not found at $FIXTURE — run create_fixture.sh first"
    exit 1
fi

echo "=== MSK Multi-Feature E2E Smoke Test ==="
echo "  STAR: $STAR"
echo "  Genome: $GENOME"
echo "  Fixture: $FIXTURE"
echo "  Output: $OUT_DIR"
echo ""

mkdir -p "$OUT_DIR"

# Write pfMultiConfig with:
#   - No global [feature] ref (all feature libs have star_feature_ref)
#   - Per-library chemistry: mRNA=TRU, PolyIII=NXT, LARRY=TRU
#   - Per-library feature ref
#   - Explicit star_library_id for output provenance
cat > "$OUT_DIR/multi_config.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_whitelist,star_feature_ref,star_library_id
$FIXTURE/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,,gex_de
$FIXTURE/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt,/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv,grna_de
$FIXTURE/LARRY,DE_30KO,Custom,Custom,TRU,/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt,/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv,larry_de
EOF

echo "Config:"
cat "$OUT_DIR/multi_config.csv"
echo ""

# Run STAR with pf-multi config
echo "Running STAR..."
# GEX reads are the mRNA library FASTQs
GEX_R2=$(ls "$FIXTURE"/mRNA/*_R2_*.fastq.gz | head -n1)
GEX_R1=$(ls "$FIXTURE"/mRNA/*_R1_*.fastq.gz | head -n1)

"$STAR" \
    --runMode alignReads \
    --runThreadN 8 \
    --genomeDir "$GENOME" \
    --readFilesIn "$GEX_R2" "$GEX_R1" \
    --readFilesCommand zcat \
    --pfMultiConfig "$OUT_DIR/multi_config.csv" \
    --defaultCrCompat yes \
    --crChemistry auto \
    --outFileNamePrefix "$OUT_DIR/" \
    --outSAMtype BAM Unsorted \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist /storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures Gene GeneFull \
    --soloBarcodeReadLength 0 \
    --readMapNumber 50000 \
    2>&1 | tee "$OUT_DIR/star_stdout.log"

STAR_EXIT=$?
echo ""
echo "STAR exit code: $STAR_EXIT"

if [ "$STAR_EXIT" -ne 0 ]; then
    echo "FATAL: STAR run failed"
    echo "Check: $OUT_DIR/Log.out"
    exit 1
fi

echo ""
echo "=== Validation ==="

PASS=0
FAIL=0

check() {
    local desc="$1"
    local cond="$2"
    if eval "$cond"; then
        echo "  PASS: $desc"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $desc"
        FAIL=$((FAIL + 1))
    fi
}

# Check that assign output directories exist for each library
check "grna_de assign output exists" \
    "[ -d '$OUT_DIR/cr_assign/CRISPR_Guide_Capture/grna_de' ]"

check "larry_de assign output exists" \
    "[ -d '$OUT_DIR/cr_assign/Custom/larry_de' ]"

# Check for chemistry logging per library in Log.out
check "grna_de chemistry logged" \
    "grep -q 'grna_de.*CRISPR Guide Capture' '$OUT_DIR/Log.out'"

check "larry_de chemistry logged" \
    "grep -q 'larry_de.*Custom' '$OUT_DIR/Log.out'"

# Check per-library star_feature_ref logged
check "grna_de per-library feature ref logged" \
    "grep -q 'star_feature_ref=.*ref_feature_geneBC.csv' '$OUT_DIR/Log.out'"

check "larry_de per-library feature ref logged" \
    "grep -q 'star_feature_ref=.*ref_feature_larryBC.csv' '$OUT_DIR/Log.out'"

# Check per-library star_whitelist logged
check "grna_de per-library whitelist logged" \
    "grep -q 'grna_de.*whitelist=.*/3M-february-2018_NXT.txt' '$OUT_DIR/Log.out'"

check "larry_de per-library whitelist logged" \
    "grep -q 'larry_de.*whitelist=.*/3M-february-2018_TRU.txt' '$OUT_DIR/Log.out'"

# Check that CRISPR is recognized as known type (no NOTICE for it)
check "CRISPR is recognized as known type (no NOTICE)" \
    "! grep -q \"NOTICE: feature_types 'CRISPR Guide Capture' is not a known\" '$OUT_DIR/Log.out'"

# Check that Custom IS flagged as unknown type
check "Custom type has NOTICE for unknown type" \
    "grep -q \"NOTICE: feature_types 'Custom' is not a known\" '$OUT_DIR/Log.out'"

# Check no global feature ref error (we intentionally omit [feature] section)
check "no global feature ref error" \
    "! grep -q 'Feature reference not provided' '$OUT_DIR/Log.out'"

echo ""
echo "=========================================="
echo "Results: $PASS passed, $FAIL failed"
echo "=========================================="
echo ""
echo "Output: $OUT_DIR"
echo "Log: $OUT_DIR/Log.out"

[ "$FAIL" -eq 0 ]

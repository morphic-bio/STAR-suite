#!/bin/bash
# Regression test comparing commit 2f89f17 (baseline) with current HEAD
# Tests MEX parity using the same dataset/parameters as run_flex_multisample_test.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
BASELINE_COMMIT="2f89f17"
BASELINE_DIR="${SCRIPT_DIR}/regress_100K_baseline"
CURRENT_DIR="${SCRIPT_DIR}/regress_100K_current"
TMP_DIR_BASE="/storage/100K/tmp/regress_100K"

# Save original commit and set up cleanup trap
ORIGINAL_COMMIT=$(cd "$REPO_DIR" && git rev-parse --short HEAD)
HAS_STASH=false
cleanup() {
    cd "$REPO_DIR"
    git checkout "$ORIGINAL_COMMIT" > /dev/null 2>&1 || true
    if [ "$HAS_STASH" = "true" ]; then
        git stash pop > /dev/null 2>&1 || true
    fi
}
trap cleanup EXIT

# Dataset and reference paths (same as run_flex_multisample_test.sh)
GENOME_DIR="/storage/flex_filtered_reference/star_index"
CB_WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
SAMPLE_WHITELIST="/storage/SC2300771_filtered_2M/sample_whitelist.tsv"
PROBE_LIST="/storage/flex_filtered_reference/filtered_reference/probe_list.txt"
SAMPLE_PROBES="/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt"
READS_R2="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz"
READS_R1="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz"

# Function to run STAR with given output directory
run_star_test() {
    local OUT_DIR="$1"
    local TMP_DIR="$2"
    local STAR_BIN="${REPO_DIR}/core/legacy/source/STAR"
    
    echo "Running STAR..."
    echo "  Output: $OUT_DIR"
    echo "  Temp: $TMP_DIR"
    
    rm -rf "$OUT_DIR"
    rm -rf "$TMP_DIR"
    mkdir -p "$OUT_DIR"
    
    "$STAR_BIN" \
      --runThreadN 8 \
      --outTmpDir "$TMP_DIR" \
      --genomeDir "$GENOME_DIR" \
      --soloType CB_UMI_Simple \
      --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
      --soloCBwhitelist "$CB_WHITELIST" \
      --flex yes \
      --soloFlexExpectedCellsPerTag 3000 \
      --soloSampleWhitelist "$SAMPLE_WHITELIST" \
      --soloProbeList "$PROBE_LIST" \
      --soloSampleProbes "$SAMPLE_PROBES" \
      --soloSampleProbeOffset 68 \
      --soloFlexAllowedTags "$SAMPLE_WHITELIST" \
      --soloFlexOutputPrefix "${OUT_DIR}/per_sample" \
      --limitIObufferSize 50000000 50000000 \
      --outSJtype None \
      --outBAMcompression 6 \
      --soloMultiMappers Rescue \
      --alignIntronMax 500000 \
      --outFilterMismatchNmax 6 \
      --outFilterMismatchNoverReadLmax 1.0 \
      --outFilterMatchNmin 25 \
      --outSAMunmapped None \
      --outFilterMatchNminOverLread 0 \
      --outFilterMultimapNmax 10000 \
      --outFilterMultimapScoreRange 4 \
      --outSAMmultNmax 10000 \
      --winAnchorMultimapNmax 200 \
      --outSAMprimaryFlag AllBestScore \
      --outFilterScoreMin 0 \
      --outFilterScoreMinOverLread 0 \
      --outSAMattributes NH HI AS nM NM GX GN ZG \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIfiltering MultiGeneUMI_CR \
      --soloUMIdedup 1MM_CR \
      --soloCellFilter None \
      --clipAdapterType CellRanger4 \
      --soloFeatures Gene \
      --alignEndsType Local \
      --soloAddTagsToUnsorted no \
      --soloStrand Unstranded \
      --chimSegmentMin 1000000 \
      --soloKeysCompat cr \
      --outSAMtype BAM Unsorted \
      --soloSampleSearchNearby no \
      --readFilesCommand zcat \
      --readFilesIn "$READS_R2" "$READS_R1" \
      --outFileNamePrefix "$OUT_DIR/"
    
    echo "  ✓ STAR run completed"
}

# Function to build STAR
build_star() {
    local COMMIT="$1"
    echo ""
    echo "=== Building STAR at commit $COMMIT ==="
    cd "$REPO_DIR"
    git checkout "$COMMIT" > /dev/null 2>&1
    cd source
    if [ ! -f "Makefile" ]; then
        echo "ERROR: Makefile not found at commit $COMMIT" >&2
        exit 1
    fi
    make clean > /dev/null 2>&1 || true
    # Build STAR - try STAR target first, fall back to default target
    if ! make -j$(nproc) STAR 2>&1 | tee /tmp/build_${COMMIT}.log | tail -10; then
        # Check if it failed due to missing target
        if grep -q "No rule to make target.*STAR" /tmp/build_${COMMIT}.log; then
            echo "  STAR target not found, using default target..."
            make clean > /dev/null 2>&1 || true
            make -j$(nproc) 2>&1 | tee -a /tmp/build_${COMMIT}.log | tail -10
        else
            echo "ERROR: Build failed at commit $COMMIT" >&2
            tail -20 /tmp/build_${COMMIT}.log >&2
            exit 1
        fi
    fi
    if [ ! -f "STAR" ]; then
        echo "ERROR: Failed to build STAR at commit $COMMIT" >&2
        echo "Build log saved to /tmp/build_${COMMIT}.log" >&2
        tail -20 /tmp/build_${COMMIT}.log >&2
        exit 1
    fi
    echo "  ✓ STAR built successfully"
}

# Get current commit
CURRENT_COMMIT="$ORIGINAL_COMMIT"

# Safety check: ensure we're in a clean state
cd "$REPO_DIR"
if [ -n "$(git status --porcelain)" ]; then
    if [ "${SKIP_CONFIRM:-}" != "yes" ]; then
        echo "WARNING: Working directory has uncommitted changes" >&2
        echo "This script will checkout different commits. Continue? (y/N)"
        read -r response
        if [ "$response" != "y" ] && [ "$response" != "Y" ]; then
            echo "Aborted."
            exit 1
        fi
    else
        echo "WARNING: Working directory has uncommitted changes (SKIP_CONFIRM=yes, stashing)"
    fi
    # Stash changes to allow clean checkout
    git stash > /dev/null 2>&1
    HAS_STASH=true
fi

# Verify baseline commit exists
if ! git rev-parse --verify "$BASELINE_COMMIT" > /dev/null 2>&1; then
    echo "ERROR: Baseline commit $BASELINE_COMMIT not found" >&2
    exit 1
fi

echo "=== 100K Regression Test ==="
echo "Baseline commit: $BASELINE_COMMIT"
echo "Current commit: $CURRENT_COMMIT"
echo "Original commit: $ORIGINAL_COMMIT"
echo ""

# Step 1: Build and run baseline
echo "=== Step 1: Running baseline ($BASELINE_COMMIT) ==="
build_star "$BASELINE_COMMIT"
run_star_test "$BASELINE_DIR" "${TMP_DIR_BASE}_baseline"

# Step 2: Build and run current HEAD
echo ""
echo "=== Step 2: Running current HEAD ($CURRENT_COMMIT) ==="
build_star "$CURRENT_COMMIT"
run_star_test "$CURRENT_DIR" "${TMP_DIR_BASE}_current"

# Step 3: Compare outputs
echo ""
echo "=== Step 3: Comparing MEX outputs ==="
"${SCRIPT_DIR}/compare_mex_outputs.sh" "$BASELINE_DIR" "$CURRENT_DIR"

# Note: cleanup trap will restore original commit automatically
echo ""
echo "=== Regression test completed ==="
echo "Baseline outputs: $BASELINE_DIR"
echo "Current outputs: $CURRENT_DIR"
echo ""
echo "To compare outputs manually:"
echo "  ${SCRIPT_DIR}/compare_mex_outputs.sh $BASELINE_DIR $CURRENT_DIR"


#!/bin/bash
# Salmon Parity Test for EC Cleaner
# Compares our EC output against Salmon's --dumpEqWeights output
#
# Usage: ./run_salmon_parity.sh [output_dir]

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FIXTURE_DIR="${SCRIPT_DIR}/../../test/fixtures/salmon_eq"
OUTPUT_DIR="${1:-/tmp/ec_parity_test}"

# Paths to fixtures
TRANSCRIPTOME="${FIXTURE_DIR}/chr22_trans.fa"
BAM_FILE="${FIXTURE_DIR}/alignments.bam"
R1="${FIXTURE_DIR}/synthetic_reads1.fq"
R2="${FIXTURE_DIR}/synthetic_reads2.fq"

# Salmon binary (use local build or system)
SALMON="${SALMON:-/mnt/pikachu/salmon/build/src/salmon}"

# Parameters (must match between Salmon and our implementation)
SCORE_EXP="1.0"                    # Salmon default: 1.0
MIN_ALN_PROB="1e-5"                # Salmon default: 1e-5
MIN_SCORE_FRACTION="0.65"          # Salmon default: 0.65
RANGE_FACTORIZATION_BINS="4"       # Salmon default: 4
HARD_FILTER="false"                # Salmon default: false
THREADS="1"                        # Single-threaded for determinism

mkdir -p "${OUTPUT_DIR}"

echo "=== Salmon Parity Test ==="
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Step 1: Build Salmon index (if not exists)
SALMON_IDX="${OUTPUT_DIR}/salmon_idx"
if [ ! -d "${SALMON_IDX}" ]; then
    echo "Step 1: Building Salmon index..."
    "${SALMON}" index -t "${TRANSCRIPTOME}" -i "${SALMON_IDX}"
else
    echo "Step 1: Using existing Salmon index"
fi

# Step 2: Run Salmon in alignment mode (if BAM exists) or mapping mode
echo ""
echo "Step 2: Running Salmon quant with --dumpEqWeights..."

SALMON_OUT="${OUTPUT_DIR}/salmon_out"

if [ -f "${BAM_FILE}" ]; then
    # Alignment mode
    echo "  Using alignment mode with BAM file"
    "${SALMON}" quant \
        -t "${TRANSCRIPTOME}" \
        -l A \
        -a "${BAM_FILE}" \
        --dumpEqWeights \
        --noLengthCorrection \
        --noEffectiveLengthCorrection \
        --noFragLengthDist \
        --validateMappings \
        --hardFilter \
        --minScoreFraction "${MIN_SCORE_FRACTION}" \
        --scoreExp "${SCORE_EXP}" \
        --minAlnProb "${MIN_ALN_PROB}" \
        --rangeFactorizationBins "${RANGE_FACTORIZATION_BINS}" \
        --threads "${THREADS}" \
        -o "${SALMON_OUT}"
else
    # Mapping mode (fallback)
    echo "  (No BAM found, using mapping mode)"
    "${SALMON}" quant \
        -i "${SALMON_IDX}" \
        -l A \
        -1 "${R1}" \
        -2 "${R2}" \
        --dumpEqWeights \
        --noLengthCorrection \
        --noEffectiveLengthCorrection \
        --noFragLengthDist \
        --validateMappings \
        --hardFilter \
        --minScoreFraction "${MIN_SCORE_FRACTION}" \
        --scoreExp "${SCORE_EXP}" \
        --minAlnProb "${MIN_ALN_PROB}" \
        --rangeFactorizationBins "${RANGE_FACTORIZATION_BINS}" \
        --threads "${THREADS}" \
        -o "${SALMON_OUT}"
fi

# Extract Salmon ECs (handle gzipped output)
SALMON_EC="${OUTPUT_DIR}/salmon_eq_classes.txt"
if [ -f "${SALMON_OUT}/aux_info/eq_classes.txt.gz" ]; then
    gunzip -c "${SALMON_OUT}/aux_info/eq_classes.txt.gz" > "${SALMON_EC}"
else
    cp "${SALMON_OUT}/aux_info/eq_classes.txt" "${SALMON_EC}"
fi
echo "  Salmon ECs: ${SALMON_EC} ($(wc -l < ${SALMON_EC}) lines)"

# Step 3: Run our EC cleaner with matching parameters
echo ""
echo "Step 3: Running our EC cleaner..."

OUR_OUT="${OUTPUT_DIR}/our_out"
OUR_EC="${OUR_OUT}/eq_classes.txt"
mkdir -p "${OUR_OUT}"

# IMPORTANT: Disable local/global pruning for Salmon parity
if [ -f "${BAM_FILE}" ]; then
    ./ec_filter_cli \
        --input "${BAM_FILE}" \
        --transcripts "${TRANSCRIPTOME}" \
        --score-exp "${SCORE_EXP}" \
        --min-aln-prob "${MIN_ALN_PROB}" \
        --min-score-fraction "${MIN_SCORE_FRACTION}" \
        --range-factorization-bins "${RANGE_FACTORIZATION_BINS}" \
        --hard-filter "${HARD_FILTER}" \
        --no-local-pruning \
        --no-global-pruning \
        --threads "${THREADS}" \
        --output-format salmon \
        -o "${OUR_EC}"
else
    ./ec_filter_cli \
        --input "${R1}" \
        --input2 "${R2}" \
        --transcripts "${TRANSCRIPTOME}" \
        --score-exp "${SCORE_EXP}" \
        --min-aln-prob "${MIN_ALN_PROB}" \
        --min-score-fraction "${MIN_SCORE_FRACTION}" \
        --range-factorization-bins "${RANGE_FACTORIZATION_BINS}" \
        --hard-filter "${HARD_FILTER}" \
        --no-local-pruning \
        --no-global-pruning \
        --threads "${THREADS}" \
        --output-format salmon \
        -o "${OUR_EC}"
fi

echo "  Our ECs: ${OUR_EC} ($(wc -l < ${OUR_EC}) lines)"

# Step 4: Compare EC files
echo ""
echo "Step 4: Comparing EC files..."

python3 "${SCRIPT_DIR}/compare_ecs.py" \
    --salmon "${SALMON_EC}" \
    --ours "${OUR_EC}" \
    --tolerance 1e-6 \
    --report "${OUTPUT_DIR}/parity_report.txt"

echo ""
echo "=== Parity Test Complete ==="
echo "Report: ${OUTPUT_DIR}/parity_report.txt"

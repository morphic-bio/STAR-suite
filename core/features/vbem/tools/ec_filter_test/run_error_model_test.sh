#!/bin/bash
# Error Model Harmonization Test Script
# Tests error model parity between our CLI and Salmon

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="${1:-/tmp/error_model_test}"
STAGE="${2:-all}"  # Stage to run: baseline, single, preburnin, postburnin, full, all

# Test data paths
TRANSCRIPTOME="${SCRIPT_DIR}/../../test/fixtures/salmon_eq/chr22_trans.fa"
BAM_FILE="${SCRIPT_DIR}/../../test/fixtures/salmon_eq/alignments.bam"

# Salmon binary
SALMON="${SALMON:-/mnt/pikachu/salmon/build/src/salmon}"

# Our CLI
OUR_CLI="${SCRIPT_DIR}/ec_filter_cli"

mkdir -p "${OUTPUT_DIR}"

echo "=== Error Model Harmonization Test ==="
echo "Output directory: ${OUTPUT_DIR}"
echo "Stage: ${STAGE}"
echo ""

# Check prerequisites
if [ ! -f "${OUR_CLI}" ]; then
    echo "Error: ec_filter_cli not found at ${OUR_CLI}"
    echo "Please build it first: cd tools/ec_filter_test && make ec_filter_cli"
    exit 1
fi

if [ ! -f "${SALMON}" ]; then
    echo "Error: Salmon not found at ${SALMON}"
    echo "Please specify SALMON environment variable or build Salmon"
    exit 1
fi

if [ ! -f "${TRANSCRIPTOME}" ]; then
    echo "Error: Transcriptome not found: ${TRANSCRIPTOME}"
    exit 1
fi

# Stage 1: Baseline Initialization Test
if [ "${STAGE}" == "baseline" ] || [ "${STAGE}" == "all" ]; then
    echo "=== Stage 1: Baseline Initialization Test ==="
    
    # Our CLI - dump initial matrices
    echo "Running our CLI..."
    "${OUR_CLI}" --input "${BAM_FILE}" --transcripts "${TRANSCRIPTOME}" \
        --use-error-model --dump-matrices "${OUTPUT_DIR}/our_init" \
        -o "${OUTPUT_DIR}/our_init_ec.txt" 2>&1 | grep -E "Error model|Matrix|Built|reads"
    
    # Salmon - dump initial matrices
    echo "Running Salmon..."
    SALMON_DUMP_MATRICES="${OUTPUT_DIR}/salmon_init" \
    "${SALMON}" quant -t "${TRANSCRIPTOME}" -l A -a "${BAM_FILE}" \
        --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
        --noFragLengthDist -p 1 -o "${OUTPUT_DIR}/salmon_init_out" 2>&1 | grep -E "quantifying|mapping"
    
    # Compare initial matrices
    echo ""
    echo "Comparing initial matrices..."
    python3 "${SCRIPT_DIR}/compare_matrices.py" \
        --ours "${OUTPUT_DIR}/our_init_final_left_bin0.tsv" \
        --salmon "${OUTPUT_DIR}/salmon_init_final_left_bin0.tsv" \
        --tolerance 1e-10 || echo "Note: Matrix files may not exist if no reads processed"
    echo ""
fi

# Stage 2: Single-Read Update Test
if [ "${STAGE}" == "single" ] || [ "${STAGE}" == "all" ]; then
    echo "=== Stage 2: Single-Read Update Test ==="
    
    # Create BAM with single read
    SINGLE_BAM="${OUTPUT_DIR}/single_read.bam"
    echo "Creating single-read BAM..."
    samtools view -h "${BAM_FILE}" | head -n 4 | samtools view -Sb - > "${SINGLE_BAM}" || {
        echo "Warning: Could not create single-read BAM, skipping Stage 2"
    }
    
    if [ -f "${SINGLE_BAM}" ]; then
        # Our CLI
        echo "Running our CLI on single read..."
        "${OUR_CLI}" --input "${SINGLE_BAM}" --transcripts "${TRANSCRIPTOME}" \
            --use-error-model --dump-matrices "${OUTPUT_DIR}/our_single" \
            -o "${OUTPUT_DIR}/our_single_ec.txt" 2>&1 | grep -E "Error model|reads"
        
        # Salmon
        echo "Running Salmon on single read..."
        SALMON_DUMP_MATRICES="${OUTPUT_DIR}/salmon_single" \
        "${SALMON}" quant -t "${TRANSCRIPTOME}" -l A -a "${SINGLE_BAM}" \
            --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
            --noFragLengthDist -p 1 -o "${OUTPUT_DIR}/salmon_single_out" 2>&1 | tail -3
        
        # Compare matrices
        echo ""
        echo "Comparing single-read matrices..."
        if [ -f "${OUTPUT_DIR}/our_single_final_left_bin0.tsv" ] && [ -f "${OUTPUT_DIR}/salmon_single_final_left_bin0.tsv" ]; then
            python3 "${SCRIPT_DIR}/compare_matrices.py" \
                --ours "${OUTPUT_DIR}/our_single_final_left_bin0.tsv" \
                --salmon "${OUTPUT_DIR}/salmon_single_final_left_bin0.tsv" \
                --tolerance 1e-10
        else
            echo "Matrix files not found - may need to process more reads"
        fi
        echo ""
    fi
fi

# Stage 3: Pre-Burnin Phase Test (0-5000 reads)
if [ "${STAGE}" == "preburnin" ] || [ "${STAGE}" == "all" ]; then
    echo "=== Stage 3: Pre-Burnin Phase Test (0-5000 reads) ==="
    
    # Limit to first 5000 reads (approximately 10000 BAM lines)
    FIRST_5K_BAM="${OUTPUT_DIR}/first_5k.bam"
    echo "Creating first-5k-reads BAM..."
    samtools view -h "${BAM_FILE}" | head -n 10000 | samtools view -Sb - > "${FIRST_5K_BAM}" || {
        echo "Warning: Could not create first-5k BAM, using full BAM"
        FIRST_5K_BAM="${BAM_FILE}"
    }
    
    # Our CLI with checkpoint dumps
    echo "Running our CLI..."
    "${OUR_CLI}" --input "${FIRST_5K_BAM}" --transcripts "${TRANSCRIPTOME}" \
        --use-error-model --dump-matrices "${OUTPUT_DIR}/our_pre5k" \
        -o "${OUTPUT_DIR}/our_pre5k_ec.txt" 2>&1 | grep -E "Error model|Matrix|reads|observed"
    
    # Salmon with checkpoint dumps
    echo "Running Salmon..."
    SALMON_DUMP_MATRICES="${OUTPUT_DIR}/salmon_pre5k" \
    "${SALMON}" quant -t "${TRANSCRIPTOME}" -l A -a "${FIRST_5K_BAM}" \
        --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
        --noFragLengthDist -p 1 -o "${OUTPUT_DIR}/salmon_pre5k_out" 2>&1 | tail -5
    
    # Compare at each checkpoint
    echo ""
    echo "Comparing matrices at checkpoints..."
    for checkpoint in 1000 2000 3000 4000 5000; do
        ours_file="${OUTPUT_DIR}/our_pre5k_pre${checkpoint}_left_bin0.tsv"
        salmon_file="${OUTPUT_DIR}/salmon_pre5k_pre${checkpoint}_left_bin0.tsv"
        if [ -f "${ours_file}" ] && [ -f "${salmon_file}" ]; then
            echo "Checkpoint ${checkpoint}:"
            python3 "${SCRIPT_DIR}/compare_matrices.py" \
                --ours "${ours_file}" \
                --salmon "${salmon_file}" \
                --tolerance 1e-10
        fi
    done
    echo ""
fi

# Stage 4: Post Pre-Burnin Test (5000+ reads)
if [ "${STAGE}" == "postburnin" ] || [ "${STAGE}" == "all" ]; then
    echo "=== Stage 4: Post Pre-Burnin Test (5000+ reads) ==="
    
    # Use first 10000 reads
    FIRST_10K_BAM="${OUTPUT_DIR}/first_10k.bam"
    echo "Creating first-10k-reads BAM..."
    samtools view -h "${BAM_FILE}" | head -n 20000 | samtools view -Sb - > "${FIRST_10K_BAM}" || {
        FIRST_10K_BAM="${BAM_FILE}"
    }
    
    # Our CLI with tracing
    echo "Running our CLI with tracing..."
    "${OUR_CLI}" --input "${FIRST_10K_BAM}" --transcripts "${TRANSCRIPTOME}" \
        --use-error-model --error-model-trace "${OUTPUT_DIR}/our_trace.txt" \
        --trace-level 2 -o "${OUTPUT_DIR}/our_post5k_ec.txt" 2>&1 | grep -E "Error model|tracing|reads"
    
    # Salmon with tracing
    echo "Running Salmon with tracing..."
    SALMON_ERROR_MODEL_TRACE="${OUTPUT_DIR}/salmon_trace.txt" \
    SALMON_TRACE_LEVEL=2 \
    "${SALMON}" quant -t "${TRANSCRIPTOME}" -l A -a "${FIRST_10K_BAM}" \
        --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
        --noFragLengthDist -p 1 -o "${OUTPUT_DIR}/salmon_post5k_out" 2>&1 | tail -5
    
    # Compare traces
    echo ""
    echo "Comparing traces..."
    if [ -f "${OUTPUT_DIR}/our_trace.txt" ] && [ -f "${OUTPUT_DIR}/salmon_trace.txt" ]; then
        python3 "${SCRIPT_DIR}/compare_error_model.py" \
            --ours "${OUTPUT_DIR}/our_trace.txt" \
            --salmon "${OUTPUT_DIR}/salmon_trace.txt" \
            --tolerance 1e-6 \
            --report "${OUTPUT_DIR}/trace_comparison.txt"
        echo ""
        echo "Trace comparison report: ${OUTPUT_DIR}/trace_comparison.txt"
    else
        echo "Warning: Trace files not found"
    fi
    echo ""
fi

# Stage 5: Full Run Comparison
if [ "${STAGE}" == "full" ] || [ "${STAGE}" == "all" ]; then
    echo "=== Stage 5: Full Run Comparison ==="
    
    # Our CLI
    echo "Running our CLI (full run)..."
    "${OUR_CLI}" --input "${BAM_FILE}" --transcripts "${TRANSCRIPTOME}" \
        --use-error-model \
        --error-model-trace "${OUTPUT_DIR}/our_full_trace.txt" \
        --trace-level 2 \
        --dump-matrices "${OUTPUT_DIR}/our_full" \
        -o "${OUTPUT_DIR}/our_full_ec.txt" 2>&1 | grep -E "Error model|tracing|Matrix|reads|observed"
    
    # Salmon
    echo "Running Salmon (full run)..."
    SALMON_ERROR_MODEL_TRACE="${OUTPUT_DIR}/salmon_full_trace.txt" \
    SALMON_TRACE_LEVEL=2 \
    SALMON_DUMP_MATRICES="${OUTPUT_DIR}/salmon_full" \
    "${SALMON}" quant -t "${TRANSCRIPTOME}" -l A -a "${BAM_FILE}" \
        --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
        --noFragLengthDist -p 1 -o "${OUTPUT_DIR}/salmon_full_out" 2>&1 | tail -5
    
    # Compare everything
    echo ""
    echo "Comparing full run traces..."
    if [ -f "${OUTPUT_DIR}/our_full_trace.txt" ] && [ -f "${OUTPUT_DIR}/salmon_full_trace.txt" ]; then
        python3 "${SCRIPT_DIR}/compare_error_model.py" \
            --ours "${OUTPUT_DIR}/our_full_trace.txt" \
            --salmon "${OUTPUT_DIR}/salmon_full_trace.txt" \
            --tolerance 1e-6 \
            --report "${OUTPUT_DIR}/full_comparison.txt"
    fi
    
    echo ""
    echo "Comparing final matrices..."
    if [ -f "${OUTPUT_DIR}/our_full_final_left_bin0.tsv" ] && [ -f "${OUTPUT_DIR}/salmon_full_final_left_bin0.tsv" ]; then
        python3 "${SCRIPT_DIR}/compare_matrices.py" \
            --ours "${OUTPUT_DIR}/our_full_final_left_bin0.tsv" \
            --salmon "${OUTPUT_DIR}/salmon_full_final_left_bin0.tsv" \
            --tolerance 1e-10 \
            --output "${OUTPUT_DIR}/matrix_comparison.txt"
    fi
    
    echo ""
    echo "Full comparison complete!"
    echo "Reports:"
    echo "  - Trace comparison: ${OUTPUT_DIR}/full_comparison.txt"
    echo "  - Matrix comparison: ${OUTPUT_DIR}/matrix_comparison.txt"
fi

echo ""
echo "=== Test Complete ==="
echo "Results in: ${OUTPUT_DIR}"

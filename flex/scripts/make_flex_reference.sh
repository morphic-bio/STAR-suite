#!/bin/bash

# make_flex_reference.sh
# Create FASTA and GTF pseudo-chromosomes from a 10x Flex probe CSV.
# Usage: ./make_flex_reference.sh <input_csv> <output_prefix> [--no-gtf-header]

set -euo pipefail

# Parse arguments
SKIP_GTF_HEADER=false
INPUT_CSV=""
OUTPUT_PREFIX=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --no-gtf-header)
            SKIP_GTF_HEADER=true
            shift
            ;;
        *)
            if [ -z "$INPUT_CSV" ]; then
                INPUT_CSV="$1"
            elif [ -z "$OUTPUT_PREFIX" ]; then
                OUTPUT_PREFIX="$1"
            else
                echo "Error: Too many arguments"
                exit 1
            fi
            shift
            ;;
    esac
done

# Check arguments
if [ -z "$INPUT_CSV" ] || [ -z "$OUTPUT_PREFIX" ]; then
    cat << EOF
Usage: $0 <input_csv> <output_prefix> [--no-gtf-header]
Example: $0 probes.csv reference
Example: $0 probes.csv reference --no-gtf-header
Outputs: reference.fa and reference.gtf
EOF
    exit 1
fi

FASTA_OUTPUT="${OUTPUT_PREFIX}.fa"
GTF_OUTPUT="${OUTPUT_PREFIX}.gtf"

if [ ! -f "$INPUT_CSV" ]; then
    echo "Error: Input file '$INPUT_CSV' not found"
    exit 1
fi

echo "Processing probe CSV: $INPUT_CSV"
echo "Output FASTA: $FASTA_OUTPUT"
echo "Output GTF: $GTF_OUTPUT"

# Initialize output files
> "$FASTA_OUTPUT"
> "$GTF_OUTPUT"

# Add GTF header (unless skipped)
if [ "$SKIP_GTF_HEADER" = false ]; then
    cat >> "$GTF_OUTPUT" << 'EOF'
##description: Flex probe reference annotation
##provider: make_flex_reference.sh
##format: gtf
EOF
fi

# Process CSV file (skip header)
{
    tail -n +2 "$INPUT_CSV"
    # Ensure there's a newline at the end if missing
    [ -n "$(tail -c1 "$INPUT_CSV")" ] && echo
} | while IFS=',' read -r gene_id probe_seq probe_id included region; do
    # Skip if not included
    if [ "$included" != "TRUE" ]; then
        continue
    fi

    # Extract gene name from probe_id (format: ENSG|GENE|<probe_key>)
    gene_name=$(echo "$probe_id" | cut -d'|' -f2)
    # Calculate sequence length
    seq_length=${#probe_seq}

    # Write FASTA entry
    echo ">${probe_id}" >> "$FASTA_OUTPUT"
    echo "$probe_seq" >> "$FASTA_OUTPUT"

    # Write GTF entry (single exon spanning entire sequence)
    # Format: seqname source feature start end score strand frame attributes
    printf "%s\tFLEX\texon\t1\t%d\t.\t+\t.\tgene_id \"%s\"; gene_name \"%s\"; transcript_id \"%s\"; probe_id \"%s\"; region \"%s\";\n" \
        "$probe_id" "$seq_length" "$gene_name" "$gene_name" "$probe_id" "$probe_id" "$region" >> "$GTF_OUTPUT"
done

echo "Generated FASTA file: $FASTA_OUTPUT"
echo "Generated GTF file: $GTF_OUTPUT"

# Count entries
fasta_count=$(grep -c "^>" "$FASTA_OUTPUT" || echo "0")
gtf_count=$(grep -c "^[^#]" "$GTF_OUTPUT" || echo "0")

echo "FASTA entries: $fasta_count"
echo "GTF entries: $gtf_count"

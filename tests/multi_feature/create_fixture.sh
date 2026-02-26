#!/usr/bin/env bash
# Create MSK 30polyKO multi-feature fixture (mRNA + PolyIII + LARRY, 100K reads each)
set -euo pipefail

SRC=/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM
OUT=${1:-/tmp/msk_multifeature_fixture_$(date +%Y%m%d_%H%M%S)}
NREADS=100000
NLINES=$((NREADS * 4))

echo "Source: $SRC"
echo "Output: $OUT"
echo "Reads per library: $NREADS"

mkdir -p "$OUT"/{mRNA,PolyIII,LARRY}

for LIB in mRNA PolyIII LARRY; do
    R1=$(ls "$SRC"/${LIB}_DE_*_L001_R1_001.fastq.gz 2>/dev/null | head -n1)
    if [ -z "$R1" ]; then
        echo "ERROR: no L001 R1 FASTQ found for $LIB"
        exit 1
    fi
    R2=${R1/_R1_/_R2_}
    if [ ! -f "$R2" ]; then
        echo "ERROR: R2 not found: $R2"
        exit 1
    fi

    echo "  $LIB: $(basename "$R1") → $NREADS reads"
    # Use set +o pipefail locally: zcat exits non-zero on SIGPIPE from head, which is expected.
    set +o pipefail
    zcat "$R1" | head -n "$NLINES" | gzip -c > "$OUT/$LIB/${LIB}_L001_R1_001.fastq.gz"
    zcat "$R2" | head -n "$NLINES" | gzip -c > "$OUT/$LIB/${LIB}_L001_R2_001.fastq.gz"
    set -o pipefail

    # Verify output is non-empty
    if [ ! -s "$OUT/$LIB/${LIB}_L001_R1_001.fastq.gz" ] || [ ! -s "$OUT/$LIB/${LIB}_L001_R2_001.fastq.gz" ]; then
        echo "ERROR: extraction failed for $LIB — output file is empty"
        exit 1
    fi
done

echo ""
echo "Fixture created at: $OUT"
echo ""
echo "Verify read counts:"
for LIB in mRNA PolyIII LARRY; do
    COUNT=$(zcat "$OUT/$LIB/${LIB}_L001_R1_001.fastq.gz" | wc -l)
    echo "  $LIB: $((COUNT / 4)) reads"
done

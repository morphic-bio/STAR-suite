#!/usr/bin/env bash
# Create a downsampled MSK 30polyKO multi-feature fixture suitable for fast
# Phase 5 validation (~2 min total runtime instead of 15+ min).
#
# Downsampling strategy:
#   - 10K reads per library (vs 100K in the full fixture)
#   - gRNA ref: used as-is (29 features, already small)
#   - LARRY ref: subset to first 500 features (vs 245K; 500x speedup in search)
set -euo pipefail

SRC=/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM
GRNA_REF=/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv
LARRY_REF=/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv
OUT=${1:-/tmp/msk_multi_fixture_ds}
NREADS=${2:-10000}
LARRY_FEATURES=${3:-500}

NLINES=$((NREADS * 4))

echo "=== Creating downsampled multi-feature fixture ==="
echo "  Source:         $SRC"
echo "  Output:         $OUT"
echo "  Reads/library:  $NREADS"
echo "  LARRY features: $LARRY_FEATURES (of $(wc -l < "$LARRY_REF") total)"
echo ""

rm -rf "$OUT"
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
    set +o pipefail
    zcat "$R1" | head -n "$NLINES" | gzip -c > "$OUT/$LIB/${LIB}_L001_R1_001.fastq.gz"
    zcat "$R2" | head -n "$NLINES" | gzip -c > "$OUT/$LIB/${LIB}_L001_R2_001.fastq.gz"
    set -o pipefail

    if [ ! -s "$OUT/$LIB/${LIB}_L001_R1_001.fastq.gz" ] || [ ! -s "$OUT/$LIB/${LIB}_L001_R2_001.fastq.gz" ]; then
        echo "ERROR: extraction failed for $LIB"
        exit 1
    fi
done

# gRNA ref is small (29 features) — copy as-is
cp "$GRNA_REF" "$OUT/ref_feature_geneBC.csv"

# LARRY ref: keep header + first N features
head -n 1 "$LARRY_REF" > "$OUT/ref_feature_larryBC.csv"
set +o pipefail
tail -n +2 "$LARRY_REF" | head -n "$LARRY_FEATURES" >> "$OUT/ref_feature_larryBC.csv"
set -o pipefail

echo ""
echo "Fixture created at: $OUT"
echo ""
echo "Read counts:"
for LIB in mRNA PolyIII LARRY; do
    COUNT=$(zcat "$OUT/$LIB/${LIB}_L001_R1_001.fastq.gz" | wc -l)
    echo "  $LIB: $((COUNT / 4)) reads"
done
echo "Feature ref sizes:"
echo "  gRNA: $(( $(wc -l < "$OUT/ref_feature_geneBC.csv") - 1 )) features"
echo "  LARRY: $(( $(wc -l < "$OUT/ref_feature_larryBC.csv") - 1 )) features"

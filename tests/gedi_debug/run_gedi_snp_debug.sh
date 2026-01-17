#!/usr/bin/env bash
set -euo pipefail

# Compile a patched GEDI class (classpath override) and run GEDI Slam SNP detection
# with targeted debug output for a locus.
#
# Usage:
#   tests/gedi_debug/run_gedi_snp_debug.sh <loc>
#
# Example:
#   tests/gedi_debug/run_gedi_snp_debug.sh 1:26429152
#
# Env overrides:
#   GEDI_SNP_DEBUG_OUT=/path/to/out.tsv
#   GEDI_SNP_DEBUG_MAX=5000
#

LOC="${1:-}"
if [[ -z "$LOC" ]]; then
  echo "Usage: $0 <loc>  (e.g. 1:26429152 or chr1:26429152)" >&2
  exit 2
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

GEDI_JAR="$ROOT/gedi.jar"
LIB_DIR="$ROOT/lib"

if [[ ! -f "$GEDI_JAR" ]]; then
  echo "Missing gedi.jar at: $GEDI_JAR" >&2
  exit 2
fi
if [[ ! -d "$LIB_DIR" ]]; then
  echo "Missing GEDI lib dir at: $LIB_DIR" >&2
  exit 2
fi

JAVA_SRC_DIR="$ROOT/tests/gedi_debug/java"
BUILD_DIR="$ROOT/tests/gedi_debug/build"
CLASSES_DIR="$BUILD_DIR/classes"
mkdir -p "$CLASSES_DIR"

echo "== Compile patched SlamDetectSnps (classpath override) =="
javac \
  -cp "$GEDI_JAR:$LIB_DIR/*" \
  -d "$CLASSES_DIR" \
  "$JAVA_SRC_DIR/gedi/slam/javapipeline/SlamDetectSnps.java"

WORK="/storage/gedi_snp_debug_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORK"

DEBUG_OUT="${GEDI_SNP_DEBUG_OUT:-$WORK/gedi_snp_debug.tsv}"
DEBUG_MAX="${GEDI_SNP_DEBUG_MAX:-5000}"

export GEDI_SNP_DEBUG_LOC="$LOC"
export GEDI_SNP_DEBUG_OUT="$DEBUG_OUT"
export GEDI_SNP_DEBUG_MAX="$DEBUG_MAX"

echo "== Run GEDI Slam (will override SlamDetectSnps) =="
echo "  debug loc: $GEDI_SNP_DEBUG_LOC"
echo "  debug out: $GEDI_SNP_DEBUG_OUT"
echo "  debug max: $GEDI_SNP_DEBUG_MAX"

# Inputs (same as our WT/no4sU run)
# GEDI expects a BAM path for -reads (see tests/run_wt_slam_gedi_comparison.sh)
READS_BAM="/storage/slam_e2e_arid1a_20260113/star/wt_Aligned.sortedByCoord.out.bam"
GENOMIC_OML="/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml"
PREFIX="$WORK/gedi_dbg"

# Match parameters used in gedi_new.param as closely as possible
java \
  --add-opens java.base/sun.nio.ch=ALL-UNNAMED \
  --add-opens java.base/sun.reflect.generics.reflectiveObjects=ALL-UNNAMED \
  --add-opens java.xml/com.sun.org.apache.xerces.internal.util=ALL-UNNAMED \
  --add-opens java.base/jdk.internal.reflect=ALL-UNNAMED \
  --add-opens java.base/jdk.internal.ref=ALL-UNNAMED \
  --add-opens java.base/java.lang=ALL-UNNAMED \
  -Xmx16G -Xms2048m \
  -cp "$CLASSES_DIR:$GEDI_JAR:$LIB_DIR/*:$ROOT/plugins/*" \
  executables.Slam \
    -reads "$READS_BAM" \
    -genomic "$GENOMIC_OML" \
    -prefix "$PREFIX" \
    -trim5p 6 -trim3p 8 \
    -strandness Sense \
    -err 0.001 \
    -snppval 0.001 \
    -snpConv 0.3 \
    -keep -D

echo "== Done =="
echo "Outputs:"
echo "  $WORK/"
echo "  debug: $DEBUG_OUT"


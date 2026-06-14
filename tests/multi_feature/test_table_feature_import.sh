#!/usr/bin/env bash
# Unit tests for PfMultiTableImport (CSV/TSV parser, duplicate collapse, rejections).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SOURCE_DIR="$REPO_ROOT/core/legacy/source"
PF_INCLUDE="$REPO_ROOT/core/features/process_features/include"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0
pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: pf-multi table feature import ==="

make -C "$SOURCE_DIR" -j8 PfMultiConfig.o PfMultiMexStub.o MexWriter.o MexWriterUtil.o >/dev/null
g++ -std=c++11 -O2 -DPF_TABLE_IMPORT_NO_PERMITS -I"$SOURCE_DIR" -I"$PF_INCLUDE" \
    -c "$SOURCE_DIR/PfMultiTableImport.cpp" -o "$WORK_DIR/PfMultiTableImport_test.o"

HARNESS="$WORK_DIR/test_table_import"
cat > "$WORK_DIR/test_table_import.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <fstream>
#include <string>
#include "PfMultiTableImport.h"

static int writeFile(const std::string& path, const std::string& body) {
    std::ofstream out(path.c_str());
    if (!out.is_open()) return 1;
    out << body;
    return 0;
}

int main() {
    const std::string work = std::getenv("TABLE_TEST_WORK") ? std::getenv("TABLE_TEST_WORK") : ".";
    const std::string whitelist = work + "/wl.txt";
    const std::string featureRef = work + "/features.csv";
    const std::string assignOut = work + "/out";
    const std::string mkdirCmd = "mkdir -p \"" + assignOut + "\"";
    if (system(mkdirCmd.c_str()) != 0) return 1;

    writeFile(whitelist, "AAACCCAAGAAACCAT\nGTATGTTCAGTAGCCT\n");
    writeFile(featureRef,
        "id,name,sequence,feature_type\n"
        "HIV_DNA,HIV_DNA,,Custom\n"
        "HIV_RNA,HIV_RNA,,Custom\n"
        "BADFEAT,BAD,,Custom\n");

    PfMultiTableImport::TableImportOptions opts;
    opts.featureTypeLabel = "Custom";
    opts.sampleName = "S1";
    opts.assignmentWhitelistNamespace = "TRU";

    // TSV happy path + duplicate collapse
    const std::string tsv = work + "/counts.tsv";
    writeFile(tsv,
        "barcode\tfeature_id\tcount\n"
        "AAACCCAAGAAACCAT\tHIV_DNA\t3\n"
        "AAACCCAAGAAACCAT\tHIV_DNA\t2\n"
        "GTATGTTCAGTAGCCT-1\tHIV_RNA\t4\n");
    auto r1 = PfMultiTableImport::runTableFeatureImport(
        whitelist, featureRef, tsv, assignOut, opts);
    if (r1.returnCode != 0) {
        std::cerr << "TSV import failed\n";
        return 1;
    }
    if (r1.stats.rowsRead != 3 || r1.stats.duplicatePairsCollapsed != 1) {
        std::cerr << "duplicate collapse mismatch\n";
        return 2;
    }
    if (r1.stats.rowsSuffixNormalized < 1) {
        std::cerr << "expected suffix normalization\n";
        return 3;
    }
    if (r1.stats.rowsRetained != 2) {
        std::cerr << "retained pair count mismatch\n";
        return 4;
    }
    if (r1.stats.barcodeNamespaceInput != "MIXED" ||
        r1.stats.barcodeNamespaceOutput != "UNSUFFIXED") {
        std::cerr << "barcode namespace provenance mismatch: input="
                  << r1.stats.barcodeNamespaceInput
                  << " output=" << r1.stats.barcodeNamespaceOutput << "\n";
        return 9;
    }

    // CSV + bad count + unknown feature
    const std::string csv = work + "/counts.csv";
    writeFile(csv,
        "barcode,feature_id,count\n"
        "AAACCCAAGAAACCAT,HIV_DNA,1\n"
        "AAACCCAAGAAACCAT,UNKNOWN,2\n"
        "AAACCCAAGAAACCAT,HIV_RNA,abc\n"
        "NOTABARCODE,HIV_DNA,1\n");
    const std::string assignOut2 = work + "/out2";
    system(("mkdir -p \"" + assignOut2 + "\"").c_str());
    auto r2 = PfMultiTableImport::runTableFeatureImport(
        whitelist, featureRef, csv, assignOut2, opts);
    if (r2.returnCode != 0) return 5;
    if (r2.stats.rowsRejectedFeature < 1 || r2.stats.rowsRejectedCount < 1) return 6;
    if (r2.stats.rowsRejectedBarcode < 1) return 7;

    // Missing columns
  bool threw = false;
  try {
    const std::string bad = work + "/bad.tsv";
    writeFile(bad, "barcode\tcount\nAAACCCAAGAAACCAT\t1\n");
    PfMultiTableImport::runTableFeatureImport(whitelist, featureRef, bad, assignOut, opts);
  } catch (...) {
    threw = true;
  }
  if (!threw) return 8;

  std::cout << "OK\n";
  return 0;
}
HARNESS_EOF

export TABLE_TEST_WORK="$WORK_DIR"
g++ -std=c++11 -O2 -I"$SOURCE_DIR" -I"$PF_INCLUDE" \
    "$WORK_DIR/test_table_import.cpp" \
    "$WORK_DIR/PfMultiTableImport_test.o" \
    "$SOURCE_DIR/PfMultiMexStub.o" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    "$SOURCE_DIR/MexWriter.o" \
    "$SOURCE_DIR/MexWriterUtil.o" \
    -lz -lpthread \
    -o "$HARNESS"

if "$HARNESS"; then
    pass "table import harness"
else
    fail "table import harness"
fi

if [[ -f "$WORK_DIR/out/matrix.mtx" && -f "$WORK_DIR/out/barcodes.tsv" && -f "$WORK_DIR/out/features.tsv" ]]; then
    pass "MEX artifacts written"
else
    fail "MEX artifacts missing"
fi

if [[ -f "$WORK_DIR/out/table_feature_import.api_run.txt" ]]; then
    pass "api_run telemetry written"
else
    fail "api_run telemetry missing"
fi

echo ""
echo "=========================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
[[ "$FAIL_COUNT" -eq 0 ]]

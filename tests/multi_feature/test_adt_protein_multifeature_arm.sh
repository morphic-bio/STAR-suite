#!/usr/bin/env bash
# Tiny multi-feature smoke: gRNA + Custom + ADT/protein arms share the pf-multi
# library routing path. ADT/protein libraries trigger assignBarcodes ADT MEX mode
# via PfMultiAssign/pf_api (not the standalone CLI --output-mode flag).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PF_DIR="$REPO_ROOT/core/features/process_features"
SOURCE_DIR="$REPO_ROOT/core/legacy/source"
PF_INCLUDE="$REPO_ROOT/core/features/process_features/include"
INC_FLAGS="-I${SOURCE_DIR} -I${REPO_ROOT}/core/features/libscrna/include -I${REPO_ROOT}/core/features/vbem/source -I${REPO_ROOT}/core/features/vbem/source/libem -I${REPO_ROOT}/slam/source -I${REPO_ROOT}/core/features/bamsort/source -I${PF_INCLUDE} -I${REPO_ROOT}/flex/source -I${REPO_ROOT}/flex/source/libflex -I${SOURCE_DIR}/htslib"
WORK_DIR="$(mktemp -d /tmp/multi_adt_arm_test.XXXXXX)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0
pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: multi-feature ADT/protein library arm ==="

make -C "$PF_DIR" -j8 libprocess_features.a >/dev/null
make -C "$SOURCE_DIR" -j8 PfMultiConfig.o PfMultiAssign.o ThreadControl.o input/CbqInputModule.o >/dev/null

HARNESS="$WORK_DIR/test_adt_arm"
cat > "$WORK_DIR/test_adt_arm.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;
#include "PfMultiConfig.h"
#include "PfMultiFeatureSpecs.h"

static string sanitizeDirName(const string& input) {
    string out = input;
    for (char& c : out) {
        if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
              (c >= '0' && c <= '9') || c == '-' || c == '_')) {
            c = '_';
        }
    }
    return out;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: test_adt_arm <config.csv> <out_prefix>" << endl;
        return 2;
    }
    const string outPrefix = argv[2];
    try {
        PfMultiConfig::Config cfg = PfMultiConfig::parseConfig(argv[1]);
        ostringstream log;
        vector<PfMultiFeatureSpecs::FeatureSpec> specs =
            PfMultiFeatureSpecs::buildFeatureSpecsFromConfig(cfg, log);

        for (size_t i = 0; i < specs.size(); ++i) {
            const auto& spec = specs[i];
            vector<PfMultiConfig::LibraryEntry> libs =
                cfg.getFeatureLibraries(spec.libraryType);
            for (size_t j = 0; j < libs.size(); ++j) {
                const auto& lib = libs[j];
                const string featureDir = sanitizeDirName(spec.libraryType);
                const string assignOut = outPrefix + "/cr_assign/" + featureDir + "/"
                    + sanitizeDirName(lib.starLibraryId);
                const bool adtMex = PfMultiFeatureSpecs::shouldEmitAdtMexOutput(
                    spec.featureRefType, spec.libraryType);
                cout << "LIBRARY"
                     << " id=" << lib.starLibraryId
                     << " feature_types=" << spec.libraryType
                     << " featureRefType=" << spec.featureRefType
                     << " assignOut=" << assignOut
                     << " adt_mex=" << (adtMex ? "yes" : "no")
                     << endl;
            }
        }
        return 0;
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
}
HARNESS_EOF

g++ -std=c++11 -O2 -I"$SOURCE_DIR" -I"$PF_INCLUDE" \
    "$WORK_DIR/test_adt_arm.cpp" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    -o "$HARNESS"

ASSIGN_HARNESS="$WORK_DIR/test_pf_multi_assign_adt"
cat > "$WORK_DIR/test_pf_multi_assign_adt.cpp" << 'ASSIGN_EOF'
#include "PfMultiAssign.h"
#include "PfMultiConfig.h"
#include "PfMultiFeatureSpecs.h"
#include "ThreadControl.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

ThreadControl g_threadChunks;

static std::string sanitizeDirName(const std::string& input) {
    std::string out = input;
    for (char& c : out) {
        if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
              (c >= '0' && c <= '9') || c == '-' || c == '_')) {
            c = '_';
        }
    }
    return out;
}

static std::string readKeyValue(const std::string& path, const std::string& key) {
    std::ifstream in(path.c_str());
    std::string line;
    const std::string prefix = key + "=";
    while (std::getline(in, line)) {
        if (line.compare(0, prefix.size(), prefix) == 0) {
            return line.substr(prefix.size());
        }
    }
    return "";
}

static bool fileExists(const std::string& path) {
    std::ifstream in(path.c_str());
    return in.good();
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: test_pf_multi_assign_adt <config.csv> <out_prefix> <whitelist.txt>\n";
        return 2;
    }

    try {
        PfMultiConfig::Config cfg = PfMultiConfig::parseConfig(argv[1]);
        const std::string outPrefix = argv[2];
        const std::string whitelist = argv[3];
        std::ostringstream log;
        std::vector<PfMultiFeatureSpecs::FeatureSpec> specs =
            PfMultiFeatureSpecs::buildFeatureSpecsFromConfig(cfg, log);

        int failures = 0;
        for (const auto& spec : specs) {
            std::vector<PfMultiConfig::LibraryEntry> libs =
                cfg.getFeatureLibraries(spec.libraryType);
            for (const auto& lib : libs) {
                const std::string assignOut = outPrefix + "/cr_assign/"
                    + sanitizeDirName(spec.libraryType) + "/"
                    + sanitizeDirName(lib.starLibraryId);

                PfMultiAssign::AssignOptions opts;
                opts.adtMexOutput = PfMultiFeatureSpecs::shouldEmitAdtMexOutput(
                    spec.featureRefType, spec.libraryType);
                opts.skipQcOutputs = true;
                opts.sampleName = lib.sample.empty() ? lib.starLibraryId : lib.sample;

                PfMultiAssign::AssignResult result = PfMultiAssign::runAssignBarcodes(
                    whitelist, lib.starFeatureRef, lib.fastqs, assignOut, opts);

                const std::string apiRun = assignOut + "/assignBarcodes.api_run.txt";
                const std::string adtFlag = readKeyValue(apiRun, "adtMexOutput");
                const std::string outputMode = readKeyValue(apiRun, "output_mode");

                std::cout << "ASSIGN_RESULT"
                          << " id=" << lib.starLibraryId
                          << " adt_mex_expected=" << (opts.adtMexOutput ? "yes" : "no")
                          << " api_adtMexOutput=" << adtFlag
                          << " api_output_mode=" << outputMode
                          << " rc=" << result.returnCode
                          << " api_run=" << apiRun
                          << "\n";

                if (result.returnCode != 0) {
                    std::cerr << "ERROR: assign failed for " << lib.starLibraryId << "\n";
                    failures++;
                    continue;
                }
                if (opts.adtMexOutput) {
                    if (adtFlag != "1" || outputMode != "adt_mex") {
                        std::cerr << "ERROR: ADT library missing pf_api ADT MEX flags in "
                                  << apiRun << "\n";
                        failures++;
                    }
                } else if (adtFlag != "0" || outputMode != "default") {
                    std::cerr << "ERROR: non-ADT library should not enable ADT MEX in "
                              << apiRun << "\n";
                    failures++;
                }
            }
        }
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
ASSIGN_EOF

g++ -std=c++11 -O2 -fopenmp ${INC_FLAGS} \
    "$WORK_DIR/test_pf_multi_assign_adt.cpp" \
    "$SOURCE_DIR/PfMultiAssign.o" \
    "$SOURCE_DIR/input/CbqInputModule.o" \
    "$SOURCE_DIR/ThreadControl.o" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    -L"$PF_DIR" -lprocess_features \
    -L"$REPO_ROOT/core/features/libscrna" -lscrna \
    -lpthread -lz -lstdc++ -lhts -lglib-2.0 \
    -o "$ASSIGN_HARNESS"

touch "$WORK_DIR/ref_grna.csv" "$WORK_DIR/ref_larry.csv" "$WORK_DIR/ref_adt.csv"

cat > "$WORK_DIR/multi_config.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref,star_library_id
$WORK_DIR/fastq_grna,PolyIII,CRISPR Guide Capture,$WORK_DIR/ref_grna.csv,grna_de
$WORK_DIR/fastq_larry,LARRY,Custom,$WORK_DIR/ref_larry.csv,larry_de
$WORK_DIR/fastq_adt,ADT,Protein,$WORK_DIR/ref_adt.csv,adt_de
EOF

ROUTING="$("$HARNESS" "$WORK_DIR/multi_config.csv" "$WORK_DIR/out" 2>&1)"
echo "$ROUTING"

if echo "$ROUTING" | grep -q 'id=grna_de .*adt_mex=no'; then
    pass "gRNA library does not request ADT MEX"
else
    fail "gRNA library should not request ADT MEX"
fi

if echo "$ROUTING" | grep -q 'id=larry_de .*adt_mex=no'; then
    pass "Custom/LARRY library does not request ADT MEX"
else
    fail "Custom library should not request ADT MEX"
fi

if echo "$ROUTING" | grep -q 'id=adt_de .*featureRefType=Antibody Capture .*adt_mex=yes'; then
    pass "Protein feature_types maps to Antibody Capture + ADT MEX"
else
    fail "Protein library should map to Antibody Capture and ADT MEX"
fi

GRNA_OUT="$WORK_DIR/out/cr_assign/CRISPR_Guide_Capture/grna_de"
LARRY_OUT="$WORK_DIR/out/cr_assign/Custom/larry_de"
ADT_OUT="$WORK_DIR/out/cr_assign/Protein/adt_de"

for path in "$GRNA_OUT" "$LARRY_OUT" "$ADT_OUT"; do
    if echo "$ROUTING" | grep -q "assignOut=$path"; then
        pass "unique assignOut path $path"
    else
        fail "missing expected assignOut path $path"
    fi
done

# Shared whitelist
cat > "$WORK_DIR/whitelist.txt" <<'EOF'
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCT
EOF

make_fastqs() {
    local dir="$1"
    local feat_seq="$2"
    mkdir -p "$dir"
    python3 - <<PY "$dir" "$feat_seq"
import os, sys
out_dir, feat = sys.argv[1], sys.argv[2]
barcodes = [
    "AAACCCAAGAAACCAT",
    "AAACCCAAGAAACCCA",
    "AAACCCAAGAAACCCT",
]
reads = [
    (0, feat, "AAAAAAAAAAAA"),
    (1, feat, "CCCCCCCCCCCC"),
]
qual = lambda n: "I" * n
with open(os.path.join(out_dir, "sample_R1_001.fastq"), "w") as r1, \
     open(os.path.join(out_dir, "sample_R2_001.fastq"), "w") as r2:
    for i, (bc_idx, seq, umi) in enumerate(reads):
        bc = barcodes[bc_idx]
        r1.write(f"@read{i}\n{bc}{umi}\n+\n{qual(len(bc)+len(umi))}\n")
        r2.write(f"@read{i}\n{seq}\n+\n{qual(len(seq))}\n")
PY
}

cat > "$WORK_DIR/ref_grna.csv" <<'EOF'
id,name,sequence,feature_type
g1,Guide-A,ATCGATCG,CRISPR Guide Capture
EOF

cat > "$WORK_DIR/ref_larry.csv" <<'EOF'
id,name,sequence,feature_type
lb1,LARRY-1,TTAATTAATTAATTAA,Custom
EOF

cat > "$WORK_DIR/ref_adt.csv" <<'EOF'
id,name,sequence,feature_type
CD29,CD29,ATCGATCGATCGATCG,ADT
EOF

make_fastqs "$WORK_DIR/fastq_grna" "ATCGATCG"
make_fastqs "$WORK_DIR/fastq_larry" "TTAATTAATTAATTAA"
make_fastqs "$WORK_DIR/fastq_adt" "ATCGATCGATCGATCG"

echo ""
echo "--- PfMultiAssign API path (feature spec -> adtMexOutput -> pf_api) ---"
ASSIGN_OUTPUT="$("$ASSIGN_HARNESS" "$WORK_DIR/multi_config.csv" "$WORK_DIR/out" "$WORK_DIR/whitelist.txt" 2>&1)" || {
    echo "$ASSIGN_OUTPUT"
    fail "PfMultiAssign harness failed"
    echo ""
    echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
    exit 1
}
echo "$ASSIGN_OUTPUT"

if echo "$ASSIGN_OUTPUT" | grep -q 'ASSIGN_RESULT id=grna_de .*adt_mex_expected=no api_adtMexOutput=0 api_output_mode=default rc=0'; then
    pass "gRNA PfMultiAssign path leaves ADT MEX off (api_run)"
else
    fail "gRNA PfMultiAssign path should not enable ADT MEX"
fi

if echo "$ASSIGN_OUTPUT" | grep -q 'ASSIGN_RESULT id=larry_de .*adt_mex_expected=no api_adtMexOutput=0 api_output_mode=default rc=0'; then
    pass "Custom PfMultiAssign path leaves ADT MEX off (api_run)"
else
    fail "Custom PfMultiAssign path should not enable ADT MEX"
fi

if echo "$ASSIGN_OUTPUT" | grep -q 'ASSIGN_RESULT id=adt_de .*adt_mex_expected=yes api_adtMexOutput=1 api_output_mode=adt_mex rc=0'; then
    pass "Protein feature spec alone enables ADT MEX via PfMultiAssign/pf_api"
else
    fail "Protein library should enable ADT MEX from feature spec without CLI flag"
fi

GRNA_API_RUN="$GRNA_OUT/assignBarcodes.api_run.txt"
ADT_API_RUN="$ADT_OUT/assignBarcodes.api_run.txt"
if [[ -f "$GRNA_API_RUN" ]] && grep -q '^output_mode=default$' "$GRNA_API_RUN" && grep -q '^adtMexOutput=0$' "$GRNA_API_RUN"; then
    pass "gRNA assignBarcodes.api_run.txt records default output mode"
else
    fail "gRNA assignBarcodes.api_run.txt missing expected pf_api flags"
fi

if [[ -f "$ADT_API_RUN" ]] && grep -q '^output_mode=adt_mex$' "$ADT_API_RUN" && grep -q '^adtMexOutput=1$' "$ADT_API_RUN"; then
    pass "ADT assignBarcodes.api_run.txt records adt_mex output mode"
else
    fail "ADT assignBarcodes.api_run.txt missing expected pf_api flags"
fi

find_sample_dir() {
    local root="$1"
    if [[ -f "$root/barcodes.tsv.gz" || -f "$root/matrix.mtx" ]]; then
        echo "$root"
        return
    fi
    for sub in "$root"/*; do
        if [[ -d "$sub" ]]; then
            find_sample_dir "$sub"
            return
        fi
    done
}

ADT_SAMPLE="$(find_sample_dir "$ADT_OUT")"
if [[ -f "$ADT_SAMPLE/protein_quant_summary.json" && -f "$ADT_SAMPLE/features.tsv.gz" ]]; then
    pass "ADT assign output has protein MEX sidecars"
else
    fail "ADT assign output missing protein MEX sidecars"
fi

if zgrep -q "Antibody Capture" "$ADT_SAMPLE/features.tsv.gz"; then
    pass "ADT features.tsv.gz uses Antibody Capture"
else
    fail "ADT features.tsv.gz missing Antibody Capture"
fi

GRNA_SAMPLE="$(find_sample_dir "$GRNA_OUT")"
LARRY_SAMPLE="$(find_sample_dir "$LARRY_OUT")"

if [[ -f "$GRNA_SAMPLE/matrix.mtx" ]]; then
    pass "gRNA legacy matrix.mtx present"
else
    fail "gRNA matrix.mtx missing"
fi

if [[ -f "$LARRY_SAMPLE/matrix.mtx" ]]; then
    pass "Custom/LARRY legacy matrix.mtx present"
else
    fail "LARRY matrix.mtx missing"
fi

if [[ ! -f "$GRNA_SAMPLE/protein_quant_summary.json" && ! -f "$LARRY_SAMPLE/protein_quant_summary.json" ]]; then
    pass "gRNA and LARRY outputs omit protein MEX sidecars"
else
    fail "non-ADT libraries should not emit protein_quant_summary.json"
fi

echo ""
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
if [[ "$FAIL_COUNT" -ne 0 ]]; then
    exit 1
fi

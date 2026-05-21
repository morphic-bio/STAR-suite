#include "OcmMultiConfig.h"
#include "OcmMultiMaterialize.h"
#include "PfMultiConfig.h"
#include "PfMultiMerge.h"
#include "VelocytoMexWriter.h"
#include "Parameters.h"
#include "InOutStreams.h"
#include "SoloFeature.h"

#include <climits>
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static int failures = 0;

static int check(bool ok, const std::string& label) {
    if (!ok) {
        std::cerr << "FAIL: " << label << "\n";
        failures++;
        return 1;
    }
    return 0;
}

static std::string fixtureRoot() {
    const char* env = std::getenv("OCM_TEST_FIXTURE_ROOT");
    if (env != nullptr && env[0] != '\0') {
        return env;
    }
    return "tests/fixtures/ocm_multi_tiny";
}

static int testBarcodeClassifier() {
    if (check(OcmMultiMaterialize::classifyBarcodeTag("AAACCCTGTAAGCGCG-1") == "OB1",
              "OB1 overhang GT") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::classifyBarcodeTag("AAACCCGCAACTAGAC-1") == "OB2",
              "OB2 overhang CA") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::classifyBarcodeTag("AAACCATTCACCTGGG-1") == "OB3",
              "OB3 overhang TC") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::classifyBarcodeTag("AAACCAAAGCATTGAT-1") == "OB4",
              "OB4 overhang AG") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::classifyBarcodeTag("AAAAAAAAAAAAAAAA-1").empty(),
              "unknown overhang rejected") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::tag8ForOcmId("OB1") == "GTGTGTGT",
              "OB1 TAG8 mapping") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::appendFlexTag8("AAACCCTGTAAGCGCG") ==
                  "AAACCCTGTAAGCGCGGTGTGTGT",
              "append OCM TAG8 from CB16 overhang") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::appendFlexTag8("AAAAAAAAAAAAAAAA").empty(),
              "do not append TAG8 to unknown OCM overhang") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::classifyBarcodeTag("AAAAAAAAAAAAAAAAGTGTGTGT-1") == "OB1",
              "classify OCM TAG8 suffix before CB16 overhang") != 0) {
        return 1;
    }
    if (check(OcmMultiMaterialize::stripFlexTag8("AAACCCTGTAAGCGCGGTGTGTGT-1") ==
                  "AAACCCTGTAAGCGCG-1",
              "strip OCM TAG8 while preserving GEM suffix") != 0) {
        return 1;
    }
    return 0;
}

static int testSampleIdValidation() {
    if (check(OcmMultiConfig::isPathSafeSampleId("GCM1-Day-4"), "JAX sample id accepted") != 0) {
        return 1;
    }
    if (check(OcmMultiConfig::isPathSafeSampleId("WT-PrS-20pct-Day-4"), "percent in sample id") != 0) {
        return 1;
    }
    if (check(!OcmMultiConfig::isPathSafeSampleId("../escape"), "reject parent traversal") != 0) {
        return 1;
    }
    if (check(!OcmMultiConfig::isPathSafeSampleId("sample/sub"), "reject slash") != 0) {
        return 1;
    }
    if (check(!OcmMultiConfig::isPathSafeSampleId(".hidden"), "reject leading dot") != 0) {
        return 1;
    }

    const std::string root = fixtureRoot();
    const std::string badConfig = root + "/config_bad_sample_id.csv";
    {
        std::ofstream out(badConfig.c_str());
        out << "[samples]\n"
            << "sample_id,ocm_barcode_ids,description\n"
            << "../bad,OB1,test\n";
    }
    std::ostringstream log;
    bool threw = false;
    try {
        (void)OcmMultiConfig::parseAndValidate(badConfig, log);
    } catch (const std::exception&) {
        threw = true;
    }
    std::remove(badConfig.c_str());
    if (check(threw, "reject path-unsafe sample_id in config") != 0) {
        return 1;
    }
    return 0;
}

static int testConfigParser() {
    const std::string root = fixtureRoot();
    const std::string configPath = root + "/config.csv";
    std::ostringstream log;
    PfMultiConfig::Config config = OcmMultiConfig::parseAndValidate(configPath, log);
    if (check(config.samples.size() == 5, "expected five configured samples") != 0) {
        return 1;
    }
    if (check(config.samples[0].sample_id == "GCM1-Day-4", "first sample id") != 0) {
        return 1;
    }
    const std::vector<std::string> unionIds = config.samples[4].resolvedOcmIds();
    if (check(unionIds.size() == 2 && unionIds[0] == "OB1" && unionIds[1] == "OB2",
              "pipe-union sample tags") != 0) {
        return 1;
    }
    return 0;
}

static bool externalMaterializeRunRequested() {
    const char* runEnv = std::getenv("OCM_TEST_RUN_DIR");
    const char* configEnv = std::getenv("OCM_TEST_CONFIG");
    return runEnv != nullptr && runEnv[0] != '\0' && configEnv != nullptr && configEnv[0] != '\0';
}

static std::string materializeConfigPath(const std::string& root) {
    if (externalMaterializeRunRequested()) {
        return std::getenv("OCM_TEST_CONFIG");
    }
    return root + "/config.csv";
}

static std::string materializeRunDir(const std::string& root) {
    if (externalMaterializeRunRequested()) {
        return std::getenv("OCM_TEST_RUN_DIR");
    }
    return root + "/run";
}

static std::string downstreamSamplesRoot(const std::string& root, const std::string& runDir) {
    if (!externalMaterializeRunRequested()) {
        std::string projectRoot = root;
        if (!projectRoot.empty() && projectRoot.back() != '/') {
            projectRoot.push_back('/');
        }
        return projectRoot;
    }
    std::string prefix = runDir;
    if (!prefix.empty() && prefix.back() != '/') {
        prefix.push_back('/');
    }
    const size_t runPos = prefix.rfind("/run/");
    if (runPos != std::string::npos) {
        const size_t samplesPos = prefix.rfind("/samples/", runPos);
        if (samplesPos != std::string::npos) {
            return prefix.substr(0, samplesPos + 1);
        }
        if (runPos == 0) {
            return "./";
        }
        return prefix.substr(0, runPos + 1);
    }
    const size_t slash = prefix.find_last_of('/');
    return (slash == std::string::npos) ? "./" : prefix.substr(0, slash + 1);
}

static int testMaterializerTiny() {
    const std::string root = fixtureRoot();
    const std::string runDir = materializeRunDir(root);
    const std::string outsDir = downstreamSamplesRoot(root, runDir) + "outs";

    Parameters P;
    P.inOut = new InOutStreams();
    const char* logPath = std::getenv("OCM_TEST_LOG");
    P.inOut->logMain.open(logPath != nullptr && logPath[0] != '\0' ? logPath : "/dev/null");
    P.runDirPerm = static_cast<mode_t>(0775);
    P.outFileNamePrefix = runDir + "/";
    P.pfMulti.ocmMultiEnable = "yes";
    P.pfMulti.ocmMultiConfig = materializeConfigPath(root);
    P.pSolo.crGexFeature = ParametersSolo::CrGexGeneFull;
    P.pSolo.featureYes.fill(false);
    P.pSolo.featureYes[SoloFeatureTypes::Velocyto] = true;

    const int ret = runOcmMultiMaterialize(P);
    delete P.inOut;
    P.inOut = nullptr;
    if (check(ret == 0, "runOcmMultiMaterialize fixture exit code") != 0) {
        return 1;
    }
    if (check(std::ifstream((outsDir + "/multi/multiplexing_analysis/cells_per_tag.json").c_str()).good(),
              "cells_per_tag.json exists") != 0) {
        return 1;
    }
    if (check(std::ifstream((outsDir + "/multi/count/raw_feature_bc_matrix/matrix.mtx.gz").c_str()).good(),
              "multi raw matrix.mtx.gz exists") != 0) {
        return 1;
    }
    if (check(std::ifstream((outsDir + "/per_sample_outs/GCM1-Day-4/count/sample_filtered_barcodes.csv").c_str()).good(),
              "sample_filtered_barcodes.csv exists") != 0) {
        return 1;
    }
    const bool haveVelocyto =
        std::ifstream((runDir + "/Solo.out/Velocyto/raw/barcodes.tsv").c_str()).good();
    if (haveVelocyto) {
        const std::string projectRoot = downstreamSamplesRoot(root, runDir);
        const std::string sampleVeloRaw = projectRoot + "samples/GCM1-Day-4/run/outs/raw_velocyto_feature_bc_matrix/spliced.mtx.gz";
        const std::string sampleVeloFiltered =
            projectRoot + "samples/GCM1-Day-4/run/outs/filtered_velocyto_feature_bc_matrix/ambiguous.mtx.gz";
        struct stat st {};
        if (check(stat(sampleVeloRaw.c_str(), &st) == 0,
                  "downstream per-sample raw_velocyto spliced.mtx.gz") != 0) {
            return 1;
        }
        if (check(stat(sampleVeloFiltered.c_str(), &st) == 0,
                  "downstream per-sample filtered_velocyto ambiguous.mtx.gz") != 0) {
            return 1;
        }
    }
    return 0;
}

static bool barcodeHasGemSuffix(const std::string& bc) {
    const size_t dashPos = bc.find_last_of('-');
    if (dashPos == std::string::npos || dashPos == bc.size() - 1) {
        return false;
    }
    for (size_t i = dashPos + 1; i < bc.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(bc[i]))) {
            return false;
        }
    }
    return true;
}

static int testCrBarcodeSuffixNormalization() {
    const std::vector<std::string> unsuffixed = {"AAACCCTGTAAGCGCG", "AAACCCGCAACTAGAC"};
    const std::vector<uint32_t> colIndices = {0, 1};
    std::ostringstream log;
    const PfMultiMerge::CrBarcodeLayout layout =
        PfMultiMerge::buildCrBarcodeLayoutForColumns(unsuffixed,
                                                     colIndices,
                                                     "1",
                                                     "TRU",
                                                     "TRU",
                                                     log);
    if (check(layout.sortedBarcodes.size() == 2, "CR layout keeps both barcodes") != 0) {
        return 1;
    }
    for (const std::string& bc : layout.sortedBarcodes) {
        if (check(barcodeHasGemSuffix(bc), "CR layout adds GEM suffix: " + bc) != 0) {
            return 1;
        }
    }
    if (check(layout.sortedBarcodes[0] < layout.sortedBarcodes[1], "CR layout sorts lexicographically") != 0) {
        return 1;
    }
    if (check(layout.sourceColToSorted[0] == 1 && layout.sourceColToSorted[1] == 0,
              "CR layout remap matches sorted order") != 0) {
        return 1;
    }

    const std::string root = fixtureRoot();
    const std::string srcDir = root + "/run/Solo.out/GeneFull/raw";
    const std::string tmpDir = root + "/run/_cr_suffix_test_raw";
    const std::string outDir = root + "/run/_cr_suffix_test_out";
    std::string cmd = "rm -rf \"" + tmpDir + "\" \"" + outDir + "\" && mkdir -p \"" + tmpDir +
                      "\" && cp \"" + srcDir + "/features.tsv\" \"" + srcDir + "/matrix.mtx\" \"" +
                      tmpDir + "/\"";
    if (system(cmd.c_str()) != 0) {
        return 1;
    }
    {
        std::ifstream in((srcDir + "/barcodes.tsv").c_str());
        std::ofstream out((tmpDir + "/barcodes.tsv").c_str());
        std::string line;
        while (std::getline(in, line)) {
            const size_t dash = line.find_last_of('-');
            if (dash != std::string::npos) {
                line = line.substr(0, dash);
            }
            out << line << "\n";
        }
    }
    const PfMultiMerge::MexAxes axes = PfMultiMerge::readMexAxes(tmpDir);
    std::vector<uint32_t> allCols(axes.barcodes.size());
    for (size_t i = 0; i < allCols.size(); ++i) {
        allCols[i] = static_cast<uint32_t>(i);
    }
    const PfMultiMerge::CrBarcodeLayout fileLayout =
        PfMultiMerge::buildCrBarcodeLayoutForColumns(axes.barcodes,
                                                     allCols,
                                                     "1",
                                                     "TRU",
                                                     "TRU",
                                                     log);
    Parameters P;
    P.runDirPerm = static_cast<mode_t>(0775);
    if (PfMultiMerge::writeColumnSubsetMexGz(tmpDir,
                                             outDir,
                                             axes,
                                             fileLayout,
                                             P,
                                             log) != 0) {
        check(false, "writeColumnSubsetMexGz on unsuffixed fixture");
        return 1;
    }
    const PfMultiMerge::MexAxes outAxes = PfMultiMerge::readMexAxes(outDir);
    for (const std::string& bc : outAxes.barcodes) {
        if (check(barcodeHasGemSuffix(bc), "streamed MEX barcodes suffixed: " + bc) != 0) {
            return 1;
        }
    }
    cmd = "rm -rf \"" + tmpDir + "\" \"" + outDir + "\"";
    (void)system(cmd.c_str());
    return 0;
}

static int testVelocytoRunWriterTiny() {
    const std::string root = fixtureRoot();
    const std::string runDir = root + "/run";
    const std::string outsDir = runDir + "/outs";
    std::ofstream log((runDir + "/velocyto_writer_test.log").c_str());

    const int ret = VelocytoMexWriter::materializeRunVelocytoMex(runDir, outsDir, log);
    if (check(ret == 0, "materializeRunVelocytoMex fixture exit code") != 0) {
        return 1;
    }
    if (check(std::ifstream((outsDir + "/velocyto_feature_bc_matrix_manifest.json").c_str()).good(),
              "Velocyto manifest exists") != 0) {
        return 1;
    }
    if (check(std::ifstream((outsDir + "/raw_velocyto_feature_bc_matrix/spliced.mtx.gz").c_str()).good(),
              "raw Velocyto spliced exists") != 0) {
        return 1;
    }
    if (check(std::ifstream((outsDir + "/filtered_velocyto_feature_bc_matrix/unspliced.mtx.gz").c_str()).good(),
              "filtered Velocyto unspliced exists") != 0) {
        return 1;
    }
    return 0;
}

int main(int argc, char** argv) {
    const std::string mode = (argc > 1) ? argv[1] : "all";
    if (mode == "classifier") {
        return testBarcodeClassifier();
    }
    if (mode == "config") {
        return testConfigParser();
    }
    if (mode == "velocyto") {
        return testVelocytoRunWriterTiny();
    }
    if (mode == "materialize") {
        return testMaterializerTiny();
    }
    if (mode == "sample_id") {
        return testSampleIdValidation();
    }
    if (mode == "cr_barcode") {
        return testCrBarcodeSuffixNormalization();
    }
    if (testBarcodeClassifier() != 0) {
        return 1;
    }
    if (testConfigParser() != 0) {
        return 1;
    }
    if (testSampleIdValidation() != 0) {
        return 1;
    }
    if (testCrBarcodeSuffixNormalization() != 0) {
        return 1;
    }
    if (testVelocytoRunWriterTiny() != 0) {
        return 1;
    }
    if (testMaterializerTiny() != 0) {
        return 1;
    }
    std::cout << "PASS: ocm multi unit tests\n";
    return failures == 0 ? 0 : 1;
}

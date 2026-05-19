#include "OcmMultiConfig.h"
#include "OcmMultiMaterialize.h"
#include "PfMultiConfig.h"
#include "VelocytoMexWriter.h"
#include "Parameters.h"
#include "InOutStreams.h"
#include "SoloFeature.h"

#include <cstdlib>
#include <fstream>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

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
    const std::string outsDir = runDir + "/outs";

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
    if (testBarcodeClassifier() != 0) {
        return 1;
    }
    if (testConfigParser() != 0) {
        return 1;
    }
    if (testSampleIdValidation() != 0) {
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

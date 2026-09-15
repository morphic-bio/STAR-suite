#include "SpatialGex.h"
#include "SpatialR1Decoder.h"
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

static void require(bool value, const std::string &message) {
    if (!value) throw std::runtime_error(message);
}

int main(int argc, char **argv) {
    require(argc == 6, "usage: test_spatial_flex_transactions CONTRACT BC1 BC2 R1_FASTQ FRESH_OUTPUT_ROOT");
    spatial_gex::PipelineConfig config;
    config.barcodeContractDirectory = argv[1];
    config.bc1OligosPath = argv[2];
    config.bc2OligosPath = argv[3];
    config.sourceRevision = STAR_SUITE_SOURCE_REVISION;
    config.starSuiteVersion = "transaction-test";
    config.expectedReads = 10000;
    config.expectedCandidates = 100000;
    config.threads = 1;
    config.featureCount = 2;
    config.flexFeatureMode = true;
    config.requirePairedCompletion = true;
    config.outputDirectory = std::string(argv[5]) + "/paired";
    config.products = spatial_gex::ProductAll;
    config.scales = spatial_gex::ScaleAll;
    std::vector<std::string> reads;
    std::ifstream input(argv[4]);
    std::string header, sequence, plus, quality;
    while (std::getline(input, header)) {
        require(bool(std::getline(input, sequence)) && bool(std::getline(input, plus))
                && bool(std::getline(input, quality)), "truncated test input");
        reads.push_back(sequence);
    }
    require(!reads.empty() && reads.size() < 9000, "small nonempty fixture required");
    // Preserve a candidate-bearing read while making its raw UMI invalid.
    std::string invalidUmi = reads.front();
    invalidUmi[0] = 'N';
    reads.push_back(invalidUmi);
    std::string error;
    auto paired = spatial_gex::Pipeline::create(config, error);
    require(bool(paired), error);
    using C = spatial_gex::FeatureEvidenceClass;
    require(!paired->completeCurrentThread(C::FlexH0, true, 0, 0, error), "completion without R1 accepted");
    quality.assign(reads.front().size(), 'I');
    require(!paired->decodeCurrentThread(reads.front().data(), reads.front().size(), quality.data(), quality.size(),
                                        uint64_t(UINT32_MAX) + 1, error), "overflowing source identity accepted");
    uint64_t assigned = 0, deny = 0, miss = 0;
    auto source = [](size_t i) {
        return i % 11 == 0 ? C::FlexUnassigned : i % 7 == 0 ? C::FlexHashDeny
            : i % 3 == 0 ? C::FlexH0 : i % 3 == 1 ? C::FlexH1 : C::FlexH1X2;
    };
    for (size_t i = 0; i < reads.size(); ++i) {
        quality.assign(reads[i].size(), 'I');
        const auto ordinal = uint64_t(i) * 3 + 1;
        require(paired->decodeCurrentThread(reads[i].data(), reads[i].size(), quality.data(), quality.size(), ordinal, error), error);
        const C cls = source(i);
        const bool keep = cls != C::FlexUnassigned && cls != C::FlexHashDeny;
        if (i == 1) {
            require(!paired->decodeCurrentThread(reads[i].data(), reads[i].size(), quality.data(), quality.size(), ordinal, error), "second decode overwrote pending R1");
            require(!paired->completeCurrentThread(cls, keep, i % 2, ordinal + 1, error), "mismatched identity accepted");
            require(!paired->completeCurrentThread(C::FlexHashDeny, true, 0, ordinal, error), "deny with feature accepted");
            require(!paired->completeCurrentThread(C::FlexH1X2, false, 0, ordinal, error), "keep without feature accepted");
            require(!paired->completeCurrentThread(C::FlexH0, true, 2, ordinal, error), "out-of-axis feature accepted");
        }
        require(paired->completeCurrentThread(cls, keep, i % 2, ordinal, error), error);
        require(!paired->completeCurrentThread(cls, keep, i % 2, ordinal, error), "duplicate completion accepted");
        assigned += keep;
        deny += cls == C::FlexHashDeny;
        miss += cls == C::FlexUnassigned;
    }
    const std::vector<std::string> genes {"ENSG_TEST_1", "ENSG_TEST_2"};
    require(paired->finalize(genes, error), error);
    const auto expected = paired->summary();
    require(expected.readsDecoded == reads.size(), "failed transactions altered decoded count");
    require(expected.featureAssignedReads == assigned && expected.flexHashDenyReads == deny
            && expected.flexUnassignedReads == miss, "terminal class accounting differs");
    require(!paired->finalize(genes, error), "double finalization accepted");
    paired.reset();

    // Feed exactly the same fixed evidence through the pre-existing indexed API.
    // This tests the new transaction lifecycle independently of read ingestion.
    config.flexFeatureMode = false;
    config.requirePairedCompletion = false;
    config.outputDirectory = std::string(argv[5]) + "/indexed";
    auto indexed = spatial_gex::Pipeline::create(config, error);
    require(bool(indexed), error);
    size_t maximumFamily = 0;
    uint64_t invalidUmis = 0;
    for (size_t i = 0; i < reads.size(); ++i) {
        spatial_r1_decoder::Result decoded;
        quality.assign(reads[i].size(), 'I');
        require(indexed->decode(0, reads[i].data(), reads[i].size(), quality.data(), quality.size(), decoded, error), error);
        maximumFamily = std::max(maximumFamily, decoded.candidates.size());
        invalidUmis += !decoded.rawUmiValid;
        if (source(i) != C::FlexUnassigned && source(i) != C::FlexHashDeny)
            require(indexed->append(0, i % 2, uint64_t(i) * 3 + 1, decoded, error), error);
    }
    require(indexed->finalize(genes, error), error);
    require(indexed->summary().joinedReads == expected.joinedReads
            && indexed->summary().candidateRows == expected.candidateRows
            && indexed->summary().hardMolecules == expected.hardMolecules,
            "paired completion altered fixed downstream evidence");
    require(invalidUmis > 0, "invalid-UMI control was not exercised");
    require(maximumFamily > 4, "fixture lacks a spatial family beyond Chromium's four-candidate cap");
    indexed.reset();

    config.flexFeatureMode = true;
    config.requirePairedCompletion = true;
    config.outputDirectory = std::string(argv[5]) + "/pending";
    auto pending = spatial_gex::Pipeline::create(config, error);
    require(bool(pending), error);
    quality.assign(reads.front().size(), 'I');
    require(pending->decodeCurrentThread(reads.front().data(), reads.front().size(), quality.data(), quality.size(), 0, error), error);
    require(!pending->finalize(genes, error) && error.find("without a terminal") != std::string::npos,
            "missing terminal completion finalized successfully");
    std::cout << "PASS: paired/indexed fixed-evidence parity; invalid, duplicate, missing and overflowing transactions; maximum family="
              << maximumFamily << "; invalid UMIs=" << invalidUmis << '\n';
}

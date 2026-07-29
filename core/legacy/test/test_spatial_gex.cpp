#include "SpatialGex.h"
#include "SpatialR1Decoder.h"

#include <cassert>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

void writeLines(const std::string &path, const std::vector<std::string> &values)
{
    std::ofstream output(path.c_str());
    assert(output.good());
    for (std::vector<std::string>::const_iterator value = values.begin();
         value != values.end(); ++value) {
        output << *value << '\n';
    }
    output.close();
    assert(output.good());
}

spatial_r1_decoder::Result decode(
    spatial_r1_decoder::Decoder &decoder, const std::string &sequence)
{
    spatial_r1_decoder::Result result;
    spatial_r1_decoder::ExactH0Counts h0;
    h0.reset(1, 1);
    std::string error;
    const std::string quality(sequence.size(), 'I');
    assert(decoder.decode(sequence.data(), sequence.size(), quality.data(),
                          quality.size(), result, &h0, error));
    assert(error.empty());
    return result;
}

void testDecoderAmbiguousBases()
{
    char directoryTemplate[] = "/tmp/star_spatial_gex_decoder_XXXXXX";
    char *directory = ::mkdtemp(directoryTemplate);
    assert(directory != NULL);
    const std::string root(directory);
    const std::string bc1Path = root + "/bc1.txt";
    const std::string bc2Path = root + "/bc2.txt";
    const std::string ambiguousBc1Path = root + "/bc1_ambiguous.txt";
    const std::string bc1 = "ACGTTGCAACGTTGC";
    const std::string bc2 = "TGCAACGTTGCAAC";
    writeLines(bc1Path, std::vector<std::string>(1, bc1));
    writeLines(bc2Path, std::vector<std::string>(1, bc2));
    writeLines(ambiguousBc1Path,
               std::vector<std::string>{bc1, "ACGGTGCAACGTTGC"});

    spatial_r1_decoder::Config config;
    config.bc1OligosPath = bc1Path;
    config.bc2OligosPath = bc2Path;
    config.gridRows = 1;
    config.gridColumns = 1;
    config.fullStartMin = 9;
    config.fullStartMax = 9;
    spatial_r1_decoder::Decoder decoder(config);

    const std::string umi = "ACGTACGTA";
    const std::string exactSequence = umi + bc1 + bc2 + "AAAAA";
    const spatial_r1_decoder::Result exact = decode(decoder, exactSequence);
    assert(exact.rawUmiValid);
    assert(!exact.rawUmiHadN);
    assert(!exact.barcodeHadN);
    assert(!exact.barcodeDpChecked);
    assert(exact.decoderAssigned);
    assert(exact.candidates.size() == 1);

    std::string bc1N = exactSequence;
    bc1N[12] = 'N';
    const spatial_r1_decoder::Result oneN = decode(decoder, bc1N);
    assert(oneN.rawUmiValid);
    assert(oneN.barcodeHadN && oneN.barcodeNCount == 1);
    assert(oneN.barcodeDpChecked);
    assert(oneN.decoderAssigned);
    assert(oneN.candidates.size() == 1);
    assert(oneN.bc1Edit == 1 && oneN.bc2Edit == 0);

    std::string bothN = exactSequence;
    bothN[12] = 'N';
    bothN[9 + bc1.size() + 3] = 'N';
    const spatial_r1_decoder::Result twoN = decode(decoder, bothN);
    assert(twoN.rawUmiValid);
    assert(twoN.barcodeHadN && twoN.barcodeNCount == 2);
    assert(twoN.barcodeDpChecked);
    assert(twoN.decoderAssigned);
    assert(twoN.candidates.size() == 1);
    assert(twoN.bc1Edit == 1 && twoN.bc2Edit == 1);

    std::string noMatch = exactSequence;
    noMatch.replace(9, bc1.size() + bc2.size(), bc1.size() + bc2.size(), 'N');
    const spatial_r1_decoder::Result unmatched = decode(decoder, noMatch);
    assert(unmatched.barcodeHadN && unmatched.barcodeDpChecked);
    assert(!unmatched.decoderAssigned);
    assert(unmatched.candidates.empty());

    std::string umiN = exactSequence;
    umiN[3] = 'N';
    const spatial_r1_decoder::Result invalidUmi = decode(decoder, umiN);
    assert(!invalidUmi.rawUmiValid);
    assert(invalidUmi.rawUmiHadN);
    assert(invalidUmi.decoderAssigned);
    assert(invalidUmi.candidates.size() == 1);

    std::string unsupported = exactSequence;
    unsupported[12] = 'X';
    const spatial_r1_decoder::Result invalidBarcode = decode(decoder, unsupported);
    assert(invalidBarcode.barcodeHadUnsupportedBase);
    assert(!invalidBarcode.barcodeDpChecked);
    assert(!invalidBarcode.decoderAssigned);
    assert(invalidBarcode.candidates.empty());

    spatial_r1_decoder::Config ambiguousConfig;
    ambiguousConfig.bc1OligosPath = ambiguousBc1Path;
    ambiguousConfig.bc2OligosPath = bc2Path;
    ambiguousConfig.gridRows = 1;
    ambiguousConfig.gridColumns = 2;
    ambiguousConfig.fullStartMin = 9;
    ambiguousConfig.fullStartMax = 9;
    spatial_r1_decoder::Decoder ambiguousDecoder(ambiguousConfig);
    const spatial_r1_decoder::Result ambiguous = decode(ambiguousDecoder, bc1N);
    assert(ambiguous.barcodeHadN && ambiguous.barcodeDpChecked);
    assert(!ambiguous.decoderAssigned);
    assert(ambiguous.candidates.size() == 2);
    assert(ambiguous.candidates[0].coordinateIndex == 0);
    assert(ambiguous.candidates[1].coordinateIndex == 1);

    assert(std::remove(bc1Path.c_str()) == 0);
    assert(std::remove(bc2Path.c_str()) == 0);
    assert(std::remove(ambiguousBc1Path.c_str()) == 0);
    assert(::rmdir(root.c_str()) == 0);
}

} // namespace

int main()
{
    std::uint32_t packedUmi = 0;
    assert(spatial_r1_decoder::encodeRawUmi9("ACGTTACGT", 9, packedUmi));
    assert(spatial_r1_decoder::decodeRawUmi9(packedUmi) == "ACGTTACGT");
    assert(!spatial_r1_decoder::encodeRawUmi9("ACGNTACGT", 9, packedUmi));
    assert(!spatial_r1_decoder::encodeRawUmi9("ACGX", 4, packedUmi));
    testDecoderAmbiguousBases();

    std::string error;
    std::uint8_t products = 0;
    assert(spatial_gex::parseProducts("all", products, error));
    assert(products == spatial_gex::ProductAll);
    assert(spatial_gex::parseProducts("strict,soft_expected,hard,gated_hard",
                                      products, error));
    assert(products == spatial_gex::ProductAll);
    assert(spatial_gex::parseProducts("hard", products, error));
    assert(products == spatial_gex::ProductHard);
    assert(!spatial_gex::parseProducts("strict,strict", products, error));
    assert(!spatial_gex::parseProducts("exact", products, error));

    const spatial_gex::PipelineConfig defaultPipelineConfig;
    assert(defaultPipelineConfig.products == spatial_gex::ProductHard);

    std::uint8_t scales = 0;
    assert(spatial_gex::parseScales("2,8,16", scales, error));
    assert(scales == spatial_gex::ScaleAll);
    assert(spatial_gex::parseScales("16,2", scales, error));
    assert(scales == (spatial_gex::Scale16um | spatial_gex::Scale2um));
    assert(!spatial_gex::parseScales("4", scales, error));

    spatial_gex::OverflowPolicy overflow = spatial_gex::OverflowPolicy::Spill;
    assert(spatial_gex::parseOverflowPolicy("Fail", overflow, error));
    assert(overflow == spatial_gex::OverflowPolicy::Fail);
    assert(spatial_gex::parseOverflowPolicy("spill", overflow, error));
    assert(overflow == spatial_gex::OverflowPolicy::Spill);
    assert(!spatial_gex::parseOverflowPolicy("sidecar", overflow, error));

    spatial_gex::Capacity small;
    small.reads = 100000;
    small.candidates = 111744;
    small.threads = 16;
    spatial_gex::MemoryModel smallModel;
    assert(spatial_gex::estimateMemory(small, smallModel, error));
    assert(smallModel.accumulationFixedBytes > 0);
    assert(smallModel.accumulationBytes > smallModel.accumulationFixedBytes);
    assert(smallModel.downstreamSpoolBytes > smallModel.accumulationFixedBytes);
    assert(smallModel.downstreamSpoolDiskBytes
           == small.reads * 256ULL + small.candidates * 64ULL
               + 1ULL * 1024ULL * 1024ULL * 1024ULL);
    assert(smallModel.peakBytes >= smallModel.accumulationBytes);
    assert(smallModel.peakBytes >= smallModel.cliqueBytes);

    spatial_gex::Capacity large = small;
    large.reads = 474131092;
    large.candidates = 529580381;
    large.threads = 64;
    spatial_gex::MemoryModel largeModel;
    assert(spatial_gex::estimateMemory(large, largeModel, error));
    assert(largeModel.peakBytes > smallModel.peakBytes);
    assert(largeModel.peakBytes == largeModel.correctionBytes);
    assert(largeModel.downstreamSpoolDiskBytes > smallModel.downstreamSpoolDiskBytes);

    std::uint64_t budget = 0;
    assert(spatial_gex::memoryFits(smallModel, smallModel.peakBytes * 2, 0.75,
                                   budget, error));
    assert(!spatial_gex::memoryFits(smallModel, smallModel.peakBytes, 0.50,
                                    budget, error));
    assert(!spatial_gex::memoryFits(smallModel, smallModel.peakBytes, 0.0,
                                    budget, error));

    error.clear();
    assert(spatial_gex::spillBudgetFits(
        smallModel, smallModel.downstreamSpoolBytes, error));
    error.clear();
    assert(!spatial_gex::spillBudgetFits(
        smallModel, smallModel.downstreamSpoolBytes - 1, error));
    assert(error.find("downstream spool estimate") != std::string::npos);
    error.clear();
    assert(!spatial_gex::spillBudgetFits(
        smallModel, smallModel.accumulationFixedBytes, error));
    assert(error.find("downstream spool estimate") != std::string::npos
           || error.find("fixed accumulation state") != std::string::npos);

    // The 100K-sized in-memory estimate and the full-slide in-memory estimate
    // differ substantially, while Spill retains the same checked downstream
    // workspace. This is the CI guard for the small-host full-slide use case.
    assert(largeModel.downstreamSpoolBytes == smallModel.downstreamSpoolBytes
           + (large.threads - small.threads) * 8ULL * 1024ULL * 1024ULL
           + (large.threads - small.threads) * 6518ULL * sizeof(std::uint64_t));
    const std::uint64_t observedHostBudget = UINT64_C(89941436006);
    error.clear();
    assert(largeModel.peakBytes > observedHostBudget);
    assert(spatial_gex::spillBudgetFits(largeModel, observedHostBudget, error));

    spatial_gex::Capacity overflowCapacity;
    overflowCapacity.reads = std::numeric_limits<std::uint64_t>::max();
    overflowCapacity.candidates = 1;
    overflowCapacity.threads = 1;
    assert(!spatial_gex::estimateMemory(overflowCapacity, largeModel, error));
    return 0;
}

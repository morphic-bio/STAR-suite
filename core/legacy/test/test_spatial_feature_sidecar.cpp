#include "SpatialFeatureSidecar.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <thread>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

namespace sfs = spatial_feature_sidecar;

namespace {

std::vector<unsigned char> readFile(const std::string &path)
{
    std::ifstream input(path.c_str(), std::ios::binary);
    assert(input.good());
    return std::vector<unsigned char>((std::istreambuf_iterator<char>(input)),
                                      std::istreambuf_iterator<char>());
}

void removePrefix(const std::string &prefix)
{
    for (const char *suffix : {".bin", ".bin.tmp", ".features.tsv",
                               ".read_name_digests.tsv", ".read_names.tmp",
                               ".summary.json"}) {
        std::remove((prefix + suffix).c_str());
    }
}

sfs::WriterConfig config(const std::string &prefix)
{
    sfs::WriterConfig value;
    value.prefix = prefix;
    value.starSuiteVersion = "test-version";
    value.sourceRevision = "test-revision";
    value.featureType = "GeneFull";
    value.policy = "GeneFull;Unique;1MM_CR;MultiGeneUMI_CR";
    value.inputManifest = "lane0\tr1.fastq.gz\nlane0\tr2.fastq.gz\n";
    value.strand = 0;
    value.crMultimapRescue = true;
    value.crIntronicFallback = true;
    return value;
}

std::vector<sfs::Record> records()
{
    std::vector<sfs::Record> result(8);
    result[0].geneIndex = 0;
    result[0].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kUniqueGene;
    result[1].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kNoGene;
    result[2].statusFlags = sfs::kRecordPresent | sfs::kMapped
        | sfs::kMultiGeneRejected;
    result[3].statusFlags = sfs::kRecordPresent | sfs::kUnmappedOrFiltered
        | sfs::kNoGene;
    result[4].geneIndex = 1;
    result[4].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kUniqueGene
        | sfs::kSameGeneGenomicMultimapper;
    result[5].geneIndex = 1;
    result[5].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kUniqueGene
        | sfs::kCrExonicRescue | (3u << sfs::kOverlapShift);
    result[6].geneIndex = 0;
    result[6].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kUniqueGene
        | sfs::kCrIntronicFallback | (2u << sfs::kOverlapShift);
    result[7].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kNoGene
        | sfs::kAlignmentEvidenceRejected | sfs::kBestScoreNaDecoy;
    return result;
}

void writeComplete(const std::string &prefix, bool threaded)
{
    sfs::Writer writer;
    std::string error;
    assert(writer.open(config(prefix), {"ENSG0", "ENSG1"}, {"zero", "one"}, error));
    const std::vector<sfs::Record> values = records();
    auto writeOne = [&](std::size_t ordinal) {
        std::string localError;
        const std::uint32_t lane = ordinal < 4 ? 0 : 1;
        assert(writer.write(ordinal, lane, "read-" + std::to_string(ordinal),
                            values[ordinal], localError));
    };
    if (threaded) {
        std::vector<std::thread> threads;
        for (std::size_t ordinal = 0; ordinal < values.size(); ++ordinal) {
            threads.push_back(std::thread(writeOne, values.size() - 1 - ordinal));
        }
        for (std::thread &thread : threads) thread.join();
    } else {
        for (std::size_t ordinal = 0; ordinal < values.size(); ++ordinal) writeOne(ordinal);
    }
    assert(writer.finalize(values.size(), error));
}

} // namespace

int main()
{
    const sfs::Record geneZero = records()[0];
    const std::vector<unsigned char> encoded = sfs::encodeRecord(geneZero);
    assert(encoded.size() == sfs::kRecordBytes);
    assert(encoded[0] == 0 && encoded[1] == 0 && encoded[2] == 0 && encoded[3] == 0);
    sfs::Record decoded;
    std::string error;
    assert(sfs::decodeRecord(encoded.data(), encoded.size(), decoded, error));
    assert(decoded.geneIndex == 0);
    assert(decoded.statusFlags == geneZero.statusFlags);
    assert(sfs::validateRecord(decoded, 2, error));
    sfs::Record missing;
    missing.geneIndex = 0;
    assert(!sfs::validateRecord(missing, 2, error));
    sfs::Record invalidRejectedUnique = geneZero;
    invalidRejectedUnique.statusFlags |= sfs::kAlignmentEvidenceRejected;
    assert(!sfs::validateRecord(invalidRejectedUnique, 2, error));
    sfs::Record invalidOrphanReason = records()[1];
    invalidOrphanReason.statusFlags |= sfs::kBestScoreNaDecoy;
    assert(!sfs::validateRecord(invalidOrphanReason, 2, error));

    char directoryTemplate[] = "/tmp/star-spatial-sidecar-test-XXXXXX";
    const char *directory = ::mkdtemp(directoryTemplate);
    assert(directory != nullptr);
    const std::string serial = std::string(directory) + "/serial";
    const std::string threaded = std::string(directory) + "/threaded";
    writeComplete(serial, false);
    writeComplete(threaded, true);
    assert(readFile(serial + ".bin") == readFile(threaded + ".bin"));
    assert(readFile(serial + ".features.tsv") == readFile(threaded + ".features.tsv"));
    assert(readFile(serial + ".read_name_digests.tsv")
           == readFile(threaded + ".read_name_digests.tsv"));

    sfs::Reader reader;
    assert(reader.open(serial + ".bin", error));
    assert(reader.header().complete);
    assert(reader.header().schemaVersion == sfs::kSchemaVersion);
    assert(reader.header().headerBytes == sfs::kHeaderBytes);
    assert(reader.header().recordBytes == sfs::kRecordBytes);
    assert(reader.header().totalReads == records().size());
    assert(reader.validateAll(error));
    for (std::size_t ordinal = 0; ordinal < records().size(); ++ordinal) {
        sfs::Record actual;
        assert(reader.read(ordinal, actual, error));
        assert(actual.geneIndex == records()[ordinal].geneIndex);
        assert(actual.statusFlags == records()[ordinal].statusFlags);
    }

    const std::string incomplete = std::string(directory) + "/incomplete";
    {
        sfs::Writer writer;
        assert(writer.open(config(incomplete), {"ENSG0"}, {"zero"}, error));
        sfs::Reader incompleteReader;
        assert(!incompleteReader.open(incomplete + ".bin.tmp", error));
    }

    const std::string missingPrefix = std::string(directory) + "/missing";
    {
        sfs::Writer writer;
        assert(writer.open(config(missingPrefix), {"ENSG0"}, {"zero"}, error));
        sfs::Record value = records()[0];
        assert(writer.write(1, 0, "read-1", value, error));
        assert(!writer.finalize(2, error));
    }

    const std::string duplicate = std::string(directory) + "/duplicate";
    {
        sfs::Writer writer;
        assert(writer.open(config(duplicate), {"ENSG0"}, {"zero"}, error));
        sfs::Record value = records()[0];
        assert(writer.write(0, 0, "read-0", value, error));
        assert(writer.write(0, 0, "read-0", value, error));
        assert(!writer.finalize(1, error));
    }

    const std::string truncated = std::string(directory) + "/truncated";
    writeComplete(truncated, false);
    struct stat info;
    assert(::stat((truncated + ".bin").c_str(), &info) == 0);
    assert(::truncate((truncated + ".bin").c_str(), info.st_size - 1) == 0);
    sfs::Reader truncatedReader;
    assert(!truncatedReader.open(truncated + ".bin", error));

    removePrefix(serial);
    removePrefix(threaded);
    removePrefix(incomplete);
    removePrefix(missingPrefix);
    removePrefix(duplicate);
    removePrefix(truncated);
    assert(::rmdir(directory) == 0);
    std::cout << "SpatialFeatureSidecar tests passed\n";
    return 0;
}

#include "SpatialFeatureSidecar.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

namespace sfs = spatial_feature_sidecar;

namespace {

void writeText(const std::string &path, const std::string &value)
{
    std::ofstream output(path.c_str());
    assert(output.good());
    output << value;
    output.close();
    assert(output.good());
}

std::string readText(const std::string &path)
{
    std::ifstream input(path.c_str());
    assert(input.good());
    return std::string((std::istreambuf_iterator<char>(input)),
                       std::istreambuf_iterator<char>());
}

void writeFastq(const std::string &path, const std::vector<std::string> &names,
                const std::vector<std::string> &sequences)
{
    assert(names.size() == sequences.size());
    gzFile output = gzopen(path.c_str(), "wb");
    assert(output != nullptr);
    for (std::size_t index = 0; index < names.size(); ++index) {
        const std::string record = "@" + names[index] + " 1:N:0:1\n" + sequences[index]
            + "\n+\n" + std::string(sequences[index].size(), 'I') + "\n";
        assert(gzwrite(output, record.data(), record.size()) == static_cast<int>(record.size()));
    }
    assert(gzclose(output) == Z_OK);
}

int run(const std::string &binary, const std::string &directory,
        const std::string &candidate, const std::string &output,
        const std::string &summary, bool swapLanes)
{
    const std::string lane0 = directory + "/lane0";
    const std::string lane1 = directory + "/lane1";
    std::string command = binary
        + " --sidecar " + directory + "/sidecar.bin"
        + " --features " + directory + "/sidecar.features.tsv"
        + " --read-name-digests " + directory + "/sidecar.read_name_digests.tsv"
        + " --candidates " + candidate
        + " --h0-prior " + directory + "/h0.tsv"
        + " --barcode-contract " + directory + "/barcode_coords.tsv";
    if (swapLanes) {
        command += " --r1-fastq " + lane1 + "_R1.fastq.gz --r1-fastq " + lane0 + "_R1.fastq.gz"
            " --r2-fastq " + lane1 + "_R2.fastq.gz --r2-fastq " + lane0 + "_R2.fastq.gz";
    } else {
        command += " --r1-fastq " + lane0 + "_R1.fastq.gz --r1-fastq " + lane1 + "_R1.fastq.gz"
            " --r2-fastq " + lane0 + "_R2.fastq.gz --r2-fastq " + lane1 + "_R2.fastq.gz";
    }
    command += " --r1-length 12 --r2-length 8 --expected-reads 4 --output " + output
        + " --summary " + summary + " >/dev/null 2>&1";
    const int status = std::system(command.c_str());
    return status == -1 ? -1 : WEXITSTATUS(status);
}

int runFused(const std::string &binary, const std::string &directory,
             const std::string &candidate, const std::string &decoded,
             const std::string &output, const std::string &summary)
{
    const std::string command = binary
        + " --sidecar " + directory + "/sidecar.bin"
        + " --features " + directory + "/sidecar.features.tsv"
        + " --read-name-digests " + directory + "/sidecar.read_name_digests.tsv"
        + " --candidates " + candidate
        + " --h0-prior " + directory + "/h0.tsv"
        + " --barcode-contract " + directory + "/barcode_coords.tsv"
        + " --decode-reads " + decoded
        + " --expected-reads 4 --output " + output
        + " --summary " + summary + " >/dev/null 2>&1";
    const int status = std::system(command.c_str());
    return status == -1 ? -1 : WEXITSTATUS(status);
}

} // namespace

int main(int argc, char **argv)
{
    assert(argc == 2);
    char directoryTemplate[] = "/tmp/star-spatial-sidecar-join-test-XXXXXX";
    const char *created = ::mkdtemp(directoryTemplate);
    assert(created != nullptr);
    const std::string directory = created;

    writeFastq(directory + "/lane0_R1.fastq.gz", {"read0", "read1"},
               {"AAAAAAAAACCC", "CCCCCCCCCAAA"});
    writeFastq(directory + "/lane0_R2.fastq.gz", {"read0", "read1"},
               {"ACGTACGT", "TGCATGCA"});
    writeFastq(directory + "/lane1_R1.fastq.gz", {"read2", "read3"},
               {"GGGGGGGGGTTT", "TTTTTTTTTGGG"});
    writeFastq(directory + "/lane1_R2.fastq.gz", {"read2", "read3"},
               {"AAAACCCC", "GGGGTTTT"});

    sfs::WriterConfig config;
    config.prefix = directory + "/sidecar";
    config.starSuiteVersion = "test";
    config.sourceRevision = "test";
    config.featureType = "GeneFull";
    config.policy = "test";
    config.inputManifest = "test";
    config.strand = 0;
    config.crMultimapRescue = true;
    config.crIntronicFallback = true;
    sfs::Writer writer;
    std::string error;
    assert(writer.open(config, {"GENE0", "GENE1"}, {"zero", "one"}, error));
    std::vector<sfs::Record> records(4);
    records[0].geneIndex = 0;
    records[0].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kUniqueGene;
    records[1].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kNoGene;
    records[2].geneIndex = 1;
    records[2].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kUniqueGene;
    records[3].statusFlags = sfs::kRecordPresent | sfs::kMapped | sfs::kMultiGeneRejected;
    for (std::size_t ordinal = 0; ordinal < records.size(); ++ordinal) {
        assert(writer.write(ordinal, ordinal < 2 ? 0 : 1,
                            "read" + std::to_string(ordinal), records[ordinal], error));
    }
    assert(writer.finalize(records.size(), error));

    writeText(directory + "/h0.tsv",
        "barcode_half\toligo_index\toligo_sequence\toligo_length\texact_h0_read_count\n"
        "BC1\t0\tAA\t2\t2\nBC1\t1\tCC\t2\t3\n"
        "BC2\t0\tGG\t2\t2\nBC2\t1\tTT\t2\t3\n");
    writeText(directory + "/barcode_coords.tsv",
        "barcode\trow2\tcol2\nAA_GG\t0\t0\nCC_GG\t0\t1\nAA_TT\t1\t0\nCC_TT\t1\t1\n");
    const std::string header =
        "global_ordinal\tread_id\traw_umi\tcandidate_count\trow2\tcol2\tmin_tier\t"
        "bc1_edit\tbc2_edit\tbc1_obs_len\tbc2_obs_len\tfull_start\t"
        "log_sequence_likelihood\n";
    const std::string candidates = directory + "/candidates.tsv";
    writeText(candidates, header
        + "0\tread0\tAAAAAAAAA\t2\t0\t0\t0\t0\t0\t2\t2\t9\t-0.1\n"
          "0\tread0\tAAAAAAAAA\t2\t1\t1\t1\t1\t0\t2\t2\t9\t-1.2\n"
          "2\tread2\tGGGGGGGGG\t1\t1\t0\t0\t0\t0\t2\t2\t9\t-0.2\n");
    const std::string output = directory + "/normalized.tsv";
    const std::string summary = directory + "/join.json";
    assert(run(argv[1], directory, candidates, output, summary, false) == 0);
    const std::string normalized = readText(output);
    assert(normalized.find("read0\tGENE0\tAAAAAAAAA\ts_002um_0_0\t")
           != std::string::npos);
    assert(normalized.find("\t8\nread0\tGENE0\tAAAAAAAAA\ts_002um_1_1\t")
           != std::string::npos);
    assert(normalized.find("read0\tGENE0\tAAAAAAAAA\ts_002um_1_1\t-1.2\t15\n")
           != std::string::npos);
    assert(normalized.find("read2\tGENE1\tGGGGGGGGG\ts_002um_1_0\t")
           != std::string::npos);
    assert(normalized.find("\t11\n") != std::string::npos);
    assert(readText(summary).find("\"candidate_less_reads\": 2") != std::string::npos);

    const std::string decodeHeader =
        "read_id\tuniverse_unique_assigned\trow2\tcol2\tfull_start\n";
    const std::string decoded = directory + "/decode_reads.tsv";
    writeText(decoded, decodeHeader
        + "read0\tTrue\t0\t0\t9\n"
          "read1\tFalse\t\t\t\n"
          "read2\tTrue\t1\t0\t9\n"
          "read3\tFalse\t\t\t\n");
    const std::string fusedOutput = directory + "/fused_normalized.tsv";
    const std::string fusedSummary = directory + "/fused_join.json";
    assert(runFused(argv[1], directory, candidates, decoded,
                    fusedOutput, fusedSummary) == 0);
    assert(readText(fusedOutput) == normalized);
    assert(readText(fusedSummary).find("\"input_mode\": \"fused_r1_tap\"")
           != std::string::npos);
    assert(readText(fusedSummary).find("\"single_star_input_stream\": true")
           != std::string::npos);

    const std::string reordered = directory + "/decode_reordered.tsv";
    writeText(reordered, decodeHeader
        + "read1\tFalse\t\t\t\nread0\tTrue\t0\t0\t9\n"
          "read2\tTrue\t1\t0\t9\nread3\tFalse\t\t\t\n");
    assert(runFused(argv[1], directory, candidates, reordered,
                    directory + "/reordered.tsv", directory + "/reordered.json") != 0);
    const std::string shortDecode = directory + "/decode_short.tsv";
    writeText(shortDecode, decodeHeader
        + "read0\tTrue\t0\t0\t9\nread1\tFalse\t\t\t\n"
          "read2\tTrue\t1\t0\t9\n");
    assert(runFused(argv[1], directory, candidates, shortDecode,
                    directory + "/short.tsv", directory + "/short.json") != 0);

    assert(run(argv[1], directory, candidates, directory + "/permuted.tsv",
               directory + "/permuted.json", true) != 0);
    const std::string badUmi = directory + "/bad_umi.tsv";
    writeText(badUmi, header
        + "0\tread0\tTTTTTTTTT\t1\t0\t0\t0\t0\t0\t2\t2\t9\t-0.1\n");
    assert(run(argv[1], directory, badUmi, directory + "/bad_umi_out.tsv",
               directory + "/bad_umi.json", false) != 0);
    const std::string forbidden = directory + "/forbidden.tsv";
    writeText(forbidden, header.substr(0, header.size() - 1) + "\tGX\n");
    assert(run(argv[1], directory, forbidden, directory + "/forbidden_out.tsv",
               directory + "/forbidden.json", false) != 0);

    const std::string cleanup = "rm -rf -- " + directory;
    assert(std::system(cleanup.c_str()) == 0);
    std::cout << "Spatial feature sidecar join tests passed\n";
    return 0;
}

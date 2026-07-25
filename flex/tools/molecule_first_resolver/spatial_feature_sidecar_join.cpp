#include "SpatialFeatureSidecar.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <openssl/sha.h>
#include <sys/stat.h>
#include <zlib.h>

namespace sfs = spatial_feature_sidecar;

namespace {

struct Arguments {
    std::string sidecar;
    std::string features;
    std::string nameDigests;
    std::string candidates;
    std::string h0Prior;
    std::string barcodeContract;
    std::string output;
    std::string summary;
    std::string decodeReads;
    std::vector<std::string> r1;
    std::vector<std::string> r2;
    std::uint64_t expectedReads = 0;
    std::size_t r1Length = 43;
    std::size_t r2Length = 75;
};

void usage(std::ostream &output)
{
    output
        << "Usage: spatial_feature_sidecar_join [options]\n"
        << "  --sidecar FILE.bin --features FILE.features.tsv\n"
        << "  --read-name-digests FILE.tsv --candidates FILE.tsv\n"
        << "  --h0-prior FILE.tsv --barcode-contract barcode_coords.tsv\n"
        << "  --r1-fastq L001_R1.fastq.gz [--r1-fastq ...]\n"
        << "  --r2-fastq L001_R2.fastq.gz [--r2-fastq ...]\n"
        << "  OR --decode-reads decode_reads.tsv (fused STAR R1-tap mode; no FASTQ rescan)\n"
        << "  --output normalized.tsv --summary summary.json\n"
        << "  [--expected-reads N] [--r1-length 43] [--r2-length 75]\n";
}

std::uint64_t unsignedInteger(const std::string &value, const std::string &label)
{
    if (value.empty() || value[0] == '-') throw std::runtime_error("invalid " + label);
    std::size_t used = 0;
    const unsigned long long parsed = std::stoull(value, &used);
    if (used != value.size()) throw std::runtime_error("invalid " + label + ": " + value);
    return static_cast<std::uint64_t>(parsed);
}

Arguments parseArguments(int argc, char **argv)
{
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--help" || option == "-h") {
            usage(std::cout);
            std::exit(0);
        }
        if (index + 1 >= argc) throw std::invalid_argument("missing value for " + option);
        const std::string value = argv[++index];
        if (option == "--sidecar") result.sidecar = value;
        else if (option == "--features") result.features = value;
        else if (option == "--read-name-digests") result.nameDigests = value;
        else if (option == "--candidates") result.candidates = value;
        else if (option == "--h0-prior") result.h0Prior = value;
        else if (option == "--barcode-contract") result.barcodeContract = value;
        else if (option == "--r1-fastq") result.r1.push_back(value);
        else if (option == "--r2-fastq") result.r2.push_back(value);
        else if (option == "--decode-reads") result.decodeReads = value;
        else if (option == "--output") result.output = value;
        else if (option == "--summary") result.summary = value;
        else if (option == "--expected-reads") result.expectedReads = unsignedInteger(value, option);
        else if (option == "--r1-length") result.r1Length = unsignedInteger(value, option);
        else if (option == "--r2-length") result.r2Length = unsignedInteger(value, option);
        else throw std::invalid_argument("unknown option: " + option);
    }
    if (result.sidecar.empty() || result.features.empty() || result.nameDigests.empty()
        || result.candidates.empty() || result.h0Prior.empty()
        || result.barcodeContract.empty() || result.output.empty() || result.summary.empty()) {
        throw std::invalid_argument("all sidecar, candidate, prior, contract, and output paths are required");
    }
    const bool fastqMode = !result.r1.empty() || !result.r2.empty();
    const bool fusedMode = !result.decodeReads.empty();
    if (fastqMode == fusedMode) {
        throw std::invalid_argument("choose exactly one input validation mode: paired FASTQs or --decode-reads");
    }
    if (fastqMode && (result.r1.empty() || result.r1.size() != result.r2.size())) {
        throw std::invalid_argument("paired R1/R2 lane lists are required in FASTQ mode");
    }
    if (result.r1Length == 0 || result.r2Length == 0) {
        throw std::invalid_argument("FASTQ lengths must be positive");
    }
    return result;
}

std::vector<std::string> tabs(const std::string &line)
{
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos) return fields;
        begin = end + 1;
    }
}

void refuseExisting(const std::string &path)
{
    if (::access(path.c_str(), F_OK) == 0) {
        throw std::runtime_error("refusing to overwrite output: " + path);
    }
}

class AtomicOutput {
  public:
    explicit AtomicOutput(const std::string &path)
        : final_(path), temporary_(path + ".tmp")
    {
        refuseExisting(final_);
        refuseExisting(temporary_);
        stream_.open(temporary_.c_str());
        if (!stream_) throw std::runtime_error("cannot create output: " + temporary_);
    }
    std::ostream &stream() { return stream_; }
    void commit()
    {
        stream_.close();
        if (!stream_ || std::rename(temporary_.c_str(), final_.c_str()) != 0) {
            throw std::runtime_error("cannot commit output: " + final_);
        }
        committed_ = true;
    }
    ~AtomicOutput() { if (!committed_) std::remove(temporary_.c_str()); }
  private:
    std::string final_, temporary_;
    std::ofstream stream_;
    bool committed_ = false;
};

std::string stripLine(std::string value)
{
    while (!value.empty() && (value.back() == '\n' || value.back() == '\r')) value.pop_back();
    return value;
}

class FastqReader {
  public:
    explicit FastqReader(const std::string &path) : path_(path), file_(gzopen(path.c_str(), "rb"))
    {
        if (file_ == nullptr) throw std::runtime_error("cannot open FASTQ: " + path);
    }
    ~FastqReader() { if (file_ != nullptr) gzclose(file_); }
    bool next(std::string &name, std::string &sequence, std::string &quality)
    {
        std::string header, plus;
        if (!line(header)) return false;
        if (!line(sequence) || !line(plus) || !line(quality)) {
            throw std::runtime_error("truncated FASTQ record: " + path_);
        }
        header = stripLine(header);
        sequence = stripLine(sequence);
        plus = stripLine(plus);
        quality = stripLine(quality);
        if (header.empty() || header[0] != '@' || plus.empty() || plus[0] != '+'
            || sequence.size() != quality.size()) {
            throw std::runtime_error("invalid FASTQ record: " + path_);
        }
        name = sfs::normalizeReadName(header);
        return true;
    }
  private:
    bool line(std::string &value)
    {
        value.clear();
        char buffer[1 << 16];
        while (gzgets(file_, buffer, sizeof(buffer)) != nullptr) {
            value += buffer;
            if (!value.empty() && value.back() == '\n') return true;
        }
        if (gzeof(file_)) return !value.empty();
        int errorNumber = 0;
        const char *message = gzerror(file_, &errorNumber);
        throw std::runtime_error("FASTQ read failed: " + path_ + ": "
                                 + (message == nullptr ? "unknown" : message));
    }
    std::string path_;
    gzFile file_;
};

std::uint64_t fnv1a64(const std::string &value)
{
    std::uint64_t hash = 1469598103934665603ULL;
    for (unsigned char byte : value) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::string digestHex(const unsigned char *digest)
{
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (std::size_t index = 0; index < SHA256_DIGEST_LENGTH; ++index) {
        output << std::setw(2) << static_cast<unsigned int>(digest[index]);
    }
    return output.str();
}

struct DigestRow {
    std::uint32_t lane = 0;
    std::uint64_t first = 0;
    std::uint64_t count = 0;
    std::string sha;
    bool operator==(const DigestRow &other) const
    {
        return lane == other.lane && first == other.first && count == other.count
            && sha == other.sha;
    }
};

std::vector<DigestRow> loadDigests(const std::string &path)
{
    std::ifstream input(path.c_str());
    std::string line;
    if (!input || !std::getline(input, line)
        || line != "lane_index\tfirst_ordinal\tread_count\tsha256") {
        throw std::runtime_error("invalid STAR read-name digest table: " + path);
    }
    std::vector<DigestRow> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = tabs(line);
        if (fields.size() != 4 || fields[3].size() != 64) {
            throw std::runtime_error("invalid STAR read-name digest row");
        }
        DigestRow row;
        row.lane = static_cast<std::uint32_t>(unsignedInteger(fields[0], "lane index"));
        row.first = unsignedInteger(fields[1], "first ordinal");
        row.count = unsignedInteger(fields[2], "digest read count");
        row.sha = fields[3];
        result.push_back(row);
    }
    if (result.empty()) throw std::runtime_error("empty STAR read-name digest table");
    return result;
}

std::vector<std::string> loadFeatures(const std::string &path, const sfs::Header &header)
{
    std::string hashError;
    const std::string digest = sfs::sha256File(path, hashError);
    if (digest.empty() || digest != header.geneDictionarySha256) {
        throw std::runtime_error(hashError.empty() ? "feature dictionary digest mismatch" : hashError);
    }
    std::ifstream input(path.c_str());
    std::string line;
    if (!input || !std::getline(input, line) || line != "gene_index\tgene_id\tgene_name") {
        throw std::runtime_error("invalid feature dictionary header");
    }
    std::vector<std::string> genes;
    while (std::getline(input, line)) {
        const std::vector<std::string> fields = tabs(line);
        if (fields.size() != 3 || fields[1].empty()
            || unsignedInteger(fields[0], "gene index") != genes.size()) {
            throw std::runtime_error("invalid feature dictionary row");
        }
        genes.push_back(fields[1]);
    }
    if (genes.size() != header.featureCount) {
        throw std::runtime_error("feature dictionary count mismatch");
    }
    return genes;
}

struct H0Prior {
    std::vector<std::uint64_t> bc1;
    std::vector<std::uint64_t> bc2;
};

H0Prior loadH0(const std::string &path)
{
    std::ifstream input(path.c_str());
    std::string line;
    if (!input || !std::getline(input, line)
        || line != "barcode_half\toligo_index\toligo_sequence\toligo_length\texact_h0_read_count") {
        throw std::runtime_error("invalid H0 prior header");
    }
    H0Prior result;
    std::map<std::string, std::vector<std::uint64_t> *> axes = {
        {"BC1", &result.bc1}, {"BC2", &result.bc2}
    };
    while (std::getline(input, line)) {
        const std::vector<std::string> fields = tabs(line);
        if (fields.size() != 5 || !axes.count(fields[0]) || fields[2].empty()
            || unsignedInteger(fields[3], "oligo length") != fields[2].size()) {
            throw std::runtime_error("invalid H0 prior row");
        }
        std::vector<std::uint64_t> &axis = *axes.at(fields[0]);
        if (unsignedInteger(fields[1], "oligo index") != axis.size()) {
            throw std::runtime_error("H0 prior axis is incomplete or out of order");
        }
        axis.push_back(unsignedInteger(fields[4], "exact H0 read count"));
    }
    if (result.bc1.empty() || result.bc2.empty()) throw std::runtime_error("empty H0 prior axis");
    std::uint64_t bc1Mass = 0, bc2Mass = 0;
    for (std::uint64_t value : result.bc1) bc1Mass += value;
    for (std::uint64_t value : result.bc2) bc2Mass += value;
    if (bc1Mass != bc2Mass) throw std::runtime_error("BC1/BC2 H0 masses do not reconcile");
    return result;
}

void validateContract(const std::string &path, std::size_t rows, std::size_t columns)
{
    std::ifstream input(path.c_str());
    std::string line;
    if (!input || !std::getline(input, line) || line != "barcode\trow2\tcol2") {
        throw std::runtime_error("invalid barcode coordinate contract header");
    }
    std::uint64_t ordinal = 0;
    while (std::getline(input, line)) {
        const std::vector<std::string> fields = tabs(line);
        if (fields.size() != 3 || fields[0].empty()) {
            throw std::runtime_error("invalid barcode coordinate contract row");
        }
        const std::uint64_t row = unsignedInteger(fields[1], "contract row2");
        const std::uint64_t column = unsignedInteger(fields[2], "contract col2");
        if (row >= rows || column >= columns || ordinal != row * columns + column) {
            throw std::runtime_error("barcode coordinate contract is not the complete row-major universe");
        }
        ++ordinal;
    }
    if (ordinal != static_cast<std::uint64_t>(rows) * columns) {
        throw std::runtime_error("barcode coordinate contract universe size mismatch");
    }
}

struct Candidate {
    std::uint64_t ordinal = 0;
    std::string readId;
    std::string rawUmi;
    std::uint64_t declaredCount = 0;
    std::uint32_t row = 0;
    std::uint32_t column = 0;
    double likelihood = 0.0;
};

class CandidateReader {
  public:
    explicit CandidateReader(const std::string &path) : input_(path.c_str())
    {
        std::string header;
        const std::string expected =
            "global_ordinal\tread_id\traw_umi\tcandidate_count\trow2\tcol2\tmin_tier\t"
            "bc1_edit\tbc2_edit\tbc1_obs_len\tbc2_obs_len\tfull_start\t"
            "log_sequence_likelihood";
        if (!input_ || !std::getline(input_, header) || header != expected) {
            throw std::runtime_error("candidate input does not match the ordinal raw-R1 schema");
        }
        advance();
    }
    bool has() const { return have_; }
    const Candidate &current() const { return current_; }
    void advance()
    {
        std::string line;
        while (std::getline(input_, line) && line.empty()) {}
        if (!input_) {
            have_ = false;
            return;
        }
        const std::vector<std::string> fields = tabs(line);
        if (fields.size() != 13 || fields[1].empty() || fields[2].empty()) {
            throw std::runtime_error("invalid candidate row at line " + std::to_string(++line_));
        }
        current_.ordinal = unsignedInteger(fields[0], "candidate ordinal");
        current_.readId = sfs::normalizeReadName(fields[1]);
        current_.rawUmi = fields[2];
        current_.declaredCount = unsignedInteger(fields[3], "candidate count");
        current_.row = static_cast<std::uint32_t>(unsignedInteger(fields[4], "row2"));
        current_.column = static_cast<std::uint32_t>(unsignedInteger(fields[5], "col2"));
        std::size_t used = 0;
        current_.likelihood = std::stod(fields[12], &used);
        if (used != fields[12].size() || !std::isfinite(current_.likelihood)) {
            throw std::runtime_error("non-finite candidate likelihood");
        }
        have_ = true;
        ++line_;
    }
  private:
    std::ifstream input_;
    Candidate current_;
    bool have_ = false;
    std::uint64_t line_ = 1;
};

class DecodeReader {
  public:
    explicit DecodeReader(const std::string &path) : input_(path.c_str())
    {
        std::string header;
        const std::string prefix =
            "read_id\tuniverse_unique_assigned\trow2\tcol2\tfull_start";
        if (!input_ || !std::getline(input_, header)
            || (header != prefix && header.compare(0, prefix.size() + 1, prefix + "\t") != 0)) {
            throw std::runtime_error("decode input does not match the complete raw-R1 schema");
        }
    }

    bool next(std::string &readId)
    {
        std::string line;
        while (std::getline(input_, line) && line.empty()) {}
        if (!input_) return false;
        const std::vector<std::string> fields = tabs(line);
        if (fields.size() < 5 || fields[0].empty()) {
            throw std::runtime_error("invalid decode row at line " + std::to_string(++line_));
        }
        readId = sfs::normalizeReadName(fields[0]);
        if (readId.empty()) {
            throw std::runtime_error("empty normalized decode read name at line "
                                     + std::to_string(line_));
        }
        ++line_;
        return true;
    }

  private:
    std::ifstream input_;
    std::uint64_t line_ = 1;
};

class DigestVerifier {
  public:
    explicit DigestVerifier(const std::vector<DigestRow> &expected) : expected_(expected) {}

    void add(std::uint64_t ordinal, const std::string &readName)
    {
        if (block_ >= expected_.size()) {
            throw std::runtime_error("decode/FASTQ stream contains more name-digest blocks than STAR");
        }
        const DigestRow &row = expected_[block_];
        if (count_ == 0) {
            if (row.first != ordinal || row.count == 0) {
                throw std::runtime_error("STAR read-name digest blocks are incomplete or out of order");
            }
            SHA256_Init(&context_);
        }
        const std::uint64_t nameHash = fnv1a64(readName);
        unsigned char nameBytes[8];
        for (int byte = 0; byte < 8; ++byte) {
            nameBytes[byte] = static_cast<unsigned char>(nameHash >> (8 * byte));
        }
        SHA256_Update(&context_, nameBytes, sizeof(nameBytes));
        ++count_;
        if (count_ == row.count) {
            unsigned char digest[SHA256_DIGEST_LENGTH];
            SHA256_Final(digest, &context_);
            if (digestHex(digest) != row.sha) {
                throw std::runtime_error("input stream and STAR blockwise read-name digests differ");
            }
            count_ = 0;
            ++block_;
        }
    }

    void finish() const
    {
        if (count_ != 0 || block_ != expected_.size()) {
            throw std::runtime_error("input stream ended before STAR read-name digest blocks");
        }
    }

    std::size_t blocks() const { return block_; }

  private:
    const std::vector<DigestRow> &expected_;
    std::size_t block_ = 0;
    std::uint64_t count_ = 0;
    SHA256_CTX context_;
};

std::string jsonEscape(const std::string &value)
{
    std::string output;
    for (char character : value) {
        if (character == '"' || character == '\\') output.push_back('\\');
        output.push_back(character);
    }
    return output;
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Arguments arguments = parseArguments(argc, argv);
        sfs::Reader sidecar;
        std::string error;
        if (!sidecar.open(arguments.sidecar, error) || !sidecar.validateAll(error)) {
            throw std::runtime_error(error);
        }
        if (arguments.expectedReads != 0
            && sidecar.header().totalReads != arguments.expectedReads) {
            throw std::runtime_error("sidecar read count differs from --expected-reads");
        }
        const std::vector<std::string> genes = loadFeatures(arguments.features, sidecar.header());
        const std::vector<DigestRow> expectedDigests = loadDigests(arguments.nameDigests);
        const H0Prior prior = loadH0(arguments.h0Prior);
        validateContract(arguments.barcodeContract, prior.bc2.size(), prior.bc1.size());
        CandidateReader candidates(arguments.candidates);

        AtomicOutput normalized(arguments.output);
        normalized.stream() << "read_id\tfeature_id\traw_umi\tcandidate\t"
                               "log_sequence_likelihood\texact_read_count\n"
                            << std::setprecision(17);

        DigestVerifier digestVerifier(expectedDigests);
        std::uint64_t ordinal = 0;
        std::uint64_t candidateRows = 0, candidateReads = 0, candidateLess = 0;
        std::uint64_t eligibleReads = 0, emittedReads = 0, emittedRows = 0;
        std::uint64_t noGene = 0, multiGene = 0, unmapped = 0;
        auto processRead = [&](const std::string &readName, const std::string *rawR1) {
            digestVerifier.add(ordinal, readName);
            std::vector<Candidate> rows;
            std::set<std::pair<std::uint32_t, std::uint32_t> > coordinates;
            while (candidates.has() && candidates.current().ordinal == ordinal) {
                const Candidate row = candidates.current();
                if (row.readId != readName) {
                    throw std::runtime_error("candidate/input-stream read-name mismatch");
                }
                if (rawR1 != nullptr && row.rawUmi != rawR1->substr(0, 9)) {
                    throw std::runtime_error("candidate raw UMI does not match raw R1");
                }
                if (row.rawUmi.size() != 9) {
                    throw std::runtime_error("candidate raw UMI is not 9 nt");
                }
                if (row.row >= prior.bc2.size() || row.column >= prior.bc1.size()
                    || !coordinates.insert(std::make_pair(row.row, row.column)).second) {
                    throw std::runtime_error("candidate is outside the contract or duplicated");
                }
                rows.push_back(row);
                candidates.advance();
            }
            if (candidates.has() && candidates.current().ordinal < ordinal) {
                throw std::runtime_error("candidate ordinals are not globally ordered");
            }
            if (rows.empty()) {
                ++candidateLess;
            } else {
                ++candidateReads;
                candidateRows += rows.size();
                if (rows.front().declaredCount != rows.size()) {
                    throw std::runtime_error("candidate_count does not equal the complete candidate set");
                }
                for (const Candidate &row : rows) {
                    if (row.declaredCount != rows.size() || row.readId != rows.front().readId
                        || row.rawUmi != rows.front().rawUmi) {
                        throw std::runtime_error("inconsistent candidate group");
                    }
                }
            }

            sfs::Record feature;
            if (!sidecar.read(ordinal, feature, error)) throw std::runtime_error(error);
            if (feature.statusFlags & sfs::kUniqueGene) {
                ++eligibleReads;
                if (!rows.empty()) {
                    ++emittedReads;
                    for (const Candidate &row : rows) {
                        const std::uint64_t left = prior.bc1[row.column];
                        const std::uint64_t right = prior.bc2[row.row];
                        if (left + 1 > std::numeric_limits<std::uint64_t>::max() / (right + 1)) {
                            throw std::runtime_error("exact H0 prior multiplication overflow");
                        }
                        const std::uint64_t exact = (left + 1) * (right + 1) - 1;
                        normalized.stream() << readName << '\t' << genes.at(feature.geneIndex)
                            << '\t' << row.rawUmi << "\ts_002um_" << row.row << '_'
                            << row.column << '\t' << row.likelihood << '\t' << exact << '\n';
                        ++emittedRows;
                    }
                }
            } else if (feature.statusFlags & sfs::kMultiGeneRejected) {
                ++multiGene;
            } else if (feature.statusFlags & sfs::kUnmappedOrFiltered) {
                ++unmapped;
            } else {
                ++noGene;
            }
            ++ordinal;
        };

        if (!arguments.decodeReads.empty()) {
            DecodeReader decoded(arguments.decodeReads);
            std::string readName;
            while (decoded.next(readName)) processRead(readName, nullptr);
        } else for (std::size_t lane = 0; lane < arguments.r1.size(); ++lane) {
            FastqReader r1(arguments.r1[lane]);
            FastqReader r2(arguments.r2[lane]);

            std::string r1Name, r1Sequence, r1Quality, r2Name, r2Sequence, r2Quality;
            std::uint64_t laneReads = 0;
            for (;;) {
                const bool haveR1 = r1.next(r1Name, r1Sequence, r1Quality);
                const bool haveR2 = r2.next(r2Name, r2Sequence, r2Quality);
                if (haveR1 != haveR2) throw std::runtime_error("R1/R2 FASTQ counts differ");
                if (!haveR1) break;
                if (r1Name != r2Name) throw std::runtime_error("R1/R2 read names differ");
                if (r1Sequence.size() != arguments.r1Length
                    || r2Sequence.size() != arguments.r2Length) {
                    throw std::runtime_error("FASTQ read length differs from the frozen fixture contract");
                }
                processRead(r1Name, &r1Sequence);
                ++laneReads;
            }
            if (laneReads == 0) throw std::runtime_error("empty FASTQ lane");
        }
        if (candidates.has()) throw std::runtime_error("candidate ordinal exceeds input stream/sidecar count");
        if (ordinal != sidecar.header().totalReads) {
            throw std::runtime_error("input stream read count differs from sidecar total");
        }
        digestVerifier.finish();
        normalized.commit();

        std::set<std::uint32_t> digestLanes;
        for (const DigestRow &row : expectedDigests) digestLanes.insert(row.lane);

        AtomicOutput summary(arguments.summary);
        summary.stream()
            << "{\n"
            << "  \"schema\": \"star_suite.spatial_feature_sidecar_join.v1\",\n"
            << "  \"status\": \"complete\",\n"
            << "  \"sidecar\": \"" << jsonEscape(arguments.sidecar) << "\",\n"
            << "  \"input_mode\": \""
            << (arguments.decodeReads.empty() ? "paired_fastq_contract" : "fused_r1_tap")
            << "\",\n"
            << "  \"total_reads\": " << ordinal << ",\n"
            << "  \"lanes\": " << digestLanes.size() << ",\n"
            << "  \"candidate_reads\": " << candidateReads << ",\n"
            << "  \"candidate_less_reads\": " << candidateLess << ",\n"
            << "  \"candidate_rows\": " << candidateRows << ",\n"
            << "  \"feature_eligible_reads\": " << eligibleReads << ",\n"
            << "  \"emitted_reads\": " << emittedReads << ",\n"
            << "  \"emitted_rows\": " << emittedRows << ",\n"
            << "  \"no_gene_reads\": " << noGene << ",\n"
            << "  \"multi_gene_rejected_reads\": " << multiGene << ",\n"
            << "  \"unmapped_or_filtered_reads\": " << unmapped << ",\n"
            << "  \"name_digest_blocks\": " << digestVerifier.blocks() << ",\n"
            << "  \"contract_coordinates\": "
            << static_cast<std::uint64_t>(prior.bc1.size()) * prior.bc2.size() << ",\n"
            << "  \"invariants\": {\n"
            << "    \"raw_umi_from_r1\": true,\n"
            << "    \"ordinal_join\": true,\n"
            << "    \"name_digests_match\": true,\n"
            << "    \"single_star_input_stream\": "
            << (arguments.decodeReads.empty() ? "false" : "true") << ",\n"
            << "    \"complete_candidate_sets\": true,\n"
            << "    \"coordinate_contract_complete\": true\n"
            << "  }\n"
            << "}\n";
        summary.commit();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "spatial_feature_sidecar_join: ERROR: " << error.what() << '\n';
        return 2;
    }
}

#include "MultiGeneUmiCr.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <unistd.h>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef STAR_SUITE_VERSION
#define STAR_SUITE_VERSION "unknown"
#endif

namespace {

const char *kHeader =
    "umi_mode\tproduct\tcandidate\tcorrected_umi\tfeature_id\tcorrected_count\t"
    "original_at_corrected_count\texpected_count\toriginal_expected_count\tmolecule_id\t"
    "member_read_count\tmember_read_ids\tread_clique_ids";

struct Arguments {
    std::string input;
    std::string resolvedSource;
    std::string outDir;
    std::string tmpDir;
    std::size_t partitions = 256;
    std::size_t partitionBufferKb = 1024;
    bool inputFeatureSorted = false;
    bool keepPartitions = false;
};

void usage(std::ostream &out)
{
    out << "Usage: molecule_first_gex_reconcile_hash --input feature_sorted_support.tsv "
        << "--input-feature-sorted --resolved-source DIR --out-dir DIR --tmp-dir DIR\n"
        << "       [--hash-partitions INT] [--partition-buffer-kb INT] "
        << "[--keep-partitions]\n\n"
        << "Uses deterministic bounded hash partitioning; no global text sort is performed.\n";
}

std::size_t parsePositiveSize(const std::string &value, const std::string &name,
                              std::size_t maximum)
{
    if (value.empty() || value[0] == '-') throw std::invalid_argument("invalid " + name);
    std::size_t used = 0;
    const unsigned long long parsed = std::stoull(value, &used);
    if (used != value.size() || parsed == 0 || parsed > maximum) {
        throw std::invalid_argument("invalid " + name + ": " + value);
    }
    return static_cast<std::size_t>(parsed);
}

Arguments parseArguments(int argc, char **argv)
{
    Arguments arguments;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--help" || option == "-h") {
            usage(std::cout);
            std::exit(0);
        }
        if (option == "--version") {
            std::cout << STAR_SUITE_VERSION << '\n';
            std::exit(0);
        }
        if (option == "--input-feature-sorted") {
            arguments.inputFeatureSorted = true;
            continue;
        }
        if (option == "--keep-partitions") {
            arguments.keepPartitions = true;
            continue;
        }
        if (index + 1 >= argc) throw std::invalid_argument("missing value for " + option);
        const std::string value = argv[++index];
        if (option == "--input") arguments.input = value;
        else if (option == "--resolved-source") arguments.resolvedSource = value;
        else if (option == "--out-dir") arguments.outDir = value;
        else if (option == "--tmp-dir") arguments.tmpDir = value;
        else if (option == "--hash-partitions") {
            arguments.partitions = parsePositiveSize(value, "hash partition count", 4096);
        } else if (option == "--partition-buffer-kb") {
            arguments.partitionBufferKb = parsePositiveSize(value, "partition buffer KiB", 16384);
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }
    if (arguments.input.empty() || arguments.resolvedSource.empty()
        || arguments.outDir.empty() || arguments.tmpDir.empty()) {
        throw std::invalid_argument(
            "--input, --resolved-source, --out-dir, and --tmp-dir are required");
    }
    if (!arguments.inputFeatureSorted) {
        throw std::invalid_argument("hash reconciliation requires --input-feature-sorted");
    }
    return arguments;
}

std::vector<std::string> splitTabs(const std::string &line)
{
    std::vector<std::string> fields;
    fields.reserve(13);
    std::size_t begin = 0;
    for (;;) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos) return fields;
        begin = end + 1;
    }
}

std::uint64_t parseCount(const std::string &value, const std::string &name)
{
    if (value.empty() || value[0] == '-') {
        throw std::invalid_argument("invalid " + name + ": " + value);
    }
    std::size_t used = 0;
    const unsigned long long parsed = std::stoull(value, &used);
    if (used != value.size()) throw std::invalid_argument("invalid " + name + ": " + value);
    return static_cast<std::uint64_t>(parsed);
}

double parseDouble(const std::string &value, const std::string &name)
{
    if (value.empty()) throw std::invalid_argument("missing " + name);
    std::size_t used = 0;
    const double parsed = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(parsed) || parsed < 0.0) {
        throw std::invalid_argument("invalid " + name + ": " + value);
    }
    return parsed;
}

bool nearlyEqual(double left, double right)
{
    return std::fabs(left - right) <= 1e-12
        + 1e-10 * std::max(std::fabs(left), std::fabs(right));
}

bool directoryEmpty(const std::string &path)
{
    DIR *directory = opendir(path.c_str());
    if (directory == nullptr) throw std::runtime_error("cannot inspect directory: " + path);
    bool empty = true;
    for (dirent *entry = readdir(directory); entry != nullptr; entry = readdir(directory)) {
        if (std::strcmp(entry->d_name, ".") != 0 && std::strcmp(entry->d_name, "..") != 0) {
            empty = false;
            break;
        }
    }
    closedir(directory);
    return empty;
}

bool isDirectory(const std::string &path)
{
    struct stat info;
    return stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
}

void prepareOutputDirectory(const std::string &path)
{
    struct stat info;
    if (stat(path.c_str(), &info) == 0) {
        if (!S_ISDIR(info.st_mode) || !directoryEmpty(path)) {
            throw std::runtime_error("output directory must be empty: " + path);
        }
        return;
    }
    if (errno != ENOENT || mkdir(path.c_str(), 0755) != 0) {
        throw std::runtime_error("cannot create output directory " + path + ": "
                                 + std::strerror(errno));
    }
}

class AtomicTable {
  public:
    AtomicTable(const std::string &directory, const std::string &name,
                const std::string &header)
        : final_(directory + "/" + name), temporary_(final_ + ".tmp"), output_(temporary_)
    {
        if (!output_) throw std::runtime_error("cannot open output: " + temporary_);
        output_ << header << '\n';
    }

    ~AtomicTable()
    {
        if (!committed_) {
            output_.close();
            std::remove(temporary_.c_str());
        }
    }

    std::ostream &stream() { return output_; }

    void commit()
    {
        output_.flush();
        if (!output_) throw std::runtime_error("failed writing output: " + temporary_);
        output_.close();
        if (std::rename(temporary_.c_str(), final_.c_str()) != 0) {
            throw std::runtime_error("cannot finalize output: " + final_);
        }
        committed_ = true;
    }

  private:
    std::string final_;
    std::string temporary_;
    std::ofstream output_;
    bool committed_ = false;
};

struct Row {
    std::string umiMode;
    std::string product;
    std::string candidate;
    std::string correctedUmi;
    std::string featureId;
    std::uint64_t correctedCount = 0;
    std::uint64_t originalCount = 0;
    double expectedCount = 0.0;
    double originalExpected = 0.0;
    std::string moleculeId;
    std::uint64_t memberReadCount = 0;
    std::string memberReadIds;
    std::string readCliqueIds;
};

Row parseRow(const std::vector<std::string> &fields)
{
    if (fields.size() != 13) throw std::invalid_argument("wrong provisional field count");
    Row row;
    row.umiMode = fields[0];
    row.product = fields[1];
    row.candidate = fields[2];
    row.correctedUmi = fields[3];
    row.featureId = fields[4];
    if (row.umiMode != "1mm_cr" && row.umiMode != "exact") {
        throw std::invalid_argument("invalid UMI mode: " + row.umiMode);
    }
    if (row.product == "soft_expected") {
        row.expectedCount = parseDouble(fields[7], "expected_count");
        row.originalExpected = parseDouble(fields[8], "original_expected_count");
    } else {
        if (row.product != "strict" && row.product != "hard"
            && row.product != "gated_hard") {
            throw std::invalid_argument("invalid product: " + row.product);
        }
        row.correctedCount = parseCount(fields[5], "corrected_count");
        row.originalCount = parseCount(fields[6], "original_at_corrected_count");
        row.moleculeId = fields[9];
        row.memberReadCount = parseCount(fields[10], "member_read_count");
        row.memberReadIds = fields[11];
        if (row.correctedCount != row.memberReadCount
            || row.originalCount > row.correctedCount || row.moleculeId.empty()
            || row.memberReadIds.empty()) {
            throw std::invalid_argument("invalid provisional integer support row");
        }
    }
    row.readCliqueIds = fields[12];
    if (row.candidate.empty() || row.correctedUmi.empty() || row.featureId.empty()
        || row.readCliqueIds.empty()) {
        throw std::invalid_argument("empty provisional key/support value");
    }
    return row;
}

struct GroupRows {
    Row first;
    bool populated = false;
    std::vector<Row> additional;

    void add(Row row)
    {
        if (!populated) {
            first = std::move(row);
            populated = true;
        } else {
            additional.push_back(std::move(row));
        }
    }

    std::size_t size() const { return populated ? additional.size() + 1 : 0; }

    const Row &at(std::size_t index) const
    {
        return index == 0 ? first : additional.at(index - 1);
    }
};

struct Stats {
    std::uint64_t groups = 0;
    std::uint64_t accepted = 0;
    std::uint64_t correctedCountTies = 0;
    std::uint64_t originalDominanceRejected = 0;
    double inputExpectedMass = 0.0;
    double outputExpectedMass = 0.0;
};

std::uint64_t stableHash(const char *data, std::size_t size)
{
    std::uint64_t value = UINT64_C(1469598103934665603);
    for (std::size_t index = 0; index < size; ++index) {
        value ^= static_cast<unsigned char>(data[index]);
        value *= UINT64_C(1099511628211);
    }
    return value;
}

std::size_t keyEnd(const std::string &line)
{
    std::size_t position = 0;
    for (int field = 0; field < 4; ++field) {
        position = line.find('\t', position);
        if (position == std::string::npos) {
            throw std::invalid_argument("provisional row has fewer than five fields");
        }
        ++position;
    }
    return position - 1;
}

std::string featureField(const std::string &line, std::size_t keyEndPosition)
{
    const std::size_t begin = keyEndPosition + 1;
    const std::size_t end = line.find('\t', begin);
    if (end == std::string::npos || end == begin) {
        throw std::invalid_argument("provisional row has an empty or missing feature_id");
    }
    return line.substr(begin, end - begin);
}

class PartitionFiles {
  public:
    PartitionFiles(const Arguments &arguments)
        : keep_(arguments.keepPartitions), bufferBytes_(arguments.partitionBufferKb * 1024)
    {
        if (!isDirectory(arguments.tmpDir)) {
            throw std::runtime_error("temporary root is not a directory: " + arguments.tmpDir);
        }
        std::ostringstream directory;
        directory << arguments.tmpDir << "/.molecule_first_gex_hash." << getpid();
        directory_ = directory.str();
        if (mkdir(directory_.c_str(), 0700) != 0) {
            throw std::runtime_error("cannot create private partition directory " + directory_
                                     + ": " + std::strerror(errno));
        }
        paths_.reserve(arguments.partitions);
        streams_.reserve(arguments.partitions);
        buffers_.reserve(arguments.partitions);
        try {
            for (std::size_t index = 0; index < arguments.partitions; ++index) {
                std::ostringstream name;
                name << directory_ << "/partition." << std::setw(4) << std::setfill('0')
                     << index << ".tsv";
                paths_.push_back(name.str());
                buffers_.emplace_back(new char[bufferBytes_]);
                streams_.emplace_back(new std::ofstream());
                streams_.back()->rdbuf()->pubsetbuf(buffers_.back().get(), bufferBytes_);
                streams_.back()->open(paths_.back().c_str(), std::ios::binary);
                if (!*streams_.back()) {
                    throw std::runtime_error("cannot create hash partition: " + paths_.back());
                }
            }
        } catch (...) {
            closeNoThrow();
            cleanupNoThrow();
            throw;
        }
    }

    ~PartitionFiles()
    {
        closeNoThrow();
        if (!keep_) cleanupNoThrow();
    }

    void write(std::size_t partition, const std::string &line)
    {
        std::ofstream &output = *streams_.at(partition);
        output.write(line.data(), static_cast<std::streamsize>(line.size()));
        output.put('\n');
        if (!output) throw std::runtime_error("failed writing hash partition");
    }

    void close()
    {
        for (std::size_t index = 0; index < streams_.size(); ++index) {
            streams_[index]->flush();
            if (!*streams_[index]) {
                throw std::runtime_error("failed flushing hash partition: " + paths_[index]);
            }
            streams_[index]->close();
            if (!*streams_[index]) {
                throw std::runtime_error("failed closing hash partition: " + paths_[index]);
            }
        }
        streams_.clear();
        buffers_.clear();
        closed_ = true;
    }

    const std::vector<std::string> &paths() const { return paths_; }

    void remove(std::size_t index)
    {
        if (keep_ || paths_.at(index).empty()) return;
        if (std::remove(paths_[index].c_str()) != 0 && errno != ENOENT) {
            throw std::runtime_error("cannot remove processed hash partition: " + paths_[index]);
        }
        paths_[index].clear();
    }

    void finish()
    {
        if (!keep_) {
            cleanupFilesNoThrow();
            if (!directory_.empty() && rmdir(directory_.c_str()) != 0 && errno != ENOENT) {
                throw std::runtime_error("cannot remove private partition directory: " + directory_);
            }
            directory_.clear();
        }
    }

  private:
    void closeNoThrow()
    {
        if (closed_) return;
        for (std::unique_ptr<std::ofstream> &stream : streams_) {
            if (stream) stream->close();
        }
        streams_.clear();
        buffers_.clear();
        closed_ = true;
    }

    void cleanupNoThrow()
    {
        cleanupFilesNoThrow();
        if (!directory_.empty()) rmdir(directory_.c_str());
        directory_.clear();
    }

    void cleanupFilesNoThrow()
    {
        for (std::string &path : paths_) {
            if (!path.empty()) std::remove(path.c_str());
            path.clear();
        }
    }

    bool keep_ = false;
    bool closed_ = false;
    std::size_t bufferBytes_ = 0;
    std::string directory_;
    std::vector<std::string> paths_;
    std::vector<std::unique_ptr<std::ofstream> > streams_;
    std::vector<std::unique_ptr<char[]> > buffers_;
};

void linkRequired(const std::string &source, const std::string &destination)
{
    if (::link(source.c_str(), destination.c_str()) != 0) {
        throw std::runtime_error("cannot hard-link " + source + " to " + destination
                                 + ": " + std::strerror(errno));
    }
}

class Reconciler {
  public:
    explicit Reconciler(const std::string &outDir)
        : strict_(outDir, "strict_molecules.tsv", moleculeHeader()),
          hard_(outDir, "hard_molecules.tsv", moleculeHeader()),
          gated_(outDir, "gated_hard_molecules.tsv", moleculeHeader()),
          soft_(outDir, "soft_expected_molecules.tsv",
              "umi_mode\tfeature_id\tcorrected_umi\tcandidate\texpected_count\tread_clique_ids")
    {
        soft_.stream() << std::setprecision(17);
    }

    void resolve(const GroupRows &group)
    {
        if (group.size() == 0) return;
        const Row &first = group.at(0);
        const std::string statsKey = first.umiMode + "." + first.product;
        Stats &summary = stats_[statsKey];
        ++summary.groups;
        for (std::size_t left = 0; left < group.size(); ++left) {
            const Row &row = group.at(left);
            if (row.umiMode != first.umiMode || row.product != first.product
                || row.candidate != first.candidate || row.correctedUmi != first.correctedUmi) {
                throw std::runtime_error("hash partition mixed reconciliation keys");
            }
            for (std::size_t right = 0; right < left; ++right) {
                if (row.featureId == group.at(right).featureId) {
                    throw std::runtime_error("duplicate gene in provisional reconciliation group");
                }
            }
        }

        std::size_t winner = 0;
        if (first.product == "soft_expected") {
            bool tied = false;
            for (std::size_t index = 0; index < group.size(); ++index) {
                const Row &row = group.at(index);
                summary.inputExpectedMass += row.expectedCount;
                if (index == 0) continue;
                if (!nearlyEqual(row.expectedCount, group.at(winner).expectedCount)
                    && row.expectedCount > group.at(winner).expectedCount) {
                    winner = index;
                    tied = false;
                } else if (nearlyEqual(row.expectedCount, group.at(winner).expectedCount)) {
                    tied = true;
                }
            }
            if (tied) {
                ++summary.correctedCountTies;
                return;
            }
            for (std::size_t index = 0; index < group.size(); ++index) {
                if (!nearlyEqual(group.at(index).originalExpected,
                                 group.at(winner).originalExpected)
                    && group.at(index).originalExpected > group.at(winner).originalExpected) {
                    ++summary.originalDominanceRejected;
                    return;
                }
            }
            ++summary.accepted;
            summary.outputExpectedMass += group.at(winner).expectedCount;
            const Row &row = group.at(winner);
            soft_.stream() << row.umiMode << '\t' << row.featureId << '\t'
                << row.correctedUmi << '\t' << row.candidate << '\t'
                << row.expectedCount << '\t' << row.readCliqueIds << '\n';
            return;
        }

        std::vector<multi_gene_umi_cr::GeneSupport> supports;
        supports.reserve(group.size());
        for (std::size_t index = 0; index < group.size(); ++index) {
            supports.push_back(multi_gene_umi_cr::GeneSupport(
                static_cast<std::uint32_t>(index), group.at(index).correctedCount,
                group.at(index).originalCount));
        }
        const multi_gene_umi_cr::Result resolution = multi_gene_umi_cr::resolve(supports);
        if (!resolution.accepted) {
            if (resolution.reason == "corrected_count_tie") ++summary.correctedCountTies;
            else if (resolution.reason == "original_umi_dominance") {
                ++summary.originalDominanceRejected;
            }
            return;
        }
        ++summary.accepted;
        const Row &row = group.at(resolution.gene);
        std::ostream *output = row.product == "strict" ? &strict_.stream()
            : row.product == "hard" ? &hard_.stream() : &gated_.stream();
        *output << row.umiMode << '\t' << row.product << '\t' << row.moleculeId << '\t'
            << row.featureId << '\t' << row.correctedUmi << '\t' << row.candidate << '\t'
            << row.memberReadCount << '\t' << row.memberReadIds << '\t'
            << row.readCliqueIds << '\n';
    }

    void commit(const Arguments &arguments, std::uint64_t inputRows)
    {
        AtomicTable config(arguments.outDir, "resolved_config.tsv", "key\tvalue");
        config.stream() << "schema\tstar_suite.molecule_first.v1\n"
            << "star_suite_version\t" << STAR_SUITE_VERSION << '\n'
            << "execution_mode\tgex_candidate_umi_hash_partitioned\n"
            << "input_order\tfeature_sorted\n"
            << "output_order\tdeterministic_hash_partition_then_first_seen\n"
            << "gex_multigene_umi_cr\t1\n"
            << "gex_reconciliation_stage\tfinal\n"
            << "gex_hash_partitions\t" << arguments.partitions << '\n';
        AtomicTable summary(arguments.outDir, "summary.tsv", "metric\tvalue");
        summary.stream() << std::setprecision(17)
            << "hash.input_rows\t" << inputRows << '\n'
            << "hash.partitions\t" << arguments.partitions << '\n';
        for (const auto &entry : stats_) {
            summary.stream() << entry.first << ".multigene_groups\t" << entry.second.groups << '\n'
                << entry.first << ".multigene_accepted\t" << entry.second.accepted << '\n'
                << entry.first << ".multigene_corrected_count_ties\t"
                << entry.second.correctedCountTies << '\n'
                << entry.first << ".multigene_original_umi_dominance_rejected\t"
                << entry.second.originalDominanceRejected << '\n';
            if (entry.first.find(".soft_expected") != std::string::npos) {
                summary.stream() << entry.first << ".multigene_input_expected_mass\t"
                    << entry.second.inputExpectedMass << '\n'
                    << entry.first << ".multigene_output_expected_mass\t"
                    << entry.second.outputExpectedMass << '\n';
            }
        }
        strict_.commit();
        hard_.commit();
        gated_.commit();
        soft_.commit();
        config.commit();
        summary.commit();
    }

  private:
    static std::string moleculeHeader()
    {
        return "umi_mode\tproduct\tmolecule_id\tfeature_id\tcorrected_umi\tcandidate\t"
               "member_read_count\tmember_read_ids\tread_clique_ids";
    }

    AtomicTable strict_;
    AtomicTable hard_;
    AtomicTable gated_;
    AtomicTable soft_;
    std::map<std::string, Stats> stats_;
};

std::uint64_t partitionInput(const Arguments &arguments, PartitionFiles &partitions)
{
    std::ifstream input(arguments.input.c_str());
    if (!input) throw std::runtime_error("cannot open input: " + arguments.input);
    std::string line;
    if (!std::getline(input, line) || line != kHeader) {
        throw std::runtime_error("provisional input header mismatch");
    }
    std::uint64_t rows = 0;
    std::string previousFeature;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t end = keyEnd(line);
        const std::string feature = featureField(line, end);
        if (!previousFeature.empty() && feature < previousFeature) {
            throw std::runtime_error("provisional input is not feature-sorted at data row "
                                     + std::to_string(rows + 1));
        }
        previousFeature = feature;
        const std::size_t partition = static_cast<std::size_t>(
            stableHash(line.data(), end) % arguments.partitions);
        partitions.write(partition, line);
        ++rows;
        if (rows % UINT64_C(100000000) == 0) {
            std::cerr << "hash partitioned rows=" << rows << '\n';
        }
    }
    if (!input.eof()) throw std::runtime_error("failed reading provisional support input");
    partitions.close();
    return rows;
}

void processPartition(const std::string &path, Reconciler &reconciler,
                      std::uint64_t &rows, std::uint64_t &groups)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        throw std::runtime_error("cannot stat hash partition: " + path);
    }
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input) throw std::runtime_error("cannot read hash partition: " + path);

    typedef std::unordered_map<std::string, GroupRows> GroupMap;
    GroupMap byKey;
    const std::uint64_t estimatedRows = static_cast<std::uint64_t>(info.st_size) / 180 + 1;
    byKey.reserve(static_cast<std::size_t>(std::min<std::uint64_t>(estimatedRows, 20000000)));
    std::vector<GroupRows *> insertionOrder;
    insertionOrder.reserve(static_cast<std::size_t>(std::min<std::uint64_t>(estimatedRows, 20000000)));

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t end = keyEnd(line);
        const std::string groupKey = line.substr(0, end);
        std::pair<GroupMap::iterator, bool> result = byKey.emplace(groupKey, GroupRows());
        if (result.second) insertionOrder.push_back(&result.first->second);
        result.first->second.add(parseRow(splitTabs(line)));
        ++rows;
    }
    if (!input.eof()) throw std::runtime_error("failed reading hash partition: " + path);
    groups += insertionOrder.size();
    for (GroupRows *group : insertionOrder) reconciler.resolve(*group);
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Arguments arguments = parseArguments(argc, argv);
        prepareOutputDirectory(arguments.outDir);
        PartitionFiles partitions(arguments);
        const std::uint64_t inputRows = partitionInput(arguments, partitions);
        Reconciler reconciler(arguments.outDir);
        std::uint64_t processedRows = 0;
        std::uint64_t groups = 0;
        for (std::size_t index = 0; index < partitions.paths().size(); ++index) {
            processPartition(partitions.paths()[index], reconciler, processedRows, groups);
            partitions.remove(index);
            if ((index + 1) % 8 == 0 || index + 1 == partitions.paths().size()) {
                std::cerr << "hash reconciled partitions=" << (index + 1)
                          << '/' << partitions.paths().size() << " rows=" << processedRows
                          << " groups=" << groups << '\n';
            }
        }
        if (processedRows != inputRows) {
            throw std::runtime_error("hash partition row conservation failure");
        }
        reconciler.commit(arguments, inputRows);
        linkRequired(arguments.resolvedSource + "/read_cliques.tsv",
                     arguments.outDir + "/read_cliques.tsv");
        linkRequired(arguments.resolvedSource + "/hard_call_audit.tsv",
                     arguments.outDir + "/hard_call_audit.tsv");
        partitions.finish();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "molecule_first_gex_reconcile_hash: " << error.what() << '\n';
        return 1;
    }
}

#include "MultiGeneUmiCr.h"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unistd.h>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef STAR_SUITE_VERSION
#define STAR_SUITE_VERSION "unknown"
#endif

namespace {

struct Arguments {
    std::string input;
    std::string resolvedSource;
    std::string outDir;
};

void usage(std::ostream &out)
{
    out << "Usage: molecule_first_gex_reconcile --input sorted_support.tsv "
        << "--resolved-source DIR --out-dir DIR\n";
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
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + option);
        }
        const std::string value = argv[++index];
        if (option == "--input") arguments.input = value;
        else if (option == "--resolved-source") arguments.resolvedSource = value;
        else if (option == "--out-dir") arguments.outDir = value;
        else throw std::invalid_argument("unknown option: " + option);
    }
    if (arguments.input.empty() || arguments.resolvedSource.empty()
        || arguments.outDir.empty()) {
        throw std::invalid_argument("--input, --resolved-source, and --out-dir are required");
    }
    return arguments;
}

std::vector<std::string> splitTabs(const std::string &line)
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
    if (directory == nullptr) throw std::runtime_error("cannot inspect output directory: " + path);
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

using GroupKey = std::tuple<std::string, std::string, std::string, std::string>;

GroupKey key(const Row &row)
{
    return std::make_tuple(row.umiMode, row.product, row.candidate, row.correctedUmi);
}

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

struct Stats {
    std::uint64_t groups = 0;
    std::uint64_t accepted = 0;
    std::uint64_t correctedCountTies = 0;
    std::uint64_t originalDominanceRejected = 0;
    double inputExpectedMass = 0.0;
    double outputExpectedMass = 0.0;
};

void linkRequired(const std::string &source, const std::string &destination)
{
    if (::link(source.c_str(), destination.c_str()) != 0) {
        throw std::runtime_error("cannot hard-link " + source + " to " + destination
                                 + ": " + std::strerror(errno));
    }
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Arguments arguments = parseArguments(argc, argv);
        prepareOutputDirectory(arguments.outDir);
        std::ifstream input(arguments.input.c_str());
        if (!input) throw std::runtime_error("cannot open input: " + arguments.input);
        std::string line;
        if (!std::getline(input, line)) throw std::runtime_error("provisional input is empty");
        const std::vector<std::string> expected = {
            "umi_mode", "product", "candidate", "corrected_umi", "feature_id",
            "corrected_count", "original_at_corrected_count", "expected_count",
            "original_expected_count", "molecule_id", "member_read_count",
            "member_read_ids", "read_clique_ids"
        };
        if (splitTabs(line) != expected) {
            throw std::runtime_error("provisional input header mismatch");
        }

        const std::string moleculeHeader =
            "umi_mode\tproduct\tmolecule_id\tfeature_id\tcorrected_umi\tcandidate\t"
            "member_read_count\tmember_read_ids\tread_clique_ids";
        AtomicTable strict(arguments.outDir, "strict_molecules.tsv", moleculeHeader);
        AtomicTable hard(arguments.outDir, "hard_molecules.tsv", moleculeHeader);
        AtomicTable gated(arguments.outDir, "gated_hard_molecules.tsv", moleculeHeader);
        AtomicTable soft(arguments.outDir, "soft_expected_molecules.tsv",
            "umi_mode\tfeature_id\tcorrected_umi\tcandidate\texpected_count\tread_clique_ids");
        soft.stream() << std::setprecision(17);
        std::map<std::string, Stats> stats;
        std::vector<Row> group;
        GroupKey current;
        bool haveCurrent = false;
        std::uint64_t lineNumber = 1;

        auto flush = [&]() {
            if (group.empty()) return;
            const std::string statsKey = group.front().umiMode + "." + group.front().product;
            Stats &summary = stats[statsKey];
            ++summary.groups;
            std::set<std::string> genes;
            for (const Row &row : group) {
                if (!genes.insert(row.featureId).second) {
                    throw std::runtime_error("duplicate gene in provisional reconciliation group");
                }
            }
            std::size_t winner = 0;
            if (group.front().product == "soft_expected") {
                bool tied = false;
                for (std::size_t index = 0; index < group.size(); ++index) {
                    summary.inputExpectedMass += group[index].expectedCount;
                    if (index == 0) continue;
                    if (!nearlyEqual(group[index].expectedCount, group[winner].expectedCount)
                        && group[index].expectedCount > group[winner].expectedCount) {
                        winner = index;
                        tied = false;
                    } else if (nearlyEqual(group[index].expectedCount,
                                           group[winner].expectedCount)) {
                        tied = true;
                    }
                }
                if (tied) {
                    ++summary.correctedCountTies;
                    group.clear();
                    return;
                }
                for (const Row &row : group) {
                    if (!nearlyEqual(row.originalExpected, group[winner].originalExpected)
                        && row.originalExpected > group[winner].originalExpected) {
                        ++summary.originalDominanceRejected;
                        group.clear();
                        return;
                    }
                }
                ++summary.accepted;
                summary.outputExpectedMass += group[winner].expectedCount;
                const Row &row = group[winner];
                soft.stream() << row.umiMode << '\t' << row.featureId << '\t'
                    << row.correctedUmi << '\t' << row.candidate << '\t'
                    << row.expectedCount << '\t' << row.readCliqueIds << '\n';
                group.clear();
                return;
            }

            std::vector<multi_gene_umi_cr::GeneSupport> supports;
            for (std::size_t index = 0; index < group.size(); ++index) {
                supports.push_back(multi_gene_umi_cr::GeneSupport(
                    static_cast<std::uint32_t>(index), group[index].correctedCount,
                    group[index].originalCount));
            }
            const multi_gene_umi_cr::Result resolution = multi_gene_umi_cr::resolve(supports);
            if (!resolution.accepted) {
                if (resolution.reason == "corrected_count_tie") ++summary.correctedCountTies;
                else if (resolution.reason == "original_umi_dominance") {
                    ++summary.originalDominanceRejected;
                }
                group.clear();
                return;
            }
            ++summary.accepted;
            const Row &row = group.at(resolution.gene);
            std::ostream *output = row.product == "strict" ? &strict.stream()
                : row.product == "hard" ? &hard.stream() : &gated.stream();
            *output << row.umiMode << '\t' << row.product << '\t' << row.moleculeId << '\t'
                << row.featureId << '\t' << row.correctedUmi << '\t' << row.candidate << '\t'
                << row.memberReadCount << '\t' << row.memberReadIds << '\t'
                << row.readCliqueIds << '\n';
            group.clear();
        };

        while (std::getline(input, line)) {
            ++lineNumber;
            if (line.empty()) continue;
            const Row row = parseRow(splitTabs(line));
            const GroupKey rowKey = key(row);
            if (haveCurrent && rowKey != current) {
                if (rowKey < current) {
                    throw std::runtime_error("provisional input is out of order at line "
                                             + std::to_string(lineNumber));
                }
                flush();
            }
            if (!haveCurrent || rowKey != current) {
                current = rowKey;
                haveCurrent = true;
            }
            group.push_back(row);
        }
        flush();

        AtomicTable config(arguments.outDir, "resolved_config.tsv", "key\tvalue");
        config.stream() << "schema\tstar_suite.molecule_first.v1\n"
            << "star_suite_version\t" << STAR_SUITE_VERSION << '\n'
            << "execution_mode\tgex_candidate_umi_streaming\n"
            << "input_order\tumi_mode_product_candidate_corrected_umi_feature_id\n"
            << "gex_multigene_umi_cr\t1\n"
            << "gex_reconciliation_stage\tfinal\n";
        AtomicTable summary(arguments.outDir, "summary.tsv", "metric\tvalue");
        summary.stream() << std::setprecision(17);
        for (const auto &entry : stats) {
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

        strict.commit();
        hard.commit();
        gated.commit();
        soft.commit();
        config.commit();
        summary.commit();
        linkRequired(arguments.resolvedSource + "/read_cliques.tsv",
                     arguments.outDir + "/read_cliques.tsv");
        linkRequired(arguments.resolvedSource + "/hard_call_audit.tsv",
                     arguments.outDir + "/hard_call_audit.tsv");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "molecule_first_gex_reconcile: " << error.what() << '\n';
        return 1;
    }
}

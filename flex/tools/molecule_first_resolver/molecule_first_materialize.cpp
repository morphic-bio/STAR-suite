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
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <tuple>
#include <vector>

#ifndef STAR_SUITE_VERSION
#define STAR_SUITE_VERSION "unknown"
#endif

namespace {

struct Arguments {
    std::string resolvedDir;
    std::string outDir;
    std::string assay;
    std::string umiMode = "1mm_cr";
};

using Key = std::pair<std::string, std::string>;
using Counts = std::map<Key, double>;

struct Product {
    std::string name;
    bool expected = false;
    Counts counts;
};

struct Scale {
    std::string name;
    int factor = 1;

    Scale(const std::string &nameValue, int factorValue)
        : name(nameValue), factor(factorValue) {}
};

void usage(std::ostream &out)
{
    out << "Usage: molecule_first_materialize --resolved-dir DIR --out-dir DIR\\n"
        << "       --assay scrna|flex|visium|visium-hd [--umi-mode 1mm_cr|exact]\\n\\n"
        << "Writes strict, soft_expected, hard, and gated_hard 10x-style MEX\\n"
        << "products on identical feature/barcode axes. Soft matrices are real.\\n";
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
        if (index + 1 >= argc) throw std::invalid_argument("missing value for " + option);
        const std::string value = argv[++index];
        if (option == "--resolved-dir") arguments.resolvedDir = value;
        else if (option == "--out-dir") arguments.outDir = value;
        else if (option == "--assay") arguments.assay = value;
        else if (option == "--umi-mode") arguments.umiMode = value;
        else throw std::invalid_argument("unknown option: " + option);
    }
    if (arguments.resolvedDir.empty() || arguments.outDir.empty() || arguments.assay.empty()) {
        throw std::invalid_argument("--resolved-dir, --out-dir, and --assay are required");
    }
    const std::set<std::string> assays = {"scrna", "flex", "visium", "visium-hd"};
    if (assays.find(arguments.assay) == assays.end()) throw std::invalid_argument("unsupported --assay");
    if (arguments.umiMode != "1mm_cr" && arguments.umiMode != "exact") {
        throw std::invalid_argument("--umi-mode must be 1mm_cr or exact");
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

double finiteDouble(const std::string &value, const std::string &path)
{
    std::size_t used = 0;
    const double parsed = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(parsed) || parsed < 0.0) {
        throw std::runtime_error("invalid non-negative count in " + path + ": " + value);
    }
    return parsed;
}

bool directoryEmpty(const std::string &path)
{
    DIR *directory = opendir(path.c_str());
    if (directory == nullptr) return false;
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

void makeDirectory(const std::string &path)
{
    struct stat info;
    if (stat(path.c_str(), &info) == 0) {
        if (!S_ISDIR(info.st_mode)) throw std::runtime_error("path is not a directory: " + path);
        return;
    }
    if (errno != ENOENT || mkdir(path.c_str(), 0755) != 0) {
        throw std::runtime_error("cannot create directory " + path + ": " + std::strerror(errno));
    }
}

void makeDirectories(const std::string &path)
{
    if (path.empty()) return;
    std::size_t begin = path[0] == '/' ? 1 : 0;
    for (;;) {
        const std::size_t end = path.find('/', begin);
        const std::string prefix = path.substr(0, end);
        if (!prefix.empty()) makeDirectory(prefix);
        if (end == std::string::npos) break;
        begin = end + 1;
    }
}

void prepareOutputDirectory(const std::string &path)
{
    struct stat info;
    if (stat(path.c_str(), &info) == 0) {
        if (!S_ISDIR(info.st_mode) || !directoryEmpty(path)) {
            throw std::runtime_error("output directory must be empty: " + path);
        }
    } else {
        makeDirectories(path);
    }
}

void requireHeader(const std::vector<std::string> &observed,
                   const std::vector<std::string> &expected, const std::string &path)
{
    if (observed != expected) throw std::runtime_error("unexpected TSV schema: " + path);
}

void loadUniverse(const std::string &path, std::set<std::string> &features,
                  std::set<std::string> &candidates)
{
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("cannot open resolver artifact: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty resolver artifact: " + path);
    requireHeader(splitTabs(line), {
        "clique_id", "feature_id", "raw_umi", "member_read_count", "member_read_ids",
        "candidate", "log_sequence_likelihood_sum", "log_exact_read_prior", "log_evidence", "posterior"
    }, path);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> row = splitTabs(line);
        if (row.size() != 10) throw std::runtime_error("malformed resolver row: " + path);
        features.insert(row[1]);
        candidates.insert(row[5]);
    }
    if (features.empty() || candidates.empty()) throw std::runtime_error("resolver universe is empty");
}

void loadMolecules(const std::string &path, const std::string &umiMode,
                   const std::string &product, Counts &counts)
{
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("cannot open resolver artifact: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty resolver artifact: " + path);
    requireHeader(splitTabs(line), {
        "umi_mode", "product", "molecule_id", "feature_id", "corrected_umi", "candidate",
        "member_read_count", "member_read_ids", "read_clique_ids"
    }, path);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> row = splitTabs(line);
        if (row.size() != 9) throw std::runtime_error("malformed resolver row: " + path);
        if (row[0] == umiMode && row[1] == product) counts[Key(row[3], row[5])] += 1.0;
    }
}

void loadSoft(const std::string &path, const std::string &umiMode, Counts &counts)
{
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("cannot open resolver artifact: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty resolver artifact: " + path);
    requireHeader(splitTabs(line), {
        "umi_mode", "feature_id", "corrected_umi", "candidate", "expected_count", "read_clique_ids"
    }, path);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> row = splitTabs(line);
        if (row.size() != 6) throw std::runtime_error("malformed resolver row: " + path);
        if (row[0] == umiMode) counts[Key(row[1], row[3])] += finiteDouble(row[4], path);
    }
}

std::string hdParent(const std::string &candidate, int factor)
{
    const std::string prefix = "s_002um_";
    if (candidate.compare(0, prefix.size(), prefix) != 0) {
        throw std::runtime_error("Visium HD candidate does not use s_002um_row_col: " + candidate);
    }
    const std::string coordinates = candidate.substr(prefix.size());
    const std::size_t separator = coordinates.find('_');
    if (separator == std::string::npos || coordinates.find('_', separator + 1) != std::string::npos) {
        throw std::runtime_error("invalid Visium HD coordinate candidate: " + candidate);
    }
    const std::string rowText = coordinates.substr(0, separator);
    const std::string columnText = coordinates.substr(separator + 1);
    std::size_t rowUsed = 0, columnUsed = 0;
    const long row = std::stol(rowText, &rowUsed);
    const long column = std::stol(columnText, &columnUsed);
    if (rowUsed != rowText.size() || columnUsed != columnText.size()) {
        throw std::runtime_error("invalid Visium HD coordinate candidate: " + candidate);
    }
    if (row < 0 || column < 0) throw std::runtime_error("negative Visium HD coordinate: " + candidate);
    std::ostringstream output;
    output << "s_" << std::setw(3) << std::setfill('0') << factor * 2 << "um_"
           << row / factor << '_' << column / factor;
    return output.str();
}

std::string scaleCandidate(const std::string &candidate, int factor, bool hd)
{
    return hd ? hdParent(candidate, factor) : candidate;
}

void atomicText(const std::string &path, const std::string &text)
{
    const std::string temporary = path + ".tmp";
    std::ofstream output(temporary.c_str());
    if (!output) throw std::runtime_error("cannot create output: " + temporary);
    output << text;
    output.close();
    if (!output || std::rename(temporary.c_str(), path.c_str()) != 0) {
        throw std::runtime_error("cannot commit output: " + path);
    }
}

std::string mexBarcode(const std::string &candidate)
{
    return candidate.size() >= 2
        && candidate.compare(candidate.size() - 2, 2, "-1") == 0
        ? candidate : candidate + "-1";
}

struct MexSummary {
    std::size_t features = 0;
    std::size_t barcodes = 0;
    std::size_t nnz = 0;
    double mass = 0.0;
};

MexSummary writeMex(const std::string &directory, const std::set<std::string> &featureUniverse,
                    const std::set<std::string> &candidateUniverse, const Counts &source,
                    const Scale &scale, bool hd, bool expected)
{
    makeDirectories(directory);
    const std::vector<std::string> features(featureUniverse.begin(), featureUniverse.end());
    std::set<std::string> scaledCandidates;
    for (const std::string &candidate : candidateUniverse) {
        scaledCandidates.insert(scaleCandidate(candidate, scale.factor, hd));
    }
    const std::vector<std::string> candidates(scaledCandidates.begin(), scaledCandidates.end());
    std::map<std::string, std::size_t> featureIndex, candidateIndex;
    for (std::size_t index = 0; index < features.size(); ++index) featureIndex[features[index]] = index + 1;
    for (std::size_t index = 0; index < candidates.size(); ++index) candidateIndex[candidates[index]] = index + 1;

    Counts aggregated;
    for (const auto &entry : source) {
        if (featureUniverse.find(entry.first.first) == featureUniverse.end()
            || candidateUniverse.find(entry.first.second) == candidateUniverse.end()) {
            throw std::runtime_error("policy contains a feature/candidate outside the read-clique universe");
        }
        aggregated[Key(entry.first.first, scaleCandidate(entry.first.second, scale.factor, hd))] += entry.second;
    }
    std::vector<std::tuple<std::size_t, std::size_t, double> > rows;
    MexSummary summary;
    summary.features = features.size();
    summary.barcodes = candidates.size();
    for (const auto &entry : aggregated) {
        if (entry.second == 0.0) continue;
        rows.push_back(std::make_tuple(candidateIndex[entry.first.second], featureIndex[entry.first.first], entry.second));
        summary.mass += entry.second;
    }
    std::sort(rows.begin(), rows.end());
    summary.nnz = rows.size();

    std::ostringstream matrix;
    matrix << "%%MatrixMarket matrix coordinate " << (expected ? "real" : "integer") << " general\n"
           << "% STAR Suite molecule-first post-collapse policy matrix\n"
           << features.size() << ' ' << candidates.size() << ' ' << rows.size() << '\n'
           << std::setprecision(17);
    for (const auto &row : rows) {
        matrix << std::get<1>(row) << ' ' << std::get<0>(row) << ' ';
        if (expected) matrix << std::get<2>(row);
        else matrix << static_cast<std::uint64_t>(std::llround(std::get<2>(row)));
        matrix << '\n';
    }
    std::ostringstream featureText, barcodeText;
    for (const std::string &feature : features) featureText << feature << '\t' << feature << "\tGene Expression\n";
    for (const std::string &candidate : candidates) barcodeText << mexBarcode(candidate) << '\n';
    atomicText(directory + "/matrix.mtx", matrix.str());
    atomicText(directory + "/features.tsv", featureText.str());
    atomicText(directory + "/barcodes.tsv", barcodeText.str());
    return summary;
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Arguments arguments = parseArguments(argc, argv);
        prepareOutputDirectory(arguments.outDir);
        std::set<std::string> features, candidates;
        loadUniverse(arguments.resolvedDir + "/read_cliques.tsv", features, candidates);

        std::vector<Product> products(4);
        products[0].name = "strict";
        products[1].name = "soft_expected";
        products[1].expected = true;
        products[2].name = "hard";
        products[3].name = "gated_hard";
        loadMolecules(arguments.resolvedDir + "/strict_molecules.tsv", arguments.umiMode,
                      "strict", products[0].counts);
        loadSoft(arguments.resolvedDir + "/soft_expected_molecules.tsv", arguments.umiMode,
                 products[1].counts);
        loadMolecules(arguments.resolvedDir + "/hard_molecules.tsv", arguments.umiMode,
                      "hard", products[2].counts);
        loadMolecules(arguments.resolvedDir + "/gated_hard_molecules.tsv", arguments.umiMode,
                      "gated_hard", products[3].counts);

        const bool hd = arguments.assay == "visium-hd";
        std::vector<Scale> scales;
        if (hd) {
            scales.push_back({"square_002um", 1});
            scales.push_back({"square_008um", 4});
            scales.push_back({"square_016um", 8});
        } else {
            scales.push_back({"raw", 1});
        }

        std::ostringstream summary;
        summary << "schema\tstar_suite.molecule_first.policy_mex.v1\n"
                << "star_suite_version\t" << STAR_SUITE_VERSION << '\n'
                << "assay\t" << arguments.assay << '\n'
                << "umi_mode\t" << arguments.umiMode << '\n'
                << "axes_source\tread_cliques.tsv\n"
                << "cell_calling_order\tafter_postcollapse_integer_matrix_only\n"
                << "soft_cell_calling_allowed\tfalse\n"
                << "product\tscale\tfeatures\tbarcodes\tnnz\tmass\tmatrix_field\n";
        for (const Product &product : products) {
            std::vector<double> masses;
            for (const Scale &scale : scales) {
                const MexSummary result = writeMex(
                    arguments.outDir + "/" + product.name + "/" + scale.name,
                    features, candidates, product.counts, scale, hd, product.expected);
                masses.push_back(result.mass);
                summary << product.name << '\t' << scale.name << '\t' << result.features << '\t'
                        << result.barcodes << '\t' << result.nnz << '\t' << std::setprecision(17)
                        << result.mass << '\t' << (product.expected ? "real" : "integer") << '\n';
            }
            if (hd) {
                const double minimum = *std::min_element(masses.begin(), masses.end());
                const double maximum = *std::max_element(masses.begin(), masses.end());
                const double tolerance = 1e-9 * std::max(1.0, maximum);
                if (maximum - minimum > tolerance) {
                    throw std::runtime_error("Visium HD scale aggregation did not conserve mass for " + product.name);
                }
            }
        }
        atomicText(arguments.outDir + "/summary.tsv", summary.str());
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "molecule_first_materialize: ERROR: " << error.what() << '\n';
        return 2;
    }
}

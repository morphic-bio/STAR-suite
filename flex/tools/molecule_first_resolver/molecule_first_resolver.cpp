#include "solo/MoleculeFirstResolver.h"

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
#include <unordered_map>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef STAR_SUITE_VERSION
#define STAR_SUITE_VERSION "unknown"
#endif

namespace {

struct Arguments {
    std::string input;
    std::string outDir;
    bool inputFeatureSorted = false;
    bool gexMultiGeneUmiCr = false;
    bool gexProvisionalFeatureSorted = false;
    molecule_first::Config config;
};

void usage(std::ostream &out)
{
    out << "Usage: molecule_first_resolver --input candidates.tsv --out-dir DIR [options]\n"
        << "\nRequired input columns:\n"
        << "  read_id feature_id raw_umi candidate log_sequence_likelihood exact_read_count\n"
        << "\nOptions:\n"
        << "  --temperature FLOAT          default 1\n"
        << "  --prior-alpha FLOAT          default 1\n"
        << "  --prior-beta FLOAT           default 1\n"
        << "  --gate-min-posterior FLOAT   default 0.95\n"
        << "  --gate-min-margin FLOAT      default 0.90\n"
        << "  --input-feature-sorted       stream feature_id/read_id/candidate sorted input\n"
        << "  --gex-multigene-umi-cr       reconcile genes after candidate-specific UMI correction\n"
        << "  --gex-provisional-feature-sorted\n"
        << "                               emit bounded per-gene supports for a later\n"
        << "                               candidate/corrected-UMI reconciliation pass\n"
        << "  --version\n";
}

double parseDouble(const std::string &text, const std::string &option)
{
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    if (used != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument("invalid finite number for " + option + ": " + text);
    }
    return value;
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
        if (option == "--gex-multigene-umi-cr") {
            arguments.gexMultiGeneUmiCr = true;
            continue;
        }
        if (option == "--gex-provisional-feature-sorted") {
            arguments.gexProvisionalFeatureSorted = true;
            continue;
        }
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + option);
        }
        const std::string value = argv[++index];
        if (option == "--input") {
            arguments.input = value;
        } else if (option == "--out-dir") {
            arguments.outDir = value;
        } else if (option == "--temperature") {
            arguments.config.temperature = parseDouble(value, option);
        } else if (option == "--prior-alpha") {
            arguments.config.priorAlpha = parseDouble(value, option);
        } else if (option == "--prior-beta") {
            arguments.config.priorBeta = parseDouble(value, option);
        } else if (option == "--gate-min-posterior") {
            arguments.config.gateMinPosterior = parseDouble(value, option);
        } else if (option == "--gate-min-margin") {
            arguments.config.gateMinMargin = parseDouble(value, option);
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }
    if (arguments.input.empty() || arguments.outDir.empty()) {
        throw std::invalid_argument("--input and --out-dir are required");
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
        if (end == std::string::npos) {
            return fields;
        }
        begin = end + 1;
    }
}

std::string join(const std::vector<std::string> &values, char separator)
{
    std::ostringstream output;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index > 0) {
            output << separator;
        }
        output << values[index];
    }
    return output.str();
}

std::uint64_t parseCount(const std::string &text)
{
    if (text.empty() || text[0] == '-') {
        throw std::invalid_argument("invalid non-negative exact_read_count: " + text);
    }
    std::size_t used = 0;
    const unsigned long long value = std::stoull(text, &used);
    if (used != text.size()) {
        throw std::invalid_argument("invalid non-negative exact_read_count: " + text);
    }
    return static_cast<std::uint64_t>(value);
}

struct LoadedInput {
    std::vector<molecule_first::CandidateRead> reads;
    molecule_first::PriorCounts priorCounts;
    std::uint64_t candidateRows = 0;
};

LoadedInput loadInput(const std::string &path)
{
    std::ifstream input(path.c_str());
    if (!input) {
        throw std::runtime_error("cannot open input: " + path);
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("input is empty: " + path);
    }
    const std::vector<std::string> expected = {
        "read_id", "feature_id", "raw_umi", "candidate",
        "log_sequence_likelihood", "exact_read_count"
    };
    if (splitTabs(line) != expected) {
        throw std::runtime_error("input header does not match the molecule-first v1 schema");
    }

    LoadedInput loaded;
    std::map<std::string, std::size_t> readPosition;
    std::map<std::string, std::set<std::string> > seenCandidates;
    std::uint64_t lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (line.empty()) {
            continue;
        }
        const std::vector<std::string> fields = splitTabs(line);
        if (fields.size() != expected.size()) {
            throw std::runtime_error("wrong field count at input line " + std::to_string(lineNumber));
        }
        for (std::size_t index = 0; index < 4; ++index) {
            if (fields[index].empty()) {
                throw std::runtime_error("empty required value at input line " + std::to_string(lineNumber));
            }
        }
        const double likelihood = parseDouble(fields[4], "log_sequence_likelihood");
        const std::uint64_t count = parseCount(fields[5]);
        const auto prior = loaded.priorCounts.find(fields[3]);
        if (prior != loaded.priorCounts.end() && prior->second != count) {
            throw std::runtime_error("inconsistent exact_read_count for candidate " + fields[3]);
        }
        loaded.priorCounts[fields[3]] = count;

        std::size_t position = 0;
        const auto found = readPosition.find(fields[0]);
        if (found == readPosition.end()) {
            position = loaded.reads.size();
            readPosition[fields[0]] = position;
            molecule_first::CandidateRead read;
            read.readId = fields[0];
            read.featureId = fields[1];
            read.rawUmi = fields[2];
            loaded.reads.push_back(read);
        } else {
            position = found->second;
            if (loaded.reads[position].featureId != fields[1]
                || loaded.reads[position].rawUmi != fields[2]) {
                throw std::runtime_error("inconsistent feature/UMI for read " + fields[0]);
            }
        }
        if (!seenCandidates[fields[0]].insert(fields[3]).second) {
            throw std::runtime_error("duplicate read/candidate row at input line "
                                     + std::to_string(lineNumber));
        }
        molecule_first::CandidateScore score;
        score.candidate = fields[3];
        score.logLikelihood = likelihood;
        loaded.reads[position].scores.push_back(score);
        ++loaded.candidateRows;
    }
    if (loaded.reads.empty()) {
        throw std::runtime_error("input contains no candidate reads");
    }
    return loaded;
}

bool directoryEmpty(const std::string &path)
{
    DIR *directory = opendir(path.c_str());
    if (directory == nullptr) {
        throw std::runtime_error("cannot inspect output directory: " + path);
    }
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
        throw std::runtime_error("cannot create output directory " + path + ": " + std::strerror(errno));
    }
}

class AtomicTable {
  public:
    AtomicTable(const std::string &directory, const std::string &name, const std::string &header)
        : finalPath_(directory + "/" + name), tempPath_(finalPath_ + ".tmp"), output_(tempPath_.c_str())
    {
        if (!output_) {
            throw std::runtime_error("cannot create output: " + tempPath_);
        }
        output_ << header << '\n';
    }

    std::ostream &stream() { return output_; }

    void commit()
    {
        output_.close();
        if (!output_ || std::rename(tempPath_.c_str(), finalPath_.c_str()) != 0) {
            throw std::runtime_error("cannot commit output: " + finalPath_);
        }
        committed_ = true;
    }

    ~AtomicTable()
    {
        if (!committed_) {
            output_.close();
        }
    }

  private:
    std::string finalPath_;
    std::string tempPath_;
    std::ofstream output_;
    bool committed_ = false;
};

void writeMolecules(AtomicTable &table, const std::vector<molecule_first::Molecule> &molecules)
{
    for (const molecule_first::Molecule &row : molecules) {
        table.stream() << row.umiMode << '\t' << row.product << '\t' << row.moleculeId << '\t'
                       << row.featureId << '\t' << row.correctedUmi << '\t' << row.candidate << '\t'
                       << row.memberReadIds.size() << '\t' << join(row.memberReadIds, ';') << '\t'
                       << join(row.readCliqueIds, ';') << '\n';
    }
}

struct StreamingSummary {
    std::uint64_t candidateRows = 0;
    std::uint64_t candidateReads = 0;
    std::uint64_t readCliques = 0;
    std::uint64_t gatedAssigned = 0;
    std::uint64_t featureChunks = 0;
    double posteriorMass = 0.0;
    std::map<std::string, std::uint64_t> productCounts;
    std::map<std::string, double> occupancyMass;
};

void processFeatureChunk(
    const std::vector<molecule_first::CandidateRead> &reads,
    const molecule_first::PriorCounts &priorCounts,
    const Arguments &arguments,
    AtomicTable &cliqueTable,
    AtomicTable &auditTable,
    AtomicTable &strictTable,
    AtomicTable &hardTable,
    AtomicTable &gatedTable,
    AtomicTable &softTable,
    AtomicTable *gexProvisionalTable,
    StreamingSummary &summary)
{
    if (reads.empty()) {
        return;
    }
    const std::vector<molecule_first::ReadClique> cliques =
        molecule_first::buildReadCliques(reads, priorCounts, arguments.config);
    const std::vector<molecule_first::HardCall> calls =
        molecule_first::gatedHardCalls(cliques, arguments.config);

    summary.candidateReads += reads.size();
    summary.readCliques += cliques.size();
    ++summary.featureChunks;
    for (const molecule_first::ReadClique &clique : cliques) {
        for (std::size_t index = 0; index < clique.candidates.size(); ++index) {
            summary.posteriorMass += clique.posterior[index];
            cliqueTable.stream() << clique.cliqueId << '\t' << clique.featureId << '\t'
                << clique.rawUmi << '\t' << clique.memberReadIds.size() << '\t'
                << join(clique.memberReadIds, ';') << '\t' << clique.candidates[index] << '\t'
                << clique.logLikelihoodSums[index] << '\t' << clique.logReadPriors[index] << '\t'
                << clique.logEvidence[index] << '\t' << clique.posterior[index] << '\n';
        }
    }
    for (const molecule_first::HardCall &call : calls) {
        summary.gatedAssigned += call.assigned ? 1 : 0;
        auditTable.stream() << call.cliqueId << '\t' << (call.assigned ? "assigned" : "deferred")
            << '\t' << call.candidate << '\t' << call.posterior << '\t' << call.margin
            << '\t' << call.reason << '\n';
    }

    for (const std::string &umiMode : {std::string("1mm_cr"), std::string("exact")}) {
        if (gexProvisionalTable != nullptr) {
            std::map<std::string, const molecule_first::ReadClique *> cliqueById;
            for (const molecule_first::ReadClique &clique : cliques) {
                cliqueById[clique.cliqueId] = &clique;
            }
            for (const std::string &product : {
                    std::string("strict"), std::string("hard"),
                    std::string("gated_hard")}) {
                const std::vector<molecule_first::Molecule> molecules =
                    molecule_first::gexPolicyMolecules(
                        cliques, umiMode, product, calls, nullptr);
                summary.productCounts[umiMode + "." + product] += molecules.size();
                for (const molecule_first::Molecule &row : molecules) {
                    std::uint64_t originalAtCorrected = 0;
                    for (const std::string &cliqueId : row.readCliqueIds) {
                        const auto found = cliqueById.find(cliqueId);
                        if (found == cliqueById.end()) {
                            throw std::logic_error("provisional molecule lost a read clique");
                        }
                        if (found->second->rawUmi == row.correctedUmi) {
                            originalAtCorrected += found->second->memberReadIds.size();
                        }
                    }
                    gexProvisionalTable->stream()
                        << row.umiMode << '\t' << row.product << '\t'
                        << row.candidate << '\t' << row.correctedUmi << '\t'
                        << row.featureId << '\t' << row.memberReadIds.size() << '\t'
                        << originalAtCorrected << "\t\t\t"
                        << row.moleculeId << '\t' << row.memberReadIds.size() << '\t'
                        << join(row.memberReadIds, ';') << '\t'
                        << join(row.readCliqueIds, ';') << '\n';
                }
            }

            const molecule_first::UmiCorrections corrections =
                molecule_first::correctedUmis(cliques, umiMode);
            const std::vector<molecule_first::Occupancy> occupancies =
                molecule_first::gexWeightedOccupancies(cliques, umiMode, nullptr);
            summary.productCounts[umiMode + ".soft_rows"] += occupancies.size();
            for (const molecule_first::Occupancy &row : occupancies) {
                double absent = 1.0;
                for (const molecule_first::ReadClique &clique : cliques) {
                    if (clique.featureId != row.featureId
                        || clique.rawUmi != row.correctedUmi) {
                        continue;
                    }
                    for (std::size_t index = 0; index < clique.candidates.size(); ++index) {
                        if (clique.candidates[index] == row.candidate
                            && corrections.at(std::make_tuple(
                                clique.featureId, row.candidate, clique.rawUmi))
                                == row.correctedUmi) {
                            absent *= 1.0 - clique.posterior[index];
                        }
                    }
                }
                const double originalExpected = 1.0 - absent;
                summary.occupancyMass[umiMode] += row.expectedCount;
                gexProvisionalTable->stream()
                    << row.umiMode << "\tsoft_expected\t" << row.candidate << '\t'
                    << row.correctedUmi << '\t' << row.featureId << "\t\t\t"
                    << row.expectedCount << '\t' << originalExpected << "\t\t\t\t"
                    << join(row.readCliqueIds, ';') << '\n';
            }
            continue;
        }
        const std::vector<molecule_first::Molecule> strict =
            molecule_first::policyMolecules(cliques, umiMode, "strict", calls);
        const std::vector<molecule_first::Molecule> hard =
            molecule_first::policyMolecules(cliques, umiMode, "hard", calls);
        const std::vector<molecule_first::Molecule> gated =
            molecule_first::policyMolecules(cliques, umiMode, "gated_hard", calls);
        const std::vector<molecule_first::Occupancy> soft =
            molecule_first::weightedOccupancies(cliques, umiMode);
        writeMolecules(strictTable, strict);
        writeMolecules(hardTable, hard);
        writeMolecules(gatedTable, gated);
        summary.productCounts[umiMode + ".strict"] += strict.size();
        summary.productCounts[umiMode + ".hard"] += hard.size();
        summary.productCounts[umiMode + ".gated_hard"] += gated.size();
        summary.productCounts[umiMode + ".soft_rows"] += soft.size();
        for (const molecule_first::Occupancy &row : soft) {
            summary.occupancyMass[umiMode] += row.expectedCount;
            softTable.stream() << row.umiMode << '\t' << row.featureId << '\t'
                << row.correctedUmi << '\t' << row.candidate << '\t' << row.expectedCount
                << '\t' << join(row.readCliqueIds, ';') << '\n';
        }
    }
}

int runFeatureSorted(const Arguments &arguments)
{
    std::ifstream input(arguments.input.c_str());
    if (!input) {
        throw std::runtime_error("cannot open input: " + arguments.input);
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("input is empty: " + arguments.input);
    }
    const std::vector<std::string> expected = {
        "read_id", "feature_id", "raw_umi", "candidate",
        "log_sequence_likelihood", "exact_read_count"
    };
    if (splitTabs(line) != expected) {
        throw std::runtime_error("input header does not match the molecule-first v1 schema");
    }

    prepareOutputDirectory(arguments.outDir);
    AtomicTable cliqueTable(arguments.outDir, "read_cliques.tsv",
        "clique_id\tfeature_id\traw_umi\tmember_read_count\tmember_read_ids\tcandidate\tlog_sequence_likelihood_sum\tlog_exact_read_prior\tlog_evidence\tposterior");
    AtomicTable auditTable(arguments.outDir, "hard_call_audit.tsv",
        "clique_id\tstatus\tcandidate\tposterior\tmargin\treason");
    const std::string moleculeHeader =
        "umi_mode\tproduct\tmolecule_id\tfeature_id\tcorrected_umi\tcandidate\tmember_read_count\tmember_read_ids\tread_clique_ids";
    AtomicTable strictTable(arguments.outDir, "strict_molecules.tsv", moleculeHeader);
    AtomicTable hardTable(arguments.outDir, "hard_molecules.tsv", moleculeHeader);
    AtomicTable gatedTable(arguments.outDir, "gated_hard_molecules.tsv", moleculeHeader);
    AtomicTable softTable(arguments.outDir, "soft_expected_molecules.tsv",
        "umi_mode\tfeature_id\tcorrected_umi\tcandidate\texpected_count\tread_clique_ids");
    AtomicTable gexProvisionalTable(arguments.outDir, "gex_provisional_support.tsv",
        "umi_mode\tproduct\tcandidate\tcorrected_umi\tfeature_id\tcorrected_count\t"
        "original_at_corrected_count\texpected_count\toriginal_expected_count\t"
        "molecule_id\tmember_read_count\tmember_read_ids\tread_clique_ids");
    cliqueTable.stream() << std::setprecision(17);
    auditTable.stream() << std::setprecision(17);
    softTable.stream() << std::setprecision(17);
    gexProvisionalTable.stream() << std::setprecision(17);

    StreamingSummary summary;
    std::unordered_map<std::string, std::uint64_t> globalPriorCounts;
    std::vector<molecule_first::CandidateRead> featureReads;
    molecule_first::PriorCounts featurePriorCounts;
    molecule_first::CandidateRead currentRead;
    std::set<std::string> currentCandidates;
    std::string currentFeature;
    std::string previousReadId;
    std::uint64_t lineNumber = 1;

    auto flushRead = [&]() {
        if (!currentRead.readId.empty()) {
            featureReads.push_back(currentRead);
            currentRead = molecule_first::CandidateRead();
            currentCandidates.clear();
        }
    };
    auto flushFeature = [&]() {
        flushRead();
        processFeatureChunk(
            featureReads, featurePriorCounts, arguments,
            cliqueTable, auditTable, strictTable, hardTable, gatedTable, softTable,
            arguments.gexProvisionalFeatureSorted ? &gexProvisionalTable : nullptr,
            summary);
        featureReads.clear();
        featurePriorCounts.clear();
        previousReadId.clear();
    };

    while (std::getline(input, line)) {
        ++lineNumber;
        if (line.empty()) {
            continue;
        }
        const std::vector<std::string> fields = splitTabs(line);
        if (fields.size() != expected.size()) {
            throw std::runtime_error("wrong field count at input line " + std::to_string(lineNumber));
        }
        for (std::size_t index = 0; index < 4; ++index) {
            if (fields[index].empty()) {
                throw std::runtime_error("empty required value at input line " + std::to_string(lineNumber));
            }
        }
        if (!currentFeature.empty() && fields[1] != currentFeature) {
            if (fields[1] < currentFeature) {
                throw std::runtime_error("feature-sorted input is out of order at line "
                                         + std::to_string(lineNumber));
            }
            flushFeature();
            currentFeature = fields[1];
        } else if (currentFeature.empty()) {
            currentFeature = fields[1];
        }
        if (!currentRead.readId.empty() && fields[0] != currentRead.readId) {
            if (fields[0] < previousReadId) {
                throw std::runtime_error("read-sorted input is out of order at line "
                                         + std::to_string(lineNumber));
            }
            flushRead();
        }
        if (currentRead.readId.empty()) {
            currentRead.readId = fields[0];
            currentRead.featureId = fields[1];
            currentRead.rawUmi = fields[2];
            previousReadId = fields[0];
        } else if (currentRead.featureId != fields[1] || currentRead.rawUmi != fields[2]) {
            throw std::runtime_error("inconsistent feature/UMI for read " + fields[0]);
        }
        if (!currentCandidates.insert(fields[3]).second) {
            throw std::runtime_error("duplicate read/candidate row at input line "
                                     + std::to_string(lineNumber));
        }
        const double likelihood = parseDouble(fields[4], "log_sequence_likelihood");
        const std::uint64_t count = parseCount(fields[5]);
        const auto globalPrior = globalPriorCounts.find(fields[3]);
        if (globalPrior != globalPriorCounts.end() && globalPrior->second != count) {
            throw std::runtime_error("inconsistent exact_read_count for candidate " + fields[3]);
        }
        globalPriorCounts[fields[3]] = count;
        featurePriorCounts[fields[3]] = count;
        molecule_first::CandidateScore score;
        score.candidate = fields[3];
        score.logLikelihood = likelihood;
        currentRead.scores.push_back(score);
        ++summary.candidateRows;
    }
    flushFeature();
    if (summary.candidateReads == 0) {
        throw std::runtime_error("input contains no candidate reads");
    }

    AtomicTable configTable(arguments.outDir, "resolved_config.tsv", "key\tvalue");
    configTable.stream() << std::setprecision(17)
        << "schema\tstar_suite.molecule_first.v1\n"
        << "star_suite_version\t" << STAR_SUITE_VERSION << '\n'
        << "execution_mode\tfeature_sorted_streaming\n"
        << "input_order\tfeature_id_read_id_candidate\n"
        << "gex_reconciliation_stage\t"
        << (arguments.gexProvisionalFeatureSorted ? "provisional" : "disabled") << '\n'
        << "temperature\t" << arguments.config.temperature << '\n'
        << "prior_alpha\t" << arguments.config.priorAlpha << '\n'
        << "prior_beta\t" << arguments.config.priorBeta << '\n'
        << "prior_units\traw_exact_reads_including_pcr_amplification\n"
        << "prior_application\tonce_per_read_clique\n"
        << "spatial_lambda\t0\n"
        << "gate_min_posterior\t" << arguments.config.gateMinPosterior << '\n'
        << "gate_min_margin\t" << arguments.config.gateMinMargin << '\n';

    AtomicTable summaryTable(arguments.outDir, "summary.tsv", "metric\tvalue");
    summaryTable.stream() << std::setprecision(17)
        << "candidate_rows\t" << summary.candidateRows << '\n'
        << "candidate_reads\t" << summary.candidateReads << '\n'
        << "feature_chunks\t" << summary.featureChunks << '\n'
        << "read_cliques\t" << summary.readCliques << '\n'
        << "posterior_mass\t" << summary.posteriorMass << '\n'
        << "gated_assigned_cliques\t" << summary.gatedAssigned << '\n'
        << "gated_deferred_cliques\t" << summary.readCliques - summary.gatedAssigned << '\n';
    for (const auto &entry : summary.productCounts) {
        summaryTable.stream() << entry.first << "_count\t" << entry.second << '\n';
    }
    for (const auto &entry : summary.occupancyMass) {
        summaryTable.stream() << entry.first << ".soft_occupancy_mass\t" << entry.second << '\n';
        summaryTable.stream() << entry.first << ".deduplicated_mass\t"
                              << summary.posteriorMass - entry.second << '\n';
    }

    cliqueTable.commit();
    auditTable.commit();
    strictTable.commit();
    hardTable.commit();
    gatedTable.commit();
    softTable.commit();
    gexProvisionalTable.commit();
    configTable.commit();
    summaryTable.commit();
    return 0;
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Arguments arguments = parseArguments(argc, argv);
        if (arguments.gexMultiGeneUmiCr && arguments.inputFeatureSorted) {
            throw std::invalid_argument(
                "--gex-multigene-umi-cr requires the global input mode so genes can be reconciled");
        }
        if (arguments.gexProvisionalFeatureSorted && !arguments.inputFeatureSorted) {
            throw std::invalid_argument(
                "--gex-provisional-feature-sorted requires --input-feature-sorted");
        }
        if (arguments.gexProvisionalFeatureSorted && arguments.gexMultiGeneUmiCr) {
            throw std::invalid_argument(
                "provisional and in-memory GEX reconciliation modes are mutually exclusive");
        }
        if (arguments.inputFeatureSorted) {
            return runFeatureSorted(arguments);
        }
        const LoadedInput input = loadInput(arguments.input);
        prepareOutputDirectory(arguments.outDir);

        const std::vector<molecule_first::ReadClique> cliques =
            molecule_first::buildReadCliques(input.reads, input.priorCounts, arguments.config);
        const std::vector<molecule_first::HardCall> calls =
            molecule_first::gatedHardCalls(cliques, arguments.config);

        AtomicTable cliqueTable(arguments.outDir, "read_cliques.tsv",
            "clique_id\tfeature_id\traw_umi\tmember_read_count\tmember_read_ids\tcandidate\tlog_sequence_likelihood_sum\tlog_exact_read_prior\tlog_evidence\tposterior");
        cliqueTable.stream() << std::setprecision(17);
        for (const molecule_first::ReadClique &clique : cliques) {
            for (std::size_t index = 0; index < clique.candidates.size(); ++index) {
                cliqueTable.stream() << clique.cliqueId << '\t' << clique.featureId << '\t'
                    << clique.rawUmi << '\t' << clique.memberReadIds.size() << '\t'
                    << join(clique.memberReadIds, ';') << '\t' << clique.candidates[index] << '\t'
                    << clique.logLikelihoodSums[index] << '\t' << clique.logReadPriors[index] << '\t'
                    << clique.logEvidence[index] << '\t' << clique.posterior[index] << '\n';
            }
        }

        AtomicTable auditTable(arguments.outDir, "hard_call_audit.tsv",
            "clique_id\tstatus\tcandidate\tposterior\tmargin\treason");
        auditTable.stream() << std::setprecision(17);
        for (const molecule_first::HardCall &call : calls) {
            auditTable.stream() << call.cliqueId << '\t' << (call.assigned ? "assigned" : "deferred")
                << '\t' << call.candidate << '\t' << call.posterior << '\t' << call.margin
                << '\t' << call.reason << '\n';
        }

        const std::string moleculeHeader =
            "umi_mode\tproduct\tmolecule_id\tfeature_id\tcorrected_umi\tcandidate\tmember_read_count\tmember_read_ids\tread_clique_ids";
        AtomicTable strictTable(arguments.outDir, "strict_molecules.tsv", moleculeHeader);
        AtomicTable hardTable(arguments.outDir, "hard_molecules.tsv", moleculeHeader);
        AtomicTable gatedTable(arguments.outDir, "gated_hard_molecules.tsv", moleculeHeader);
        AtomicTable softTable(arguments.outDir, "soft_expected_molecules.tsv",
            "umi_mode\tfeature_id\tcorrected_umi\tcandidate\texpected_count\tread_clique_ids");
        softTable.stream() << std::setprecision(17);

        std::map<std::string, std::uint64_t> productCounts;
        std::map<std::string, double> occupancyMass;
        std::map<std::string, molecule_first::GexReconciliationStats> reconciliation;
        for (const std::string &umiMode : {std::string("1mm_cr"), std::string("exact")}) {
            molecule_first::GexReconciliationStats strictStats, hardStats, gatedStats, softStats;
            const std::vector<molecule_first::Molecule> strict =
                arguments.gexMultiGeneUmiCr
                ? molecule_first::gexPolicyMolecules(
                    cliques, umiMode, "strict", calls, &strictStats)
                : molecule_first::policyMolecules(cliques, umiMode, "strict", calls);
            const std::vector<molecule_first::Molecule> hard =
                arguments.gexMultiGeneUmiCr
                ? molecule_first::gexPolicyMolecules(
                    cliques, umiMode, "hard", calls, &hardStats)
                : molecule_first::policyMolecules(cliques, umiMode, "hard", calls);
            const std::vector<molecule_first::Molecule> gated =
                arguments.gexMultiGeneUmiCr
                ? molecule_first::gexPolicyMolecules(
                    cliques, umiMode, "gated_hard", calls, &gatedStats)
                : molecule_first::policyMolecules(cliques, umiMode, "gated_hard", calls);
            const std::vector<molecule_first::Occupancy> soft =
                arguments.gexMultiGeneUmiCr
                ? molecule_first::gexWeightedOccupancies(cliques, umiMode, &softStats)
                : molecule_first::weightedOccupancies(cliques, umiMode);
            if (arguments.gexMultiGeneUmiCr) {
                reconciliation[umiMode + ".strict"] = strictStats;
                reconciliation[umiMode + ".hard"] = hardStats;
                reconciliation[umiMode + ".gated_hard"] = gatedStats;
                reconciliation[umiMode + ".soft_expected"] = softStats;
            }
            writeMolecules(strictTable, strict);
            writeMolecules(hardTable, hard);
            writeMolecules(gatedTable, gated);
            productCounts[umiMode + ".strict"] = strict.size();
            productCounts[umiMode + ".hard"] = hard.size();
            productCounts[umiMode + ".gated_hard"] = gated.size();
            productCounts[umiMode + ".soft_rows"] = soft.size();
            for (const molecule_first::Occupancy &row : soft) {
                occupancyMass[umiMode] += row.expectedCount;
                softTable.stream() << row.umiMode << '\t' << row.featureId << '\t'
                    << row.correctedUmi << '\t' << row.candidate << '\t' << row.expectedCount
                    << '\t' << join(row.readCliqueIds, ';') << '\n';
            }
        }

        AtomicTable configTable(arguments.outDir, "resolved_config.tsv", "key\tvalue");
        configTable.stream() << std::setprecision(17)
            << "schema\tstar_suite.molecule_first.v1\n"
            << "star_suite_version\t" << STAR_SUITE_VERSION << '\n'
            << "gex_multigene_umi_cr\t" << (arguments.gexMultiGeneUmiCr ? 1 : 0) << '\n'
            << "soft_expected_semantics\t"
            << (arguments.gexMultiGeneUmiCr
                ? "candidate_weighted_expected_occupancy_then_multigene_reconciliation"
                : "candidate_weighted_expected_occupancy") << '\n'
            << "temperature\t" << arguments.config.temperature << '\n'
            << "prior_alpha\t" << arguments.config.priorAlpha << '\n'
            << "prior_beta\t" << arguments.config.priorBeta << '\n'
            << "prior_units\traw_exact_reads_including_pcr_amplification\n"
            << "prior_application\tonce_per_read_clique\n"
            << "spatial_lambda\t0\n"
            << "gate_min_posterior\t" << arguments.config.gateMinPosterior << '\n'
            << "gate_min_margin\t" << arguments.config.gateMinMargin << '\n';

        double posteriorMass = 0.0;
        for (const molecule_first::ReadClique &clique : cliques) {
            for (double value : clique.posterior) {
                posteriorMass += value;
            }
        }
        std::uint64_t gatedAssigned = 0;
        for (const molecule_first::HardCall &call : calls) {
            gatedAssigned += call.assigned ? 1 : 0;
        }
        AtomicTable summaryTable(arguments.outDir, "summary.tsv", "metric\tvalue");
        summaryTable.stream() << std::setprecision(17)
            << "candidate_rows\t" << input.candidateRows << '\n'
            << "candidate_reads\t" << input.reads.size() << '\n'
            << "read_cliques\t" << cliques.size() << '\n'
            << "posterior_mass\t" << posteriorMass << '\n'
            << "gated_assigned_cliques\t" << gatedAssigned << '\n'
            << "gated_deferred_cliques\t" << calls.size() - gatedAssigned << '\n';
        for (const auto &entry : productCounts) {
            summaryTable.stream() << entry.first << "_count\t" << entry.second << '\n';
        }
        for (const auto &entry : occupancyMass) {
            summaryTable.stream() << entry.first << ".soft_occupancy_mass\t" << entry.second << '\n';
            summaryTable.stream() << entry.first << ".deduplicated_mass\t"
                                  << posteriorMass - entry.second << '\n';
        }
        for (const auto &entry : reconciliation) {
            summaryTable.stream()
                << entry.first << ".multigene_groups\t" << entry.second.groups << '\n'
                << entry.first << ".multigene_accepted\t" << entry.second.accepted << '\n'
                << entry.first << ".multigene_corrected_count_ties\t"
                << entry.second.correctedCountTies << '\n'
                << entry.first << ".multigene_original_umi_dominance_rejected\t"
                << entry.second.originalUmiDominanceRejected << '\n';
            if (entry.first.find("soft_expected") != std::string::npos) {
                summaryTable.stream()
                    << entry.first << ".multigene_input_expected_mass\t"
                    << entry.second.inputExpectedMass << '\n'
                    << entry.first << ".multigene_output_expected_mass\t"
                    << entry.second.outputExpectedMass << '\n';
            }
        }

        cliqueTable.commit();
        auditTable.commit();
        strictTable.commit();
        hardTable.commit();
        gatedTable.commit();
        softTable.commit();
        configTable.commit();
        summaryTable.commit();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "molecule_first_resolver: ERROR: " << error.what() << '\n';
        return 2;
    }
}

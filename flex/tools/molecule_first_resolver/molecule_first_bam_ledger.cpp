#include <htslib/sam.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef STAR_SUITE_VERSION
#define STAR_SUITE_VERSION "unknown"
#endif

namespace {

struct Arguments {
    std::string input;
    std::string whitelist;
    std::string output;
    std::string summary;
    std::string assay;
    std::string sampleId;
    std::string featureTag = "GX";
};

struct Counters {
    std::uint64_t records = 0;
    std::uint64_t primary = 0;
    std::uint64_t eligible = 0;
    std::uint64_t exactEligible = 0;
    std::uint64_t missingRawTags = 0;
    std::uint64_t ambiguousFeature = 0;
    std::uint64_t invalidQuality = 0;
    std::uint64_t unsupported = 0;
    std::uint64_t candidateReads = 0;
    std::uint64_t candidateRows = 0;
};

void usage(std::ostream &out)
{
    out << "Usage: molecule_first_bam_ledger --input STAR.bam --whitelist barcodes.txt\\n"
        << "       --output candidate_reads.tsv --summary summary.tsv\\n"
        << "       --assay scrna|flex|visium [--sample-id ID] [--feature-tag GX]\\n\\n"
        << "The input may be SAM, BAM, or CRAM. Only raw CR/CY/UR and the declared\\n"
        << "feature tag are read; corrected CB/UB and cell calls are never used.\\n";
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
        else if (option == "--whitelist") arguments.whitelist = value;
        else if (option == "--output") arguments.output = value;
        else if (option == "--summary") arguments.summary = value;
        else if (option == "--assay") arguments.assay = value;
        else if (option == "--sample-id") arguments.sampleId = value;
        else if (option == "--feature-tag") arguments.featureTag = value;
        else throw std::invalid_argument("unknown option: " + option);
    }
    if (arguments.input.empty() || arguments.whitelist.empty() || arguments.output.empty()
        || arguments.summary.empty() || arguments.assay.empty()) {
        throw std::invalid_argument("--input, --whitelist, --output, --summary, and --assay are required");
    }
    if (arguments.assay != "scrna" && arguments.assay != "flex"
        && arguments.assay != "visium") {
        throw std::invalid_argument("--assay must be scrna, flex, or visium");
    }
    if (arguments.assay == "flex" && arguments.sampleId.empty()) {
        throw std::invalid_argument("--sample-id is required for Flex; run one demultiplexed sample per ledger");
    }
    if (arguments.featureTag.size() != 2) {
        throw std::invalid_argument("--feature-tag must be a two-character SAM tag");
    }
    if (arguments.output == arguments.summary) {
        throw std::invalid_argument("--output and --summary must differ");
    }
    return arguments;
}

bool isDna(const std::string &value)
{
    if (value.empty()) return false;
    for (char base : value) {
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') return false;
    }
    return true;
}

std::vector<std::string> loadWhitelist(const std::string &path)
{
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("cannot open whitelist: " + path);
    std::set<std::string> unique;
    std::string line;
    std::size_t length = 0;
    while (std::getline(input, line)) {
        const std::size_t fieldEnd = line.find_first_of("\t ,\r");
        std::string barcode = line.substr(0, fieldEnd);
        if (barcode.size() > 2 && barcode.compare(barcode.size() - 2, 2, "-1") == 0) {
            barcode.resize(barcode.size() - 2);
        }
        if (barcode.empty()) continue;
        if (!isDna(barcode)) throw std::runtime_error("whitelist contains a non-ACGT barcode: " + barcode);
        if (length == 0) length = barcode.size();
        if (barcode.size() != length) throw std::runtime_error("whitelist barcodes have inconsistent lengths");
        if (!unique.insert(barcode).second) throw std::runtime_error("duplicate whitelist barcode: " + barcode);
    }
    if (unique.empty()) throw std::runtime_error("whitelist is empty: " + path);
    return std::vector<std::string>(unique.begin(), unique.end());
}

std::string auxString(const bam1_t *record, const std::string &tag)
{
    const std::uint8_t *value = bam_aux_get(record, tag.c_str());
    if (value == nullptr) return std::string();
    const char *text = bam_aux2Z(value);
    return text == nullptr ? std::string() : std::string(text);
}

bool primaryRecord(const bam1_t *record)
{
    return (record->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0;
}

bool ambiguousFeature(const std::string &feature)
{
    return feature.empty() || feature == "-" || feature.find(';') != std::string::npos
        || feature.find(',') != std::string::npos;
}

bool starHeader(const bam_hdr_t *header)
{
    const char *text = header == nullptr ? nullptr : header->text;
    if (text == nullptr) return false;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        if (line.compare(0, 3, "@PG") != 0) continue;
        std::istringstream fields(line);
        std::string field;
        while (std::getline(fields, field, '\t')) {
            if (field == "PN:STAR" || field == "ID:STAR") return true;
        }
    }
    return false;
}

struct InputHandle {
    samFile *file = nullptr;
    bam_hdr_t *header = nullptr;
    bam1_t *record = nullptr;

    explicit InputHandle(const std::string &path)
    {
        file = sam_open(path.c_str(), "r");
        if (file == nullptr) throw std::runtime_error("cannot open alignment input: " + path);
        header = sam_hdr_read(file);
        if (header == nullptr) throw std::runtime_error("cannot read alignment header: " + path);
        record = bam_init1();
        if (record == nullptr) throw std::runtime_error("cannot allocate BAM record");
    }

    ~InputHandle()
    {
        if (record != nullptr) bam_destroy1(record);
        if (header != nullptr) bam_hdr_destroy(header);
        if (file != nullptr) sam_close(file);
    }
};

bool eligible(const bam1_t *record, const Arguments &arguments, Counters &counters,
              std::string &barcode, std::string &quality, std::string &umi,
              std::string &feature)
{
    ++counters.records;
    if (!primaryRecord(record)) return false;
    ++counters.primary;
    barcode = auxString(record, "CR");
    quality = auxString(record, "CY");
    umi = auxString(record, "UR");
    feature = auxString(record, arguments.featureTag);
    if (barcode.empty() || quality.empty() || umi.empty()) {
        ++counters.missingRawTags;
        return false;
    }
    if (ambiguousFeature(feature)) {
        ++counters.ambiguousFeature;
        return false;
    }
    if (quality.size() != barcode.size()) {
        ++counters.invalidQuality;
        return false;
    }
    ++counters.eligible;
    return true;
}

std::vector<std::string> candidatesFor(const std::string &observed,
                                       const std::unordered_set<std::string> &whitelist)
{
    const auto exact = whitelist.find(observed);
    if (exact != whitelist.end()) return std::vector<std::string>(1, observed);
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::vector<std::string> candidates;
    std::string candidate = observed;
    for (std::size_t position = 0; position < observed.size(); ++position) {
        for (char base : bases) {
            if (base == observed[position]) continue;
            candidate[position] = base;
            if (whitelist.find(candidate) != whitelist.end()) candidates.push_back(candidate);
        }
        candidate[position] = observed[position];
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    return candidates;
}

double logSequenceLikelihood(const std::string &observed, const std::string &quality,
                             const std::string &candidate)
{
    if (observed.size() != candidate.size() || observed.size() != quality.size()) {
        throw std::runtime_error("barcode/quality/candidate length mismatch");
    }
    double output = 0.0;
    for (std::size_t index = 0; index < observed.size(); ++index) {
        if (observed[index] != 'A' && observed[index] != 'C'
            && observed[index] != 'G' && observed[index] != 'T') {
            output += std::log(0.25);
            continue;
        }
        int phred = static_cast<unsigned char>(quality[index]) - 33;
        phred = std::max(0, std::min(60, phred));
        const double error = std::pow(10.0, -static_cast<double>(phred) / 10.0);
        if (observed[index] == candidate[index]) {
            output += std::log(std::max(1e-300, 1.0 - error));
        } else {
            output += std::log(std::max(1e-300, error / 3.0));
        }
    }
    return output;
}

void requireFresh(const std::string &path)
{
    std::ifstream existing(path.c_str());
    if (existing.good()) throw std::runtime_error("refusing to overwrite existing output: " + path);
}

void writeSummary(const Arguments &arguments, const Counters &priorPass,
                  const Counters &outputPass, std::uint64_t exactPriorMass)
{
    const std::string temporary = arguments.summary + ".tmp";
    std::ofstream output(temporary.c_str());
    if (!output) throw std::runtime_error("cannot create summary: " + temporary);
    output << "key\tvalue\n"
        << "schema\tstar_suite.molecule_first.bam_ledger.v1\n"
        << "star_suite_version\t" << STAR_SUITE_VERSION << '\n'
        << "assay\t" << arguments.assay << '\n'
        << "sample_id\t" << arguments.sampleId << '\n'
        << "feature_tag\t" << arguments.featureTag << '\n'
        << "star_header_verified\ttrue\n"
        << "barcode_source\tCR\nquality_source\tCY\numi_source\tUR\n"
        << "corrected_tags_used\tfalse\ncalled_cell_fields_used\tfalse\n"
        << "candidate_rule\texact_or_all_hamming1_substitutions\n"
        << "prior_units\traw_exact_reads_including_pcr_amplification\n"
        << "prior_eligible_records\t" << priorPass.eligible << '\n'
        << "prior_exact_reads\t" << priorPass.exactEligible << '\n'
        << "prior_exact_mass\t" << exactPriorMass << '\n'
        << "input_records\t" << outputPass.records << '\n'
        << "primary_records\t" << outputPass.primary << '\n'
        << "eligible_records\t" << outputPass.eligible << '\n'
        << "missing_raw_tags\t" << outputPass.missingRawTags << '\n'
        << "ambiguous_or_missing_feature\t" << outputPass.ambiguousFeature << '\n'
        << "invalid_quality_length\t" << outputPass.invalidQuality << '\n'
        << "unsupported_records\t" << outputPass.unsupported << '\n'
        << "candidate_reads\t" << outputPass.candidateReads << '\n'
        << "candidate_rows\t" << outputPass.candidateRows << '\n';
    output.close();
    if (!output || std::rename(temporary.c_str(), arguments.summary.c_str()) != 0) {
        throw std::runtime_error("cannot commit summary: " + arguments.summary);
    }
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Arguments arguments = parseArguments(argc, argv);
        requireFresh(arguments.output);
        requireFresh(arguments.summary);
        const std::vector<std::string> orderedWhitelist = loadWhitelist(arguments.whitelist);
        const std::unordered_set<std::string> whitelist(orderedWhitelist.begin(), orderedWhitelist.end());
        std::unordered_map<std::string, std::uint64_t> exactCounts;

        Counters priorCounters;
        {
            InputHandle input(arguments.input);
            if (!starHeader(input.header)) {
                throw std::runtime_error("alignment input is not verifiably STAR-produced (@PG PN:STAR/ID:STAR required)");
            }
            std::string barcode, quality, umi, feature;
            int readStatus = 0;
            while ((readStatus = sam_read1(input.file, input.header, input.record)) >= 0) {
                if (!eligible(input.record, arguments, priorCounters, barcode, quality, umi, feature)) continue;
                if (whitelist.find(barcode) != whitelist.end()) {
                    ++exactCounts[barcode];
                    ++priorCounters.exactEligible;
                }
            }
            if (readStatus < -1) throw std::runtime_error("alignment read error during prior pass");
        }

        const std::string temporary = arguments.output + ".tmp";
        std::ofstream output(temporary.c_str());
        if (!output) throw std::runtime_error("cannot create output: " + temporary);
        output << "read_id\tfeature_id\traw_umi\tcandidate\tlog_sequence_likelihood\texact_read_count\n";
        output << std::setprecision(17);
        Counters outputCounters;
        {
            InputHandle input(arguments.input);
            if (!starHeader(input.header)) throw std::runtime_error("STAR provenance changed between passes");
            std::string barcode, quality, umi, feature;
            int readStatus = 0;
            while ((readStatus = sam_read1(input.file, input.header, input.record)) >= 0) {
                if (!eligible(input.record, arguments, outputCounters, barcode, quality, umi, feature)) continue;
                if (barcode.size() != orderedWhitelist.front().size()) {
                    ++outputCounters.unsupported;
                    continue;
                }
                const std::vector<std::string> candidates = candidatesFor(barcode, whitelist);
                if (candidates.empty()) {
                    ++outputCounters.unsupported;
                    continue;
                }
                ++outputCounters.candidateReads;
                const std::string readId = bam_get_qname(input.record);
                for (const std::string &candidate : candidates) {
                    output << readId << '\t' << feature << '\t' << umi << '\t' << candidate << '\t'
                           << logSequenceLikelihood(barcode, quality, candidate) << '\t'
                           << exactCounts[candidate] << '\n';
                    ++outputCounters.candidateRows;
                }
            }
            if (readStatus < -1) throw std::runtime_error("alignment read error during candidate pass");
        }
        output.close();
        if (!output || std::rename(temporary.c_str(), arguments.output.c_str()) != 0) {
            throw std::runtime_error("cannot commit output: " + arguments.output);
        }
        std::uint64_t exactPriorMass = 0;
        for (const auto &entry : exactCounts) exactPriorMass += entry.second;
        writeSummary(arguments, priorCounters, outputCounters, exactPriorMass);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "molecule_first_bam_ledger: ERROR: " << error.what() << '\n';
        return 2;
    }
}

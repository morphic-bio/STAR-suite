// Stream one raw MatrixMarket MEX and report per-cell statistics for explicit,
// tag-aware cohorts.  The sparse matrix itself is never materialized.

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <zlib.h>

namespace fs = std::filesystem;

namespace {

[[noreturn]] void fail(const std::string& message)
{
    throw std::runtime_error(message);
}

class TextReader {
public:
    explicit TextReader(const fs::path& path) : path_(path), gzip_(path.extension() == ".gz")
    {
        if (gzip_) {
            inputGzip_ = gzopen(path.c_str(), "rb");
            if (inputGzip_ == nullptr) fail("cannot open " + path.string());
            gzbuffer(inputGzip_, 4U << 20);
        } else {
            inputPlain_.open(path);
            if (!inputPlain_) fail("cannot open " + path.string());
        }
    }

    TextReader(const TextReader&) = delete;
    TextReader& operator=(const TextReader&) = delete;

    ~TextReader()
    {
        if (inputGzip_ != nullptr) gzclose(inputGzip_);
    }

    bool getline(std::string& line)
    {
        if (!gzip_) {
            if (!std::getline(inputPlain_, line)) return false;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            return true;
        }

        line.clear();
        constexpr int chunkBytes = 1 << 16;
        char chunk[chunkBytes];
        while (true) {
            char* result = gzgets(inputGzip_, chunk, chunkBytes);
            if (result == nullptr) {
                if (!line.empty()) return true;
                int errorNumber = Z_OK;
                const char* errorMessage = gzerror(inputGzip_, &errorNumber);
                if (errorNumber != Z_OK && errorNumber != Z_STREAM_END) {
                    fail("gzip read failed for " + path_.string() + ": " + errorMessage);
                }
                return false;
            }
            const std::size_t length = std::char_traits<char>::length(chunk);
            if (length != 0 && chunk[length - 1] == '\n') {
                line.append(chunk, length - 1);
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            line.append(chunk, length);
        }
    }

private:
    fs::path path_;
    bool gzip_ = false;
    gzFile inputGzip_ = nullptr;
    std::ifstream inputPlain_;
};

fs::path plainOrGzip(const fs::path& directory, const std::string& name)
{
    const fs::path plain = directory / name;
    if (fs::is_regular_file(plain)) return plain;
    const fs::path compressed = directory / (name + ".gz");
    if (fs::is_regular_file(compressed)) return compressed;
    fail("missing " + plain.string() + "[.gz]");
}

fs::path featurePath(const fs::path& directory)
{
    const fs::path features = directory / "features.tsv";
    if (fs::is_regular_file(features)) return features;
    if (fs::is_regular_file(features.string() + ".gz")) return features.string() + ".gz";
    const fs::path genes = directory / "genes.tsv";
    if (fs::is_regular_file(genes)) return genes;
    if (fs::is_regular_file(genes.string() + ".gz")) return genes.string() + ".gz";
    fail("missing " + features.string() + "[.gz] (or genes.tsv[.gz])");
}

std::string firstField(std::string line)
{
    const std::size_t separator = line.find_first_of("\t ");
    if (separator != std::string::npos) line.resize(separator);
    return line;
}

std::string normalizeComposite(std::string barcode)
{
    barcode = firstField(std::move(barcode));
    if (barcode.size() >= 2 && barcode.compare(barcode.size() - 2, 2, "-1") == 0) {
        barcode.resize(barcode.size() - 2);
    }
    return barcode;
}

bool parseUnsigned(const char*& cursor, const char* end, std::uint64_t& value)
{
    while (cursor != end && (*cursor == ' ' || *cursor == '\t')) ++cursor;
    if (cursor == end) return false;
    const auto parsed = std::from_chars(cursor, end, value);
    if (parsed.ec != std::errc()) return false;
    cursor = parsed.ptr;
    return true;
}

void checkedAdd(std::uint64_t& total, std::uint64_t value, const std::string& what)
{
    if (value > std::numeric_limits<std::uint64_t>::max() - total) fail(what + " overflow");
    total += value;
}

struct CellStats {
    std::string barcode;
    bool found = false;
    std::uint64_t mexColumn = 0;
    std::uint64_t totalCounts = 0;
    std::uint64_t detectedGenes = 0;
    std::uint64_t nnz = 0;
};

struct Cohort {
    std::string label;
    fs::path path;
    std::vector<std::size_t> cells;
};

struct Options {
    fs::path mexDirectory;
    std::vector<std::pair<std::string, fs::path>> cohortSpecs;
    fs::path outputPrefix;
    std::string inputLabel;
};

void printHelp(std::ostream& output)
{
    output
        << "Usage: flex_mex_cell_stats --mex-dir DIR --cohort LABEL=FILE [...]\n"
        << "                           --out-prefix PATH [--input-label LABEL]\n\n"
        << "Streams matrix.mtx[.gz], barcodes.tsv[.gz], and features.tsv[.gz]\n"
        << "(or genes.tsv[.gz]) without materializing the sparse matrix. Cell IDs\n"
        << "are compared exactly after removing only a trailing -1; CB16+TAG8 IDs\n"
        << "are never truncated. Outputs PATH.cells.tsv and PATH.summary.tsv.\n";
}

Options parseOptions(int argc, char** argv)
{
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        auto value = [&]() -> std::string {
            if (++index >= argc) fail("missing value after " + argument);
            return argv[index];
        };
        if (argument == "--mex-dir") {
            options.mexDirectory = value();
        } else if (argument == "--cohort") {
            const std::string specification = value();
            const std::size_t equals = specification.find('=');
            if (equals == std::string::npos || equals == 0 || equals + 1 == specification.size()) {
                fail("--cohort must be LABEL=FILE: " + specification);
            }
            options.cohortSpecs.emplace_back(specification.substr(0, equals),
                                             specification.substr(equals + 1));
        } else if (argument == "--out-prefix") {
            options.outputPrefix = value();
        } else if (argument == "--input-label") {
            options.inputLabel = value();
        } else if (argument == "--help" || argument == "-h") {
            printHelp(std::cout);
            std::exit(0);
        } else {
            fail("unknown argument: " + argument);
        }
    }
    if (options.mexDirectory.empty()) fail("--mex-dir is required");
    if (!fs::is_directory(options.mexDirectory)) {
        fail("not a MEX directory: " + options.mexDirectory.string());
    }
    if (options.cohortSpecs.empty()) fail("at least one --cohort is required");
    if (options.outputPrefix.empty()) fail("--out-prefix is required");
    if (options.inputLabel.empty()) options.inputLabel = options.mexDirectory.filename().string();
    return options;
}

void readCohorts(const Options& options,
                 std::vector<Cohort>& cohorts,
                 std::vector<CellStats>& cells,
                 std::unordered_map<std::string, std::size_t>& cellByBarcode)
{
    std::unordered_set<std::string> labels;
    for (const auto& specification : options.cohortSpecs) {
        if (!labels.insert(specification.first).second) {
            fail("duplicate cohort label: " + specification.first);
        }
        Cohort cohort;
        cohort.label = specification.first;
        cohort.path = specification.second;
        TextReader input(cohort.path);
        std::unordered_set<std::string> seenInCohort;
        std::string line;
        std::uint64_t lineNumber = 0;
        while (input.getline(line)) {
            ++lineNumber;
            std::string barcode = normalizeComposite(std::move(line));
            if (barcode.empty() || barcode[0] == '#') continue;
            if (!seenInCohort.insert(barcode).second) {
                fail("duplicate barcode " + barcode + " in " + cohort.path.string() + ":"
                     + std::to_string(lineNumber));
            }
            auto inserted = cellByBarcode.emplace(barcode, cells.size());
            if (inserted.second) cells.push_back(CellStats{barcode});
            cohort.cells.push_back(inserted.first->second);
        }
        if (cohort.cells.empty()) fail("empty cohort: " + cohort.path.string());
        cohorts.push_back(std::move(cohort));
    }
}

std::uint64_t countAxisRows(const fs::path& path)
{
    TextReader input(path);
    std::string line;
    std::uint64_t count = 0;
    while (input.getline(line)) {
        if (!line.empty()) ++count;
    }
    return count;
}

std::vector<std::int64_t> readBarcodes(
    const fs::path& path,
    const std::unordered_map<std::string, std::size_t>& cellByBarcode,
    std::vector<CellStats>& cells)
{
    TextReader input(path);
    std::vector<std::int64_t> columnToCell;
    std::string line;
    while (input.getline(line)) {
        std::string barcode = normalizeComposite(std::move(line));
        if (barcode.empty()) fail("empty barcode in " + path.string());
        const auto selected = cellByBarcode.find(barcode);
        if (selected == cellByBarcode.end()) {
            columnToCell.push_back(-1);
            continue;
        }
        CellStats& stats = cells[selected->second];
        if (stats.found) fail("selected barcode occurs more than once in " + path.string() + ": " + barcode);
        stats.found = true;
        stats.mexColumn = columnToCell.size() + 1;
        columnToCell.push_back(static_cast<std::int64_t>(selected->second));
    }
    return columnToCell;
}

void streamMatrix(const fs::path& path,
                  std::uint64_t featureCount,
                  const std::vector<std::int64_t>& columnToCell,
                  std::vector<CellStats>& cells)
{
    TextReader input(path);
    std::string line;
    if (!input.getline(line) || line.rfind("%%MatrixMarket matrix coordinate", 0) != 0) {
        fail("unsupported MatrixMarket header in " + path.string());
    }
    do {
        if (!input.getline(line)) fail("missing MatrixMarket dimensions in " + path.string());
    } while (line.empty() || line[0] == '%');

    const char* cursor = line.data();
    const char* end = cursor + line.size();
    std::uint64_t matrixRows = 0;
    std::uint64_t matrixColumns = 0;
    std::uint64_t expectedNnz = 0;
    if (!parseUnsigned(cursor, end, matrixRows) || !parseUnsigned(cursor, end, matrixColumns)
        || !parseUnsigned(cursor, end, expectedNnz)) {
        fail("invalid MatrixMarket dimensions in " + path.string());
    }
    if (matrixRows != featureCount) {
        fail("MatrixMarket row count does not match features in " + path.string());
    }
    if (matrixColumns != columnToCell.size()) {
        fail("MatrixMarket column count does not match barcodes in " + path.string());
    }

    std::uint64_t observedNnz = 0;
    while (input.getline(line)) {
        if (line.empty() || line[0] == '%') continue;
        cursor = line.data();
        end = cursor + line.size();
        std::uint64_t row = 0;
        std::uint64_t column = 0;
        std::uint64_t value = 0;
        if (!parseUnsigned(cursor, end, row) || !parseUnsigned(cursor, end, column)
            || !parseUnsigned(cursor, end, value) || row == 0 || row > matrixRows || column == 0
            || column > matrixColumns) {
            fail("invalid MatrixMarket entry " + std::to_string(observedNnz + 1) + " in "
                 + path.string());
        }
        ++observedNnz;
        const std::int64_t selected = columnToCell[column - 1];
        if (selected < 0) continue;
        CellStats& stats = cells[static_cast<std::size_t>(selected)];
        checkedAdd(stats.totalCounts, value, "per-cell count");
        checkedAdd(stats.nnz, 1, "per-cell nnz");
        if (value != 0) checkedAdd(stats.detectedGenes, 1, "per-cell detected genes");
    }
    if (observedNnz != expectedNnz) {
        fail("MatrixMarket nnz header says " + std::to_string(expectedNnz) + " but read "
             + std::to_string(observedNnz) + " in " + path.string());
    }
}

std::uint64_t nearestRank(const std::vector<std::uint64_t>& sortedValues, unsigned percentile)
{
    if (sortedValues.empty()) return 0;
    if (percentile == 0) return sortedValues.front();
    const std::size_t rank = static_cast<std::size_t>(
        (static_cast<std::uint64_t>(percentile) * sortedValues.size() + 99) / 100);
    return sortedValues[std::min(rank, sortedValues.size()) - 1];
}

void ensureOutputParent(const fs::path& prefix)
{
    if (!prefix.parent_path().empty()) fs::create_directories(prefix.parent_path());
}

void writeCells(const Options& options,
                const std::vector<Cohort>& cohorts,
                const std::vector<CellStats>& cells)
{
    const fs::path path = options.outputPrefix.string() + ".cells.tsv";
    std::ofstream output(path);
    if (!output) fail("cannot write " + path.string());
    output << "input\tcohort\tbarcode\tstatus\tmex_column\ttotal_counts\tdetected_genes\tnnz\n";
    for (const Cohort& cohort : cohorts) {
        for (const std::size_t cellIndex : cohort.cells) {
            const CellStats& stats = cells[cellIndex];
            output << options.inputLabel << '\t' << cohort.label << '\t' << stats.barcode << '\t';
            if (!stats.found) {
                output << "missing\tNA\tNA\tNA\tNA\n";
            } else {
                output << "found\t" << stats.mexColumn << '\t' << stats.totalCounts << '\t'
                       << stats.detectedGenes << '\t' << stats.nnz << '\n';
            }
        }
    }
}

void writeSummary(const Options& options,
                  const std::vector<Cohort>& cohorts,
                  const std::vector<CellStats>& cells)
{
    const fs::path path = options.outputPrefix.string() + ".summary.tsv";
    std::ofstream output(path);
    if (!output) fail("cannot write " + path.string());
    output
        << "input\tcohort\texpected_cells\tfound_cells\tmissing_cells\tcells_with_counts"
        << "\ttotal_counts\ttotal_detected_genes\ttotal_nnz"
        << "\tcounts_min\tcounts_q25\tcounts_median\tcounts_q75\tcounts_p90"
        << "\tcounts_p95\tcounts_p99\tcounts_max"
        << "\tgenes_min\tgenes_q25\tgenes_median\tgenes_q75\tgenes_p90"
        << "\tgenes_p95\tgenes_p99\tgenes_max\n";

    for (const Cohort& cohort : cohorts) {
        std::uint64_t found = 0;
        std::uint64_t cellsWithCounts = 0;
        std::uint64_t totalCounts = 0;
        std::uint64_t totalGenes = 0;
        std::uint64_t totalNnz = 0;
        std::vector<std::uint64_t> counts;
        std::vector<std::uint64_t> genes;
        counts.reserve(cohort.cells.size());
        genes.reserve(cohort.cells.size());
        for (const std::size_t cellIndex : cohort.cells) {
            const CellStats& stats = cells[cellIndex];
            if (!stats.found) continue;
            ++found;
            if (stats.totalCounts != 0) ++cellsWithCounts;
            checkedAdd(totalCounts, stats.totalCounts, "cohort total counts");
            checkedAdd(totalGenes, stats.detectedGenes, "cohort total detected genes");
            checkedAdd(totalNnz, stats.nnz, "cohort total nnz");
            counts.push_back(stats.totalCounts);
            genes.push_back(stats.detectedGenes);
        }
        std::sort(counts.begin(), counts.end());
        std::sort(genes.begin(), genes.end());
        output << options.inputLabel << '\t' << cohort.label << '\t' << cohort.cells.size() << '\t'
               << found << '\t' << cohort.cells.size() - found << '\t' << cellsWithCounts << '\t'
               << totalCounts << '\t' << totalGenes << '\t' << totalNnz;
        for (const unsigned percentile : {0U, 25U, 50U, 75U, 90U, 95U, 99U, 100U}) {
            output << '\t' << nearestRank(counts, percentile);
        }
        for (const unsigned percentile : {0U, 25U, 50U, 75U, 90U, 95U, 99U, 100U}) {
            output << '\t' << nearestRank(genes, percentile);
        }
        output << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const Options options = parseOptions(argc, argv);
        std::vector<Cohort> cohorts;
        std::vector<CellStats> cells;
        std::unordered_map<std::string, std::size_t> cellByBarcode;
        readCohorts(options, cohorts, cells, cellByBarcode);

        const fs::path barcodes = plainOrGzip(options.mexDirectory, "barcodes.tsv");
        const fs::path features = featurePath(options.mexDirectory);
        const fs::path matrix = plainOrGzip(options.mexDirectory, "matrix.mtx");
        const std::uint64_t featureCount = countAxisRows(features);
        if (featureCount == 0) fail("empty feature axis: " + features.string());
        const std::vector<std::int64_t> columnToCell =
            readBarcodes(barcodes, cellByBarcode, cells);
        streamMatrix(matrix, featureCount, columnToCell, cells);

        ensureOutputParent(options.outputPrefix);
        writeCells(options, cohorts, cells);
        writeSummary(options, cohorts, cells);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "flex_mex_cell_stats: " << error.what() << '\n';
        return 1;
    }
}

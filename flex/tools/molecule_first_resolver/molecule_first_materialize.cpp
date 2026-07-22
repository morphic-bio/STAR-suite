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
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <unistd.h>

#ifndef STAR_SUITE_VERSION
#define STAR_SUITE_VERSION "unknown"
#endif

namespace {

const char *const kCliqueHeader =
    "clique_id\tfeature_id\traw_umi\tmember_read_count\tmember_read_ids\tcandidate\t"
    "log_sequence_likelihood_sum\tlog_exact_read_prior\tlog_evidence\tposterior";
const char *const kMoleculeHeader =
    "umi_mode\tproduct\tmolecule_id\tfeature_id\tcorrected_umi\tcandidate\t"
    "member_read_count\tmember_read_ids\tread_clique_ids";
const char *const kSoftHeader =
    "umi_mode\tfeature_id\tcorrected_umi\tcandidate\texpected_count\tread_clique_ids";

struct Arguments {
    std::string resolvedDir;
    std::string outDir;
    std::string assay;
    std::string umiMode = "1mm_cr";
    std::string tmpDir;
    std::uint64_t sortMemoryMb = 1024;
};

struct ProductSpec {
    std::string name;
    bool expected;
    std::string fileName;

    ProductSpec(const std::string &nameValue, bool expectedValue,
                const std::string &fileNameValue)
        : name(nameValue), expected(expectedValue), fileName(fileNameValue) {}
};

struct Scale {
    std::string name;
    int factor = 1;

    Scale(const std::string &nameValue, int factorValue)
        : name(nameValue), factor(factorValue) {}
};

struct MexSummary {
    std::size_t features = 0;
    std::size_t barcodes = 0;
    std::uint64_t nnz = 0;
    double mass = 0.0;
};

const std::vector<ProductSpec> &products()
{
    static const std::vector<ProductSpec> values = {
        {"strict", false, "strict_molecules.tsv"},
        {"soft_expected", true, "soft_expected_molecules.tsv"},
        {"hard", false, "hard_molecules.tsv"},
        {"gated_hard", false, "gated_hard_molecules.tsv"}
    };
    return values;
}

void usage(std::ostream &out)
{
    out << "Usage: molecule_first_materialize --resolved-dir DIR --out-dir DIR\n"
        << "       --assay scrna|flex|visium|visium-hd [--umi-mode 1mm_cr|exact]\n"
        << "       [--sort-memory-mb INT] [--tmp-dir DIR]\n\n"
        << "Writes strict, soft_expected, hard, and gated_hard 10x-style MEX\n"
        << "products on identical feature/barcode axes. Soft matrices are real.\n"
        << "Visium HD uses bounded external aggregation; --sort-memory-mb is the\n"
        << "combined in-memory record budget (default 1024 MiB).\n";
}

std::uint64_t parsePositiveInteger(const std::string &text, const std::string &option)
{
    if (text.empty() || text[0] == '-') {
        throw std::invalid_argument("invalid positive integer for " + option + ": " + text);
    }
    std::size_t used = 0;
    const unsigned long long value = std::stoull(text, &used);
    if (used != text.size() || value == 0) {
        throw std::invalid_argument("invalid positive integer for " + option + ": " + text);
    }
    return static_cast<std::uint64_t>(value);
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
        else if (option == "--tmp-dir") arguments.tmpDir = value;
        else if (option == "--sort-memory-mb") {
            arguments.sortMemoryMb = parsePositiveInteger(value, option);
        } else throw std::invalid_argument("unknown option: " + option);
    }
    if (arguments.resolvedDir.empty() || arguments.outDir.empty() || arguments.assay.empty()) {
        throw std::invalid_argument("--resolved-dir, --out-dir, and --assay are required");
    }
    const std::set<std::string> assays = {"scrna", "flex", "visium", "visium-hd"};
    if (assays.find(arguments.assay) == assays.end()) throw std::invalid_argument("unsupported --assay");
    if (arguments.umiMode != "1mm_cr" && arguments.umiMode != "exact") {
        throw std::invalid_argument("--umi-mode must be 1mm_cr or exact");
    }
    const std::uint64_t maxMb = std::numeric_limits<std::size_t>::max() / (1024ULL * 1024ULL);
    if (arguments.sortMemoryMb > maxMb) throw std::invalid_argument("--sort-memory-mb is too large");
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

void requireDirectory(const std::string &path)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !S_ISDIR(info.st_mode)) {
        throw std::runtime_error("temporary root is not a directory: " + path);
    }
}

void requireHeader(const std::vector<std::string> &observed,
                   const std::vector<std::string> &expected, const std::string &path)
{
    if (observed != expected) throw std::runtime_error("unexpected TSV schema: " + path);
}

class AtomicOutput {
  public:
    explicit AtomicOutput(const std::string &path)
        : finalPath_(path), tempPath_(path + ".tmp"), output_(tempPath_.c_str(), std::ios::binary)
    {
        if (!output_) throw std::runtime_error("cannot create output: " + tempPath_);
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

    ~AtomicOutput()
    {
        if (!committed_) {
            output_.close();
            std::remove(tempPath_.c_str());
        }
    }

  private:
    std::string finalPath_;
    std::string tempPath_;
    std::ofstream output_;
    bool committed_ = false;
};

void atomicText(const std::string &path, const std::string &text)
{
    AtomicOutput output(path);
    output.stream() << text;
    output.commit();
}

void atomicLinkOrCopy(const std::string &source, const std::string &destination)
{
    const std::string temporary = destination + ".tmp";
    std::remove(temporary.c_str());
    if (::link(source.c_str(), temporary.c_str()) == 0) {
        if (std::rename(temporary.c_str(), destination.c_str()) != 0) {
            std::remove(temporary.c_str());
            throw std::runtime_error("cannot commit linked axis: " + destination);
        }
        return;
    }
    std::ifstream input(source.c_str(), std::ios::binary);
    if (!input) throw std::runtime_error("cannot read shared axis: " + source);
    AtomicOutput output(destination);
    output.stream() << input.rdbuf();
    if (!input && !input.eof()) throw std::runtime_error("cannot copy shared axis: " + source);
    output.commit();
}

class PrivateTempDirectory {
  public:
    explicit PrivateTempDirectory(const std::string &root)
    {
        requireDirectory(root);
        std::string pattern = root + "/.molecule_first_materialize.XXXXXX";
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        char *created = mkdtemp(buffer.data());
        if (created == nullptr) {
            throw std::runtime_error("cannot create private temporary directory in " + root
                                     + ": " + std::strerror(errno));
        }
        path_ = created;
    }

    const std::string &path() const { return path_; }

    ~PrivateTempDirectory()
    {
        if (!path_.empty()) rmdir(path_.c_str());
    }

  private:
    std::string path_;
};

struct FieldSpan {
    std::size_t begin = 0;
    std::size_t end = 0;
};

bool scanFields(const std::string &line, FieldSpan *fields, std::size_t expected)
{
    std::size_t field = 0;
    std::size_t begin = 0;
    for (std::size_t position = 0; position <= line.size(); ++position) {
        if (position == line.size() || line[position] == '\t') {
            if (field >= expected) return false;
            fields[field].begin = begin;
            fields[field].end = position;
            ++field;
            begin = position + 1;
        }
    }
    return field == expected;
}

bool spanEquals(const std::string &line, const FieldSpan &span, const std::string &value)
{
    return span.end - span.begin == value.size()
        && line.compare(span.begin, value.size(), value) == 0;
}

std::string spanString(const std::string &line, const FieldSpan &span)
{
    return line.substr(span.begin, span.end - span.begin);
}

double spanFiniteDouble(const std::string &line, const FieldSpan &span, const std::string &path)
{
    const char *begin = line.c_str() + span.begin;
    char *end = nullptr;
    errno = 0;
    const double value = std::strtod(begin, &end);
    if (end != line.c_str() + span.end || errno == ERANGE
        || !std::isfinite(value) || value < 0.0) {
        throw std::runtime_error("invalid non-negative count in " + path + ": "
                                 + spanString(line, span));
    }
    return value;
}

std::uint32_t spanUnsigned(const std::string &line, std::size_t begin, std::size_t end,
                           const std::string &candidate)
{
    if (begin == end) throw std::runtime_error("invalid Visium HD coordinate candidate: " + candidate);
    if (line[begin] == '+') ++begin;
    if (begin == end) throw std::runtime_error("invalid Visium HD coordinate candidate: " + candidate);
    std::uint64_t value = 0;
    for (std::size_t position = begin; position < end; ++position) {
        const unsigned char character = static_cast<unsigned char>(line[position]);
        if (character < '0' || character > '9') {
            throw std::runtime_error("invalid Visium HD coordinate candidate: " + candidate);
        }
        value = value * 10 + static_cast<unsigned>(character - '0');
        if (value > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("Visium HD coordinate exceeds uint32 range: " + candidate);
        }
    }
    return static_cast<std::uint32_t>(value);
}

struct Coordinate {
    std::uint32_t row = 0;
    std::uint32_t column = 0;
};

std::uint64_t coordinateKey(const Coordinate &coordinate)
{
    return (static_cast<std::uint64_t>(coordinate.row) << 32) | coordinate.column;
}

Coordinate coordinateFromKey(std::uint64_t key)
{
    Coordinate coordinate;
    coordinate.row = static_cast<std::uint32_t>(key >> 32);
    coordinate.column = static_cast<std::uint32_t>(key);
    return coordinate;
}

Coordinate parseHdCandidate(const std::string &line, const FieldSpan &span)
{
    const std::string prefix = "s_002um_";
    const std::string candidate = spanString(line, span);
    if (span.end - span.begin < prefix.size()
        || line.compare(span.begin, prefix.size(), prefix) != 0) {
        throw std::runtime_error("Visium HD candidate does not use s_002um_row_col: " + candidate);
    }
    const std::size_t coordinatesBegin = span.begin + prefix.size();
    const std::size_t separator = line.find('_', coordinatesBegin);
    if (separator == std::string::npos || separator >= span.end
        || line.find('_', separator + 1) < span.end) {
        throw std::runtime_error("invalid Visium HD coordinate candidate: " + candidate);
    }
    Coordinate coordinate;
    coordinate.row = spanUnsigned(line, coordinatesBegin, separator, candidate);
    coordinate.column = spanUnsigned(line, separator + 1, span.end, candidate);
    return coordinate;
}

class CoordinateSet {
  public:
    void insert(const Coordinate &coordinate)
    {
        std::vector<std::uint64_t> &words = rows_[coordinate.row];
        const std::size_t word = coordinate.column / 64;
        if (words.size() <= word) words.resize(word + 1, 0);
        const std::uint64_t mask = 1ULL << (coordinate.column % 64);
        if ((words[word] & mask) == 0) {
            words[word] |= mask;
            ++size_;
            maxRow_ = std::max(maxRow_, coordinate.row);
            maxColumn_ = std::max(maxColumn_, coordinate.column);
        }
    }

    std::vector<Coordinate> coordinates() const
    {
        std::vector<Coordinate> values;
        values.reserve(size_);
        for (const auto &rowEntry : rows_) {
            for (std::size_t word = 0; word < rowEntry.second.size(); ++word) {
                std::uint64_t bits = rowEntry.second[word];
                while (bits != 0) {
#if defined(__GNUC__) || defined(__clang__)
                    const unsigned offset = static_cast<unsigned>(__builtin_ctzll(bits));
#else
                    unsigned offset = 0;
                    while (((bits >> offset) & 1ULL) == 0) ++offset;
#endif
                    Coordinate coordinate;
                    coordinate.row = rowEntry.first;
                    coordinate.column = static_cast<std::uint32_t>(word * 64 + offset);
                    values.push_back(coordinate);
                    bits &= bits - 1;
                }
            }
        }
        return values;
    }

    std::size_t size() const { return size_; }
    std::uint32_t maxRow() const { return maxRow_; }
    std::uint32_t maxColumn() const { return maxColumn_; }

  private:
    std::unordered_map<std::uint32_t, std::vector<std::uint64_t> > rows_;
    std::size_t size_ = 0;
    std::uint32_t maxRow_ = 0;
    std::uint32_t maxColumn_ = 0;
};

bool decimalComponentLess(std::uint32_t left, std::uint32_t right, bool followedByUnderscore)
{
    char leftBuffer[11], rightBuffer[11];
    char *leftBegin = leftBuffer + sizeof(leftBuffer);
    char *rightBegin = rightBuffer + sizeof(rightBuffer);
    do {
        *--leftBegin = static_cast<char>('0' + left % 10);
        left /= 10;
    } while (left != 0);
    do {
        *--rightBegin = static_cast<char>('0' + right % 10);
        right /= 10;
    } while (right != 0);
    const std::size_t leftSize = static_cast<std::size_t>(leftBuffer + sizeof(leftBuffer) - leftBegin);
    const std::size_t rightSize = static_cast<std::size_t>(rightBuffer + sizeof(rightBuffer) - rightBegin);
    const std::size_t shared = std::min(leftSize, rightSize);
    const int comparison = std::memcmp(leftBegin, rightBegin, shared);
    if (comparison != 0) return comparison < 0;
    if (leftSize == rightSize) return false;
    // A row component is followed by '_'. Digits sort before that delimiter,
    // so "1000_" precedes "1_" in the canonical candidate string. A final
    // column has no delimiter and therefore uses ordinary prefix ordering.
    return followedByUnderscore ? leftSize > rightSize : leftSize < rightSize;
}

std::vector<std::uint32_t> decimalRanks(std::uint32_t maximum, bool followedByUnderscore)
{
    std::vector<std::uint32_t> values(static_cast<std::size_t>(maximum) + 1);
    for (std::size_t index = 0; index < values.size(); ++index) {
        values[index] = static_cast<std::uint32_t>(index);
    }
    std::sort(values.begin(), values.end(), [&](std::uint32_t left, std::uint32_t right) {
        return decimalComponentLess(left, right, followedByUnderscore);
    });
    std::vector<std::uint32_t> ranks(values.size());
    for (std::size_t index = 0; index < values.size(); ++index) ranks[values[index]] = index;
    return ranks;
}

class CoordinateAxis {
  public:
    static CoordinateAxis build(std::vector<Coordinate> coordinates, int factor,
                                const std::string &name)
    {
        if (coordinates.empty()) throw std::runtime_error("Visium HD coordinate axis is empty");
        CoordinateAxis axis;
        axis.factor_ = factor;
        axis.name_ = name;
        std::uint32_t maxRow = 0, maxColumn = 0;
        for (const Coordinate &coordinate : coordinates) {
            maxRow = std::max(maxRow, coordinate.row);
            maxColumn = std::max(maxColumn, coordinate.column);
        }
        if (maxRow <= 1000000 && maxColumn <= 1000000) {
            const std::vector<std::uint32_t> rowRanks = decimalRanks(maxRow, true);
            const std::vector<std::uint32_t> columnRanks = decimalRanks(maxColumn, false);
            std::sort(coordinates.begin(), coordinates.end(),
                [&](const Coordinate &left, const Coordinate &right) {
                    if (rowRanks[left.row] != rowRanks[right.row]) {
                        return rowRanks[left.row] < rowRanks[right.row];
                    }
                    return columnRanks[left.column] < columnRanks[right.column];
                });
        } else {
            std::sort(coordinates.begin(), coordinates.end(),
                [](const Coordinate &left, const Coordinate &right) {
                    if (left.row != right.row) return decimalComponentLess(left.row, right.row, true);
                    return decimalComponentLess(left.column, right.column, false);
                });
        }
        axis.coordinates_ = std::move(coordinates);

        const std::uint64_t width = static_cast<std::uint64_t>(maxColumn) + 1;
        const std::uint64_t cells = (static_cast<std::uint64_t>(maxRow) + 1) * width;
        const std::uint64_t densityLimit = std::max<std::uint64_t>(1000000, axis.coordinates_.size() * 4ULL);
        if (cells <= densityLimit && cells <= 50000000ULL) {
            axis.denseWidth_ = static_cast<std::size_t>(width);
            axis.denseIndex_.assign(static_cast<std::size_t>(cells), 0);
            for (std::size_t index = 0; index < axis.coordinates_.size(); ++index) {
                const Coordinate &coordinate = axis.coordinates_[index];
                axis.denseIndex_[static_cast<std::size_t>(coordinate.row) * axis.denseWidth_
                                 + coordinate.column] = static_cast<std::uint32_t>(index + 1);
            }
        } else {
            axis.sparseIndex_.reserve(axis.coordinates_.size() * 2);
            for (std::size_t index = 0; index < axis.coordinates_.size(); ++index) {
                axis.sparseIndex_[coordinateKey(axis.coordinates_[index])] =
                    static_cast<std::uint32_t>(index + 1);
            }
        }
        return axis;
    }

    std::uint32_t index(const Coordinate &coordinate) const
    {
        if (!denseIndex_.empty()) {
            if (coordinate.column >= denseWidth_) return 0;
            const std::uint64_t position = static_cast<std::uint64_t>(coordinate.row) * denseWidth_
                + coordinate.column;
            if (position >= denseIndex_.size()) return 0;
            return denseIndex_[static_cast<std::size_t>(position)];
        }
        const auto found = sparseIndex_.find(coordinateKey(coordinate));
        return found == sparseIndex_.end() ? 0 : found->second;
    }

    const Coordinate &coordinate(std::uint32_t oneBasedIndex) const
    {
        if (oneBasedIndex == 0 || oneBasedIndex > coordinates_.size()) {
            throw std::runtime_error("internal Visium HD axis index is out of range");
        }
        return coordinates_[oneBasedIndex - 1];
    }

    std::size_t size() const { return coordinates_.size(); }
    int factor() const { return factor_; }
    const std::string &name() const { return name_; }
    const std::vector<Coordinate> &coordinates() const { return coordinates_; }

  private:
    int factor_ = 1;
    std::string name_;
    std::vector<Coordinate> coordinates_;
    std::size_t denseWidth_ = 0;
    std::vector<std::uint32_t> denseIndex_;
    std::unordered_map<std::uint64_t, std::uint32_t> sparseIndex_;
};

struct HdAxes {
    std::vector<std::string> features;
    std::unordered_map<std::string, std::uint32_t> featureIndex;
    CoordinateAxis fine;
    CoordinateAxis eight;
    CoordinateAxis sixteen;
};

HdAxes loadHdAxes(const std::string &path)
{
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("cannot open resolver artifact: " + path);
    std::string line;
    if (!std::getline(input, line) || line != kCliqueHeader) {
        throw std::runtime_error("unexpected TSV schema: " + path);
    }

    std::set<std::string> featureSet;
    CoordinateSet candidateSet;
    std::string previousFeature;
    FieldSpan fields[10];
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        if (!scanFields(line, fields, 10)) throw std::runtime_error("malformed resolver row: " + path);
        if (!spanEquals(line, fields[1], previousFeature)) {
            previousFeature = spanString(line, fields[1]);
            featureSet.insert(previousFeature);
        }
        candidateSet.insert(parseHdCandidate(line, fields[5]));
    }
    if (featureSet.empty() || candidateSet.size() == 0) {
        throw std::runtime_error("resolver universe is empty");
    }

    HdAxes axes;
    axes.features.assign(featureSet.begin(), featureSet.end());
    axes.featureIndex.reserve(axes.features.size() * 2);
    for (std::size_t index = 0; index < axes.features.size(); ++index) {
        axes.featureIndex[axes.features[index]] = static_cast<std::uint32_t>(index + 1);
    }

    std::vector<Coordinate> fineCoordinates = candidateSet.coordinates();
    CoordinateSet eightSet, sixteenSet;
    for (const Coordinate &coordinate : fineCoordinates) {
        Coordinate parent;
        parent.row = coordinate.row / 4;
        parent.column = coordinate.column / 4;
        eightSet.insert(parent);
        parent.row = coordinate.row / 8;
        parent.column = coordinate.column / 8;
        sixteenSet.insert(parent);
    }
    axes.fine = CoordinateAxis::build(std::move(fineCoordinates), 1, "square_002um");
    axes.eight = CoordinateAxis::build(eightSet.coordinates(), 4, "square_008um");
    axes.sixteen = CoordinateAxis::build(sixteenSet.coordinates(), 8, "square_016um");
    return axes;
}

std::string mexBarcode(const std::string &candidate)
{
    return candidate.size() >= 2
        && candidate.compare(candidate.size() - 2, 2, "-1") == 0
        ? candidate : candidate + "-1";
}

void writeFeatureAxis(const std::string &path, const std::vector<std::string> &features)
{
    AtomicOutput output(path);
    for (const std::string &feature : features) {
        output.stream() << feature << '\t' << feature << "\tGene Expression\n";
    }
    output.commit();
}

void writeHdBarcodeAxis(const std::string &path, const CoordinateAxis &axis)
{
    AtomicOutput output(path);
    const int microns = axis.factor() * 2;
    for (const Coordinate &coordinate : axis.coordinates()) {
        output.stream() << "s_" << std::setw(3) << std::setfill('0') << microns << "um_"
                        << std::setfill(' ') << coordinate.row << '_' << coordinate.column << "-1\n";
    }
    output.commit();
}

void prepareHdAxisFiles(const std::string &outDir, const HdAxes &axes)
{
    const CoordinateAxis *scaleAxes[] = {&axes.fine, &axes.eight, &axes.sixteen};
    for (const ProductSpec &product : products()) {
        for (const CoordinateAxis *axis : scaleAxes) {
            makeDirectories(outDir + "/" + product.name + "/" + axis->name());
        }
    }

    const std::string canonicalFeatures = outDir + "/strict/square_002um/features.tsv";
    writeFeatureAxis(canonicalFeatures, axes.features);
    for (const CoordinateAxis *axis : scaleAxes) {
        const std::string canonicalBarcode = outDir + "/strict/" + axis->name() + "/barcodes.tsv";
        writeHdBarcodeAxis(canonicalBarcode, *axis);
        for (const ProductSpec &product : products()) {
            const std::string directory = outDir + "/" + product.name + "/" + axis->name();
            const std::string featurePath = directory + "/features.tsv";
            const std::string barcodePath = directory + "/barcodes.tsv";
            if (featurePath != canonicalFeatures) atomicLinkOrCopy(canonicalFeatures, featurePath);
            if (barcodePath != canonicalBarcode) atomicLinkOrCopy(canonicalBarcode, barcodePath);
        }
    }
}

struct SortRecord {
    std::uint64_t key = 0;
    std::uint64_t order = 0;
    double value = 0.0;
};

static_assert(sizeof(SortRecord) == 24, "private spill record must remain compact");

std::uint64_t matrixKey(std::uint32_t column, std::uint32_t feature)
{
    return (static_cast<std::uint64_t>(column) << 32) | feature;
}

std::uint32_t keyColumn(std::uint64_t key) { return static_cast<std::uint32_t>(key >> 32); }
std::uint32_t keyFeature(std::uint64_t key) { return static_cast<std::uint32_t>(key); }

bool recordLess(const SortRecord &left, const SortRecord &right)
{
    if (left.key != right.key) return left.key < right.key;
    return left.order < right.order;
}

class ExternalSorter {
  public:
    ExternalSorter(const std::string &temporaryDirectory, const std::string &label,
                   std::size_t memoryBytes)
        : temporaryDirectory_(temporaryDirectory), label_(label)
    {
        capacity_ = std::max<std::size_t>(1, memoryBytes / sizeof(SortRecord));
        records_.reserve(capacity_);
    }

    void add(std::uint32_t column, std::uint32_t feature, double value, std::uint64_t order)
    {
        SortRecord record;
        record.key = matrixKey(column, feature);
        record.order = order;
        record.value = value;
        records_.push_back(record);
        if (records_.size() >= capacity_) spill();
    }

    void finish()
    {
        if (finished_) return;
        if (!records_.empty()) spill();
        std::vector<SortRecord>().swap(records_);
        finished_ = true;
    }

    template <typename Callback>
    void forEachAggregated(Callback callback) const
    {
        if (!finished_) throw std::runtime_error("internal external sorter was not finalized");
        struct Reader {
            explicit Reader(const std::string &path) : input(path.c_str(), std::ios::binary)
            {
                if (!input) throw std::runtime_error("cannot read private spill: " + path);
            }
            bool next(SortRecord &record)
            {
                input.read(reinterpret_cast<char *>(&record), sizeof(record));
                if (input) return true;
                if (input.eof() && input.gcount() == 0) return false;
                throw std::runtime_error("truncated private materializer spill");
            }
            std::ifstream input;
        };
        struct Node {
            SortRecord record;
            std::size_t reader = 0;
        };
        struct Later {
            bool operator()(const Node &left, const Node &right) const
            {
                if (recordLess(right.record, left.record)) return true;
                if (recordLess(left.record, right.record)) return false;
                return left.reader > right.reader;
            }
        };

        std::vector<std::unique_ptr<Reader> > readers;
        std::priority_queue<Node, std::vector<Node>, Later> queue;
        for (std::size_t index = 0; index < chunks_.size(); ++index) {
            readers.emplace_back(new Reader(chunks_[index]));
            Node node;
            node.reader = index;
            if (readers.back()->next(node.record)) queue.push(node);
        }

        bool haveGroup = false;
        std::uint64_t groupKey = 0;
        double groupValue = 0.0;
        while (!queue.empty()) {
            const Node node = queue.top();
            queue.pop();
            if (!haveGroup || node.record.key != groupKey) {
                if (haveGroup && groupValue != 0.0) callback(groupKey, groupValue);
                haveGroup = true;
                groupKey = node.record.key;
                groupValue = node.record.value;
            } else {
                groupValue += node.record.value;
            }
            Node next;
            next.reader = node.reader;
            if (readers[node.reader]->next(next.record)) queue.push(next);
        }
        if (haveGroup && groupValue != 0.0) callback(groupKey, groupValue);
    }

    void cleanup()
    {
        for (const std::string &path : chunks_) std::remove(path.c_str());
        chunks_.clear();
    }

    ~ExternalSorter() { cleanup(); }

  private:
    void spill()
    {
        std::sort(records_.begin(), records_.end(), recordLess);
        std::ostringstream name;
        name << temporaryDirectory_ << '/' << label_ << '.' << std::setw(6) << std::setfill('0')
             << chunks_.size() << ".bin";
        const std::string path = name.str();
        std::ofstream output(path.c_str(), std::ios::binary);
        if (!output) throw std::runtime_error("cannot create private spill: " + path);
        output.write(reinterpret_cast<const char *>(records_.data()),
                     static_cast<std::streamsize>(records_.size() * sizeof(SortRecord)));
        output.close();
        if (!output) throw std::runtime_error("cannot write private spill: " + path);
        chunks_.push_back(path);
        records_.clear();
    }

    std::string temporaryDirectory_;
    std::string label_;
    std::size_t capacity_ = 1;
    std::vector<SortRecord> records_;
    std::vector<std::string> chunks_;
    bool finished_ = false;
};

bool resolvedIsFeatureSorted(const std::string &resolvedDir)
{
    std::ifstream input((resolvedDir + "/resolved_config.tsv").c_str());
    if (!input) return false;
    std::string line;
    while (std::getline(input, line)) {
        if (line == "execution_mode\tfeature_sorted_streaming") return true;
    }
    return false;
}

void streamHdProduct(const std::string &path, const ProductSpec &product,
                     const Arguments &arguments, const HdAxes &axes,
                     bool featureSorted, ExternalSorter &sorter)
{
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("cannot open resolver artifact: " + path);
    std::string line;
    const std::string expectedHeader = product.expected ? kSoftHeader : kMoleculeHeader;
    if (!std::getline(input, line) || line != expectedHeader) {
        throw std::runtime_error("unexpected TSV schema: " + path);
    }

    const std::size_t fieldCount = product.expected ? 6 : 9;
    FieldSpan fields[9];
    std::string currentFeature;
    std::uint32_t currentFeatureIndex = 0;
    std::unordered_map<std::uint64_t, double> featureCounts;
    std::uint64_t lineNumber = 1;

    const auto flushFeature = [&]() {
        if (!featureSorted || currentFeatureIndex == 0) return;
        for (const auto &entry : featureCounts) {
            const Coordinate coordinate = coordinateFromKey(entry.first);
            const std::uint32_t column = axes.fine.index(coordinate);
            if (column == 0) {
                throw std::runtime_error("policy contains a candidate outside the read-clique universe");
            }
            sorter.add(column, currentFeatureIndex, entry.second, 0);
        }
        featureCounts.clear();
    };

    while (std::getline(input, line)) {
        ++lineNumber;
        if (line.empty()) continue;
        if (!scanFields(line, fields, fieldCount)) throw std::runtime_error("malformed resolver row: " + path);
        const std::size_t featureField = product.expected ? 1 : 3;
        const std::size_t candidateField = product.expected ? 3 : 5;
        const bool featureChanged = !spanEquals(line, fields[featureField], currentFeature);
        if (featureChanged) {
            const std::string nextFeature = spanString(line, fields[featureField]);
            if (featureSorted && !currentFeature.empty() && nextFeature < currentFeature) {
                throw std::runtime_error("feature-sorted resolver artifact is out of order: " + path);
            }
            flushFeature();
            currentFeature = nextFeature;
            currentFeatureIndex = 0;
        }

        const bool selected = spanEquals(line, fields[0], arguments.umiMode)
            && (product.expected || spanEquals(line, fields[1], product.name));
        if (!selected) continue;
        if (currentFeatureIndex == 0) {
            const auto feature = axes.featureIndex.find(currentFeature);
            if (feature == axes.featureIndex.end()) {
                throw std::runtime_error("policy contains a feature outside the read-clique universe");
            }
            currentFeatureIndex = feature->second;
        }
        const Coordinate coordinate = parseHdCandidate(line, fields[candidateField]);
        const std::uint32_t column = axes.fine.index(coordinate);
        if (column == 0) {
            throw std::runtime_error("policy contains a candidate outside the read-clique universe");
        }
        const double value = product.expected ? spanFiniteDouble(line, fields[4], path) : 1.0;
        if (featureSorted) {
            featureCounts[coordinateKey(coordinate)] += value;
        } else {
            sorter.add(column, currentFeatureIndex, value, lineNumber);
        }
    }
    flushFeature();
}

MexSummary summarizeSorted(const ExternalSorter &sorter, std::size_t features,
                           std::size_t barcodes)
{
    MexSummary summary;
    summary.features = features;
    summary.barcodes = barcodes;
    sorter.forEachAggregated([&](std::uint64_t, double value) {
        ++summary.nnz;
        summary.mass += value;
    });
    return summary;
}

void writeSortedMatrix(const std::string &path, const ExternalSorter &sorter,
                       const MexSummary &summary, bool expected)
{
    AtomicOutput output(path);
    output.stream() << "%%MatrixMarket matrix coordinate " << (expected ? "real" : "integer")
                    << " general\n"
                    << "% STAR Suite molecule-first post-collapse policy matrix\n"
                    << summary.features << ' ' << summary.barcodes << ' ' << summary.nnz << '\n'
                    << std::setprecision(17);
    std::uint64_t written = 0;
    sorter.forEachAggregated([&](std::uint64_t key, double value) {
        output.stream() << keyFeature(key) << ' ' << keyColumn(key) << ' ';
        if (expected) output.stream() << value;
        else output.stream() << static_cast<std::uint64_t>(std::llround(value));
        output.stream() << '\n';
        ++written;
    });
    if (written != summary.nnz) throw std::runtime_error("matrix NNZ changed between deterministic merge passes");
    output.commit();
}

std::vector<MexSummary> materializeHdProduct(const Arguments &arguments, const HdAxes &axes,
                                             const ProductSpec &product, bool featureSorted,
                                             const std::string &temporaryDirectory)
{
    const std::size_t totalBytes = static_cast<std::size_t>(arguments.sortMemoryMb * 1024ULL * 1024ULL);
    const std::size_t coarseBytes = std::max<std::size_t>(sizeof(SortRecord), totalBytes / 2);
    ExternalSorter fine(temporaryDirectory, product.name + ".002", totalBytes);
    streamHdProduct(arguments.resolvedDir + "/" + product.fileName, product, arguments,
                    axes, featureSorted, fine);
    fine.finish();

    ExternalSorter eight(temporaryDirectory, product.name + ".008", coarseBytes);
    ExternalSorter sixteen(temporaryDirectory, product.name + ".016", coarseBytes);
    MexSummary fineSummary;
    fineSummary.features = axes.features.size();
    fineSummary.barcodes = axes.fine.size();
    fine.forEachAggregated([&](std::uint64_t key, double value) {
        ++fineSummary.nnz;
        fineSummary.mass += value;
        const std::uint32_t fineColumn = keyColumn(key);
        const Coordinate &coordinate = axes.fine.coordinate(fineColumn);
        Coordinate parent;
        parent.row = coordinate.row / 4;
        parent.column = coordinate.column / 4;
        const std::uint32_t eightColumn = axes.eight.index(parent);
        parent.row = coordinate.row / 8;
        parent.column = coordinate.column / 8;
        const std::uint32_t sixteenColumn = axes.sixteen.index(parent);
        if (eightColumn == 0 || sixteenColumn == 0) {
            throw std::runtime_error("internal Visium HD barcode hierarchy is incomplete");
        }
        eight.add(eightColumn, keyFeature(key), value, fineColumn);
        sixteen.add(sixteenColumn, keyFeature(key), value, fineColumn);
    });
    eight.finish();
    sixteen.finish();

    writeSortedMatrix(arguments.outDir + "/" + product.name + "/square_002um/matrix.mtx",
                      fine, fineSummary, product.expected);
    fine.cleanup();

    MexSummary eightSummary = summarizeSorted(eight, axes.features.size(), axes.eight.size());
    writeSortedMatrix(arguments.outDir + "/" + product.name + "/square_008um/matrix.mtx",
                      eight, eightSummary, product.expected);
    eight.cleanup();

    MexSummary sixteenSummary = summarizeSorted(
        sixteen, axes.features.size(), axes.sixteen.size());
    writeSortedMatrix(arguments.outDir + "/" + product.name + "/square_016um/matrix.mtx",
                      sixteen, sixteenSummary, product.expected);
    sixteen.cleanup();

    return {fineSummary, eightSummary, sixteenSummary};
}

using Key = std::pair<std::string, std::string>;
using Counts = std::map<Key, double>;

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

void prepareGenericAxisFiles(const std::string &outDir, const std::vector<std::string> &features,
                             const std::vector<std::string> &candidates)
{
    for (const ProductSpec &product : products()) makeDirectories(outDir + "/" + product.name + "/raw");
    const std::string canonicalFeatures = outDir + "/strict/raw/features.tsv";
    const std::string canonicalBarcodes = outDir + "/strict/raw/barcodes.tsv";
    writeFeatureAxis(canonicalFeatures, features);
    {
        AtomicOutput output(canonicalBarcodes);
        for (const std::string &candidate : candidates) output.stream() << mexBarcode(candidate) << '\n';
        output.commit();
    }
    for (const ProductSpec &product : products()) {
        const std::string directory = outDir + "/" + product.name + "/raw";
        if (directory + "/features.tsv" != canonicalFeatures) {
            atomicLinkOrCopy(canonicalFeatures, directory + "/features.tsv");
        }
        if (directory + "/barcodes.tsv" != canonicalBarcodes) {
            atomicLinkOrCopy(canonicalBarcodes, directory + "/barcodes.tsv");
        }
    }
}

MexSummary writeGenericMatrix(const std::string &path, const std::set<std::string> &featureUniverse,
                              const std::set<std::string> &candidateUniverse, const Counts &source,
                              bool expected)
{
    const std::vector<std::string> features(featureUniverse.begin(), featureUniverse.end());
    const std::vector<std::string> candidates(candidateUniverse.begin(), candidateUniverse.end());
    std::map<std::string, std::size_t> featureIndex, candidateIndex;
    for (std::size_t index = 0; index < features.size(); ++index) featureIndex[features[index]] = index + 1;
    for (std::size_t index = 0; index < candidates.size(); ++index) candidateIndex[candidates[index]] = index + 1;

    std::vector<std::tuple<std::size_t, std::size_t, double> > rows;
    MexSummary summary;
    summary.features = features.size();
    summary.barcodes = candidates.size();
    for (const auto &entry : source) {
        if (featureUniverse.find(entry.first.first) == featureUniverse.end()
            || candidateUniverse.find(entry.first.second) == candidateUniverse.end()) {
            throw std::runtime_error("policy contains a feature/candidate outside the read-clique universe");
        }
        if (entry.second == 0.0) continue;
        rows.push_back(std::make_tuple(candidateIndex[entry.first.second],
                                       featureIndex[entry.first.first], entry.second));
        summary.mass += entry.second;
    }
    std::sort(rows.begin(), rows.end());
    summary.nnz = rows.size();

    AtomicOutput output(path);
    output.stream() << "%%MatrixMarket matrix coordinate " << (expected ? "real" : "integer")
                    << " general\n"
                    << "% STAR Suite molecule-first post-collapse policy matrix\n"
                    << features.size() << ' ' << candidates.size() << ' ' << rows.size() << '\n'
                    << std::setprecision(17);
    for (const auto &row : rows) {
        output.stream() << std::get<1>(row) << ' ' << std::get<0>(row) << ' ';
        if (expected) output.stream() << std::get<2>(row);
        else output.stream() << static_cast<std::uint64_t>(std::llround(std::get<2>(row)));
        output.stream() << '\n';
    }
    output.commit();
    return summary;
}

void appendSummaryRow(std::ostringstream &summary, const ProductSpec &product,
                      const std::string &scale, const MexSummary &result)
{
    summary << product.name << '\t' << scale << '\t' << result.features << '\t'
            << result.barcodes << '\t' << result.nnz << '\t' << std::setprecision(17)
            << result.mass << '\t' << (product.expected ? "real" : "integer") << '\n';
}

std::ostringstream summaryHeader(const Arguments &arguments)
{
    std::ostringstream summary;
    summary << "schema\tstar_suite.molecule_first.policy_mex.v1\n"
            << "star_suite_version\t" << STAR_SUITE_VERSION << '\n'
            << "assay\t" << arguments.assay << '\n'
            << "umi_mode\t" << arguments.umiMode << '\n'
            << "axes_source\tread_cliques.tsv\n"
            << "cell_calling_order\tafter_postcollapse_integer_matrix_only\n"
            << "soft_cell_calling_allowed\tfalse\n"
            << "product\tscale\tfeatures\tbarcodes\tnnz\tmass\tmatrix_field\n";
    return summary;
}

void runHd(const Arguments &arguments)
{
    HdAxes axes = loadHdAxes(arguments.resolvedDir + "/read_cliques.tsv");
    prepareHdAxisFiles(arguments.outDir, axes);
    const std::string temporaryRoot = arguments.tmpDir.empty() ? arguments.outDir : arguments.tmpDir;
    PrivateTempDirectory temporary(temporaryRoot);
    const bool featureSorted = resolvedIsFeatureSorted(arguments.resolvedDir);
    std::ostringstream summary = summaryHeader(arguments);
    for (const ProductSpec &product : products()) {
        const std::vector<MexSummary> results = materializeHdProduct(
            arguments, axes, product, featureSorted, temporary.path());
        appendSummaryRow(summary, product, "square_002um", results[0]);
        appendSummaryRow(summary, product, "square_008um", results[1]);
        appendSummaryRow(summary, product, "square_016um", results[2]);
        const double minimum = std::min(results[0].mass, std::min(results[1].mass, results[2].mass));
        const double maximum = std::max(results[0].mass, std::max(results[1].mass, results[2].mass));
        const double tolerance = 1e-9 * std::max(1.0, maximum);
        if (maximum - minimum > tolerance) {
            throw std::runtime_error("Visium HD scale aggregation did not conserve mass for " + product.name);
        }
    }
    atomicText(arguments.outDir + "/summary.tsv", summary.str());
}

void runGeneric(const Arguments &arguments)
{
    std::set<std::string> featureUniverse, candidateUniverse;
    loadUniverse(arguments.resolvedDir + "/read_cliques.tsv", featureUniverse, candidateUniverse);
    const std::vector<std::string> features(featureUniverse.begin(), featureUniverse.end());
    const std::vector<std::string> candidates(candidateUniverse.begin(), candidateUniverse.end());
    prepareGenericAxisFiles(arguments.outDir, features, candidates);
    std::ostringstream summary = summaryHeader(arguments);
    for (const ProductSpec &product : products()) {
        Counts counts;
        if (product.expected) {
            loadSoft(arguments.resolvedDir + "/" + product.fileName, arguments.umiMode, counts);
        } else {
            loadMolecules(arguments.resolvedDir + "/" + product.fileName,
                          arguments.umiMode, product.name, counts);
        }
        const MexSummary result = writeGenericMatrix(
            arguments.outDir + "/" + product.name + "/raw/matrix.mtx",
            featureUniverse, candidateUniverse, counts, product.expected);
        appendSummaryRow(summary, product, "raw", result);
    }
    atomicText(arguments.outDir + "/summary.tsv", summary.str());
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Arguments arguments = parseArguments(argc, argv);
        prepareOutputDirectory(arguments.outDir);
        if (arguments.assay == "visium-hd") runHd(arguments);
        else runGeneric(arguments);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "molecule_first_materialize: ERROR: " << error.what() << '\n';
        return 2;
    }
}

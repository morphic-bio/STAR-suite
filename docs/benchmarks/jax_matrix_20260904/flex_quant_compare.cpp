// One-pass STAR-Flex quantification comparison against Cell Ranger MEX output.
//
// This is a native replacement for the repeated Python MEX passes in:
//   concordance_vs_cr.py
//   paper_protocol_concordance.py
//   umi_totals_vs_cr_320k.py
//   cell_calling_pr.py (STAR operating point and AUPRC table only)
//
// Each query and Cell Ranger matrix is parsed once per sample group, then the
// reports are calculated from the shared sparse representation. Cell identity
// is CB16|tag: the same CB16 under two tags represents two independent cells.
// If --star-run and --sample-whitelist are supplied, the raw MEX is streamed
// once to add the STAR cell-calling report.
//
// Build:
//   g++ -O3 -DNDEBUG -std=c++17 -fopenmp flex_quant_compare.cpp -lz
//       -o flex_quant_compare
// For cyto H5AD input, additionally use -DFLEX_QUANT_COMPARE_HDF5 and the
// compiler/linker flags printed by `pkg-config --cflags --libs hdf5`.
//
// Run:
//   ./flex_quant_compare --query-root STAR_RUN/per_sample --label star-bgzf
//       --cr-root CR_RUN/outs/per_sample_outs --cr-config CR_RUN_CONFIG.csv
//       --tag-map sample_whitelist.tsv
//       --out-prefix analysis/star-bgzf --threads 48

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <omp.h>
#include <zlib.h>
#ifdef FLEX_QUANT_COMPARE_HDF5
#include <hdf5.h>
#endif

namespace fs = std::filesystem;

namespace {

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

class TextReader {
public:
    explicit TextReader(const fs::path &path) : path_(path), gzip_(path.extension() == ".gz") {
        if (gzip_) {
            gz_ = gzopen(path.c_str(), "rb");
            if (gz_ == nullptr) fail("cannot open " + path.string());
            gzbuffer(gz_, 4U << 20);
        } else {
            in_.open(path);
            if (!in_) fail("cannot open " + path.string());
        }
    }

    TextReader(const TextReader &) = delete;
    TextReader &operator=(const TextReader &) = delete;

    ~TextReader() {
        if (gz_ != nullptr) gzclose(gz_);
    }

    bool getline(std::string &line) {
        if (!gzip_) return static_cast<bool>(std::getline(in_, line));
        line.clear();
        constexpr int chunk_size = 1 << 20;
        thread_local std::vector<char> chunk(chunk_size);
        while (true) {
            char *result = gzgets(gz_, chunk.data(), chunk_size);
            if (result == nullptr) {
                if (!line.empty()) return true;
                int error_number = Z_OK;
                const char *message = gzerror(gz_, &error_number);
                if (error_number != Z_OK && error_number != Z_STREAM_END) {
                    fail("gzip read failed for " + path_.string() + ": " + message);
                }
                return false;
            }
            const std::size_t length = std::char_traits<char>::length(chunk.data());
            if (length != 0 && chunk[length - 1] == '\n') {
                line.append(chunk.data(), length - 1);
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            line.append(chunk.data(), length);
        }
    }

private:
    fs::path path_;
    bool gzip_ = false;
    gzFile gz_ = nullptr;
    std::ifstream in_;
};

fs::path plain_or_gz(const fs::path &directory, const std::string &stem) {
    const fs::path plain = directory / stem;
    if (fs::exists(plain)) return plain;
    const fs::path compressed = directory / (stem + ".gz");
    if (fs::exists(compressed)) return compressed;
    fail("missing " + plain.string() + "[.gz]");
}

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string &value, char delimiter) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = value.find(delimiter, begin);
        fields.push_back(value.substr(begin, end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

std::vector<std::string> csv_fields(const std::string &line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
        if (ch == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                field.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (ch == ',' && !quoted) {
            fields.push_back(field);
            field.clear();
        } else {
            field.push_back(ch);
        }
    }
    if (quoted) fail("unterminated quoted CSV field: " + line);
    fields.push_back(field);
    return fields;
}

using TagMap = std::unordered_map<std::string, std::string>;

std::string strip_barcode_suffix(std::string barcode) {
    const std::size_t tab = barcode.find('\t');
    if (tab != std::string::npos) barcode.resize(tab);
    const std::size_t dash = barcode.find('-');
    if (dash != std::string::npos) barcode.resize(dash);
    return barcode;
}

std::string normalize_barcode(std::string barcode) {
    barcode = strip_barcode_suffix(std::move(barcode));
    if (barcode.size() > 16) barcode.resize(16);
    return barcode;
}

std::string tagged_cell_key(const std::string &raw_barcode,
                            const std::string &tag_id) {
    const std::string cb16 = normalize_barcode(raw_barcode);
    if (cb16.size() != 16) fail("cell barcode is not 16 bases: " + raw_barcode);
    if (tag_id.empty()) fail("empty tag id for cell barcode " + raw_barcode);
    return cb16 + "|" + tag_id;
}

std::string cr_cell_key(const std::string &raw_barcode,
                        const std::vector<std::string> &group_tag_ids,
                        const TagMap &tag8_to_id) {
    const std::string composite = strip_barcode_suffix(raw_barcode);
    if (composite.size() < 16) {
        fail("Cell Ranger cell barcode is shorter than 16 bases: " + raw_barcode);
    }
    if (composite.size() >= 24) {
        const std::string tag8 = composite.substr(16, 8);
        const auto found = tag8_to_id.find(tag8);
        if (found != tag8_to_id.end()) {
            if (std::find(group_tag_ids.begin(), group_tag_ids.end(), found->second) ==
                group_tag_ids.end()) {
                fail("Cell Ranger tag " + tag8 + " maps to " + found->second +
                     ", which is outside this sample group");
            }
            return tagged_cell_key(composite, found->second);
        }
        if (!tag8_to_id.empty() || group_tag_ids.size() > 1) {
            fail("Cell Ranger tag " + tag8 + " is absent from --tag-map");
        }
    }
    if (group_tag_ids.size() == 1) {
        return tagged_cell_key(composite, group_tag_ids.front());
    }
    fail("grouped Cell Ranger comparison requires CB16+TAG8 barcodes and --tag-map");
}

std::vector<std::string> read_barcodes(const fs::path &directory, bool normalize = true) {
    TextReader input(plain_or_gz(directory, "barcodes.tsv"));
    std::vector<std::string> result;
    std::string line;
    while (input.getline(line)) {
        line = trim(line);
        if (!line.empty()) result.push_back(normalize ? normalize_barcode(line) : line);
    }
    return result;
}

std::vector<std::string> read_genes(const fs::path &directory) {
    TextReader input(plain_or_gz(directory, "features.tsv"));
    std::vector<std::string> result;
    std::string line;
    while (input.getline(line)) {
        if (!line.empty()) result.push_back(line.substr(0, line.find('\t')));
    }
    return result;
}

bool parse_unsigned(const char *&cursor, const char *end, std::uint64_t &value) {
    while (cursor != end && (*cursor == ' ' || *cursor == '\t')) ++cursor;
    if (cursor == end) return false;
    const auto parsed = std::from_chars(cursor, end, value);
    if (parsed.ec != std::errc()) return false;
    cursor = parsed.ptr;
    return true;
}

struct SparseMatrix {
    std::vector<std::string> barcodes;
    std::vector<std::string> genes;
    std::vector<std::uint64_t> offsets;
    std::vector<std::uint32_t> columns;
    std::vector<std::uint64_t> values;
    long double total = 0;

    std::size_t rows() const { return barcodes.size(); }
};

void normalize_sparse_rows(SparseMatrix &matrix) {
    const std::size_t cell_count = matrix.rows();
    bool needs_normalization = false;
    for (std::size_t row = 0; row < cell_count && !needs_normalization; ++row) {
        for (std::uint64_t i = matrix.offsets[row] + 1; i < matrix.offsets[row + 1]; ++i) {
            if (matrix.columns[i - 1] >= matrix.columns[i]) {
                needs_normalization = true;
                break;
            }
        }
    }
    if (!needs_normalization) return;

    std::vector<std::uint64_t> offsets(cell_count + 1, 0);
    std::vector<std::uint32_t> columns;
    std::vector<std::uint64_t> values;
    columns.reserve(matrix.columns.size());
    values.reserve(matrix.values.size());
    std::vector<std::pair<std::uint32_t, std::uint64_t>> row_values;
    for (std::size_t row = 0; row < cell_count; ++row) {
        row_values.clear();
        for (std::uint64_t i = matrix.offsets[row]; i < matrix.offsets[row + 1]; ++i) {
            row_values.emplace_back(matrix.columns[i], matrix.values[i]);
        }
        std::sort(row_values.begin(), row_values.end());
        for (const auto &[column, value] : row_values) {
            if (!columns.empty() && offsets[row] < columns.size() && columns.back() == column) {
                values.back() += value;
            } else {
                columns.push_back(column);
                values.push_back(value);
            }
        }
        offsets[row + 1] = columns.size();
    }
    matrix.offsets.swap(offsets);
    matrix.columns.swap(columns);
    matrix.values.swap(values);
}

SparseMatrix read_mex(const fs::path &directory, bool normalize_barcodes = true) {
    SparseMatrix matrix;
    matrix.barcodes = read_barcodes(directory, normalize_barcodes);
    matrix.genes = read_genes(directory);

    TextReader input(plain_or_gz(directory, "matrix.mtx"));
    std::string line;
    if (!input.getline(line) || line.rfind("%%MatrixMarket", 0) != 0) {
        fail("unexpected MatrixMarket header in " + directory.string());
    }
    do {
        if (!input.getline(line)) fail("missing MatrixMarket dimensions in " + directory.string());
    } while (line.empty() || line[0] == '%');

    const char *cursor = line.data();
    const char *end = cursor + line.size();
    std::uint64_t feature_count = 0, cell_count = 0, expected_nnz = 0;
    if (!parse_unsigned(cursor, end, feature_count) ||
        !parse_unsigned(cursor, end, cell_count) ||
        !parse_unsigned(cursor, end, expected_nnz)) {
        fail("invalid MatrixMarket dimensions in " + directory.string());
    }
    if (feature_count != matrix.genes.size() || cell_count != matrix.barcodes.size()) {
        fail("MatrixMarket axes disagree with features/barcodes in " + directory.string());
    }

    std::vector<std::uint32_t> input_cells;
    std::vector<std::uint32_t> input_genes;
    std::vector<std::uint64_t> input_values;
    input_cells.reserve(expected_nnz);
    input_genes.reserve(expected_nnz);
    input_values.reserve(expected_nnz);
    std::vector<std::uint64_t> row_counts(cell_count, 0);

    while (input.getline(line)) {
        if (line.empty() || line[0] == '%') continue;
        cursor = line.data();
        end = cursor + line.size();
        std::uint64_t gene = 0, cell = 0, value = 0;
        if (!parse_unsigned(cursor, end, gene) || !parse_unsigned(cursor, end, cell) ||
            !parse_unsigned(cursor, end, value) || gene == 0 || cell == 0 ||
            gene > feature_count || cell > cell_count) {
            fail("invalid MatrixMarket entry in " + directory.string() + ": " + line);
        }
        input_genes.push_back(static_cast<std::uint32_t>(gene - 1));
        input_cells.push_back(static_cast<std::uint32_t>(cell - 1));
        input_values.push_back(value);
        ++row_counts[cell - 1];
        matrix.total += static_cast<long double>(value);
    }
    if (input_cells.size() != expected_nnz) {
        fail("MatrixMarket nnz mismatch in " + directory.string());
    }

    matrix.offsets.resize(cell_count + 1, 0);
    for (std::size_t row = 0; row < cell_count; ++row) {
        matrix.offsets[row + 1] = matrix.offsets[row] + row_counts[row];
    }
    matrix.columns.resize(expected_nnz);
    matrix.values.resize(expected_nnz);
    std::vector<std::uint64_t> cursor_by_row = matrix.offsets;
    for (std::size_t i = 0; i < input_cells.size(); ++i) {
        const std::uint64_t destination = cursor_by_row[input_cells[i]]++;
        matrix.columns[destination] = input_genes[i];
        matrix.values[destination] = input_values[i];
    }

    input_cells.clear();
    input_genes.clear();
    input_values.clear();

    normalize_sparse_rows(matrix);
    return matrix;
}

#ifdef FLEX_QUANT_COMPARE_HDF5
template <herr_t (*Closer)(hid_t)>
class H5Object {
public:
    explicit H5Object(hid_t value = -1) : value_(value) {}
    H5Object(const H5Object &) = delete;
    H5Object &operator=(const H5Object &) = delete;
    H5Object(H5Object &&other) noexcept : value_(other.value_) { other.value_ = -1; }
    ~H5Object() { if (value_ >= 0) Closer(value_); }
    operator hid_t() const { return value_; }
    bool valid() const { return value_ >= 0; }
private:
    hid_t value_;
};

using H5File = H5Object<H5Fclose>;
using H5Dataset = H5Object<H5Dclose>;
using H5Space = H5Object<H5Sclose>;
using H5Type = H5Object<H5Tclose>;

std::size_t h5_vector_size(hid_t dataset, const std::string &path) {
    H5Space space(H5Dget_space(dataset));
    if (!space.valid() || H5Sget_simple_extent_ndims(space) != 1) {
        fail("expected one-dimensional HDF5 dataset " + path);
    }
    hsize_t size = 0;
    if (H5Sget_simple_extent_dims(space, &size, nullptr) < 0) {
        fail("cannot read HDF5 dimensions for " + path);
    }
    return static_cast<std::size_t>(size);
}

template <typename T>
std::vector<T> read_h5_vector(hid_t file, const std::string &path, hid_t memory_type) {
    H5Dataset dataset(H5Dopen2(file, path.c_str(), H5P_DEFAULT));
    if (!dataset.valid()) fail("cannot open HDF5 dataset " + path);
    std::vector<T> result(h5_vector_size(dataset, path));
    if (!result.empty() && H5Dread(dataset, memory_type, H5S_ALL, H5S_ALL,
                                   H5P_DEFAULT, result.data()) < 0) {
        fail("cannot read HDF5 dataset " + path);
    }
    return result;
}

std::vector<std::string> read_h5_strings(hid_t file, const std::string &path) {
    H5Dataset dataset(H5Dopen2(file, path.c_str(), H5P_DEFAULT));
    if (!dataset.valid()) fail("cannot open HDF5 string dataset " + path);
    H5Space space(H5Dget_space(dataset));
    const std::size_t size = h5_vector_size(dataset, path);
    H5Type memory_type(H5Dget_type(dataset));
    if (!space.valid() || !memory_type.valid() || H5Tget_class(memory_type) != H5T_STRING ||
        H5Tis_variable_str(memory_type) <= 0) {
        fail("cannot prepare HDF5 string reader for " + path);
    }
    std::vector<char *> pointers(size, nullptr);
    if (!pointers.empty() && H5Dread(dataset, memory_type, H5S_ALL, H5S_ALL,
                                     H5P_DEFAULT, pointers.data()) < 0) {
        fail("cannot read HDF5 strings from " + path);
    }
    std::vector<std::string> result;
    result.reserve(size);
    for (const char *value : pointers) result.emplace_back(value == nullptr ? "" : value);
    if (!pointers.empty() && H5Dvlen_reclaim(memory_type, space, H5P_DEFAULT, pointers.data()) < 0) {
        fail("cannot reclaim HDF5 strings from " + path);
    }
    return result;
}

SparseMatrix read_cyto_h5ad(const fs::path &path) {
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    H5File file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    if (!file.valid()) fail("cannot open " + path.string());
    SparseMatrix matrix;
    matrix.barcodes = read_h5_strings(file, "/obs/barcode/values");
    for (auto &barcode : matrix.barcodes) barcode = normalize_barcode(barcode);
    matrix.genes = read_h5_strings(file, "/var/feature/values");
    const auto input_offsets = read_h5_vector<std::int64_t>(file, "/X/indptr", H5T_NATIVE_LLONG);
    const auto input_columns = read_h5_vector<std::int64_t>(file, "/X/indices", H5T_NATIVE_LLONG);
    const auto input_values = read_h5_vector<double>(file, "/X/data", H5T_NATIVE_DOUBLE);
    if (input_offsets.size() != matrix.barcodes.size() + 1 ||
        input_columns.size() != input_values.size() || input_offsets.empty() ||
        input_offsets.back() < 0 || static_cast<std::size_t>(input_offsets.back()) != input_values.size()) {
        fail("inconsistent CSR arrays in " + path.string());
    }
    matrix.offsets.reserve(input_offsets.size());
    for (const std::int64_t value : input_offsets) {
        if (value < 0) fail("negative CSR offset in " + path.string());
        matrix.offsets.push_back(static_cast<std::uint64_t>(value));
    }
    matrix.columns.reserve(input_columns.size());
    matrix.values.reserve(input_values.size());
    for (std::size_t i = 0; i < input_values.size(); ++i) {
        if (input_columns[i] < 0 || static_cast<std::size_t>(input_columns[i]) >= matrix.genes.size()) {
            fail("invalid CSR column in " + path.string());
        }
        const double rounded = std::round(input_values[i]);
        if (!std::isfinite(input_values[i]) || input_values[i] < 0 ||
            std::abs(input_values[i] - rounded) > 1e-6) {
            fail("non-integral count in " + path.string());
        }
        matrix.columns.push_back(static_cast<std::uint32_t>(input_columns[i]));
        matrix.values.push_back(static_cast<std::uint64_t>(rounded));
        matrix.total += rounded;
    }
    normalize_sparse_rows(matrix);
    return matrix;
}
#endif

void append_row(const SparseMatrix &source, std::size_t row, SparseMatrix &destination) {
    for (std::uint64_t i = source.offsets[row]; i < source.offsets[row + 1]; ++i) {
        destination.columns.push_back(source.columns[i]);
        destination.values.push_back(source.values[i]);
    }
}

void append_merged_rows(const std::vector<std::pair<const SparseMatrix *, std::size_t>> &rows,
                        SparseMatrix &destination) {
    std::unordered_map<std::uint32_t, std::uint64_t> sums;
    std::size_t expected = 0;
    for (const auto &[matrix, row] : rows) expected += matrix->offsets[row + 1] - matrix->offsets[row];
    sums.reserve(expected);
    for (const auto &[matrix, row] : rows) {
        for (std::uint64_t i = matrix->offsets[row]; i < matrix->offsets[row + 1]; ++i) {
            sums[matrix->columns[i]] += matrix->values[i];
        }
    }
    std::vector<std::pair<std::uint32_t, std::uint64_t>> ordered(sums.begin(), sums.end());
    std::sort(ordered.begin(), ordered.end());
    for (const auto &[column, value] : ordered) {
        destination.columns.push_back(column);
        destination.values.push_back(value);
    }
}

SparseMatrix combine_and_collapse(std::vector<SparseMatrix> parts) {
    if (parts.empty()) fail("no matrices to combine");
    for (std::size_t i = 1; i < parts.size(); ++i) {
        if (parts[i].genes != parts[0].genes) fail("gene axes differ across grouped outputs");
    }
    struct RowReference {
        std::string barcode;
        std::size_t part;
        std::size_t row;
    };
    std::vector<RowReference> references;
    for (std::size_t part = 0; part < parts.size(); ++part) {
        for (std::size_t row = 0; row < parts[part].rows(); ++row) {
            references.push_back({parts[part].barcodes[row], part, row});
        }
    }
    std::sort(references.begin(), references.end(), [](const auto &a, const auto &b) {
        if (a.barcode != b.barcode) return a.barcode < b.barcode;
        if (a.part != b.part) return a.part < b.part;
        return a.row < b.row;
    });

    SparseMatrix result;
    result.genes = parts[0].genes;
    result.offsets.push_back(0);
    for (const auto &part : parts) result.total += part.total;
    for (std::size_t begin = 0; begin < references.size();) {
        std::size_t end = begin + 1;
        while (end < references.size() && references[end].barcode == references[begin].barcode) ++end;
        result.barcodes.push_back(references[begin].barcode);
        if (end == begin + 1) {
            const auto &ref = references[begin];
            append_row(parts[ref.part], ref.row, result);
        } else {
            std::vector<std::pair<const SparseMatrix *, std::size_t>> rows;
            rows.reserve(end - begin);
            for (std::size_t i = begin; i < end; ++i) {
                rows.emplace_back(&parts[references[i].part], references[i].row);
            }
            append_merged_rows(rows, result);
        }
        result.offsets.push_back(result.columns.size());
        begin = end;
    }
    return result;
}

struct Group {
    std::string sample;
    std::vector<std::string> barcodes;
};

std::vector<Group> read_groups(const fs::path &config) {
    std::ifstream input(config);
    if (!input) fail("cannot open " + config.string());
    std::vector<Group> groups;
    std::string line;
    bool samples = false;
    bool header = false;
    while (std::getline(input, line)) {
        const auto fields = csv_fields(line);
        const std::string first = fields.empty() ? std::string() : trim(fields[0]);
        if (!samples) {
            if (first == "[samples]") samples = true;
            continue;
        }
        if (!header) {
            if (fields.size() < 2 || trim(fields[0]) != "sample_id" ||
                trim(fields[1]) != "probe_barcode_ids") {
                fail("unexpected [samples] header in " + config.string());
            }
            header = true;
            continue;
        }
        if (first.empty() || first[0] == '[') break;
        if (fields.size() < 2) fail("invalid sample row in " + config.string());
        Group group;
        group.sample = first;
        for (std::string barcode : split(fields[1], '|')) {
            barcode = trim(barcode);
            if (!barcode.empty()) group.barcodes.push_back(barcode);
        }
        if (group.barcodes.empty()) fail("sample has no probe barcode IDs: " + first);
        groups.push_back(std::move(group));
    }
    if (groups.empty()) fail("no [samples] rows in " + config.string());
    return groups;
}

std::string join(const std::vector<std::string> &values, const std::string &delimiter) {
    std::ostringstream out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << delimiter;
        out << values[i];
    }
    return out.str();
}

SparseMatrix read_star_group(const fs::path &root, const Group &group) {
    std::vector<SparseMatrix> parts;
    parts.reserve(group.barcodes.size());
    for (const auto &tag_id : group.barcodes) {
        SparseMatrix part = read_mex(root / tag_id / "Gene" / "filtered");
        for (auto &barcode : part.barcodes) barcode = tagged_cell_key(barcode, tag_id);
        parts.push_back(std::move(part));
    }
    return combine_and_collapse(std::move(parts));
}

SparseMatrix read_cr_group(const fs::path &root, const Group &group,
                           const TagMap &tag8_to_id) {
    std::vector<SparseMatrix> one;
    SparseMatrix matrix = read_mex(
        root / group.sample / "count" / "sample_filtered_feature_bc_matrix", false);
    for (auto &barcode : matrix.barcodes) {
        barcode = cr_cell_key(barcode, group.barcodes, tag8_to_id);
    }
    one.push_back(std::move(matrix));
    return combine_and_collapse(std::move(one));
}

std::unordered_map<std::string, std::string>
read_unique_cr_symbol_map(const fs::path &cr_root, const Group &first_group) {
    const fs::path directory = cr_root / first_group.sample / "count" /
                               "sample_filtered_feature_bc_matrix";
    TextReader input(plain_or_gz(directory, "features.tsv"));
    std::vector<std::pair<std::string, std::string>> rows;
    std::unordered_map<std::string, std::size_t> counts;
    std::string line;
    while (input.getline(line)) {
        const auto fields = split(line, '\t');
        if (fields.size() < 2) continue;
        rows.emplace_back(fields[0], fields[1]);
        ++counts[fields[1]];
    }
    std::unordered_map<std::string, std::string> result;
    result.reserve(rows.size());
    for (const auto &[id, symbol] : rows) {
        if (counts[symbol] == 1) result.emplace(symbol, id);
    }
    return result;
}

SparseMatrix read_cyto_group(
    const fs::path &root, const Group &group,
    const std::unordered_map<std::string, std::string> &symbol_to_id) {
#ifdef FLEX_QUANT_COMPARE_HDF5
    std::vector<SparseMatrix> parts;
    parts.reserve(group.barcodes.size());
    for (const auto &tag_id : group.barcodes) {
        SparseMatrix part = read_cyto_h5ad(root / (tag_id + ".filt.h5ad"));
        for (auto &barcode : part.barcodes) barcode = tagged_cell_key(barcode, tag_id);
        parts.push_back(std::move(part));
    }
    for (auto &part : parts) {
        for (auto &gene : part.genes) {
            const auto translated = symbol_to_id.find(gene);
            if (translated != symbol_to_id.end()) gene = translated->second;
        }
    }
    return combine_and_collapse(std::move(parts));
#else
    static_cast<void>(root);
    static_cast<void>(group);
    static_cast<void>(symbol_to_id);
    fail("cyto input requires a build with -DFLEX_QUANT_COMPARE_HDF5");
#endif
}

struct Alignment {
    std::vector<std::int32_t> query_gene_to_common;
    std::vector<std::int32_t> cr_gene_to_common;
    std::vector<std::pair<std::size_t, std::size_t>> cells;
    std::size_t common_gene_count = 0;
    double jaccard = std::numeric_limits<double>::quiet_NaN();
};

std::unordered_map<std::string, std::size_t> unique_index(const std::vector<std::string> &values) {
    std::unordered_map<std::string, std::size_t> index;
    std::unordered_map<std::string, std::size_t> counts;
    counts.reserve(values.size());
    for (const auto &value : values) ++counts[value];
    index.reserve(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (counts[values[i]] == 1) index.emplace(values[i], i);
    }
    return index;
}

Alignment align_matrices(const SparseMatrix &query, const SparseMatrix &cr) {
    Alignment result;
    const auto query_index = unique_index(query.genes);
    const auto cr_index = unique_index(cr.genes);
    std::vector<std::string> common_genes;
    common_genes.reserve(std::min(query_index.size(), cr_index.size()));
    for (const auto &[gene, unused] : query_index) {
        if (cr_index.find(gene) != cr_index.end()) common_genes.push_back(gene);
    }
    std::sort(common_genes.begin(), common_genes.end());
    result.common_gene_count = common_genes.size();
    result.query_gene_to_common.assign(query.genes.size(), -1);
    result.cr_gene_to_common.assign(cr.genes.size(), -1);
    for (std::size_t i = 0; i < common_genes.size(); ++i) {
        result.query_gene_to_common[query_index.at(common_genes[i])] = static_cast<std::int32_t>(i);
        result.cr_gene_to_common[cr_index.at(common_genes[i])] = static_cast<std::int32_t>(i);
    }

    std::size_t qi = 0, ci = 0, union_count = 0;
    while (qi < query.barcodes.size() || ci < cr.barcodes.size()) {
        if (ci == cr.barcodes.size() ||
            (qi < query.barcodes.size() && query.barcodes[qi] < cr.barcodes[ci])) {
            ++qi;
        } else if (qi == query.barcodes.size() || cr.barcodes[ci] < query.barcodes[qi]) {
            ++ci;
        } else {
            result.cells.emplace_back(qi++, ci++);
        }
        ++union_count;
    }
    result.jaccard = union_count == 0 ? std::numeric_limits<double>::quiet_NaN()
                                      : static_cast<double>(result.cells.size()) / union_count;
    return result;
}

void fill_row(const SparseMatrix &matrix, std::size_t row,
              const std::vector<std::int32_t> &gene_map,
              std::vector<std::uint64_t> &dense) {
    std::fill(dense.begin(), dense.end(), 0);
    for (std::uint64_t i = matrix.offsets[row]; i < matrix.offsets[row + 1]; ++i) {
        const std::int32_t common = gene_map[matrix.columns[i]];
        if (common >= 0) dense[common] += matrix.values[i];
    }
}

double pearson(const std::vector<std::uint64_t> &x, const std::vector<std::uint64_t> &y,
               const std::vector<std::size_t> *selection = nullptr, bool log_transform = false) {
    const std::size_t n = selection == nullptr ? x.size() : selection->size();
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    long double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    for (std::size_t p = 0; p < n; ++p) {
        const std::size_t i = selection == nullptr ? p : (*selection)[p];
        const long double a = log_transform ? std::log1p(static_cast<long double>(x[i])) : x[i];
        const long double b = log_transform ? std::log1p(static_cast<long double>(y[i])) : y[i];
        sx += a;
        sy += b;
        sxx += a * a;
        syy += b * b;
        sxy += a * b;
    }
    const long double vx = sxx - sx * sx / n;
    const long double vy = syy - sy * sy / n;
    if (vx <= 0 || vy <= 0) return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>((sxy - sx * sy / n) / std::sqrt(vx * vy));
}

std::unordered_map<std::uint64_t, double>
average_ranks(const std::vector<std::uint64_t> &values,
              const std::vector<std::size_t> *selection) {
    const std::size_t n = selection == nullptr ? values.size() : selection->size();
    std::unordered_map<std::uint64_t, std::size_t> counts;
    counts.reserve(128);
    for (std::size_t p = 0; p < n; ++p) {
        const std::size_t i = selection == nullptr ? p : (*selection)[p];
        ++counts[values[i]];
    }
    std::vector<std::pair<std::uint64_t, std::size_t>> ordered(counts.begin(), counts.end());
    std::sort(ordered.begin(), ordered.end());
    std::unordered_map<std::uint64_t, double> ranks;
    ranks.reserve(ordered.size());
    std::size_t below = 0;
    for (const auto &[value, count] : ordered) {
        ranks[value] = 0.5 * (static_cast<double>(below + 1) + static_cast<double>(below + count));
        below += count;
    }
    return ranks;
}

double spearman(const std::vector<std::uint64_t> &x, const std::vector<std::uint64_t> &y,
                const std::vector<std::size_t> *selection = nullptr) {
    const std::size_t n = selection == nullptr ? x.size() : selection->size();
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    const auto xr = average_ranks(x, selection);
    const auto yr = average_ranks(y, selection);
    const long double mean = 0.5L * (n + 1.0L);
    long double xx = 0, yy = 0, xy = 0;
    for (std::size_t p = 0; p < n; ++p) {
        const std::size_t i = selection == nullptr ? p : (*selection)[p];
        const long double a = xr.at(x[i]) - mean;
        const long double b = yr.at(y[i]) - mean;
        xx += a * a;
        yy += b * b;
        xy += a * b;
    }
    if (xx <= 0 || yy <= 0) return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(xy / std::sqrt(xx * yy));
}

double mean(std::vector<double> values) {
    values.erase(std::remove_if(values.begin(), values.end(), [](double value) { return std::isnan(value); }), values.end());
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double median(std::vector<double> values) {
    values.erase(std::remove_if(values.begin(), values.end(), [](double value) { return std::isnan(value); }), values.end());
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    if (values.size() % 2) return values[middle];
    return 0.5 * (values[middle - 1] + values[middle]);
}

double percentile(std::vector<double> values, double p) {
    values.erase(std::remove_if(values.begin(), values.end(), [](double value) { return std::isnan(value); }), values.end());
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const double position = (values.size() - 1) * p;
    const std::size_t low = static_cast<std::size_t>(std::floor(position));
    const std::size_t high = static_cast<std::size_t>(std::ceil(position));
    const double fraction = position - low;
    return values[low] * (1.0 - fraction) + values[high] * fraction;
}

double minimum(const std::vector<double> &values) {
    double result = std::numeric_limits<double>::infinity();
    for (double value : values) if (!std::isnan(value)) result = std::min(result, value);
    return std::isinf(result) ? std::numeric_limits<double>::quiet_NaN() : result;
}

struct GroupMetrics {
    std::string sample;
    std::string tags;
    std::size_t query_cells = 0;
    std::size_t cr_cells = 0;
    std::size_t common_cells = 0;
    std::size_t paper_genes = 0;
    double jaccard = 0;
    double raw_cell_pear_median = 0;
    double raw_fraction_exact = 0;
    double raw_p01 = 0;
    double raw_worst = 0;
    double raw_cell_spear_median = 0;
    double raw_gene_pear = 0;
    double paper_cell_pear_mean = 0;
    double paper_cell_pear_median = 0;
    double paper_cell_spear_mean = 0;
    double paper_cell_spear_median = 0;
    double paper_gene_pear = 0;
    double own_umi_ratio = 0;
    double common_umi_ratio = 0;
};

GroupMetrics compare_group(
    const fs::path &query_root, const fs::path &cr_root, const Group &group,
    const std::string &input_kind,
    const std::unordered_map<std::string, std::string> &symbol_to_id,
    const TagMap &tag8_to_id) {
    std::cerr << "reading " << group.sample << " (" << join(group.barcodes, "+") << ")\n";
    SparseMatrix query = input_kind == "star"
                             ? read_star_group(query_root, group)
                             : read_cyto_group(query_root, group, symbol_to_id);
    SparseMatrix cr = read_cr_group(cr_root, group, tag8_to_id);
    const Alignment alignment = align_matrices(query, cr);
    if (alignment.common_gene_count == 0 || alignment.cells.empty()) {
        fail("no common genes or cells for " + group.sample);
    }

    const std::size_t genes = alignment.common_gene_count;
    const std::size_t cells = alignment.cells.size();
    std::vector<std::uint64_t> query_totals(genes, 0), cr_totals(genes, 0);
    std::vector<std::uint32_t> query_detected(genes, 0), cr_detected(genes, 0);
    for (const auto &[query_row, cr_row] : alignment.cells) {
        for (std::uint64_t i = query.offsets[query_row]; i < query.offsets[query_row + 1]; ++i) {
            const std::int32_t common = alignment.query_gene_to_common[query.columns[i]];
            if (common >= 0) {
                query_totals[common] += query.values[i];
                if (query.values[i] != 0) ++query_detected[common];
            }
        }
        for (std::uint64_t i = cr.offsets[cr_row]; i < cr.offsets[cr_row + 1]; ++i) {
            const std::int32_t common = alignment.cr_gene_to_common[cr.columns[i]];
            if (common >= 0) {
                cr_totals[common] += cr.values[i];
                if (cr.values[i] != 0) ++cr_detected[common];
            }
        }
    }
    std::vector<std::size_t> paper_genes;
    for (std::size_t gene = 0; gene < genes; ++gene) {
        if (query_totals[gene] >= 20 && cr_totals[gene] >= 20 &&
            query_detected[gene] >= 0.01 * cells && cr_detected[gene] >= 0.01 * cells) {
            paper_genes.push_back(gene);
        }
    }

    std::vector<double> raw_pear(cells), raw_spear(cells), paper_pear(cells), paper_spear(cells);
#pragma omp parallel
    {
        std::vector<std::uint64_t> query_dense(genes), cr_dense(genes);
#pragma omp for schedule(dynamic, 32)
        for (std::int64_t cell = 0; cell < static_cast<std::int64_t>(cells); ++cell) {
            fill_row(query, alignment.cells[cell].first, alignment.query_gene_to_common, query_dense);
            fill_row(cr, alignment.cells[cell].second, alignment.cr_gene_to_common, cr_dense);
            raw_pear[cell] = pearson(query_dense, cr_dense);
            raw_spear[cell] = spearman(query_dense, cr_dense);
            paper_pear[cell] = pearson(query_dense, cr_dense, &paper_genes, true);
            paper_spear[cell] = spearman(query_dense, cr_dense, &paper_genes);
        }
    }

    std::size_t exact = 0, finite = 0;
    for (double value : raw_pear) {
        if (!std::isnan(value)) {
            ++finite;
            if (value >= 0.999999) ++exact;
        }
    }
    const long double query_common_sum = std::accumulate(
        query_totals.begin(), query_totals.end(), static_cast<long double>(0));
    const long double cr_common_sum = std::accumulate(
        cr_totals.begin(), cr_totals.end(), static_cast<long double>(0));

    GroupMetrics metrics;
    metrics.sample = group.sample;
    metrics.tags = join(group.barcodes, "+");
    metrics.query_cells = query.rows();
    metrics.cr_cells = cr.rows();
    metrics.common_cells = cells;
    metrics.paper_genes = paper_genes.size();
    metrics.jaccard = alignment.jaccard;
    metrics.raw_cell_pear_median = median(raw_pear);
    metrics.raw_fraction_exact = finite == 0 ? std::numeric_limits<double>::quiet_NaN()
                                             : static_cast<double>(exact) / finite;
    metrics.raw_p01 = percentile(raw_pear, 0.01);
    metrics.raw_worst = minimum(raw_pear);
    metrics.raw_cell_spear_median = median(raw_spear);
    metrics.raw_gene_pear = pearson(query_totals, cr_totals);
    metrics.paper_cell_pear_mean = mean(paper_pear);
    metrics.paper_cell_pear_median = median(paper_pear);
    metrics.paper_cell_spear_mean = mean(paper_spear);
    metrics.paper_cell_spear_median = median(paper_spear);
    metrics.paper_gene_pear = pearson(query_totals, cr_totals, &paper_genes, true);
    metrics.own_umi_ratio = static_cast<double>(query.total / cr.total);
    metrics.common_umi_ratio = static_cast<double>(query_common_sum / cr_common_sum);
    return metrics;
}

void write_raw(const fs::path &path, const std::string &label, const std::string &input_kind,
               const fs::path &query_root,
               const fs::path &cr_root, const fs::path &config,
               const std::vector<GroupMetrics> &metrics) {
    std::ofstream out(path);
    if (!out) fail("cannot write " + path.string());
    out << "# concordance vs Cell Ranger 9.0.1 - " << label << '\n'
        << "# query_kind=" << input_kind << " query_outputs=" << query_root.string() << '\n'
        << "# genes matched by Ensembl ID via Cell Ranger's feature table\n"
        << "# cr_root=" << cr_root.string() << " cr_config=" << config.string() << '\n';
    out << std::left << std::setw(29) << "sample" << std::setw(14) << "tags"
        << std::right << std::setw(18) << "cells query/CR" << std::setw(9) << "Jaccard"
        << std::setw(10) << "cellPear" << std::setw(10) << "r>=1-1e6"
        << std::setw(9) << "p01" << std::setw(9) << "worst" << std::setw(11) << "cellSpear"
        << std::setw(10) << "genePear" << '\n';
    std::vector<double> j, cp, exact, p01, worst, cs, gp;
    out << std::fixed;
    for (const auto &m : metrics) {
        std::ostringstream cells;
        cells << m.query_cells << '/' << m.cr_cells;
        out << std::left << std::setw(29) << m.sample << std::setw(14) << m.tags
            << std::right << std::setw(18) << cells.str() << std::setprecision(4)
            << std::setw(9) << m.jaccard << std::setprecision(6) << std::setw(10)
            << m.raw_cell_pear_median << std::setprecision(4) << std::setw(10)
            << m.raw_fraction_exact << std::setw(9) << m.raw_p01 << std::setw(9)
            << m.raw_worst << std::setprecision(6) << std::setw(11)
            << m.raw_cell_spear_median << std::setw(10) << m.raw_gene_pear << '\n';
        j.push_back(m.jaccard); cp.push_back(m.raw_cell_pear_median);
        exact.push_back(m.raw_fraction_exact); p01.push_back(m.raw_p01);
        worst.push_back(m.raw_worst); cs.push_back(m.raw_cell_spear_median);
        gp.push_back(m.raw_gene_pear);
    }
    out << std::left << std::setw(29) << "median" << std::setw(14) << "" << std::setw(18) << ""
        << std::right << std::setprecision(4) << std::setw(9) << median(j)
        << std::setprecision(6) << std::setw(10) << median(cp)
        << std::setprecision(4) << std::setw(10) << median(exact)
        << std::setw(9) << median(p01) << std::setw(9) << minimum(worst)
        << std::setprecision(6) << std::setw(11) << median(cs) << std::setw(10) << median(gp) << '\n'
        << "# r>=1-1e6 is the fraction of common cells at Pearson >= 0.999999; p01 is the 1st percentile and worst is the minimum.\n"
        << "# Cells called by only one tool are excluded from correlation columns; the Jaccard column reports that disagreement.\n";
}

void write_paper(const fs::path &path, const std::string &label, const fs::path &cr_root,
                 const fs::path &config, const std::vector<GroupMetrics> &metrics) {
    std::ofstream out(path);
    if (!out) fail("cannot write " + path.string());
    out << "# paper-protocol concordance vs Cell Ranger 9.0.1 - " << label << '\n'
        << "# log1p; genes >=20 counts in both and detected in >=1% of common cells; Ensembl-ID matched\n"
        << "# cr_root=" << cr_root.string() << " cr_config=" << config.string() << '\n';
    out << std::left << std::setw(29) << "sample" << std::setw(14) << "tags"
        << std::right << std::setw(12) << "cells" << std::setw(7) << "genes"
        << std::setw(9) << "Jaccard" << std::setw(15) << "cellPear mean"
        << std::setw(9) << "median" << std::setw(16) << "cellSpear mean"
        << std::setw(9) << "median" << std::setw(10) << "genePear" << '\n';
    std::vector<double> j, cpm, cpmed, csm, csmed, gp;
    out << std::fixed;
    for (const auto &m : metrics) {
        out << std::left << std::setw(29) << m.sample << std::setw(14) << m.tags
            << std::right << std::setw(12) << m.common_cells << std::setw(7) << m.paper_genes
            << std::setprecision(4) << std::setw(9) << m.jaccard << std::setprecision(6)
            << std::setw(15) << m.paper_cell_pear_mean << std::setw(9)
            << m.paper_cell_pear_median << std::setw(16) << m.paper_cell_spear_mean
            << std::setw(9) << m.paper_cell_spear_median << std::setw(10)
            << m.paper_gene_pear << '\n';
        j.push_back(m.jaccard); cpm.push_back(m.paper_cell_pear_mean);
        cpmed.push_back(m.paper_cell_pear_median); csm.push_back(m.paper_cell_spear_mean);
        csmed.push_back(m.paper_cell_spear_median); gp.push_back(m.paper_gene_pear);
    }
    for (const auto &[name, use_mean] :
         std::vector<std::pair<std::string, bool>>{{"mean over samples", true}, {"median over samples", false}}) {
        const auto reduce = [use_mean](const std::vector<double> &values) {
            return use_mean ? mean(values) : median(values);
        };
        out << std::left << std::setw(29) << name << std::setw(14) << "" << std::setw(12) << ""
            << std::setw(7) << "" << std::right << std::setprecision(4) << std::setw(9) << reduce(j)
            << std::setprecision(6) << std::setw(15) << reduce(cpm) << std::setw(9) << reduce(cpmed)
            << std::setw(16) << reduce(csm) << std::setw(9) << reduce(csmed)
            << std::setw(10) << reduce(gp) << '\n';
    }
}

void write_umi(const fs::path &path, const std::string &label,
               const std::vector<GroupMetrics> &metrics) {
    std::ofstream out(path);
    if (!out) fail("cannot write " + path.string());
    out << "# UMI totals vs Cell Ranger - " << label << '\n';
    out << std::left << std::setw(29) << "sample" << std::setw(14) << "tags"
        << std::right << std::setw(15) << "query/CR own" << std::setw(17)
        << "query/CR common" << std::setw(13) << "common cells" << '\n';
    std::vector<double> own, common;
    out << std::fixed << std::setprecision(6);
    for (const auto &m : metrics) {
        out << std::left << std::setw(29) << m.sample << std::setw(14) << m.tags
            << std::right << std::setw(15) << m.own_umi_ratio << std::setw(17)
            << m.common_umi_ratio << std::setw(13) << m.common_cells << '\n';
        own.push_back(m.own_umi_ratio);
        common.push_back(m.common_umi_ratio);
    }
    out << std::left << std::setw(29) << "median" << std::setw(14) << "" << std::right
        << std::setw(15) << median(own) << std::setw(17) << median(common) << '\n'
        << std::left << std::setw(29) << "mean" << std::setw(14) << "" << std::right
        << std::setw(15) << mean(own) << std::setw(17) << mean(common) << '\n';
}

struct RawTotals {
    std::vector<std::string> barcodes;
    std::vector<std::uint64_t> totals;
};

RawTotals read_raw_totals(const fs::path &directory) {
    RawTotals result;
    result.barcodes = read_barcodes(directory, false);
    TextReader input(plain_or_gz(directory, "matrix.mtx"));
    std::string line;
    if (!input.getline(line) || line.rfind("%%MatrixMarket", 0) != 0) {
        fail("unexpected MatrixMarket header in " + directory.string());
    }
    do {
        if (!input.getline(line)) fail("missing MatrixMarket dimensions in " + directory.string());
    } while (line.empty() || line[0] == '%');
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    std::uint64_t feature_count = 0, cell_count = 0, expected_nnz = 0;
    if (!parse_unsigned(cursor, end, feature_count) ||
        !parse_unsigned(cursor, end, cell_count) ||
        !parse_unsigned(cursor, end, expected_nnz) || cell_count != result.barcodes.size()) {
        fail("invalid MatrixMarket dimensions in " + directory.string());
    }
    result.totals.assign(cell_count, 0);
    std::uint64_t observed_nnz = 0;
    while (input.getline(line)) {
        if (line.empty() || line[0] == '%') continue;
        cursor = line.data();
        end = cursor + line.size();
        std::uint64_t feature = 0, cell = 0, value = 0;
        if (!parse_unsigned(cursor, end, feature) || !parse_unsigned(cursor, end, cell) ||
            !parse_unsigned(cursor, end, value) || feature == 0 || feature > feature_count ||
            cell == 0 || cell > cell_count) {
            fail("invalid MatrixMarket entry in " + directory.string() + ": " + line);
        }
        result.totals[cell - 1] += value;
        ++observed_nnz;
    }
    if (observed_nnz != expected_nnz) {
        fail("MatrixMarket nnz mismatch in " + directory.string());
    }
    return result;
}

TagMap read_sample_tags(const fs::path &path) {
    std::ifstream input(path);
    if (!input) fail("cannot open " + path.string());
    TagMap result;
    std::string line;
    std::size_t number = 0;
    while (std::getline(input, line)) {
        ++number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto fields = split(line, '\t');
        if (fields.size() < 2 || fields[0].empty() || fields[1].empty()) {
            fail("invalid sample whitelist row " + std::to_string(number) + " in " + path.string());
        }
        if (!result.emplace(fields[0], fields[1]).second) {
            fail("duplicate tag id " + fields[0] + " in " + path.string());
        }
    }
    return result;
}

TagMap invert_sample_tags(const TagMap &id_to_tag8, const fs::path &path) {
    TagMap result;
    result.reserve(id_to_tag8.size());
    for (const auto &[tag_id, tag8] : id_to_tag8) {
        if (!result.emplace(tag8, tag_id).second) {
            fail("duplicate tag sequence " + tag8 + " in " + path.string());
        }
    }
    return result;
}

std::unordered_set<std::string> tagged_barcode_set(const fs::path &directory,
                                                   const std::string &tag_id) {
    const auto values = read_barcodes(directory);
    std::unordered_set<std::string> result;
    result.reserve(values.size());
    for (const auto &barcode : values) result.insert(tagged_cell_key(barcode, tag_id));
    return result;
}

std::unordered_set<std::string> cr_barcode_set(const fs::path &directory,
                                               const Group &group,
                                               const TagMap &tag8_to_id) {
    const auto values = read_barcodes(directory, false);
    std::unordered_set<std::string> result;
    result.reserve(values.size());
    for (const auto &barcode : values) {
        result.insert(cr_cell_key(barcode, group.barcodes, tag8_to_id));
    }
    return result;
}

std::unordered_set<std::string> cyto_barcode_set(const fs::path &root,
                                                 const Group &group) {
#ifdef FLEX_QUANT_COMPARE_HDF5
    std::unordered_set<std::string> result;
    for (const auto &barcode_id : group.barcodes) {
        H5File file(H5Fopen((root / (barcode_id + ".filt.h5ad")).c_str(),
                           H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file.valid()) fail("cannot open cyto output for " + barcode_id);
        auto barcodes = read_h5_strings(file, "/obs/barcode/values");
        for (const auto &barcode : barcodes) {
            result.insert(tagged_cell_key(barcode, barcode_id));
        }
    }
    return result;
#else
    static_cast<void>(root);
    static_cast<void>(group);
    fail("cyto input requires a build with -DFLEX_QUANT_COMPARE_HDF5");
#endif
}

void write_cell_calling(const fs::path &path, const fs::path &star_run,
                        const fs::path &cr_root, const fs::path &config,
                        const fs::path &sample_whitelist,
                        const fs::path &cyto_counts,
                        const std::vector<Group> &groups,
                        const TagMap &tag8_to_id) {
    std::cerr << "streaming STAR raw matrix for cell-calling curves\n";
    const RawTotals raw = read_raw_totals(star_run / "Solo.out" / "Gene" / "raw");
    const auto tags = read_sample_tags(sample_whitelist);
    const auto raw_tag8_to_id = invert_sample_tags(tags, sample_whitelist);
    std::ofstream out(path);
    if (!out) fail("cannot write " + path.string());
    out << "# Cell calling vs Cell Ranger 9.0.1 filtered barcodes\n"
        << "# Frontier from STAR-Flex raw UMI ranking; operating points marked.\n"
        << "# cr_root=" << cr_root.string() << " cr_config=" << config.string() << '\n';
    out << std::left << std::setw(29) << "sample" << std::setw(14) << "tags"
        << std::setw(11) << "tool" << std::right << std::setw(8) << "called"
        << std::setw(8) << "CR" << std::setw(8) << "TP" << std::setw(11)
        << "precision" << std::setw(9) << "recall" << std::setw(8) << "F1"
        << std::setw(8) << "AUPRC" << '\n';

    for (const auto &group : groups) {
        std::unordered_set<std::string> star_calls;
        for (const auto &barcode_id : group.barcodes) {
            const auto tag = tags.find(barcode_id);
            if (tag == tags.end()) fail(barcode_id + " is absent from " + sample_whitelist.string());
            const auto calls = tagged_barcode_set(
                star_run / "per_sample" / barcode_id / "Gene" / "filtered", barcode_id);
            star_calls.insert(calls.begin(), calls.end());
        }
        const auto cr_calls = cr_barcode_set(
            cr_root / group.sample / "count" / "sample_filtered_feature_bc_matrix",
            group, tag8_to_id);

        // std::map reproduces np.unique's lexicographic CB order.  The stable
        // score sort below therefore also reproduces NumPy's tie behavior.
        std::map<std::string, std::uint64_t> score_by_barcode;
        for (std::size_t row = 0; row < raw.barcodes.size(); ++row) {
            const std::string &composite = raw.barcodes[row];
            if (composite.size() < 24) continue;
            const auto tag = raw_tag8_to_id.find(composite.substr(16, 8));
            if (tag == raw_tag8_to_id.end() ||
                std::find(group.barcodes.begin(), group.barcodes.end(), tag->second) ==
                    group.barcodes.end()) {
                continue;
            }
            score_by_barcode[tagged_cell_key(composite, tag->second)] += raw.totals[row];
        }
        std::vector<std::pair<std::string, std::uint64_t>> candidates(
            score_by_barcode.begin(), score_by_barcode.end());
        std::vector<std::size_t> order(candidates.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return candidates[a].second > candidates[b].second;
        });
        std::size_t positives = 0;
        for (const auto &[barcode, unused] : candidates) {
            positives += cr_calls.find(barcode) != cr_calls.end();
        }
        std::vector<double> precision(order.size()), recall(order.size());
        std::size_t hits = 0;
        for (std::size_t rank = 0; rank < order.size(); ++rank) {
            hits += cr_calls.find(candidates[order[rank]].first) != cr_calls.end();
            precision[rank] = static_cast<double>(hits) / (rank + 1);
            recall[rank] = static_cast<double>(hits) / std::max<std::size_t>(positives, 1);
        }
        double auprc = std::numeric_limits<double>::quiet_NaN();
        if (recall.size() > 1) {
            auprc = 0;
            for (std::size_t i = 1; i < recall.size(); ++i) {
                auprc += 0.5 * (precision[i - 1] + precision[i]) * (recall[i] - recall[i - 1]);
            }
        }
        const auto operating_point = [&](const std::unordered_set<std::string> &calls,
                                         const std::string &tool) {
            std::size_t true_positives = 0;
            for (const auto &barcode : calls) {
                true_positives += cr_calls.find(barcode) != cr_calls.end();
            }
            const double p = calls.empty() ? std::numeric_limits<double>::quiet_NaN()
                                            : static_cast<double>(true_positives) / calls.size();
            const double r = cr_calls.empty() ? std::numeric_limits<double>::quiet_NaN()
                                               : static_cast<double>(true_positives) / cr_calls.size();
            const double f1 = p + r == 0 ? std::numeric_limits<double>::quiet_NaN()
                                         : 2 * p * r / (p + r);
            out << std::left << std::setw(29) << group.sample << std::setw(14)
                << join(group.barcodes, "+") << std::setw(11) << tool << std::right
                << std::setw(8) << calls.size() << std::setw(8) << cr_calls.size()
                << std::setw(8) << true_positives << std::fixed << std::setprecision(4)
                << std::setw(11) << p << std::setw(9) << r << std::setw(8) << f1
                << std::setw(8) << auprc << '\n';
        };
        operating_point(star_calls, "STAR-Flex");
        if (!cyto_counts.empty()) operating_point(cyto_barcode_set(cyto_counts, group), "cyto");
    }
    out << "# Native report omits the optional PNG; use cell_calling_pr.py only when a figure is required.\n";
}

struct Options {
    fs::path query_root;
    fs::path star_run;
    fs::path cr_root;
    fs::path cr_config;
    fs::path tag_map;
    fs::path sample_whitelist;
    fs::path cyto_counts;
    fs::path out_prefix;
    std::string label;
    std::string input_kind = "star";
    int threads = omp_get_max_threads();
};

Options parse_options(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) fail("missing value for " + argument);
            return argv[i];
        };
        if (argument == "--query-root") options.query_root = value();
        else if (argument == "--star-run") options.star_run = value();
        else if (argument == "--cr-root") options.cr_root = value();
        else if (argument == "--cr-config") options.cr_config = value();
        else if (argument == "--tag-map") options.tag_map = value();
        else if (argument == "--sample-whitelist") options.sample_whitelist = value();
        else if (argument == "--cyto-counts") options.cyto_counts = value();
        else if (argument == "--out-prefix") options.out_prefix = value();
        else if (argument == "--label") options.label = value();
        else if (argument == "--input-kind") options.input_kind = value();
        else if (argument == "--threads") options.threads = std::stoi(value());
        else if (argument == "--help" || argument == "-h") {
            std::cout << "Usage: flex_quant_compare --query-root DIR --input-kind star|cyto "
                         "--label NAME --cr-root DIR "
                         "--cr-config FILE --out-prefix PATH [--tag-map FILE] [--threads N] "
                         "[--star-run DIR --sample-whitelist FILE [--cyto-counts DIR]]\n";
            std::exit(0);
        } else fail("unknown argument: " + argument);
    }
    if (options.query_root.empty() || options.cr_root.empty() || options.cr_config.empty() ||
        options.out_prefix.empty() || options.label.empty()) {
        fail("--query-root, --label, --cr-root, --cr-config, and --out-prefix are required");
    }
    if (options.threads < 1) fail("--threads must be positive");
    if (options.input_kind != "star" && options.input_kind != "cyto") {
        fail("--input-kind must be star or cyto");
    }
    if (options.star_run.empty() != options.sample_whitelist.empty()) {
        fail("--star-run and --sample-whitelist must be supplied together");
    }
    if (!options.cyto_counts.empty() && options.star_run.empty()) {
        fail("--cyto-counts requires --star-run and --sample-whitelist");
    }
#ifndef FLEX_QUANT_COMPARE_HDF5
    if (options.input_kind == "cyto" || !options.cyto_counts.empty()) {
        fail("cyto input requires a build with -DFLEX_QUANT_COMPARE_HDF5");
    }
#endif
    return options;
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        omp_set_num_threads(options.threads);
        const auto groups = read_groups(options.cr_config);
        const fs::path tag_map_path = options.tag_map.empty()
                                          ? options.sample_whitelist
                                          : options.tag_map;
        const TagMap id_to_tag8 = tag_map_path.empty()
                                      ? TagMap()
                                      : read_sample_tags(tag_map_path);
        const TagMap tag8_to_id = tag_map_path.empty()
                                      ? TagMap()
                                      : invert_sample_tags(id_to_tag8, tag_map_path);
        for (const auto &group : groups) {
            if (group.barcodes.size() > 1 && tag_map_path.empty()) {
                fail("grouped comparison for " + group.sample + " requires --tag-map");
            }
            if (!tag_map_path.empty()) {
                for (const auto &tag_id : group.barcodes) {
                    if (id_to_tag8.find(tag_id) == id_to_tag8.end()) {
                        fail(tag_id + " is absent from " + tag_map_path.string());
                    }
                }
            }
        }
        const auto symbol_to_id = options.input_kind == "cyto"
                                      ? read_unique_cr_symbol_map(options.cr_root, groups.front())
                                      : std::unordered_map<std::string, std::string>();
        std::vector<GroupMetrics> metrics;
        metrics.reserve(groups.size());
        for (const auto &group : groups) {
            metrics.push_back(compare_group(options.query_root, options.cr_root, group,
                                            options.input_kind, symbol_to_id, tag8_to_id));
        }
        if (!options.out_prefix.parent_path().empty()) fs::create_directories(options.out_prefix.parent_path());
        write_raw(options.out_prefix.string() + ".raw-concordance.txt", options.label,
                  options.input_kind, options.query_root, options.cr_root,
                  options.cr_config, metrics);
        write_paper(options.out_prefix.string() + ".paper-protocol.txt", options.label,
                    options.cr_root, options.cr_config, metrics);
        write_umi(options.out_prefix.string() + ".umi-totals.txt", options.label, metrics);
        if (!options.star_run.empty()) {
            write_cell_calling(options.out_prefix.string() + ".cell-calling-pr.txt",
                               options.star_run, options.cr_root, options.cr_config,
                               options.sample_whitelist, options.cyto_counts, groups,
                               tag8_to_id);
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "flex_quant_compare: " << error.what() << '\n';
        return 1;
    }
}

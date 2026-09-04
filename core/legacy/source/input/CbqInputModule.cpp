#include "input/CbqInputModule.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cstring>
#include <dlfcn.h>
#include <fstream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace star {
namespace input {

InputSourcePlan make_cbq_input_source_plan(
    const std::vector<std::vector<std::string>>& read_files_names,
    const std::vector<std::string>& read_groups,
    uint32_t mate_count) {
    InputSourcePlan plan;
    plan.format = SourceFormat::Binseq;
    plan.module_name = "Cbq";
    plan.mate_files = read_files_names;
    plan.read_groups = read_groups;
    plan.read_files_n = read_files_names.empty()
        ? 0
        : static_cast<uint32_t>(read_files_names.front().size());
    plan.mate_count = mate_count;
    plan.preserves_source_order = true;
    plan.uses_helper_stream = false;
    plan.uses_internal_gzip = false;
    return plan;
}

namespace {

const uint64_t PRESENCE_PAIRED = 1ULL << 0;
const uint64_t PRESENCE_QUALITIES = 1ULL << 1;
const uint64_t PRESENCE_HEADERS = 1ULL << 2;
const uint64_t PRESENCE_FLAGS = 1ULL << 3;

bool set_error(std::string* error, const std::string& message) {
    if (error != nullptr) {
        *error = message;
    }
    return false;
}

uint64_t read_le64(const uint8_t* bytes) {
    uint64_t value = 0;
    for (unsigned ii = 0; ii < 8; ++ii) {
        value |= static_cast<uint64_t>(bytes[ii]) << (8 * ii);
    }
    return value;
}

uint16_t read_le16(const uint8_t* bytes) {
    return static_cast<uint16_t>(bytes[0]) |
           static_cast<uint16_t>(static_cast<uint16_t>(bytes[1]) << 8);
}

bool checked_size(uint64_t value, size_t* out) {
    if (value > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        return false;
    }
    *out = static_cast<size_t>(value);
    return true;
}

bool is_space_char(char value) {
    return std::isspace(static_cast<unsigned char>(value)) != 0;
}

CbqByteSpan make_span(const std::string& value, size_t begin, size_t size) {
    CbqByteSpan span;
    span.data = value.data() + begin;
    span.size = size;
    return span;
}

CbqByteSpan make_empty_span() {
    CbqByteSpan span;
    span.data = "";
    span.size = 0;
    return span;
}

void assign_from_span(std::string* out, CbqByteSpan span) {
    if (span.data == nullptr || span.size == 0) {
        out->clear();
        return;
    }
    out->assign(span.data, span.size);
}

CbqSegmentRole role_for_segment(uint32_t segment_index) {
    if (segment_index == 0) {
        return CbqSegmentRole::Read1;
    }
    if (segment_index == 1) {
        return CbqSegmentRole::Read2;
    }
    return CbqSegmentRole::Unknown;
}

char parse_read_filter_span(char record_type, CbqByteSpan read_name_extra) {
    if (record_type != '@') {
        return 'N';
    }

    size_t field_size = 0;
    while (field_size < read_name_extra.size && !is_space_char(read_name_extra.data[field_size])) {
        ++field_size;
    }
    if (field_size > 3 &&
        read_name_extra.data[1] == ':' &&
        read_name_extra.data[2] == 'Y' &&
        read_name_extra.data[3] == ':') {
        return 'Y';
    }
    return 'N';
}

void apply_quality_conversion(std::string& quality, int add) {
    if (add == 0) {
        return;
    }
    for (char& qchar : quality) {
        int qs = static_cast<int>(qchar) + add;
        if (qs < 33) {
            qs = 33;
        } else if (qs > 126) {
            qs = 126;
        }
        qchar = static_cast<char>(qs);
    }
}

struct ParsedHeaderSpan {
    CbqByteSpan read_name;
    CbqByteSpan read_name_extra;
    char read_filter = 'N';
};

ParsedHeaderSpan parse_header_payload_span(
    const std::string& backing,
    size_t payload_begin,
    size_t payload_size,
    bool has_qualities,
    const std::vector<char>& read_name_separator_chars) {

    char record_type = has_qualities ? '@' : '>';
    size_t cursor = payload_begin;
    const size_t payload_end = payload_begin + payload_size;
    if (cursor < payload_end && (backing[cursor] == '@' || backing[cursor] == '>')) {
        record_type = backing[cursor];
        ++cursor;
    }

    while (cursor < payload_end && is_space_char(backing[cursor])) {
        ++cursor;
    }

    const size_t name_begin = cursor;
    while (cursor < payload_end && !is_space_char(backing[cursor])) {
        ++cursor;
    }
    size_t name_end = cursor;
    for (char separator : read_name_separator_chars) {
        for (size_t pos = name_begin; pos < name_end; ++pos) {
            if (backing[pos] == separator) {
                name_end = pos;
                break;
            }
        }
    }

    while (cursor < payload_end && is_space_char(backing[cursor])) {
        ++cursor;
    }

    ParsedHeaderSpan parsed;
    parsed.read_name = make_span(backing, name_begin, name_end - name_begin);
    parsed.read_name_extra = make_span(backing, cursor, payload_end - cursor);
    parsed.read_filter = parse_read_filter_span(record_type, parsed.read_name_extra);
    return parsed;
}

struct FileHeaderFields {
    uint8_t version = 0;
    uint64_t presence_flags = 0;
    uint64_t compression_level = 0;
    uint64_t block_size = 0;

    bool is_paired() const {
        return (presence_flags & PRESENCE_PAIRED) != 0;
    }
    bool has_qualities() const {
        return (presence_flags & PRESENCE_QUALITIES) != 0;
    }
    bool has_headers() const {
        return (presence_flags & PRESENCE_HEADERS) != 0;
    }
    bool has_flags() const {
        return (presence_flags & PRESENCE_FLAGS) != 0;
    }
};

bool parse_file_header(const std::array<uint8_t, 64>& bytes,
                       FileHeaderFields* header,
                       std::string* error) {
    if (std::memcmp(bytes.data(), "CBQFILE", 7) != 0) {
        return set_error(error, "CBQ file header magic is not CBQFILE");
    }
    header->version = bytes[7];
    if (header->version != 1) {
        std::ostringstream msg;
        msg << "unsupported CBQ file version: " << static_cast<unsigned>(header->version);
        return set_error(error, msg.str());
    }
    header->presence_flags = read_le64(bytes.data() + 8);
    header->compression_level = read_le64(bytes.data() + 16);
    header->block_size = read_le64(bytes.data() + 24);
    return true;
}

struct BlockHeaderFields {
    uint8_t version = 0;
    uint64_t len_z_seq_len = 0;
    uint64_t len_z_header_len = 0;
    uint64_t len_z_npos = 0;
    uint64_t len_z_seq = 0;
    uint64_t len_z_flags = 0;
    uint64_t len_z_headers = 0;
    uint64_t len_z_qual = 0;
    uint64_t nuclen = 0;
    uint64_t len_nef = 0;
    uint64_t num_records = 0;
    uint64_t num_sequences = 0;
};

bool parse_block_header(const std::array<uint8_t, 96>& bytes,
                        BlockHeaderFields* header,
                        std::string* error) {
    if (std::memcmp(bytes.data(), "BLK", 3) != 0) {
        return set_error(error, "CBQ block header magic is not BLK");
    }
    header->version = bytes[3];
    if (header->version != 1) {
        std::ostringstream msg;
        msg << "unsupported CBQ block version: " << static_cast<unsigned>(header->version);
        return set_error(error, msg.str());
    }
    header->len_z_seq_len = read_le64(bytes.data() + 8);
    header->len_z_header_len = read_le64(bytes.data() + 16);
    header->len_z_npos = read_le64(bytes.data() + 24);
    header->len_z_seq = read_le64(bytes.data() + 32);
    header->len_z_flags = read_le64(bytes.data() + 40);
    header->len_z_headers = read_le64(bytes.data() + 48);
    header->len_z_qual = read_le64(bytes.data() + 56);
    header->nuclen = read_le64(bytes.data() + 64);
    header->len_nef = read_le64(bytes.data() + 72);
    header->num_records = read_le64(bytes.data() + 80);
    header->num_sequences = read_le64(bytes.data() + 88);
    return true;
}

class ZstdRuntime {
public:
    ZstdRuntime() : handle_(nullptr), decompress_(nullptr), is_error_(nullptr), get_error_name_(nullptr) {}

    ~ZstdRuntime() {
        if (handle_ != nullptr) {
            dlclose(handle_);
        }
    }

    bool load(std::string* error) {
        std::lock_guard<std::mutex> lock(load_mutex_);
        if (decompress_ != nullptr) {
            return true;
        }
        handle_ = dlopen("libzstd.so.1", RTLD_LAZY | RTLD_LOCAL);
        if (handle_ == nullptr) {
            handle_ = dlopen("libzstd.so", RTLD_LAZY | RTLD_LOCAL);
        }
        if (handle_ == nullptr) {
            return set_error(error, "could not load libzstd.so.1 or libzstd.so");
        }

        decompress_ = reinterpret_cast<DecompressFn>(dlsym(handle_, "ZSTD_decompress"));
        is_error_ = reinterpret_cast<IsErrorFn>(dlsym(handle_, "ZSTD_isError"));
        get_error_name_ = reinterpret_cast<GetErrorNameFn>(dlsym(handle_, "ZSTD_getErrorName"));
        if (decompress_ == nullptr || is_error_ == nullptr || get_error_name_ == nullptr) {
            return set_error(error, "libzstd is missing required ZSTD_decompress symbols");
        }
        return true;
    }

    bool decompress(const std::vector<uint8_t>& compressed,
                    size_t expected_size,
                    std::vector<uint8_t>* out,
                    const std::string& column_name,
                    std::string* error) {
        if (!load(error)) {
            return false;
        }
        out->assign(expected_size, 0);
        if (expected_size == 0 && compressed.empty()) {
            return true;
        }
        const size_t nbytes = decompress_(
            out->data(),
            out->size(),
            compressed.data(),
            compressed.size());
        if (is_error_(nbytes)) {
            std::ostringstream msg;
            msg << "zstd decompression failed for CBQ column " << column_name
                << ": " << get_error_name_(nbytes);
            return set_error(error, msg.str());
        }
        if (nbytes != expected_size) {
            std::ostringstream msg;
            msg << "zstd decompressed CBQ column " << column_name << " to "
                << nbytes << " bytes, expected " << expected_size;
            return set_error(error, msg.str());
        }
        return true;
    }

private:
    using DecompressFn = size_t (*)(void*, size_t, const void*, size_t);
    using IsErrorFn = unsigned (*)(size_t);
    using GetErrorNameFn = const char* (*)(size_t);

    void* handle_;
    DecompressFn decompress_;
    IsErrorFn is_error_;
    GetErrorNameFn get_error_name_;
    std::mutex load_mutex_;
};

ZstdRuntime& zstd_runtime() {
    static ZstdRuntime runtime;
    return runtime;
}

struct BitVectorLite {
    std::vector<uint64_t> words;
    uint64_t len = 0;

    bool bit(uint64_t pos) const {
        if (pos >= len) {
            return false;
        }
        return ((words[static_cast<size_t>(pos / 64)] >> (pos % 64)) & 1ULL) != 0;
    }

    uint64_t bits(uint64_t pos, uint64_t nbits) const {
        uint64_t value = 0;
        for (uint64_t ii = 0; ii < nbits; ++ii) {
            if (bit(pos + ii)) {
                value |= 1ULL << ii;
            }
        }
        return value;
    }
};

class SliceReader {
public:
    explicit SliceReader(const std::vector<uint8_t>& bytes) : bytes_(bytes), offset_(0) {}

    bool read_u8(uint8_t* value, std::string* error) {
        if (offset_ + 1 > bytes_.size()) {
            return set_error(error, "truncated sucds payload while reading u8");
        }
        *value = bytes_[offset_++];
        return true;
    }

    bool read_bool(bool* value, std::string* error) {
        uint8_t byte = 0;
        if (!read_u8(&byte, error)) {
            return false;
        }
        *value = byte != 0;
        return true;
    }

    bool read_u16(uint16_t* value, std::string* error) {
        if (offset_ + 2 > bytes_.size()) {
            return set_error(error, "truncated sucds payload while reading u16");
        }
        *value = read_le16(bytes_.data() + offset_);
        offset_ += 2;
        return true;
    }

    bool read_u64(uint64_t* value, std::string* error) {
        if (offset_ + 8 > bytes_.size()) {
            return set_error(error, "truncated sucds payload while reading u64");
        }
        *value = read_le64(bytes_.data() + offset_);
        offset_ += 8;
        return true;
    }

    bool skip_bytes(uint64_t nbytes, std::string* error) {
        size_t n = 0;
        if (!checked_size(nbytes, &n) || offset_ + n > bytes_.size()) {
            return set_error(error, "truncated sucds payload while skipping bytes");
        }
        offset_ += n;
        return true;
    }

    bool read_bitvector(BitVectorLite* bv, std::string* error) {
        if (!read_vec_u64(&bv->words, error)) {
            return false;
        }
        if (!read_u64(&bv->len, error)) {
            return false;
        }
        const uint64_t expected_words = (bv->len + 63) / 64;
        if (expected_words > bv->words.size()) {
            return set_error(error, "sucds BitVector has fewer words than its bit length requires");
        }
        return true;
    }

    bool read_darray_bitvector(BitVectorLite* high_bits, std::string* error) {
        if (!read_bitvector(high_bits, error)) {
            return false;
        }
        if (!skip_darray_index(error)) {
            return false;
        }
        if (!skip_optional_darray_index(error)) {
            return false;
        }
        return skip_optional_rank9sel_index(error);
    }

    bool at_end() const {
        return offset_ == bytes_.size();
    }

private:
    bool read_vec_u64(std::vector<uint64_t>* values, std::string* error) {
        uint64_t len64 = 0;
        if (!read_u64(&len64, error)) {
            return false;
        }
        size_t len = 0;
        if (!checked_size(len64, &len)) {
            return set_error(error, "sucds Vec length exceeds platform size");
        }
        values->clear();
        values->reserve(len);
        for (size_t ii = 0; ii < len; ++ii) {
            uint64_t value = 0;
            if (!read_u64(&value, error)) {
                return false;
            }
            values->push_back(value);
        }
        return true;
    }

    bool skip_vec_i64(std::string* error) {
        uint64_t len = 0;
        if (!read_u64(&len, error)) {
            return false;
        }
        return skip_bytes(len * 8, error);
    }

    bool skip_vec_u16(std::string* error) {
        uint64_t len = 0;
        if (!read_u64(&len, error)) {
            return false;
        }
        return skip_bytes(len * 2, error);
    }

    bool skip_vec_u64(std::string* error) {
        uint64_t len = 0;
        if (!read_u64(&len, error)) {
            return false;
        }
        return skip_bytes(len * 8, error);
    }

    bool skip_optional_vec_u64(std::string* error) {
        bool present = false;
        if (!read_bool(&present, error)) {
            return false;
        }
        return present ? skip_vec_u64(error) : true;
    }

    bool skip_darray_index(std::string* error) {
        uint64_t unused = 0;
        bool unused_bool = false;
        return skip_vec_i64(error) &&
               skip_vec_u16(error) &&
               skip_vec_u64(error) &&
               read_u64(&unused, error) &&
               read_bool(&unused_bool, error);
    }

    bool skip_optional_darray_index(std::string* error) {
        bool present = false;
        if (!read_bool(&present, error)) {
            return false;
        }
        return present ? skip_darray_index(error) : true;
    }

    bool skip_rank9sel_index(std::string* error) {
        uint64_t unused = 0;
        return read_u64(&unused, error) &&
               skip_vec_u64(error) &&
               skip_optional_vec_u64(error) &&
               skip_optional_vec_u64(error);
    }

    bool skip_optional_rank9sel_index(std::string* error) {
        bool present = false;
        if (!read_bool(&present, error)) {
            return false;
        }
        return present ? skip_rank9sel_index(error) : true;
    }

    const std::vector<uint8_t>& bytes_;
    size_t offset_;
};

bool decode_elias_fano_positions(const std::vector<uint8_t>& bytes,
                                 std::vector<uint64_t>* positions,
                                 std::string* error) {
    positions->clear();
    if (bytes.empty()) {
        return true;
    }

    SliceReader reader(bytes);
    BitVectorLite high_bits;
    BitVectorLite low_bits;
    uint64_t low_len = 0;
    uint64_t universe = 0;

    if (!reader.read_darray_bitvector(&high_bits, error) ||
        !reader.read_bitvector(&low_bits, error) ||
        !reader.read_u64(&low_len, error) ||
        !reader.read_u64(&universe, error)) {
        return false;
    }
    if (low_len > 63) {
        return set_error(error, "unsupported sucds Elias-Fano low-bit length > 63");
    }

    uint64_t ordinal = 0;
    for (uint64_t pos = 0; pos < high_bits.len; ++pos) {
        if (!high_bits.bit(pos)) {
            continue;
        }
        if (low_len != 0 && ordinal * low_len + low_len > low_bits.len) {
            return set_error(error, "sucds Elias-Fano low-bit vector is too short");
        }
        const uint64_t high = pos - ordinal;
        const uint64_t low = low_bits.bits(ordinal * low_len, low_len);
        const uint64_t value = (high << low_len) | low;
        if (value >= universe) {
            return set_error(error, "decoded Elias-Fano value exceeds universe");
        }
        positions->push_back(value);
        ++ordinal;
    }

    if (!reader.at_end()) {
        return set_error(error, "trailing bytes after sucds Elias-Fano payload");
    }
    return true;
}

bool read_exact(std::ifstream& stream,
                uint8_t* dest,
                size_t nbytes,
                const std::string& description,
                std::string* error) {
    stream.read(reinterpret_cast<char*>(dest), static_cast<std::streamsize>(nbytes));
    if (stream.gcount() != static_cast<std::streamsize>(nbytes)) {
        return set_error(error, "truncated CBQ stream while reading " + description);
    }
    return true;
}

bool read_vector(std::ifstream& stream,
                 uint64_t nbytes64,
                 std::vector<uint8_t>* dest,
                 const std::string& description,
                 std::string* error) {
    size_t nbytes = 0;
    if (!checked_size(nbytes64, &nbytes)) {
        return set_error(error, "CBQ column length exceeds platform size for " + description);
    }
    dest->assign(nbytes, 0);
    if (nbytes == 0) {
        return true;
    }
    return read_exact(stream, dest->data(), nbytes, description, error);
}

bool bytes_to_u64_vector(const std::vector<uint8_t>& bytes,
                         std::vector<uint64_t>* values,
                         const std::string& description,
                         std::string* error) {
    if (bytes.size() % 8 != 0) {
        return set_error(error, "CBQ column " + description + " byte length is not a multiple of 8");
    }
    values->clear();
    values->reserve(bytes.size() / 8);
    for (size_t offset = 0; offset < bytes.size(); offset += 8) {
        values->push_back(read_le64(bytes.data() + offset));
    }
    return true;
}

struct DecodedBlock {
    uint64_t num_records = 0;
    uint64_t num_sequences = 0;
    std::vector<uint64_t> seq_lengths;
    std::vector<uint64_t> seq_offsets;
    std::vector<uint64_t> header_lengths;
    std::vector<uint64_t> header_offsets;
    std::vector<uint8_t> seq_word_bytes;
    std::vector<uint64_t> n_positions;
    std::string sequences;
    std::string headers;
    std::string qualities;
    size_t record_index = 0;

    void clear() {
        num_records = 0;
        num_sequences = 0;
        seq_lengths.clear();
        seq_offsets.clear();
        header_lengths.clear();
        header_offsets.clear();
        seq_word_bytes.clear();
        n_positions.clear();
        sequences.clear();
        headers.clear();
        qualities.clear();
        record_index = 0;
    }
};

bool make_offsets(const std::vector<uint64_t>& lengths,
                  uint64_t expected_total,
                  std::vector<uint64_t>* offsets,
                  const std::string& description,
                  std::string* error) {
    offsets->assign(lengths.size() + 1, 0);
    uint64_t total = 0;
    for (size_t ii = 0; ii < lengths.size(); ++ii) {
        total += lengths[ii];
        (*offsets)[ii + 1] = total;
    }
    if (total != expected_total) {
        std::ostringstream msg;
        msg << "CBQ " << description << " lengths sum to " << total
            << ", expected " << expected_total;
        return set_error(error, msg.str());
    }
    return true;
}

enum class BlockLoadStatus {
    Block,
    End,
    Error
};

} // namespace

struct CbqBatchBacking {
    DecodedBlock block;
    std::vector<std::string> synthetic_read_names;
    std::vector<CbqSegmentView> segment_views;
    std::vector<CbqReadView> read_views;
};

size_t cbq_segment_sequence_length(const CbqSegmentView& segment) {
    if (segment.packed_sequence.available) {
        return segment.packed_sequence.length;
    }
    return segment.sequence.size;
}

namespace {

bool cbq_segment_base_is_n(const CbqPackedSequenceView& packed, uint64_t global_offset) {
    if (packed.n_positions == nullptr || packed.n_positions_count == 0) {
        return false;
    }
    const uint64_t* begin = packed.n_positions;
    const uint64_t* end = begin + packed.n_positions_count;
    return std::binary_search(begin, end, global_offset);
}

unsigned char cbq_ascii_base_to_star_number(char base) {
    switch (base) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return 4;
    }
}

unsigned char cbq_packed_base_number(const CbqPackedSequenceView& packed, size_t index) {
    const uint64_t global_offset = packed.base_offset + static_cast<uint64_t>(index);
    if (cbq_segment_base_is_n(packed, global_offset)) {
        return 4;
    }
    const uint64_t word_index = global_offset / 32;
    const uint64_t offset_in_word = global_offset % 32;
    const size_t byte_offset = static_cast<size_t>(word_index * 8);
    if (packed.words == nullptr || byte_offset + 8 > packed.word_bytes) {
        return 4;
    }
    const uint64_t word = read_le64(packed.words + byte_offset);
    return static_cast<unsigned char>((word >> (offset_in_word * 2)) & 0x3ULL);
}

char cbq_number_to_ascii(unsigned char base) {
    static const char lookup[5] = {'A', 'C', 'G', 'T', 'N'};
    return lookup[base < 5 ? base : 4];
}

struct CbqPackedByteDecode {
    char ascii[4];
};

const std::array<CbqPackedByteDecode, 256>& cbq_packed_ascii_lookup() {
    static const std::array<CbqPackedByteDecode, 256> lookup = [] {
        std::array<CbqPackedByteDecode, 256> table = {};
        static const char ascii_lookup[4] = {'A', 'C', 'G', 'T'};
        for (unsigned byte = 0; byte < table.size(); ++byte) {
            for (unsigned base = 0; base < 4; ++base) {
                const unsigned code = (byte >> (base * 2)) & 0x3U;
                table[byte].ascii[base] = ascii_lookup[code];
            }
        }
        return table;
    }();
    return lookup;
}

} // namespace

char cbq_segment_base_ascii(const CbqSegmentView& segment, size_t index) {
    if (segment.sequence.data != nullptr && index < segment.sequence.size) {
        return segment.sequence.data[index];
    }
    if (segment.packed_sequence.available && index < segment.packed_sequence.length) {
        return cbq_number_to_ascii(cbq_packed_base_number(segment.packed_sequence, index));
    }
    return 'N';
}

unsigned char cbq_segment_base_star_number(const CbqSegmentView& segment, size_t index) {
    if (segment.packed_sequence.available && index < segment.packed_sequence.length) {
        return cbq_packed_base_number(segment.packed_sequence, index);
    }
    if (segment.sequence.data != nullptr && index < segment.sequence.size) {
        return cbq_ascii_base_to_star_number(segment.sequence.data[index]);
    }
    return 4;
}

void materialize_cbq_segment_sequence(const CbqSegmentView& segment, std::string* sequence) {
    if (sequence == nullptr) {
        return;
    }
    const size_t length = cbq_segment_sequence_length(segment);
    sequence->resize(length);
    for (size_t ii = 0; ii < length; ++ii) {
        (*sequence)[ii] = cbq_segment_base_ascii(segment, ii);
    }
}

bool materialize_cbq_segment_sequence_to_buffer(const CbqSegmentView& segment,
                                                char* dest,
                                                size_t capacity,
                                                size_t* length_out,
                                                std::string* error) {
    if (dest == nullptr) {
        return set_error(error, "CBQ sequence materialization destination is null");
    }
    const size_t length = cbq_segment_sequence_length(segment);
    if (length >= capacity) {
        std::ostringstream msg;
        msg << "CBQ sequence length " << length
            << " exceeds destination capacity " << capacity;
        return set_error(error, msg.str());
    }

    if (segment.sequence.data != nullptr && segment.sequence.size == length) {
        std::memcpy(dest, segment.sequence.data, length);
        dest[length] = '\0';
        if (length_out != nullptr) {
            *length_out = length;
        }
        return true;
    }

    if (!segment.packed_sequence.available) {
        std::fill(dest, dest + length, 'N');
        dest[length] = '\0';
        if (length_out != nullptr) {
            *length_out = length;
        }
        return true;
    }

    const CbqPackedSequenceView& packed = segment.packed_sequence;
    const auto& lookup = cbq_packed_ascii_lookup();
    size_t ii = 0;
    size_t byte_offset = static_cast<size_t>(packed.base_offset >> 2);
    const unsigned first_base = static_cast<unsigned>(packed.base_offset & 0x3ULL);

    if (first_base != 0 && ii < length) {
        const size_t take = std::min<size_t>(4U - first_base, length);
        if (packed.words == nullptr || byte_offset >= packed.word_bytes) {
            std::fill(dest, dest + take, 'N');
        } else {
            const CbqPackedByteDecode& decoded = lookup[packed.words[byte_offset]];
            std::memcpy(dest, decoded.ascii + first_base, take);
        }
        ii = take;
        ++byte_offset;
    }

    while (length - ii >= 4) {
        if (packed.words == nullptr || byte_offset >= packed.word_bytes) {
            std::memset(dest + ii, 'N', 4);
        } else {
            const CbqPackedByteDecode& decoded = lookup[packed.words[byte_offset]];
            std::memcpy(dest + ii, decoded.ascii, 4);
        }
        ii += 4;
        ++byte_offset;
    }

    if (ii < length) {
        const size_t take = length - ii;
        if (packed.words == nullptr || byte_offset >= packed.word_bytes) {
            std::fill(dest + ii, dest + length, 'N');
        } else {
            const CbqPackedByteDecode& decoded = lookup[packed.words[byte_offset]];
            std::memcpy(dest + ii, decoded.ascii, take);
        }
    }

    if (packed.n_positions != nullptr) {
        for (size_t inpos = 0; inpos < packed.n_positions_count; ++inpos) {
            const uint64_t npos = packed.n_positions[inpos];
            if (npos < packed.base_offset) {
                continue;
            }
            const uint64_t local = npos - packed.base_offset;
            if (local >= length) {
                break;
            }
            dest[local] = 'N';
        }
    }

    dest[length] = '\0';
    if (length_out != nullptr) {
        *length_out = length;
    }
    return true;
}

struct CbqInputModule::Impl {
    struct BlockIndexEntry {
        uint64_t offset = 0;
        uint64_t cumulative_records = 0;
    };

    std::ifstream stream;
    uint32_t lane_index = 0;
    uint64_t read_ordinal = 0;
    uint64_t lane_record_index = 0;
    bool opened = false;
    FileHeaderFields file_header;
    std::shared_ptr<CbqBatchBacking> batch;
    std::vector<BlockIndexEntry> block_index;
    uint64_t current_lane_records = 0;
    uint64_t batch_first_record = 0;
    bool range_mode = false;
    uint64_t range_first_record = 0;
    uint64_t range_end_record = 0;

    void clear_batch_views() {
        if (!batch) {
            return;
        }
        batch->synthetic_read_names.clear();
        batch->segment_views.clear();
        batch->read_views.clear();
    }

    bool build_batch_views(const InputSourcePlan& plan, std::string* error) {
        if (!batch) {
            return set_error(error, "CBQ batch storage is not initialized");
        }
        DecodedBlock& block = batch->block;
        std::vector<std::string>& synthetic_read_names = batch->synthetic_read_names;
        std::vector<CbqSegmentView>& segment_views = batch->segment_views;
        std::vector<CbqReadView>& read_views = batch->read_views;

        clear_batch_views();

        size_t num_records = 0;
        if (!checked_size(block.num_records, &num_records)) {
            return set_error(error, "CBQ block record count exceeds platform size");
        }
        if (num_records > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
            return set_error(error, "CBQ block record count exceeds uint32_t");
        }
        if (plan.mate_count == 0 ||
            num_records > std::numeric_limits<size_t>::max() / plan.mate_count) {
            return set_error(error, "CBQ block segment count exceeds platform size");
        }

        const size_t segment_count = num_records * plan.mate_count;
        read_views.assign(num_records, CbqReadView());
        segment_views.assign(segment_count, CbqSegmentView());
        synthetic_read_names.reserve((file_header.has_headers() && !block.headers.empty()) ? 0 : num_records);

        for (size_t irecord = 0; irecord < num_records; ++irecord) {
            const uint64_t first_sequence_index = static_cast<uint64_t>(irecord) * plan.mate_count;

            ParsedHeaderSpan parsed;
            if (file_header.has_headers() && !block.headers.empty()) {
                const size_t begin = static_cast<size_t>(block.header_offsets[static_cast<size_t>(first_sequence_index)]);
                const size_t end = static_cast<size_t>(block.header_offsets[static_cast<size_t>(first_sequence_index + 1)]);
                parsed = parse_header_payload_span(
                    block.headers,
                    begin,
                    end - begin,
                    file_header.has_qualities(),
                    plan.read_name_separator_chars);
            } else {
                synthetic_read_names.push_back(std::to_string(lane_record_index + irecord + 1));
                parsed.read_name = make_span(synthetic_read_names.back(), 0, synthetic_read_names.back().size());
                parsed.read_name_extra = make_empty_span();
                parsed.read_filter = 'N';
            }

            CbqReadView& record = read_views[irecord];
            record.read_name = parsed.read_name;
            record.read_name_extra = parsed.read_name_extra;
            record.lane_index = lane_index;
            record.read_ordinal = read_ordinal + irecord + 1;
            record.lane_read_ordinal = lane_record_index + irecord + 1;
            record.read_filter = parsed.read_filter;
            record.segment_count = plan.mate_count;
            record.segments = segment_views.data() + (irecord * plan.mate_count);

            for (uint32_t isegment = 0; isegment < plan.mate_count; ++isegment) {
                const uint64_t seq_index = first_sequence_index + isegment;
                const size_t begin = static_cast<size_t>(block.seq_offsets[static_cast<size_t>(seq_index)]);
                const size_t end = static_cast<size_t>(block.seq_offsets[static_cast<size_t>(seq_index + 1)]);
                const size_t length = end - begin;
                if (length > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
                    return set_error(error, "CBQ sequence length exceeds uint32_t");
                }

                CbqSegmentView& segment = segment_views[irecord * plan.mate_count + isegment];
                segment.role = role_for_segment(isegment);
                segment.source_index = isegment;
                segment.sequence = make_empty_span();
                segment.quality = make_span(block.qualities, begin, length);
                segment.packed_sequence.words = block.seq_word_bytes.empty() ? nullptr : block.seq_word_bytes.data();
                segment.packed_sequence.word_bytes = block.seq_word_bytes.size();
                segment.packed_sequence.base_offset = begin;
                segment.packed_sequence.length = static_cast<uint32_t>(length);
                const auto n_begin = std::lower_bound(block.n_positions.begin(), block.n_positions.end(), begin);
                const auto n_end = std::lower_bound(block.n_positions.begin(), block.n_positions.end(), end);
                segment.packed_sequence.n_positions_count =
                    static_cast<size_t>(std::distance(n_begin, n_end));
                segment.packed_sequence.n_positions = segment.packed_sequence.n_positions_count == 0
                    ? nullptr
                    : &(*n_begin);
                segment.packed_sequence.available = !block.seq_word_bytes.empty();
                segment.original_length = static_cast<uint32_t>(length);
                segment.has_quality = file_header.has_qualities();
            }
        }

        read_ordinal += static_cast<uint64_t>(num_records);
        lane_record_index += static_cast<uint64_t>(num_records);
        return true;
    }

    bool read_current_lane_index(const InputSourcePlan& plan, std::string* error) {
        const std::string& path = plan.mate_files.front().at(lane_index);
        stream.clear();
        stream.seekg(0, std::ios::end);
        const std::streampos end_pos = stream.tellg();
        if (end_pos < 0) {
            return set_error(error, "could not seek CBQ file for index: " + path);
        }
        const uint64_t file_size = static_cast<uint64_t>(end_pos);
        const uint64_t index_header_size = 24U;
        const uint64_t index_footer_size = 16U;
        if (file_size < 64U + index_header_size + index_footer_size) {
            return set_error(error, "CBQ range mode requires a CBQINDEX footer: " + path);
        }

        std::array<uint8_t, 16> footer{};
        stream.seekg(static_cast<std::streamoff>(file_size - index_footer_size), std::ios::beg);
        if (!read_exact(stream, footer.data(), footer.size(), "index footer", error)) {
            return false;
        }
        if (std::memcmp(footer.data() + 8, "CBQINDEX", 8) != 0) {
            return set_error(error, "CBQ range mode requires a CBQINDEX footer: " + path);
        }

        const uint64_t z_index_size = read_le64(footer.data());
        if (file_size < 64U + index_header_size + index_footer_size + z_index_size) {
            return set_error(error, "CBQ index footer points before file header: " + path);
        }
        const uint64_t index_header_offset =
            file_size - index_footer_size - z_index_size - index_header_size;

        std::array<uint8_t, 24> index_header{};
        stream.seekg(static_cast<std::streamoff>(index_header_offset), std::ios::beg);
        if (!read_exact(stream, index_header.data(), index_header.size(), "index header", error)) {
            return false;
        }
        if (std::memcmp(index_header.data(), "CBQINDEX", 8) != 0) {
            return set_error(error, "CBQ index header magic mismatch: " + path);
        }

        const uint64_t index_size = read_le64(index_header.data() + 8);
        const uint64_t compressed_index_size = read_le64(index_header.data() + 16);
        if (compressed_index_size != z_index_size) {
            return set_error(error, "CBQ index footer/header compressed-size mismatch: " + path);
        }
        if (index_size % 16U != 0) {
            return set_error(error, "CBQ index size is not a multiple of 16 bytes: " + path);
        }
        size_t index_size_size_t = 0;
        if (!checked_size(index_size, &index_size_size_t)) {
            return set_error(error, "CBQ index size exceeds platform size: " + path);
        }

        std::vector<uint8_t> z_index;
        if (!read_vector(stream, compressed_index_size, &z_index, "z_index", error)) {
            return false;
        }

        std::vector<uint8_t> index_bytes;
        if (!zstd_runtime().decompress(z_index, index_size_size_t, &index_bytes, "index", error)) {
            return false;
        }
        if (index_bytes.size() != index_size_size_t) {
            return set_error(error, "CBQ index decompressed to an unexpected size: " + path);
        }

        block_index.clear();
        block_index.reserve(index_bytes.size() / 16U);
        uint64_t previous_cumulative_records = 0;
        for (size_t offset = 0; offset < index_bytes.size(); offset += 16U) {
            BlockIndexEntry entry;
            entry.offset = read_le64(index_bytes.data() + offset);
            entry.cumulative_records = read_le64(index_bytes.data() + offset + 8U);
            if (entry.offset < 64U || entry.offset >= index_header_offset) {
                return set_error(error, "CBQ index contains an invalid block offset: " + path);
            }
            if (entry.cumulative_records < previous_cumulative_records) {
                return set_error(error, "CBQ index cumulative records are not monotonic: " + path);
            }
            previous_cumulative_records = entry.cumulative_records;
            block_index.push_back(entry);
        }

        current_lane_records = block_index.empty() ? 0U : block_index.back().cumulative_records;
        stream.clear();
        return true;
    }

    bool seek_to_range_start(uint64_t first_record, std::string* error) {
        if (first_record >= current_lane_records || block_index.empty()) {
            lane_record_index = current_lane_records;
            read_ordinal = current_lane_records;
            batch.reset();
            return true;
        }

        std::vector<BlockIndexEntry>::const_iterator it =
            std::upper_bound(block_index.begin(),
                             block_index.end(),
                             first_record,
                             [](uint64_t value, const BlockIndexEntry& entry) {
                                 return value < entry.cumulative_records;
                             });
        if (it == block_index.end()) {
            lane_record_index = current_lane_records;
            read_ordinal = current_lane_records;
            batch.reset();
            return true;
        }

        std::vector<BlockIndexEntry>::const_iterator index_begin = block_index.begin();
        const size_t block_index_position =
            static_cast<size_t>(std::distance(index_begin, it));
        const uint64_t block_first_record = block_index_position == 0
            ? 0U
            : block_index[block_index_position - 1U].cumulative_records;

        stream.clear();
        stream.seekg(static_cast<std::streamoff>(it->offset), std::ios::beg);
        if (!stream.good()) {
            return set_error(error, "could not seek CBQ stream to indexed block");
        }
        lane_record_index = block_first_record;
        read_ordinal = block_first_record;
        batch.reset();
        return true;
    }

    bool open_lane(const InputSourcePlan& plan, uint32_t lane, std::string* error) {
        close_lane();
        const std::string& path = plan.mate_files.front().at(lane);
        stream.open(path.c_str(), std::ios::binary);
        if (!stream.good()) {
            return set_error(error, "could not open CBQ file: " + path);
        }

        std::array<uint8_t, 64> header_bytes{};
        if (!read_exact(stream, header_bytes.data(), header_bytes.size(), "file header", error)) {
            close_lane();
            return false;
        }
        if (!parse_file_header(header_bytes, &file_header, error)) {
            close_lane();
            return false;
        }

        const bool file_paired = file_header.is_paired();
        if ((plan.mate_count == 2) != file_paired) {
            std::ostringstream msg;
            msg << "CBQ mate-count mismatch for " << path << ": file paired="
                << (file_paired ? "true" : "false")
                << " but plan mate_count=" << plan.mate_count;
            close_lane();
            return set_error(error, msg.str());
        }

        lane_index = lane;
        lane_record_index = 0;
        batch_first_record = 0;
        batch.reset();
        opened = true;
        return true;
    }

    void close_lane() {
        if (stream.is_open()) {
            stream.close();
        }
        batch.reset();
        block_index.clear();
        current_lane_records = 0;
        batch_first_record = 0;
        range_mode = false;
        range_first_record = 0;
        range_end_record = 0;
    }

    BlockLoadStatus load_next_block(const InputSourcePlan& plan, std::string* error) {
        batch.reset();

        std::array<uint8_t, 8> first{};
        stream.read(reinterpret_cast<char*>(first.data()), static_cast<std::streamsize>(first.size()));
        if (stream.gcount() == 0 && stream.eof()) {
            return BlockLoadStatus::End;
        }
        if (stream.gcount() != static_cast<std::streamsize>(first.size())) {
            set_error(error, "truncated CBQ stream while reading block/index magic");
            return BlockLoadStatus::Error;
        }
        if (std::memcmp(first.data(), "CBQINDEX", 8) == 0) {
            return BlockLoadStatus::End;
        }

        std::array<uint8_t, 96> header_bytes{};
        std::copy(first.begin(), first.end(), header_bytes.begin());
        if (!read_exact(stream,
                        header_bytes.data() + first.size(),
                        header_bytes.size() - first.size(),
                        "block header",
                        error)) {
            return BlockLoadStatus::Error;
        }

        BlockHeaderFields header;
        if (!parse_block_header(header_bytes, &header, error)) {
            return BlockLoadStatus::Error;
        }

        const uint64_t expected_sequences = header.num_records * plan.mate_count;
        if (header.num_sequences != expected_sequences) {
            std::ostringstream msg;
            msg << "CBQ block has num_sequences=" << header.num_sequences
                << " but expected " << expected_sequences
                << " for num_records=" << header.num_records
                << " mate_count=" << plan.mate_count;
            set_error(error, msg.str());
            return BlockLoadStatus::Error;
        }

        batch = std::make_shared<CbqBatchBacking>();
        DecodedBlock& block = batch->block;

        std::vector<uint8_t> z_seq_len;
        std::vector<uint8_t> z_header_len;
        std::vector<uint8_t> z_npos;
        std::vector<uint8_t> z_seq;
        std::vector<uint8_t> z_flags;
        std::vector<uint8_t> z_headers;
        std::vector<uint8_t> z_qual;

        if (!read_vector(stream, header.len_z_seq_len, &z_seq_len, "z_seq_len", error) ||
            !read_vector(stream, header.len_z_header_len, &z_header_len, "z_header_len", error) ||
            !read_vector(stream, header.len_z_npos, &z_npos, "z_npos", error) ||
            !read_vector(stream, header.len_z_seq, &z_seq, "z_seq", error) ||
            !read_vector(stream, header.len_z_flags, &z_flags, "z_flags", error) ||
            !read_vector(stream, header.len_z_headers, &z_headers, "z_headers", error) ||
            !read_vector(stream, header.len_z_qual, &z_qual, "z_qual", error)) {
            return BlockLoadStatus::Error;
        }

        size_t num_sequences = 0;
        size_t nuclen = 0;
        size_t ef_len = 0;
        if (!checked_size(header.num_sequences, &num_sequences) ||
            !checked_size(header.nuclen, &nuclen) ||
            !checked_size(header.len_nef, &ef_len)) {
            set_error(error, "CBQ block metadata exceeds platform size");
            return BlockLoadStatus::Error;
        }

        std::vector<uint8_t> seq_len_bytes;
        std::vector<uint8_t> header_len_bytes;
        std::vector<uint8_t> npos_bytes;
        std::vector<uint8_t> flags_bytes;
        std::vector<uint8_t> header_bytes_column;
        std::vector<uint8_t> qual_bytes;

        const size_t seq_len_bytes_expected = num_sequences * 8;
        const size_t seq_word_bytes_expected = ((nuclen + 31) / 32) * 8;
        if (!zstd_runtime().decompress(z_seq_len, seq_len_bytes_expected, &seq_len_bytes, "seq_len", error)) {
            return BlockLoadStatus::Error;
        }
        if (!bytes_to_u64_vector(seq_len_bytes, &block.seq_lengths, "seq_len", error)) {
            return BlockLoadStatus::Error;
        }
        if (!make_offsets(block.seq_lengths, header.nuclen, &block.seq_offsets, "sequence", error)) {
            return BlockLoadStatus::Error;
        }

        if (file_header.has_headers()) {
            if (!zstd_runtime().decompress(z_header_len, seq_len_bytes_expected, &header_len_bytes, "header_len", error) ||
                !bytes_to_u64_vector(header_len_bytes, &block.header_lengths, "header_len", error)) {
                return BlockLoadStatus::Error;
            }
        } else {
            block.header_lengths.assign(num_sequences, 0);
        }

        if (!z_npos.empty()) {
            if (!zstd_runtime().decompress(z_npos, ef_len, &npos_bytes, "npos", error)) {
                return BlockLoadStatus::Error;
            }
        }

        if (!decode_elias_fano_positions(npos_bytes, &block.n_positions, error)) {
            return BlockLoadStatus::Error;
        }

        if (!zstd_runtime().decompress(z_seq, seq_word_bytes_expected, &block.seq_word_bytes, "seq", error)) {
            return BlockLoadStatus::Error;
        }
        for (uint64_t npos : block.n_positions) {
            if (npos >= header.nuclen) {
                set_error(error, "CBQ N-position exceeds sequence length");
                return BlockLoadStatus::Error;
            }
        }

        if (file_header.has_flags()) {
            const size_t flags_expected = static_cast<size_t>(header.num_records) * 8;
            if (!zstd_runtime().decompress(z_flags, flags_expected, &flags_bytes, "flags", error)) {
                return BlockLoadStatus::Error;
            }
        }

        uint64_t total_header_len = 0;
        for (uint64_t length : block.header_lengths) {
            total_header_len += length;
        }
        if (!make_offsets(block.header_lengths, total_header_len, &block.header_offsets, "header", error)) {
            return BlockLoadStatus::Error;
        }
        if (file_header.has_headers()) {
            size_t total_header_size = 0;
            if (!checked_size(total_header_len, &total_header_size)) {
                set_error(error, "CBQ header payload exceeds platform size");
                return BlockLoadStatus::Error;
            }
            if (!zstd_runtime().decompress(z_headers, total_header_size, &header_bytes_column, "headers", error)) {
                return BlockLoadStatus::Error;
            }
            block.headers.assign(reinterpret_cast<const char*>(header_bytes_column.data()), header_bytes_column.size());
        }

        if (file_header.has_qualities()) {
            if (!zstd_runtime().decompress(z_qual, nuclen, &qual_bytes, "qual", error)) {
                return BlockLoadStatus::Error;
            }
            block.qualities.assign(reinterpret_cast<const char*>(qual_bytes.data()), qual_bytes.size());
        } else {
            block.qualities.assign(nuclen, 'A');
        }

        block.num_records = header.num_records;
        block.num_sequences = header.num_sequences;
        block.record_index = 0;
        batch_first_record = lane_record_index;
        if (!build_batch_views(plan, error)) {
            return BlockLoadStatus::Error;
        }
        return BlockLoadStatus::Block;
    }

    BlockLoadStatus load_next_available_block(const InputSourcePlan& plan, std::string* error) {
        while (lane_index < plan.read_files_n) {
            if (range_mode &&
                range_first_record >= range_end_record) {
                opened = false;
                return BlockLoadStatus::End;
            }
            if (batch &&
                batch->block.record_index < static_cast<size_t>(batch->block.num_records)) {
                return BlockLoadStatus::Block;
            }
            if (range_mode && lane_record_index >= range_end_record) {
                opened = false;
                return BlockLoadStatus::End;
            }

            const BlockLoadStatus status = load_next_block(plan, error);
            if (status == BlockLoadStatus::Error) {
                return BlockLoadStatus::Error;
            }
            if (status == BlockLoadStatus::Block) {
                if (!batch || batch->block.num_records == 0) {
                    continue;
                }
                return BlockLoadStatus::Block;
            }

            if (range_mode) {
                opened = false;
                return BlockLoadStatus::End;
            }
            close_lane();
            ++lane_index;
            if (lane_index >= plan.read_files_n) {
                opened = false;
                return BlockLoadStatus::End;
            }
            if (!open_lane(plan, lane_index, error)) {
                return BlockLoadStatus::Error;
            }
        }

        opened = false;
        return BlockLoadStatus::End;
    }
};

CbqInputModule::CbqInputModule() : configured_(false), impl_(new Impl()) {
    plan_.format = SourceFormat::Binseq;
    plan_.module_name = "Cbq";
}

CbqInputModule::~CbqInputModule() = default;

const char* CbqInputModule::name() const {
    return "CbqInputModule";
}

bool CbqInputModule::configure(const InputSourcePlan& plan, std::string* error) {
    if (plan.format != SourceFormat::Binseq) {
        return set_error(error, "CbqInputModule received a non-BINSEQ input source plan");
    }
    if (plan.mate_files.size() != 1) {
        return set_error(error, "CBQ input expects one CBQ path per lane, not split external mate files");
    }
    if (plan.mate_files.front().empty()) {
        return set_error(error, "CBQ input requires at least one CBQ input path");
    }
    if (plan.mate_count != 1 && plan.mate_count != 2) {
        return set_error(error, "CBQ input currently supports mate_count 1 or 2");
    }
    if (plan.read_files_n != 0 && plan.read_files_n != plan.mate_files.front().size()) {
        std::ostringstream msg;
        msg << "CBQ input read_files_n=" << plan.read_files_n
            << " does not match CBQ lane count=" << plan.mate_files.front().size();
        return set_error(error, msg.str());
    }
    if (!plan.read_groups.empty() && plan.read_groups.size() != plan.mate_files.front().size()) {
        std::ostringstream msg;
        msg << "CBQ input has " << plan.mate_files.front().size()
            << " lanes but " << plan.read_groups.size() << " read groups";
        return set_error(error, msg.str());
    }

    plan_ = plan;
    plan_.format = SourceFormat::Binseq;
    plan_.module_name = "Cbq";
    plan_.preserves_source_order = true;
    plan_.uses_helper_stream = false;
    plan_.read_files_n = static_cast<uint32_t>(plan.mate_files.front().size());
    if (plan_.read_name_separator_chars.empty()) {
        plan_.read_name_separator_chars.push_back(' ');
    }
    configured_ = true;
    return true;
}

const InputSourcePlan& CbqInputModule::plan() const {
    return plan_;
}

bool CbqInputModule::supports_record_iteration() const {
    return true;
}

bool CbqInputModule::open(std::string* error) {
    if (!configured_) {
        return set_error(error, "CbqInputModule must be configured before open()");
    }
    close();
    if (plan_.read_files_n == 0) {
        return set_error(error, "CBQ input source plan has no lanes");
    }

    impl_->lane_index = 0;
    impl_->read_ordinal = 0;
    impl_->opened = false;
    if (!impl_->open_lane(plan_, impl_->lane_index, error)) {
        return false;
    }
    impl_->opened = true;
    return true;
}

bool CbqInputModule::open_range(uint32_t lane_index,
                                uint64_t first_record,
                                uint64_t record_count,
                                std::string* error) {
    if (!configured_) {
        return set_error(error, "CbqInputModule must be configured before open_range()");
    }
    close();
    if (plan_.read_files_n == 0) {
        return set_error(error, "CBQ input source plan has no lanes");
    }
    if (lane_index >= plan_.read_files_n) {
        std::ostringstream msg;
        msg << "CBQ range lane " << lane_index
            << " is outside lane count " << plan_.read_files_n;
        return set_error(error, msg.str());
    }

    impl_->lane_index = lane_index;
    impl_->read_ordinal = 0;
    impl_->opened = false;
    if (!impl_->open_lane(plan_, lane_index, error)) {
        return false;
    }
    if (!impl_->read_current_lane_index(plan_, error)) {
        close();
        return false;
    }
    if (first_record > impl_->current_lane_records) {
        std::ostringstream msg;
        msg << "CBQ range starts at record " << first_record
            << " but lane has only " << impl_->current_lane_records << " records";
        close();
        return set_error(error, msg.str());
    }

    uint64_t end_record = impl_->current_lane_records;
    if (record_count != std::numeric_limits<uint64_t>::max()) {
        if (record_count <= impl_->current_lane_records - first_record) {
            end_record = first_record + record_count;
        }
    }

    impl_->range_mode = true;
    impl_->range_first_record = first_record;
    impl_->range_end_record = end_record;
    if (!impl_->seek_to_range_start(first_record, error)) {
        close();
        return false;
    }
    impl_->opened = true;
    return true;
}

uint64_t CbqInputModule::current_lane_record_count() const {
    return impl_ ? impl_->current_lane_records : 0U;
}

InputStatus CbqInputModule::next_batch(CbqReadBatchView* batch, std::string* error) {
    if (batch == nullptr) {
        if (error != nullptr) {
            *error = "next_batch() requires a non-null CbqReadBatchView";
        }
        return InputStatus::Error;
    }
    if (!impl_->opened) {
        if (error != nullptr) {
            *error = "CbqInputModule is not open";
        }
        return InputStatus::Error;
    }

    for (;;) {
        const BlockLoadStatus status = impl_->load_next_available_block(plan_, error);
        if (status == BlockLoadStatus::Error) {
            return InputStatus::Error;
        }
        if (status == BlockLoadStatus::End) {
            return InputStatus::End;
        }
        if (!impl_->batch) {
            if (error != nullptr) {
                *error = "CBQ batch storage is not available";
            }
            return InputStatus::Error;
        }

        DecodedBlock& block = impl_->batch->block;
        size_t start = block.record_index;
        size_t remaining = static_cast<size_t>(block.num_records) - start;
        if (impl_->range_mode) {
            const uint64_t block_first_record = impl_->batch_first_record;
            const uint64_t current_record = block_first_record + static_cast<uint64_t>(start);
            if (current_record < impl_->range_first_record) {
                const uint64_t local_start = impl_->range_first_record - block_first_record;
                if (local_start >= block.num_records) {
                    block.record_index = static_cast<size_t>(block.num_records);
                    continue;
                }
                start = static_cast<size_t>(local_start);
                remaining = static_cast<size_t>(block.num_records) - start;
            }

            const uint64_t absolute_start = block_first_record + static_cast<uint64_t>(start);
            if (absolute_start >= impl_->range_end_record) {
                block.record_index = static_cast<size_t>(block.num_records);
                continue;
            }
            const uint64_t range_remaining = impl_->range_end_record - absolute_start;
            if (range_remaining < static_cast<uint64_t>(remaining)) {
                remaining = static_cast<size_t>(range_remaining);
            }
        }

        if (remaining == 0) {
            block.record_index = static_cast<size_t>(block.num_records);
            continue;
        }
        if (remaining > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
            if (error != nullptr) {
                *error = "CBQ unread batch size exceeds uint32_t";
            }
            return InputStatus::Error;
        }

        batch->records = impl_->batch->read_views.data() + start;
        batch->record_count = static_cast<uint32_t>(remaining);
        batch->lane_index = impl_->lane_index;
        batch->preserves_source_order = true;
        batch->backing_storage_owned_by_reader = true;
        batch->backing = impl_->batch;

        block.record_index = static_cast<size_t>(block.num_records);
        return InputStatus::Record;
    }
}

InputStatus CbqInputModule::next_record(InputRecord* record, std::string* error) {
    if (record == nullptr) {
        if (error != nullptr) {
            *error = "next_record() requires a non-null InputRecord";
        }
        return InputStatus::Error;
    }
    if (!impl_->opened) {
        if (error != nullptr) {
            *error = "CbqInputModule is not open";
        }
        return InputStatus::Error;
    }

    const CbqReadView* selected = nullptr;
    for (;;) {
        const BlockLoadStatus status = impl_->load_next_available_block(plan_, error);
        if (status == BlockLoadStatus::Error) {
            return InputStatus::Error;
        }
        if (status == BlockLoadStatus::End) {
            return InputStatus::End;
        }
        if (!impl_->batch) {
            if (error != nullptr) {
                *error = "CBQ batch storage is not available";
            }
            return InputStatus::Error;
        }

        DecodedBlock& block = impl_->batch->block;
        if (impl_->range_mode) {
            const uint64_t block_first_record = impl_->batch_first_record;
            const uint64_t current_record =
                block_first_record + static_cast<uint64_t>(block.record_index);
            if (current_record < impl_->range_first_record) {
                const uint64_t local_start = impl_->range_first_record - block_first_record;
                if (local_start >= block.num_records) {
                    block.record_index = static_cast<size_t>(block.num_records);
                    continue;
                }
                block.record_index = static_cast<size_t>(local_start);
            }

            const uint64_t absolute_record =
                block_first_record + static_cast<uint64_t>(block.record_index);
            if (absolute_record >= impl_->range_end_record) {
                block.record_index = static_cast<size_t>(block.num_records);
                continue;
            }
        }

        selected = &impl_->batch->read_views[block.record_index++];
        break;
    }

    const CbqReadView& view = *selected;

    assign_from_span(&record->read_name, view.read_name);
    assign_from_span(&record->read_name_extra, view.read_name_extra);
    record->lane_index = view.lane_index;
    record->read_ordinal = view.read_ordinal;
    record->read_filter = view.read_filter;
    record->mate_count = view.segment_count;
    record->mates.clear();
    record->mates.reserve(view.segment_count);

    for (uint32_t isegment = 0; isegment < view.segment_count; ++isegment) {
        const CbqSegmentView& segment = view.segments[isegment];
        InputMateRecord mate;
        materialize_cbq_segment_sequence(segment, &mate.sequence);
        assign_from_span(&mate.quality, segment.quality);
        apply_quality_conversion(mate.quality, plan_.out_qs_conversion_add);
        mate.original_length = segment.original_length;
        mate.has_quality = segment.has_quality;
        record->mates.push_back(mate);
    }

    return InputStatus::Record;
}

void CbqInputModule::close() {
    if (impl_) {
        impl_->close_lane();
        impl_->lane_index = 0;
        impl_->opened = false;
    }
}

} // namespace input
} // namespace star

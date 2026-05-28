#include "input/CbqInputModule.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstring>
#include <dlfcn.h>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace star {
namespace input {
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

std::string first_field(const std::string& value) {
    std::istringstream in(value);
    std::string field;
    in >> field;
    return field;
}

std::string canonical_read_name(const std::string& raw_name, const std::vector<char>& separators) {
    std::string name = raw_name;
    for (char separator : separators) {
        const size_t pos = name.find(separator);
        if (pos != std::string::npos) {
            name.resize(pos);
        }
    }
    return name;
}

char parse_read_filter(char record_type, const std::string& read_name_extra) {
    if (record_type != '@') {
        return 'N';
    }

    const std::string field2 = first_field(read_name_extra);
    if (field2.size() > 3 && field2[1] == ':' && field2[2] == 'Y' && field2[3] == ':') {
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

struct ParsedHeader {
    std::string read_name;
    std::string read_name_extra;
    char read_filter = 'N';
};

ParsedHeader parse_header_payload(
    std::string payload,
    bool has_qualities,
    const std::vector<char>& read_name_separator_chars) {

    char record_type = has_qualities ? '@' : '>';
    if (!payload.empty() && (payload[0] == '@' || payload[0] == '>')) {
        record_type = payload[0];
        payload.erase(payload.begin());
    }

    std::istringstream payload_stream(payload);
    std::string raw_name;
    payload_stream >> raw_name;
    payload_stream >> std::ws;
    std::string extra;
    std::getline(payload_stream, extra);

    ParsedHeader parsed;
    parsed.read_name = canonical_read_name(raw_name, read_name_separator_chars);
    parsed.read_name_extra = extra;
    parsed.read_filter = parse_read_filter(record_type, extra);
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

std::string decode_twobit_sequence(const std::vector<uint8_t>& bytes, uint64_t nuclen) {
    static const char lookup[4] = {'A', 'C', 'G', 'T'};
    std::string sequence;
    sequence.resize(static_cast<size_t>(nuclen));
    for (uint64_t ibase = 0; ibase < nuclen; ++ibase) {
        const uint64_t word_index = ibase / 32;
        const uint64_t offset_in_word = ibase % 32;
        const uint64_t word = read_le64(bytes.data() + static_cast<size_t>(word_index * 8));
        sequence[static_cast<size_t>(ibase)] = lookup[(word >> (offset_in_word * 2)) & 0x3ULL];
    }
    return sequence;
}

struct DecodedBlock {
    uint64_t num_records = 0;
    uint64_t num_sequences = 0;
    std::vector<uint64_t> seq_lengths;
    std::vector<uint64_t> seq_offsets;
    std::vector<uint64_t> header_lengths;
    std::vector<uint64_t> header_offsets;
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

struct CbqInputModule::Impl {
    std::ifstream stream;
    uint32_t lane_index = 0;
    uint64_t read_ordinal = 0;
    uint64_t lane_record_index = 0;
    bool opened = false;
    FileHeaderFields file_header;
    DecodedBlock block;

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
        block.clear();
        opened = true;
        return true;
    }

    void close_lane() {
        if (stream.is_open()) {
            stream.close();
        }
        block.clear();
    }

    BlockLoadStatus load_next_block(const InputSourcePlan& plan, std::string* error) {
        block.clear();

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
        std::vector<uint8_t> seq_word_bytes;
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

        std::vector<uint64_t> n_positions;
        if (!decode_elias_fano_positions(npos_bytes, &n_positions, error)) {
            return BlockLoadStatus::Error;
        }

        if (!zstd_runtime().decompress(z_seq, seq_word_bytes_expected, &seq_word_bytes, "seq", error)) {
            return BlockLoadStatus::Error;
        }
        block.sequences = decode_twobit_sequence(seq_word_bytes, header.nuclen);
        for (uint64_t npos : n_positions) {
            if (npos >= header.nuclen) {
                set_error(error, "CBQ N-position exceeds decoded sequence length");
                return BlockLoadStatus::Error;
            }
            block.sequences[static_cast<size_t>(npos)] = 'N';
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
        return BlockLoadStatus::Block;
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

    while (impl_->lane_index < plan_.read_files_n) {
        if (impl_->block.record_index >= static_cast<size_t>(impl_->block.num_records)) {
            const BlockLoadStatus status = impl_->load_next_block(plan_, error);
            if (status == BlockLoadStatus::Error) {
                return InputStatus::Error;
            }
            if (status == BlockLoadStatus::End) {
                impl_->close_lane();
                ++impl_->lane_index;
                if (impl_->lane_index >= plan_.read_files_n) {
                    impl_->opened = false;
                    return InputStatus::End;
                }
                if (!impl_->open_lane(plan_, impl_->lane_index, error)) {
                    return InputStatus::Error;
                }
                continue;
            }
        }

        const uint64_t record_index = static_cast<uint64_t>(impl_->block.record_index);
        const uint64_t first_sequence_index = record_index * plan_.mate_count;
        std::string header_payload;
        if (impl_->file_header.has_headers() && !impl_->block.headers.empty()) {
            const uint64_t begin = impl_->block.header_offsets[static_cast<size_t>(first_sequence_index)];
            const uint64_t end = impl_->block.header_offsets[static_cast<size_t>(first_sequence_index + 1)];
            header_payload = impl_->block.headers.substr(static_cast<size_t>(begin), static_cast<size_t>(end - begin));
        } else {
            header_payload = std::to_string(impl_->lane_record_index + 1);
        }
        ParsedHeader parsed = parse_header_payload(
            header_payload,
            impl_->file_header.has_qualities(),
            plan_.read_name_separator_chars);

        record->read_name = parsed.read_name;
        record->read_name_extra = parsed.read_name_extra;
        record->lane_index = impl_->lane_index;
        record->read_ordinal = ++impl_->read_ordinal;
        record->read_filter = parsed.read_filter;
        record->mate_count = plan_.mate_count;
        record->mates.clear();
        record->mates.reserve(plan_.mate_count);

        for (uint32_t imate = 0; imate < plan_.mate_count; ++imate) {
            const uint64_t seq_index = first_sequence_index + imate;
            const uint64_t begin = impl_->block.seq_offsets[static_cast<size_t>(seq_index)];
            const uint64_t end = impl_->block.seq_offsets[static_cast<size_t>(seq_index + 1)];
            InputMateRecord mate;
            mate.sequence = impl_->block.sequences.substr(static_cast<size_t>(begin), static_cast<size_t>(end - begin));
            mate.quality = impl_->block.qualities.substr(static_cast<size_t>(begin), static_cast<size_t>(end - begin));
            apply_quality_conversion(mate.quality, plan_.out_qs_conversion_add);
            mate.original_length = static_cast<uint32_t>(mate.sequence.size());
            mate.has_quality = impl_->file_header.has_qualities();
            record->mates.push_back(mate);
        }

        ++impl_->block.record_index;
        ++impl_->lane_record_index;
        return InputStatus::Record;
    }

    impl_->opened = false;
    return InputStatus::End;
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

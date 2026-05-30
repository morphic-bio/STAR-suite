#include "input/CbqWriter.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <dlfcn.h>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

const uint64_t PRESENCE_PAIRED = 1ULL << 0;
const uint64_t PRESENCE_QUALITIES = 1ULL << 1;
const uint64_t PRESENCE_HEADERS = 1ULL << 2;

void append_le16(std::vector<uint8_t>* out, uint16_t value) {
    out->push_back(static_cast<uint8_t>(value & 0xffU));
    out->push_back(static_cast<uint8_t>((value >> 8) & 0xffU));
}

void append_le64(std::vector<uint8_t>* out, uint64_t value) {
    for (unsigned ii = 0; ii < 8; ++ii) {
        out->push_back(static_cast<uint8_t>((value >> (8 * ii)) & 0xffU));
    }
}

void append_i64(std::vector<uint8_t>* out, int64_t value) {
    append_le64(out, static_cast<uint64_t>(value));
}

void append_bool(std::vector<uint8_t>* out, bool value) {
    out->push_back(value ? 1U : 0U);
}

void append_vec_u64(std::vector<uint8_t>* out, const std::vector<uint64_t>& values) {
    append_le64(out, static_cast<uint64_t>(values.size()));
    for (size_t ii = 0; ii < values.size(); ++ii) {
        append_le64(out, values[ii]);
    }
}

void append_vec_i64(std::vector<uint8_t>* out, const std::vector<int64_t>& values) {
    append_le64(out, static_cast<uint64_t>(values.size()));
    for (size_t ii = 0; ii < values.size(); ++ii) {
        append_i64(out, values[ii]);
    }
}

void append_vec_u16(std::vector<uint8_t>* out, const std::vector<uint16_t>& values) {
    append_le64(out, static_cast<uint64_t>(values.size()));
    for (size_t ii = 0; ii < values.size(); ++ii) {
        append_le16(out, values[ii]);
    }
}

void append_bytes(std::vector<uint8_t>* out, const std::string& value) {
    out->insert(out->end(), value.begin(), value.end());
}

uint64_t floor_log2_u64(uint64_t value) {
    if (value == 0) {
        return 0;
    }
    uint64_t result = 0;
    while (value >>= 1U) {
        ++result;
    }
    return result;
}

std::string normalize_header_payload(const std::string& header_payload) {
    if (!header_payload.empty() && (header_payload[0] == '@' || header_payload[0] == '>')) {
        return header_payload.substr(1);
    }
    return header_payload;
}

class ZstdRuntime {
public:
    ZstdRuntime()
        : handle_(NULL),
          compress_bound_(NULL),
          compress_(NULL),
          is_error_(NULL),
          get_error_name_(NULL) {}

    ~ZstdRuntime() {
        if (handle_ != NULL) {
            dlclose(handle_);
        }
    }

    bool load(std::string* error) {
        if (compress_ != NULL) {
            return true;
        }
        handle_ = dlopen("libzstd.so.1", RTLD_LAZY | RTLD_LOCAL);
        if (handle_ == NULL) {
            handle_ = dlopen("libzstd.so", RTLD_LAZY | RTLD_LOCAL);
        }
        if (handle_ == NULL) {
            if (error != NULL) {
                *error = "could not load libzstd.so.1 or libzstd.so";
            }
            return false;
        }

        compress_bound_ = reinterpret_cast<CompressBoundFn>(dlsym(handle_, "ZSTD_compressBound"));
        compress_ = reinterpret_cast<CompressFn>(dlsym(handle_, "ZSTD_compress"));
        is_error_ = reinterpret_cast<IsErrorFn>(dlsym(handle_, "ZSTD_isError"));
        get_error_name_ = reinterpret_cast<GetErrorNameFn>(dlsym(handle_, "ZSTD_getErrorName"));
        if (compress_bound_ == NULL || compress_ == NULL ||
            is_error_ == NULL || get_error_name_ == NULL) {
            if (error != NULL) {
                *error = "libzstd is missing required compression symbols";
            }
            return false;
        }
        return true;
    }

    bool compress(const std::vector<uint8_t>& input,
                  int compression_level,
                  std::vector<uint8_t>* output,
                  const std::string& column_name,
                  std::string* error) {
        if (input.empty()) {
            output->clear();
            return true;
        }
        if (!load(error)) {
            return false;
        }
        const size_t bound = compress_bound_(input.size());
        output->assign(bound, 0);
        const size_t nbytes = compress_(
            output->data(),
            output->size(),
            input.data(),
            input.size(),
            compression_level);
        if (is_error_(nbytes)) {
            if (error != NULL) {
                std::ostringstream msg;
                msg << "zstd compression failed for CBQ column " << column_name
                    << ": " << get_error_name_(nbytes);
                *error = msg.str();
            }
            return false;
        }
        output->resize(nbytes);
        return true;
    }

private:
    typedef size_t (*CompressBoundFn)(size_t);
    typedef size_t (*CompressFn)(void*, size_t, const void*, size_t, int);
    typedef unsigned (*IsErrorFn)(size_t);
    typedef const char* (*GetErrorNameFn)(size_t);

    void* handle_;
    CompressBoundFn compress_bound_;
    CompressFn compress_;
    IsErrorFn is_error_;
    GetErrorNameFn get_error_name_;
};

struct BitVector {
    std::vector<uint64_t> words;
    uint64_t len;

    BitVector() : len(0) {}

    void assign_false(uint64_t nbits) {
        len = nbits;
        words.assign(static_cast<size_t>((nbits + 63U) / 64U), 0);
    }

    void set_bit(uint64_t pos) {
        if (pos >= len) {
            throw std::runtime_error("bit position exceeds BitVector length");
        }
        words[static_cast<size_t>(pos / 64U)] |= 1ULL << (pos % 64U);
    }

    void push_bits(uint64_t bits, uint64_t nbits) {
        uint64_t consumed = 0;
        while (consumed < nbits) {
            const uint64_t pos_in_word = len % 64U;
            if (pos_in_word == 0) {
                words.push_back(0);
            }
            const uint64_t take = std::min<uint64_t>(64U - pos_in_word, nbits - consumed);
            const uint64_t mask = take == 64U ? std::numeric_limits<uint64_t>::max()
                                               : ((1ULL << take) - 1ULL);
            words.back() |= ((bits >> consumed) & mask) << pos_in_word;
            consumed += take;
            len += take;
        }
    }

    void serialize(std::vector<uint8_t>* out) const {
        append_vec_u64(out, words);
        append_le64(out, len);
    }
};

struct DArrayIndex {
    std::vector<int64_t> block_inventory;
    std::vector<uint16_t> subblock_inventory;
    std::vector<uint64_t> overflow_positions;
    uint64_t num_positions;
    bool over_one;

    DArrayIndex() : num_positions(0), over_one(true) {}

    static void flush_block(std::vector<uint64_t>* positions,
                            std::vector<int64_t>* block_inventory,
                            std::vector<uint16_t>* subblock_inventory,
                            std::vector<uint64_t>* overflow_positions) {
        if (positions->empty()) {
            return;
        }
        const uint64_t first = positions->front();
        const uint64_t last = positions->back();
        if (last - first < 65536ULL) {
            block_inventory->push_back(static_cast<int64_t>(first));
            for (size_t ii = 0; ii < positions->size(); ii += 32) {
                subblock_inventory->push_back(static_cast<uint16_t>((*positions)[ii] - first));
            }
        } else {
            block_inventory->push_back(-static_cast<int64_t>(overflow_positions->size() + 1));
            for (size_t ii = 0; ii < positions->size(); ++ii) {
                overflow_positions->push_back((*positions)[ii]);
            }
            for (size_t ii = 0; ii < positions->size(); ii += 32) {
                subblock_inventory->push_back(std::numeric_limits<uint16_t>::max());
            }
        }
        positions->clear();
    }

    static DArrayIndex build_select1(const BitVector& bit_vector) {
        DArrayIndex index;
        std::vector<uint64_t> current_block;
        current_block.reserve(1024);
        for (size_t word_idx = 0; word_idx < bit_vector.words.size(); ++word_idx) {
            uint64_t word = bit_vector.words[word_idx];
            while (word != 0) {
                const unsigned bit = static_cast<unsigned>(__builtin_ctzll(word));
                const uint64_t pos = static_cast<uint64_t>(word_idx) * 64ULL + bit;
                if (pos >= bit_vector.len) {
                    break;
                }
                current_block.push_back(pos);
                ++index.num_positions;
                if (current_block.size() == 1024) {
                    flush_block(&current_block,
                                &index.block_inventory,
                                &index.subblock_inventory,
                                &index.overflow_positions);
                }
                word &= word - 1ULL;
            }
        }
        flush_block(&current_block,
                    &index.block_inventory,
                    &index.subblock_inventory,
                    &index.overflow_positions);
        return index;
    }

    void serialize(std::vector<uint8_t>* out) const {
        append_vec_i64(out, block_inventory);
        append_vec_u16(out, subblock_inventory);
        append_vec_u64(out, overflow_positions);
        append_le64(out, num_positions);
        append_bool(out, over_one);
    }
};

void serialize_darray_select1(const BitVector& bit_vector, std::vector<uint8_t>* out) {
    bit_vector.serialize(out);
    DArrayIndex::build_select1(bit_vector).serialize(out);
    append_bool(out, false);
    append_bool(out, false);
}

std::vector<uint8_t> encode_elias_fano(const std::vector<uint64_t>& positions, uint64_t universe) {
    std::vector<uint8_t> out;
    if (positions.empty()) {
        return out;
    }
    if (universe == 0) {
        throw std::runtime_error("cannot encode N positions with empty universe");
    }
    const uint64_t num_values = static_cast<uint64_t>(positions.size());
    const uint64_t low_len = floor_log2_u64(universe / num_values);
    if (low_len > 63U) {
        throw std::runtime_error("unsupported Elias-Fano low-bit width > 63");
    }

    const uint64_t high_len = (num_values + 1U) + (universe >> low_len) + 1U;
    BitVector high_bits;
    high_bits.assign_false(high_len);
    BitVector low_bits;

    uint64_t last = 0;
    for (size_t ii = 0; ii < positions.size(); ++ii) {
        const uint64_t value = positions[ii];
        if (ii != 0 && value < last) {
            throw std::runtime_error("N positions must be sorted");
        }
        if (value >= universe) {
            throw std::runtime_error("N position exceeds Elias-Fano universe");
        }
        last = value;
        if (low_len != 0) {
            const uint64_t low_mask = (1ULL << low_len) - 1ULL;
            low_bits.push_bits(value & low_mask, low_len);
        }
        high_bits.set_bit((value >> low_len) + static_cast<uint64_t>(ii));
    }

    serialize_darray_select1(high_bits, &out);
    low_bits.serialize(&out);
    append_le64(&out, low_len);
    append_le64(&out, universe);
    return out;
}

uint8_t base_code(char base, bool* is_ambiguous) {
    switch (base) {
        case 'A':
        case 'a':
            *is_ambiguous = false;
            return 0;
        case 'C':
        case 'c':
            *is_ambiguous = false;
            return 1;
        case 'G':
        case 'g':
            *is_ambiguous = false;
            return 2;
        case 'T':
        case 't':
            *is_ambiguous = false;
            return 3;
        default:
            *is_ambiguous = true;
            return 0;
    }
}

struct CbqBlock {
    std::vector<uint64_t> seq_lengths;
    std::vector<uint64_t> header_lengths;
    std::vector<uint64_t> n_positions;
    std::vector<uint64_t> seq_words;
    std::vector<uint8_t> headers;
    std::vector<uint8_t> qualities;
    uint64_t num_records;
    uint64_t num_sequences;
    uint64_t nuclen;
    uint64_t virtual_size;

    CbqBlock()
        : num_records(0), num_sequences(0), nuclen(0), virtual_size(0) {}

    void clear() {
        seq_lengths.clear();
        header_lengths.clear();
        n_positions.clear();
        seq_words.clear();
        headers.clear();
        qualities.clear();
        num_records = 0;
        num_sequences = 0;
        nuclen = 0;
        virtual_size = 0;
    }

    bool empty() const {
        return num_records == 0;
    }

    uint64_t estimate_record_size(const std::vector<star::input::CbqWriterSegment>& segments) const {
        uint64_t size = 0;
        for (size_t ii = 0; ii < segments.size(); ++ii) {
            size += ((static_cast<uint64_t>(segments[ii].sequence.size()) + 31U) / 32U) * 8U;
            size += static_cast<uint64_t>(normalize_header_payload(segments[ii].header_payload).size());
            size += static_cast<uint64_t>(segments[ii].quality.size());
        }
        return size;
    }

    void add_sequence(const std::string& sequence,
                      const std::string& quality,
                      const std::string& header_payload) {
        if (sequence.size() != quality.size()) {
            throw std::runtime_error("internal CBQ sequence/quality length mismatch");
        }
        const std::string normalized_header = normalize_header_payload(header_payload);
        seq_lengths.push_back(static_cast<uint64_t>(sequence.size()));
        header_lengths.push_back(static_cast<uint64_t>(normalized_header.size()));
        append_bytes(&headers, normalized_header);
        append_bytes(&qualities, quality);

        for (size_t ii = 0; ii < sequence.size(); ++ii) {
            const uint64_t global_offset = nuclen + static_cast<uint64_t>(ii);
            const size_t word_index = static_cast<size_t>(global_offset / 32U);
            if (seq_words.size() <= word_index) {
                seq_words.resize(word_index + 1, 0);
            }
            bool ambiguous = false;
            const uint8_t code = base_code(sequence[ii], &ambiguous);
            if (ambiguous) {
                n_positions.push_back(global_offset);
            }
            seq_words[word_index] |= static_cast<uint64_t>(code) << ((global_offset % 32U) * 2U);
        }
        nuclen += static_cast<uint64_t>(sequence.size());
        ++num_sequences;
    }

    void add_record(const std::vector<star::input::CbqWriterSegment>& segments) {
        for (size_t ii = 0; ii < segments.size(); ++ii) {
            add_sequence(segments[ii].sequence, segments[ii].quality, segments[ii].header_payload);
        }
        ++num_records;
        virtual_size += estimate_record_size(segments);
    }
};

struct EncodedBlock {
    std::vector<uint8_t> z_seq_len;
    std::vector<uint8_t> z_header_len;
    std::vector<uint8_t> z_npos;
    std::vector<uint8_t> z_seq;
    std::vector<uint8_t> z_flags;
    std::vector<uint8_t> z_headers;
    std::vector<uint8_t> z_qual;
    uint64_t nuclen;
    uint64_t len_nef;
    uint64_t num_records;
    uint64_t num_sequences;

    uint64_t payload_size() const {
        return static_cast<uint64_t>(z_seq_len.size() + z_header_len.size() +
                                     z_npos.size() + z_seq.size() + z_flags.size() +
                                     z_headers.size() + z_qual.size());
    }
};

std::vector<uint8_t> u64_column_bytes(const std::vector<uint64_t>& values) {
    std::vector<uint8_t> bytes;
    bytes.reserve(values.size() * 8U);
    for (size_t ii = 0; ii < values.size(); ++ii) {
        append_le64(&bytes, values[ii]);
    }
    return bytes;
}

std::vector<uint8_t> seq_word_bytes(const std::vector<uint64_t>& words, uint64_t nuclen) {
    const size_t expected_words = static_cast<size_t>((nuclen + 31U) / 32U);
    std::vector<uint8_t> bytes;
    bytes.reserve(expected_words * 8U);
    for (size_t ii = 0; ii < expected_words; ++ii) {
        append_le64(&bytes, ii < words.size() ? words[ii] : 0);
    }
    return bytes;
}

bool encode_block(const CbqBlock& block,
                  ZstdRuntime* zstd,
                  int compression_level,
                  EncodedBlock* encoded,
                  std::string* error) {
    encoded->nuclen = block.nuclen;
    encoded->num_records = block.num_records;
    encoded->num_sequences = block.num_sequences;

    const std::vector<uint8_t> seq_len_bytes = u64_column_bytes(block.seq_lengths);
    const std::vector<uint8_t> header_len_bytes = u64_column_bytes(block.header_lengths);
    const std::vector<uint8_t> npos_bytes = encode_elias_fano(block.n_positions, block.nuclen);
    const std::vector<uint8_t> seq_bytes = seq_word_bytes(block.seq_words, block.nuclen);
    encoded->len_nef = static_cast<uint64_t>(npos_bytes.size());

    return zstd->compress(seq_len_bytes, compression_level, &encoded->z_seq_len, "seq_len", error) &&
           zstd->compress(header_len_bytes, compression_level, &encoded->z_header_len, "header_len", error) &&
           zstd->compress(npos_bytes, compression_level, &encoded->z_npos, "npos", error) &&
           zstd->compress(seq_bytes, compression_level, &encoded->z_seq, "seq", error) &&
           zstd->compress(block.headers, compression_level, &encoded->z_headers, "headers", error) &&
           zstd->compress(block.qualities, compression_level, &encoded->z_qual, "qual", error);
}

void write_file_header(std::ostream& out,
                       bool paired,
                       uint64_t compression_level,
                       uint64_t block_size) {
    std::vector<uint8_t> header;
    header.insert(header.end(), "CBQFILE", "CBQFILE" + 7);
    header.push_back(1U);
    uint64_t flags = PRESENCE_QUALITIES | PRESENCE_HEADERS;
    if (paired) {
        flags |= PRESENCE_PAIRED;
    }
    append_le64(&header, flags);
    append_le64(&header, compression_level);
    append_le64(&header, block_size);
    header.resize(64, 0);
    out.write(reinterpret_cast<const char*>(header.data()), static_cast<std::streamsize>(header.size()));
}

void write_block_header(std::ostream& out, const EncodedBlock& block) {
    std::vector<uint8_t> header;
    header.insert(header.end(), "BLK", "BLK" + 3);
    header.push_back(1U);
    header.insert(header.end(), 4, 42U);
    append_le64(&header, static_cast<uint64_t>(block.z_seq_len.size()));
    append_le64(&header, static_cast<uint64_t>(block.z_header_len.size()));
    append_le64(&header, static_cast<uint64_t>(block.z_npos.size()));
    append_le64(&header, static_cast<uint64_t>(block.z_seq.size()));
    append_le64(&header, static_cast<uint64_t>(block.z_flags.size()));
    append_le64(&header, static_cast<uint64_t>(block.z_headers.size()));
    append_le64(&header, static_cast<uint64_t>(block.z_qual.size()));
    append_le64(&header, block.nuclen);
    append_le64(&header, block.len_nef);
    append_le64(&header, block.num_records);
    append_le64(&header, block.num_sequences);
    if (header.size() != 96U) {
        throw std::runtime_error("internal CBQ block header size mismatch");
    }
    out.write(reinterpret_cast<const char*>(header.data()), static_cast<std::streamsize>(header.size()));
}

void write_column(std::ostream& out, const std::vector<uint8_t>& bytes) {
    if (!bytes.empty()) {
        out.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    }
}

struct BlockRange {
    uint64_t offset;
    uint64_t cumulative_records;
};

} // namespace

namespace star {
namespace input {

struct CbqWriter::Impl {
    Impl()
        : mate_count(0),
          paired(false),
          compression_level(0),
          block_size(0),
          bytes_written(0),
          cumulative_records(0),
          opened(false),
          finished(false) {}

    bool flush(std::string* error) {
        if (block.empty()) {
            return true;
        }
        EncodedBlock encoded;
        try {
            if (!encode_block(block, &zstd, compression_level, &encoded, error)) {
                return false;
            }
        } catch (const std::exception& ex) {
            if (error != NULL) {
                *error = ex.what();
            }
            return false;
        }

        const uint64_t block_offset = bytes_written;
        write_block_header(out, encoded);
        write_column(out, encoded.z_seq_len);
        write_column(out, encoded.z_header_len);
        write_column(out, encoded.z_npos);
        write_column(out, encoded.z_seq);
        write_column(out, encoded.z_flags);
        write_column(out, encoded.z_headers);
        write_column(out, encoded.z_qual);
        if (!out.good()) {
            if (error != NULL) {
                *error = "failed while writing CBQ block to: " + path;
            }
            return false;
        }

        cumulative_records += encoded.num_records;
        BlockRange range;
        range.offset = block_offset;
        range.cumulative_records = cumulative_records;
        ranges.push_back(range);
        bytes_written += 96U + encoded.payload_size();
        block.clear();
        return true;
    }

    bool write_index(std::string* error) {
        std::vector<uint8_t> index_bytes;
        index_bytes.reserve(ranges.size() * 16U);
        for (size_t ii = 0; ii < ranges.size(); ++ii) {
            append_le64(&index_bytes, ranges[ii].offset);
            append_le64(&index_bytes, ranges[ii].cumulative_records);
        }

        std::vector<uint8_t> z_index;
        if (!zstd.compress(index_bytes, 0, &z_index, "index", error)) {
            return false;
        }

        std::vector<uint8_t> index_header;
        index_header.insert(index_header.end(), "CBQINDEX", "CBQINDEX" + 8);
        append_le64(&index_header, static_cast<uint64_t>(index_bytes.size()));
        append_le64(&index_header, static_cast<uint64_t>(z_index.size()));
        out.write(reinterpret_cast<const char*>(index_header.data()),
                  static_cast<std::streamsize>(index_header.size()));
        write_column(out, z_index);

        std::vector<uint8_t> footer;
        append_le64(&footer, static_cast<uint64_t>(z_index.size()));
        footer.insert(footer.end(), "CBQINDEX", "CBQINDEX" + 8);
        out.write(reinterpret_cast<const char*>(footer.data()), static_cast<std::streamsize>(footer.size()));
        if (!out.good()) {
            if (error != NULL) {
                *error = "failed while writing CBQ index to: " + path;
            }
            return false;
        }
        out.flush();
        if (!out.good()) {
            if (error != NULL) {
                *error = "failed while flushing CBQ output: " + path;
            }
            return false;
        }
        return true;
    }

    std::string path;
    uint32_t mate_count;
    bool paired;
    int compression_level;
    uint64_t block_size;
    std::ofstream out;
    ZstdRuntime zstd;
    CbqBlock block;
    std::vector<BlockRange> ranges;
    uint64_t bytes_written;
    uint64_t cumulative_records;
    bool opened;
    bool finished;
};

CbqWriter::CbqWriter()
    : impl_(new Impl()) {}

CbqWriter::~CbqWriter() {}

bool CbqWriter::open(const std::string& path,
                     uint32_t mate_count,
                     int compression_level,
                     uint64_t block_size,
                     std::string* error) {
    if (mate_count != 1 && mate_count != 2) {
        if (error != NULL) {
            *error = "CBQ writer supports mate_count 1 or 2";
        }
        return false;
    }
    if (block_size == 0) {
        if (error != NULL) {
            *error = "CBQ writer block_size must be positive";
        }
        return false;
    }
    impl_.reset(new Impl());
    impl_->path = path;
    impl_->mate_count = mate_count;
    impl_->paired = (mate_count == 2);
    impl_->compression_level = compression_level;
    impl_->block_size = block_size;
    impl_->out.open(path.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!impl_->out.good()) {
        if (error != NULL) {
            *error = "could not open CBQ output: " + path;
        }
        return false;
    }
    write_file_header(impl_->out, impl_->paired, static_cast<uint64_t>(compression_level), block_size);
    if (!impl_->out.good()) {
        if (error != NULL) {
            *error = "failed while writing CBQ header to: " + path;
        }
        return false;
    }
    impl_->bytes_written = 64U;
    impl_->opened = true;
    return true;
}

bool CbqWriter::add_record(const std::vector<CbqWriterSegment>& segments, std::string* error) {
    if (!impl_->opened || impl_->finished) {
        if (error != NULL) {
            *error = "CBQ writer is not open";
        }
        return false;
    }
    if (segments.size() != impl_->mate_count) {
        if (error != NULL) {
            std::ostringstream msg;
            msg << "CBQ writer expected " << impl_->mate_count
                << " segment(s), received " << segments.size();
            *error = msg.str();
        }
        return false;
    }
    for (size_t ii = 0; ii < segments.size(); ++ii) {
        if (segments[ii].sequence.size() != segments[ii].quality.size()) {
            if (error != NULL) {
                std::ostringstream msg;
                msg << "CBQ writer sequence/quality length mismatch for segment " << (ii + 1)
                    << ": sequence=" << segments[ii].sequence.size()
                    << " quality=" << segments[ii].quality.size();
                *error = msg.str();
            }
            return false;
        }
    }

    const uint64_t estimated = impl_->block.estimate_record_size(segments);
    if (!impl_->block.empty() && impl_->block.virtual_size + estimated > impl_->block_size) {
        if (!impl_->flush(error)) {
            return false;
        }
    }

    try {
        impl_->block.add_record(segments);
    } catch (const std::exception& ex) {
        if (error != NULL) {
            *error = ex.what();
        }
        return false;
    }

    if (impl_->block.virtual_size >= impl_->block_size) {
        return impl_->flush(error);
    }
    return true;
}

bool CbqWriter::finish(std::string* error) {
    if (!impl_->opened) {
        return true;
    }
    if (impl_->finished) {
        return true;
    }
    if (!impl_->flush(error)) {
        return false;
    }
    if (!impl_->write_index(error)) {
        return false;
    }
    impl_->out.close();
    impl_->finished = true;
    return true;
}

bool CbqWriter::is_open() const {
    return impl_->opened && !impl_->finished;
}

} // namespace input
} // namespace star

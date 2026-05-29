#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

namespace {

const uint64_t PRESENCE_PAIRED = 1ULL << 0;
const uint64_t PRESENCE_QUALITIES = 1ULL << 1;
const uint64_t PRESENCE_HEADERS = 1ULL << 2;
const uint64_t DEFAULT_BLOCK_SIZE = 1024ULL * 1024ULL;

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

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

std::string trim_fastq_header_payload(const std::string& header_line) {
    if (!header_line.empty() && header_line[0] == '@') {
        return header_line.substr(1);
    }
    return header_line;
}

std::string canonical_read_id(const std::string& header_payload) {
    size_t end = header_payload.find_first_of(" \t\r\n");
    std::string id = header_payload.substr(0, end == std::string::npos ? header_payload.size() : end);
    if (id.size() > 2 && id[id.size() - 2] == '/' &&
        (id[id.size() - 1] == '1' || id[id.size() - 1] == '2')) {
        id.resize(id.size() - 2);
    }
    return id;
}

std::string parse_error_context(const std::string& path, uint64_t record_ordinal) {
    std::ostringstream msg;
    msg << path << " record " << record_ordinal;
    return msg.str();
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
            *error = "could not load libzstd.so.1 or libzstd.so";
            return false;
        }

        compress_bound_ = reinterpret_cast<CompressBoundFn>(dlsym(handle_, "ZSTD_compressBound"));
        compress_ = reinterpret_cast<CompressFn>(dlsym(handle_, "ZSTD_compress"));
        is_error_ = reinterpret_cast<IsErrorFn>(dlsym(handle_, "ZSTD_isError"));
        get_error_name_ = reinterpret_cast<GetErrorNameFn>(dlsym(handle_, "ZSTD_getErrorName"));
        if (compress_bound_ == NULL || compress_ == NULL ||
            is_error_ == NULL || get_error_name_ == NULL) {
            *error = "libzstd is missing required compression symbols";
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
            std::ostringstream msg;
            msg << "zstd compression failed for CBQ column " << column_name
                << ": " << get_error_name_(nbytes);
            *error = msg.str();
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

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

class FastqReader {
public:
    explicit FastqReader(const std::string& path)
        : path_(path),
          gzip_(ends_with(path, ".gz")),
          gz_(NULL),
          line_number_(0),
          record_ordinal_(0) {}

    ~FastqReader() {
        close();
    }

    bool open(std::string* error) {
        close();
        if (gzip_) {
            gz_ = gzopen(path_.c_str(), "rb");
            if (gz_ == NULL) {
                *error = "could not open gzip FASTQ: " + path_;
                return false;
            }
        } else {
            stream_.open(path_.c_str(), std::ios::in | std::ios::binary);
            if (!stream_.good()) {
                *error = "could not open FASTQ: " + path_;
                return false;
            }
        }
        return true;
    }

    void close() {
        if (gz_ != NULL) {
            gzclose(gz_);
            gz_ = NULL;
        }
        if (stream_.is_open()) {
            stream_.close();
        }
    }

    bool read_record(FastqRecord* record, bool* got_record, std::string* error) {
        got_record = got_record == NULL ? &local_got_record_ : got_record;
        *got_record = false;
        std::string plus;

        if (!getline(&record->header, error)) {
            return error->empty();
        }
        *got_record = true;
        ++record_ordinal_;
        if (!getline(&record->sequence, error) ||
            !getline(&plus, error) ||
            !getline(&record->quality, error)) {
            if (error->empty()) {
                *error = "truncated FASTQ at " + parse_error_context(path_, record_ordinal_);
            }
            return false;
        }
        if (record->header.empty() || record->header[0] != '@') {
            *error = "FASTQ header does not start with @ at " +
                     parse_error_context(path_, record_ordinal_);
            return false;
        }
        if (plus.empty() || plus[0] != '+') {
            *error = "FASTQ plus line does not start with + at " +
                     parse_error_context(path_, record_ordinal_);
            return false;
        }
        if (record->sequence.size() != record->quality.size()) {
            std::ostringstream msg;
            msg << "FASTQ sequence/quality length mismatch at "
                << parse_error_context(path_, record_ordinal_)
                << ": sequence=" << record->sequence.size()
                << " quality=" << record->quality.size();
            *error = msg.str();
            return false;
        }
        return true;
    }

    uint64_t record_ordinal() const {
        return record_ordinal_;
    }

private:
    bool getline(std::string* line, std::string* error) {
        line->clear();
        if (gzip_) {
            return getline_gzip(line, error);
        }

        if (!std::getline(stream_, *line)) {
            if (stream_.eof()) {
                error->clear();
                return false;
            }
            *error = "failed while reading FASTQ: " + path_;
            return false;
        }
        ++line_number_;
        if (!line->empty() && (*line)[line->size() - 1] == '\r') {
            line->resize(line->size() - 1);
        }
        return true;
    }

    bool getline_gzip(std::string* line, std::string* error) {
        char buffer[65536];
        while (true) {
            char* result = gzgets(gz_, buffer, static_cast<int>(sizeof(buffer)));
            if (result == NULL) {
                if (!line->empty() && gzeof(gz_)) {
                    ++line_number_;
                    if (!line->empty() && (*line)[line->size() - 1] == '\r') {
                        line->resize(line->size() - 1);
                    }
                    return true;
                }
                if (gzeof(gz_)) {
                    error->clear();
                    return false;
                }
                int zerr = Z_OK;
                const char* zmsg = gzerror(gz_, &zerr);
                std::ostringstream msg;
                msg << "failed while reading gzip FASTQ " << path_;
                if (zmsg != NULL) {
                    msg << ": " << zmsg;
                }
                *error = msg.str();
                return false;
            }

            const size_t len = std::strlen(buffer);
            line->append(buffer, len);
            if (!line->empty() && (*line)[line->size() - 1] == '\n') {
                line->resize(line->size() - 1);
                if (!line->empty() && (*line)[line->size() - 1] == '\r') {
                    line->resize(line->size() - 1);
                }
                ++line_number_;
                return true;
            }
            if (len + 1 < sizeof(buffer) && gzeof(gz_)) {
                ++line_number_;
                if (!line->empty() && (*line)[line->size() - 1] == '\r') {
                    line->resize(line->size() - 1);
                }
                return true;
            }
        }
    }

    std::string path_;
    bool gzip_;
    gzFile gz_;
    std::ifstream stream_;
    uint64_t line_number_;
    uint64_t record_ordinal_;
    bool local_got_record_;
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
    append_bool(out, false); // no select0 index
    append_bool(out, false); // no rank9 index
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

    uint64_t estimate_record_size(const FastqRecord& first, const FastqRecord* second) const {
        uint64_t size = 0;
        size += ((static_cast<uint64_t>(first.sequence.size()) + 31U) / 32U) * 8U;
        size += static_cast<uint64_t>(trim_fastq_header_payload(first.header).size());
        size += static_cast<uint64_t>(first.quality.size());
        if (second != NULL) {
            size += ((static_cast<uint64_t>(second->sequence.size()) + 31U) / 32U) * 8U;
            size += static_cast<uint64_t>(trim_fastq_header_payload(second->header).size());
            size += static_cast<uint64_t>(second->quality.size());
        }
        return size;
    }

    void add_sequence(const std::string& sequence,
                      const std::string& quality,
                      const std::string& header_payload) {
        if (sequence.size() != quality.size()) {
            throw std::runtime_error("internal sequence/quality length mismatch");
        }
        seq_lengths.push_back(static_cast<uint64_t>(sequence.size()));
        header_lengths.push_back(static_cast<uint64_t>(header_payload.size()));
        append_bytes(&headers, header_payload);
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

    void add_record(const FastqRecord& first, const FastqRecord* second) {
        add_sequence(first.sequence, first.quality, trim_fastq_header_payload(first.header));
        if (second != NULL) {
            add_sequence(second->sequence, second->quality, trim_fastq_header_payload(second->header));
        }
        ++num_records;
        virtual_size += estimate_record_size(first, second);
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

class CbqWriter {
public:
    CbqWriter(const std::string& path,
              bool paired,
              int compression_level,
              uint64_t block_size)
        : path_(path),
          paired_(paired),
          compression_level_(compression_level),
          block_size_(block_size),
          bytes_written_(0),
          cumulative_records_(0) {}

    bool open(std::string* error) {
        out_.open(path_.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!out_.good()) {
            *error = "could not open CBQ output: " + path_;
            return false;
        }
        write_file_header(out_, paired_, static_cast<uint64_t>(compression_level_), block_size_);
        bytes_written_ = 64U;
        return true;
    }

    bool add_record(const FastqRecord& first, const FastqRecord* second, std::string* error) {
        const uint64_t estimated = block_.estimate_record_size(first, second);
        if (!block_.empty() && block_.virtual_size + estimated > block_size_) {
            if (!flush(error)) {
                return false;
            }
        }
        block_.add_record(first, second);
        if (block_.virtual_size >= block_size_) {
            return flush(error);
        }
        return true;
    }

    bool finish(std::string* error) {
        if (!flush(error)) {
            return false;
        }
        return write_index(error);
    }

private:
    bool flush(std::string* error) {
        if (block_.empty()) {
            return true;
        }
        EncodedBlock encoded;
        try {
            if (!encode_block(block_, &zstd_, compression_level_, &encoded, error)) {
                return false;
            }
        } catch (const std::exception& ex) {
            *error = ex.what();
            return false;
        }

        const uint64_t block_offset = bytes_written_;
        write_block_header(out_, encoded);
        write_column(out_, encoded.z_seq_len);
        write_column(out_, encoded.z_header_len);
        write_column(out_, encoded.z_npos);
        write_column(out_, encoded.z_seq);
        write_column(out_, encoded.z_flags);
        write_column(out_, encoded.z_headers);
        write_column(out_, encoded.z_qual);
        if (!out_.good()) {
            *error = "failed while writing CBQ block to: " + path_;
            return false;
        }

        cumulative_records_ += encoded.num_records;
        BlockRange range;
        range.offset = block_offset;
        range.cumulative_records = cumulative_records_;
        ranges_.push_back(range);
        bytes_written_ += 96U + encoded.payload_size();
        block_.clear();
        return true;
    }

    bool write_index(std::string* error) {
        std::vector<uint8_t> index_bytes;
        index_bytes.reserve(ranges_.size() * 16U);
        for (size_t ii = 0; ii < ranges_.size(); ++ii) {
            append_le64(&index_bytes, ranges_[ii].offset);
            append_le64(&index_bytes, ranges_[ii].cumulative_records);
        }

        std::vector<uint8_t> z_index;
        if (!zstd_.compress(index_bytes, 0, &z_index, "index", error)) {
            return false;
        }

        std::vector<uint8_t> index_header;
        index_header.insert(index_header.end(), "CBQINDEX", "CBQINDEX" + 8);
        append_le64(&index_header, static_cast<uint64_t>(index_bytes.size()));
        append_le64(&index_header, static_cast<uint64_t>(z_index.size()));
        out_.write(reinterpret_cast<const char*>(index_header.data()),
                   static_cast<std::streamsize>(index_header.size()));
        write_column(out_, z_index);

        std::vector<uint8_t> footer;
        append_le64(&footer, static_cast<uint64_t>(z_index.size()));
        footer.insert(footer.end(), "CBQINDEX", "CBQINDEX" + 8);
        out_.write(reinterpret_cast<const char*>(footer.data()), static_cast<std::streamsize>(footer.size()));
        if (!out_.good()) {
            *error = "failed while writing CBQ index to: " + path_;
            return false;
        }
        return true;
    }

    std::string path_;
    bool paired_;
    int compression_level_;
    uint64_t block_size_;
    std::ofstream out_;
    ZstdRuntime zstd_;
    CbqBlock block_;
    std::vector<BlockRange> ranges_;
    uint64_t bytes_written_;
    uint64_t cumulative_records_;
};

struct Options {
    std::vector<std::string> read_files;
    std::string output;
    int compression_level;
    uint64_t block_size;
    bool validate_read_names;

    Options()
        : compression_level(0),
          block_size(DEFAULT_BLOCK_SIZE),
          validate_read_names(true) {}
};

void usage(std::ostream& out, const char* prog) {
    out << "Usage: " << prog << " --readFilesIn R1.fastq[.gz] [R2.fastq[.gz]] --outFile out.cbq [options]\n"
        << "Options:\n"
        << "  --compressionLevel N, -l N   Zstandard compression level (default: 0)\n"
        << "  --blockSize N                CBQ virtual block size in bytes (default: 1048576)\n"
        << "  --no-validate-read-names     Do not require paired FASTQ read ids to match\n";
}

bool parse_u64(const std::string& text, uint64_t* value) {
    char* end = NULL;
    errno = 0;
    unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0') {
        return false;
    }
    *value = static_cast<uint64_t>(parsed);
    return true;
}

bool parse_i32(const std::string& text, int* value) {
    char* end = NULL;
    errno = 0;
    long parsed = std::strtol(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0' ||
        parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    *value = static_cast<int>(parsed);
    return true;
}

bool require_value(int argc, char** argv, int* index, std::string* value, std::string* error) {
    if (*index + 1 >= argc) {
        *error = std::string("missing value for ") + argv[*index];
        return false;
    }
    ++(*index);
    *value = argv[*index];
    return true;
}

bool parse_args(int argc, char** argv, Options* opts, std::string* error) {
    for (int ii = 1; ii < argc; ++ii) {
        const std::string arg = argv[ii];
        std::string value;
        if (arg == "--help" || arg == "-h") {
            usage(std::cout, argv[0]);
            std::exit(0);
        } else if (arg == "--readFilesIn") {
            while (ii + 1 < argc) {
                const std::string next = argv[ii + 1];
                if (!next.empty() && next[0] == '-') {
                    break;
                }
                ++ii;
                opts->read_files.push_back(argv[ii]);
                if (opts->read_files.size() == 2) {
                    break;
                }
            }
        } else if (arg == "--outFile" || arg == "-o") {
            if (!require_value(argc, argv, &ii, &opts->output, error)) {
                return false;
            }
        } else if (arg == "--compressionLevel" || arg == "-l") {
            if (!require_value(argc, argv, &ii, &value, error) ||
                !parse_i32(value, &opts->compression_level)) {
                *error = "invalid compression level: " + value;
                return false;
            }
        } else if (arg == "--blockSize") {
            if (!require_value(argc, argv, &ii, &value, error) ||
                !parse_u64(value, &opts->block_size) ||
                opts->block_size == 0) {
                *error = "invalid block size: " + value;
                return false;
            }
        } else if (arg == "--no-validate-read-names") {
            opts->validate_read_names = false;
        } else {
            *error = "unknown argument: " + arg;
            return false;
        }
    }

    if (opts->read_files.empty() || opts->read_files.size() > 2) {
        *error = "--readFilesIn requires one or two FASTQ paths";
        return false;
    }
    if (opts->output.empty()) {
        *error = "--outFile is required";
        return false;
    }
    return true;
}

bool encode_fastq_to_cbq(const Options& opts, std::string* error) {
    const bool paired = opts.read_files.size() == 2;
    FastqReader r1(opts.read_files[0]);
    FastqReader r2(paired ? opts.read_files[1] : std::string());
    if (!r1.open(error)) {
        return false;
    }
    if (paired && !r2.open(error)) {
        return false;
    }

    CbqWriter writer(opts.output, paired, opts.compression_level, opts.block_size);
    if (!writer.open(error)) {
        return false;
    }

    uint64_t records = 0;
    while (true) {
        FastqRecord rec1;
        FastqRecord rec2;
        bool got1 = false;
        bool got2 = false;
        if (!r1.read_record(&rec1, &got1, error)) {
            return false;
        }
        if (paired) {
            if (!r2.read_record(&rec2, &got2, error)) {
                return false;
            }
            if (got1 != got2) {
                *error = "paired FASTQ files have different record counts";
                return false;
            }
        }
        if (!got1) {
            break;
        }

        if (paired && opts.validate_read_names) {
            const std::string id1 = canonical_read_id(trim_fastq_header_payload(rec1.header));
            const std::string id2 = canonical_read_id(trim_fastq_header_payload(rec2.header));
            if (id1 != id2) {
                std::ostringstream msg;
                msg << "paired FASTQ read-name mismatch at record " << (records + 1)
                    << ": " << id1 << " vs " << id2;
                *error = msg.str();
                return false;
            }
        }

        if (!writer.add_record(rec1, paired ? &rec2 : NULL, error)) {
            return false;
        }
        ++records;
        if (records % 1000000ULL == 0) {
            std::cerr << "encoded " << records << " records\n";
        }
    }

    if (!writer.finish(error)) {
        return false;
    }
    std::cerr << "encoded " << records << " records to " << opts.output << "\n";
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Options opts;
    std::string error;
    if (!parse_args(argc, argv, &opts, &error)) {
        std::cerr << "ERROR: " << error << "\n\n";
        usage(std::cerr, argv[0]);
        return 2;
    }
    if (!encode_fastq_to_cbq(opts, &error)) {
        std::cerr << "ERROR: " << error << "\n";
        return 1;
    }
    return 0;
}

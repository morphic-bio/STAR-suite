#include "input/CbqStarAdapter.h"

#include "SequenceFuns.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>

namespace star {
namespace input {
namespace {

bool set_error(std::string* error, const std::string& message) {
    if (error != nullptr) {
        *error = message;
    }
    return false;
}

std::string span_string(CbqByteSpan span) {
    if (span.data == nullptr || span.size == 0) {
        return std::string();
    }
    return std::string(span.data, span.size);
}

struct PackedByteDecode {
    char ascii[4];
    char number[4];
};

const std::array<PackedByteDecode, 256>& packed_byte_decode_lookup() {
    static const std::array<PackedByteDecode, 256> lookup = [] {
        std::array<PackedByteDecode, 256> table = {};
        static const char ascii_lookup[4] = {'A', 'C', 'G', 'T'};
        for (unsigned byte = 0; byte < table.size(); ++byte) {
            for (unsigned base = 0; base < 4; ++base) {
                const unsigned code = (byte >> (base * 2)) & 0x3U;
                table[byte].ascii[base] = ascii_lookup[code];
                table[byte].number[base] = static_cast<char>(code);
            }
        }
        return table;
    }();
    return lookup;
}

void copy_packed_sequence_to_star_buffers(const CbqPackedSequenceView& packed,
                                          char* read0,
                                          char* read1) {
    const auto& lookup = packed_byte_decode_lookup();

    uint32_t ii = 0;
    while (ii < packed.length) {
        const uint64_t global_offset = packed.base_offset + static_cast<uint64_t>(ii);
        const size_t byte_offset = static_cast<size_t>(global_offset >> 2);
        const unsigned base_in_byte = static_cast<unsigned>(global_offset & 0x3ULL);
        const uint32_t take = std::min<uint32_t>(4U - base_in_byte, packed.length - ii);

        if (packed.words == nullptr || byte_offset >= packed.word_bytes) {
            for (uint32_t jj = 0; jj < take; ++jj) {
                read0[ii + jj] = 'N';
                read1[ii + jj] = static_cast<char>(4);
            }
        } else if (base_in_byte == 0 && take == 4) {
            const PackedByteDecode& decoded = lookup[packed.words[byte_offset]];
            std::memcpy(read0 + ii, decoded.ascii, 4);
            std::memcpy(read1 + ii, decoded.number, 4);
        } else {
            const PackedByteDecode& decoded = lookup[packed.words[byte_offset]];
            for (uint32_t jj = 0; jj < take; ++jj) {
                read0[ii + jj] = decoded.ascii[base_in_byte + jj];
                read1[ii + jj] = decoded.number[base_in_byte + jj];
            }
        }

        ii += take;
    }

    if (packed.n_positions != nullptr) {
        for (size_t inpos = 0; inpos < packed.n_positions_count; ++inpos) {
            const uint64_t npos = packed.n_positions[inpos];
            if (npos < packed.base_offset) {
                continue;
            }
            const uint64_t local = npos - packed.base_offset;
            if (local >= packed.length) {
                break;
            }
            read0[local] = 'N';
            read1[local] = static_cast<char>(4);
        }
    }
    read0[packed.length] = '\0';
}

void copy_packed_sequence_to_ascii(const CbqPackedSequenceView& packed,
                                   char* dest) {
    const auto& lookup = packed_byte_decode_lookup();

    uint32_t ii = 0;
    while (ii < packed.length) {
        const uint64_t global_offset = packed.base_offset + static_cast<uint64_t>(ii);
        const size_t byte_offset = static_cast<size_t>(global_offset >> 2);
        const unsigned base_in_byte = static_cast<unsigned>(global_offset & 0x3ULL);
        const uint32_t take = std::min<uint32_t>(4U - base_in_byte, packed.length - ii);

        if (packed.words == nullptr || byte_offset >= packed.word_bytes) {
            std::fill(dest + ii, dest + ii + take, 'N');
        } else if (base_in_byte == 0 && take == 4) {
            const PackedByteDecode& decoded = lookup[packed.words[byte_offset]];
            std::memcpy(dest + ii, decoded.ascii, 4);
        } else {
            const PackedByteDecode& decoded = lookup[packed.words[byte_offset]];
            for (uint32_t jj = 0; jj < take; ++jj) {
                dest[ii + jj] = decoded.ascii[base_in_byte + jj];
            }
        }

        ii += take;
    }

    if (packed.n_positions != nullptr) {
        for (size_t inpos = 0; inpos < packed.n_positions_count; ++inpos) {
            const uint64_t npos = packed.n_positions[inpos];
            if (npos < packed.base_offset) {
                continue;
            }
            const uint64_t local = npos - packed.base_offset;
            if (local >= packed.length) {
                break;
            }
            dest[local] = 'N';
        }
    }
    dest[packed.length] = '\0';
}

bool copy_mate_sequence_to_ascii(const CbqStarChunkMate& mate,
                                 char* dest,
                                 std::string* error) {
    if (dest == nullptr) {
        return set_error(error, "CBQ STAR adapter ASCII sequence destination is null");
    }
    if (mate.length > DEF_readSeqLengthMax) {
        std::ostringstream msg;
        msg << "CBQ STAR adapter sequence length " << mate.length
            << " exceeds DEF_readSeqLengthMax=" << DEF_readSeqLengthMax;
        return set_error(error, msg.str());
    }
    if (mate.packed_sequence.available) {
        copy_packed_sequence_to_ascii(mate.packed_sequence, dest);
    } else {
        if (mate.sequence.size != mate.length) {
            return set_error(error, "CBQ STAR adapter fallback sequence length is inconsistent");
        }
        if (mate.length != 0) {
            std::memcpy(dest, mate.sequence.data, mate.length);
        }
        dest[mate.length] = '\0';
    }
    return true;
}

bool copy_span_to_cstr(CbqByteSpan span,
                       char* dest,
                       size_t dest_size,
                       const std::string& field_name,
                       std::string* error) {
    if (dest == nullptr) {
        return set_error(error, "CBQ STAR adapter destination is null for " + field_name);
    }
    if (span.size + 1 > dest_size) {
        std::ostringstream msg;
        msg << "CBQ STAR adapter " << field_name << " length " << span.size
            << " exceeds destination capacity " << dest_size;
        return set_error(error, msg.str());
    }
    if (span.size != 0) {
        std::memcpy(dest, span.data, span.size);
    }
    dest[span.size] = '\0';
    return true;
}

bool copy_string_to_cstr(const std::string& value,
                         char* dest,
                         size_t dest_size,
                         const std::string& field_name,
                         std::string* error) {
    CbqByteSpan span;
    span.data = value.data();
    span.size = value.size();
    return copy_span_to_cstr(span, dest, dest_size, field_name, error);
}

bool copy_prefixed_span_to_cstr(char prefix,
                                CbqByteSpan span,
                                char* dest,
                                size_t dest_size,
                                const std::string& field_name,
                                size_t* copied_size,
                                std::string* error) {
    if (dest == nullptr) {
        return set_error(error, "CBQ STAR adapter destination is null for " + field_name);
    }
    if (span.size + 2 > dest_size) {
        std::ostringstream msg;
        msg << "CBQ STAR adapter " << field_name << " length " << (span.size + 1)
            << " exceeds destination capacity " << dest_size;
        return set_error(error, msg.str());
    }
    dest[0] = prefix;
    if (span.size != 0) {
        std::memcpy(dest + 1, span.data, span.size);
    }
    dest[span.size + 1] = '\0';
    if (copied_size != nullptr) {
        *copied_size = span.size + 2;
    }
    return true;
}

void apply_quality_conversion(char* quality, uint length, int add) {
    if (add == 0) {
        return;
    }
    for (uint ii = 0; ii < length; ++ii) {
        int qs = static_cast<int>(quality[ii]) + add;
        if (qs < 33) {
            qs = 33;
        } else if (qs > 126) {
            qs = 126;
        }
        quality[ii] = static_cast<char>(qs);
    }
}

void set_clip_info_from_fastq_plus_line(ClipMate& clip_mate,
                                        char* sequence,
                                        uint length) {
    clip_mate.clippedInfo = '+';

    if (clip_mate.type != 10 || clip_mate.cr4 == nullptr || clip_mate.adSeqNum.empty()) {
        return;
    }

    clip_mate.cr4->opalFillOneSeq(0, sequence, length);
    clip_mate.cr4->opalAlign(reinterpret_cast<uint8_t*>(clip_mate.adSeqNum.data()),
                             clip_mate.adSeqNum.size(),
                             1);

    const int clip_length = clip_mate.cr4->opalRes[0].endLocationTarget + 1;
    const int score = clip_mate.cr4->opalRes[0].score;
    const bool no_clip = score < 20 || (score == 20 && clip_length > 26) ||
                         (score == 21 && clip_length > 30);
    clip_mate.clippedInfo = static_cast<char>(no_clip ? 0 : clip_length);
}

bool validate_buffers(const CbqStarReadBuffers& buffers, std::string* error) {
    if (buffers.read_name_mates == nullptr ||
        buffers.read0 == nullptr ||
        buffers.read1 == nullptr ||
        buffers.qual0 == nullptr ||
        buffers.read_name_extra == nullptr ||
        buffers.read_length == nullptr ||
        buffers.read_length_original == nullptr ||
        buffers.i_read_all == nullptr ||
        buffers.read_files_index == nullptr ||
        buffers.read_filter == nullptr ||
        buffers.read_file_type == nullptr) {
        return set_error(error, "CBQ STAR adapter received incomplete STAR buffer pointers");
    }
    return true;
}

bool chunk_mate_from_segment(const CbqSegmentView& segment,
                             CbqStarChunkMate* mate,
                             std::string* error) {
    if (mate == nullptr) {
        return set_error(error, "CBQ STAR chunk mate destination is null");
    }
    const size_t sequence_length = cbq_segment_sequence_length(segment);
    if (sequence_length > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
        std::ostringstream msg;
        msg << "CBQ STAR chunk sequence length " << sequence_length
            << " exceeds uint32_t";
        return set_error(error, msg.str());
    }
    mate->sequence = segment.sequence;
    mate->quality = segment.quality;
    mate->packed_sequence = segment.packed_sequence;
    mate->length = static_cast<uint32_t>(sequence_length);
    mate->original_length = segment.original_length;
    mate->has_quality = segment.has_quality;
    return true;
}

bool load_cbq_star_read_into_buffers(const CbqStarChunkRead& read,
                                     const CbqStarChunkMate* mates,
                                     const CbqStarAdapterOptions& options,
                                     CbqStarReadBuffers* buffers,
                                     std::vector<std::vector<ClipMate>>* clip_mates,
                                     std::string* error) {
    if (buffers == nullptr) {
        return set_error(error, "CBQ STAR adapter requires non-null buffers");
    }
    if (!validate_buffers(*buffers, error)) {
        return false;
    }
    if (options.read_nends == 0) {
        return set_error(error, "CBQ STAR adapter read_nends is zero");
    }
    if (mates == nullptr) {
        return set_error(error, "CBQ STAR adapter received null mate records");
    }
    if (read.mate_count < options.read_nends) {
        std::ostringstream msg;
        msg << "CBQ STAR adapter record has " << read.mate_count
            << " mates but STAR expects " << options.read_nends;
        return set_error(error, msg.str());
    }
    if (clip_mates != nullptr && clip_mates->size() < options.read_nends) {
        return set_error(error, "CBQ STAR adapter clip_mates has fewer mates than read_nends");
    }
    if (buffers->read_name_extra->size() < options.read_nends) {
        buffers->read_name_extra->resize(options.read_nends);
    }

    const char header_prefix = read.has_quality ? '@' : '>';

    size_t read_name_bytes = 0;
    if (options.out_sam_read_id_number) {
        const std::string read_name_core = std::to_string(read.read_ordinal);
        const std::string star_read_name = std::string(1, header_prefix) + read_name_core;
        if (!copy_string_to_cstr(star_read_name,
                                 buffers->read_name_mates[0],
                                 DEF_readNameLengthMax,
                                 "read name",
                                 error)) {
            return false;
        }
        read_name_bytes = star_read_name.size() + 1;
    } else {
        if (!copy_prefixed_span_to_cstr(header_prefix,
                                        read.read_name,
                                        buffers->read_name_mates[0],
                                        DEF_readNameLengthMax,
                                        "read name",
                                        &read_name_bytes,
                                        error)) {
            return false;
        }
    }
    for (uint32 imate = 1; imate < options.read_nends; ++imate) {
        std::memcpy(buffers->read_name_mates[imate], buffers->read_name_mates[0], read_name_bytes);
    }

    *buffers->i_read_all = read.read_ordinal;
    *buffers->read_files_index = read.lane_index;
    *buffers->read_filter = read.read_filter;
    *buffers->read_file_type = read.has_quality ? 2 : 1;

    for (uint32 imate = 0; imate < options.read_nends; ++imate) {
        const CbqStarChunkMate& mate = mates[imate];
        const size_t sequence_length = mate.length;
        if (sequence_length > DEF_readSeqLengthMax) {
            std::ostringstream msg;
            msg << "CBQ STAR adapter mate " << (imate + 1)
                << " sequence length " << sequence_length
                << " exceeds DEF_readSeqLengthMax=" << DEF_readSeqLengthMax;
            return set_error(error, msg.str());
        }
        if (mate.has_quality && mate.quality.size != sequence_length) {
            std::ostringstream msg;
            msg << "CBQ STAR adapter mate " << (imate + 1)
                << " quality length " << mate.quality.size
                << " does not match sequence length " << sequence_length;
            return set_error(error, msg.str());
        }

        const uint length = static_cast<uint>(sequence_length);
        if (mate.packed_sequence.available) {
            copy_packed_sequence_to_star_buffers(mate.packed_sequence,
                                                 buffers->read0[imate],
                                                 buffers->read1[imate]);
        } else {
            if (!copy_span_to_cstr(mate.sequence,
                                   buffers->read0[imate],
                                   DEF_readSeqLengthMax + 1,
                                   "sequence",
                                   error)) {
                return false;
            }
            convertNucleotidesToNumbers(buffers->read0[imate], buffers->read1[imate], length);
        }

        if (mate.has_quality) {
            if (!copy_span_to_cstr(mate.quality,
                                   buffers->qual0[imate],
                                   DEF_readSeqLengthMax + 1,
                                   "quality",
                                   error)) {
                return false;
            }
            apply_quality_conversion(buffers->qual0[imate], length, options.out_qs_conversion_add);
        } else {
            std::fill(buffers->qual0[imate], buffers->qual0[imate] + length, 'A');
            buffers->qual0[imate][length] = '\0';
        }

        buffers->read_length[imate] = length;
        buffers->read_length_original[imate] = length;

        if (options.preserve_read_name_extra) {
            (*buffers->read_name_extra)[imate] = span_string(read.read_name_extra);
        } else {
            (*buffers->read_name_extra)[imate].clear();
        }

        if (clip_mates != nullptr && mate.has_quality) {
            if (!mate.clip_info_5p_prepared) {
                set_clip_info_from_fastq_plus_line((*clip_mates)[imate][0],
                                                  buffers->read0[imate],
                                                  length);
            } else {
                (*clip_mates)[imate][0].clippedInfo = mate.clip_info_5p;
            }
        }
        if (!options.trim_cutadapt_enabled && clip_mates != nullptr) {
            (*clip_mates)[imate][0].clip(buffers->read_length[imate], buffers->read1[imate]);
            (*clip_mates)[imate][1].clip(buffers->read_length[imate], buffers->read1[imate]);
        }
    }

    return true;
}

} // namespace

void CbqStarChunk::clear() {
    backings.clear();
    reads.clear();
    mates.clear();
}

size_t CbqStarChunk::read_count() const {
    return reads.size();
}

bool append_cbq_batch_to_star_chunk(const CbqReadBatchView& batch,
                                    uint32_t record_count,
                                    uint32_t read_nends,
                                    CbqStarChunk* chunk,
                                    std::string* error) {
    if (chunk == nullptr) {
        return set_error(error, "CBQ STAR chunk destination is null");
    }
    if (batch.records == nullptr && record_count != 0) {
        return set_error(error, "CBQ STAR chunk source records are null");
    }
    if (record_count > batch.record_count) {
        std::ostringstream msg;
        msg << "CBQ STAR chunk requested " << record_count
            << " records from a batch with " << batch.record_count;
        return set_error(error, msg.str());
    }
    if (read_nends == 0) {
        return set_error(error, "CBQ STAR chunk read_nends is zero");
    }
    if (record_count == 0) {
        return true;
    }

    if (batch.backing) {
        chunk->backings.push_back(batch.backing);
    }
    chunk->reads.reserve(chunk->reads.size() + record_count);
    chunk->mates.reserve(chunk->mates.size() + static_cast<size_t>(record_count) * read_nends);

    for (uint32_t irecord = 0; irecord < record_count; ++irecord) {
        const CbqReadView& source = batch.records[irecord];
        if (source.segment_count < read_nends || source.segments == nullptr) {
            std::ostringstream msg;
            msg << "CBQ STAR chunk source record has " << source.segment_count
                << " segments but STAR expects " << read_nends;
            return set_error(error, msg.str());
        }
        if (chunk->mates.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
            return set_error(error, "CBQ STAR chunk mate offset exceeds uint32_t");
        }

        CbqStarChunkRead read;
        read.read_name = source.read_name;
        read.read_name_extra = source.read_name_extra;
        read.lane_index = source.lane_index;
        read.read_ordinal = source.read_ordinal;
        read.read_filter = source.read_filter;
        read.mate_offset = static_cast<uint32_t>(chunk->mates.size());
        read.mate_count = read_nends;
        read.has_quality = source.segments[0].has_quality;

        chunk->reads.push_back(read);
        for (uint32_t imate = 0; imate < read_nends; ++imate) {
            CbqStarChunkMate mate;
            if (!chunk_mate_from_segment(source.segments[imate], &mate, error)) {
                return false;
            }
            mate.clip_info_5p = '+';
            mate.clip_info_5p_prepared = true;
            chunk->mates.push_back(mate);
        }
    }

    return true;
}

bool prepare_cbq_star_chunk_clip_info(CbqStarChunk* chunk,
                                      uint32_t mate_index,
                                      ClipMate* clip_mate,
                                      std::string* error) {
    if (chunk == nullptr) {
        return set_error(error, "CBQ STAR chunk is null for clip-info preparation");
    }
    if (clip_mate == nullptr ||
        clip_mate->type != 10 ||
        clip_mate->cr4 == nullptr ||
        clip_mate->adSeqNum.empty()) {
        return true;
    }
    if (clip_mate->cr4->dbN <= 0) {
        return set_error(error, "CBQ STAR adapter CR4 clip database size is not positive");
    }

    std::vector<size_t> record_indices(static_cast<size_t>(clip_mate->cr4->dbN), 0);
    std::vector<char> sequence(DEF_readSeqLengthMax + 1, '\0');
    size_t irecord = 0;

    while (irecord < chunk->reads.size()) {
        int dbN1 = 0;
        while (dbN1 < clip_mate->cr4->dbN && irecord < chunk->reads.size()) {
            const CbqStarChunkRead& read = chunk->reads[irecord];
            if (mate_index >= read.mate_count ||
                read.mate_offset > chunk->mates.size() ||
                read.mate_count > chunk->mates.size() - read.mate_offset) {
                return set_error(error, "CBQ STAR chunk mate offsets are inconsistent during clip-info preparation");
            }

            CbqStarChunkMate& mate = chunk->mates[read.mate_offset + mate_index];
            mate.clip_info_5p = '+';
            mate.clip_info_5p_prepared = true;
            if (!copy_mate_sequence_to_ascii(mate, sequence.data(), error)) {
                return false;
            }

            clip_mate->cr4->opalFillOneSeq(static_cast<uint32>(dbN1),
                                           sequence.data(),
                                           mate.length);
            record_indices[static_cast<size_t>(dbN1)] = irecord;
            ++dbN1;
            ++irecord;
        }

        clip_mate->cr4->opalAlign(reinterpret_cast<uint8_t*>(clip_mate->adSeqNum.data()),
                                  clip_mate->adSeqNum.size(),
                                  dbN1);

        for (int idb = 0; idb < dbN1; ++idb) {
            const size_t clipped_record = record_indices[static_cast<size_t>(idb)];
            const CbqStarChunkRead& read = chunk->reads[clipped_record];
            CbqStarChunkMate& mate = chunk->mates[read.mate_offset + mate_index];
            const int clip_length = clip_mate->cr4->opalRes[idb].endLocationTarget + 1;
            const int score = clip_mate->cr4->opalRes[idb].score;
            const bool no_clip = score < 20 || (score == 20 && clip_length > 26) ||
                                 (score == 21 && clip_length > 30);
            mate.clip_info_5p = static_cast<char>(no_clip ? 0 : clip_length);
            mate.clip_info_5p_prepared = true;
        }
    }

    return true;
}

bool load_cbq_read_view_into_star_mates(const CbqReadView& view,
                                        const CbqStarAdapterOptions& options,
                                        CbqStarReadBuffers* buffers,
                                        std::vector<std::vector<ClipMate>>* clip_mates,
                                        std::string* error) {
    if (options.read_nends == 0) {
        return set_error(error, "CBQ STAR adapter read_nends is zero");
    }
    if (view.segment_count < options.read_nends || view.segments == nullptr) {
        std::ostringstream msg;
        msg << "CBQ STAR adapter record has " << view.segment_count
            << " segments but STAR expects " << options.read_nends;
        return set_error(error, msg.str());
    }

    CbqStarChunkRead read;
    read.read_name = view.read_name;
    read.read_name_extra = view.read_name_extra;
    read.lane_index = view.lane_index;
    read.read_ordinal = view.read_ordinal;
    read.read_filter = view.read_filter;
    read.mate_offset = 0;
    read.mate_count = options.read_nends;
    read.has_quality = view.segments[0].has_quality;

    std::vector<CbqStarChunkMate> mates(options.read_nends);
    for (uint32 imate = 0; imate < options.read_nends; ++imate) {
        if (!chunk_mate_from_segment(view.segments[imate], &mates[imate], error)) {
            return false;
        }
    }

    return load_cbq_star_read_into_buffers(read, mates.data(), options, buffers, clip_mates, error);
}

bool load_cbq_star_chunk_read_into_star_mates(const CbqStarChunk& chunk,
                                              size_t read_index,
                                              const CbqStarAdapterOptions& options,
                                              CbqStarReadBuffers* buffers,
                                              std::vector<std::vector<ClipMate>>* clip_mates,
                                              std::string* error) {
    if (read_index >= chunk.reads.size()) {
        std::ostringstream msg;
        msg << "CBQ STAR adapter read index " << read_index
            << " exceeds chunk read count " << chunk.reads.size();
        return set_error(error, msg.str());
    }
    const CbqStarChunkRead& read = chunk.reads[read_index];
    if (read.mate_offset > chunk.mates.size() ||
        read.mate_count > chunk.mates.size() - read.mate_offset) {
        return set_error(error, "CBQ STAR adapter chunk mate offsets are inconsistent");
    }
    return load_cbq_star_read_into_buffers(read,
                                           chunk.mates.data() + read.mate_offset,
                                           options,
                                           buffers,
                                           clip_mates,
                                           error);
}

} // namespace input
} // namespace star

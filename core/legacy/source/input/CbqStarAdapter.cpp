#include "input/CbqStarAdapter.h"

#include "SequenceFuns.h"

#include <algorithm>
#include <cstring>
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

} // namespace

bool load_cbq_read_view_into_star_mates(const CbqReadView& view,
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
    if (view.segment_count < options.read_nends || view.segments == nullptr) {
        std::ostringstream msg;
        msg << "CBQ STAR adapter record has " << view.segment_count
            << " segments but STAR expects " << options.read_nends;
        return set_error(error, msg.str());
    }
    if (clip_mates != nullptr && clip_mates->size() < options.read_nends) {
        return set_error(error, "CBQ STAR adapter clip_mates has fewer mates than read_nends");
    }
    if (buffers->read_name_extra->size() < options.read_nends) {
        buffers->read_name_extra->resize(options.read_nends);
    }

    const bool has_quality = view.segments[0].has_quality;
    const char header_prefix = has_quality ? '@' : '>';
    const std::string read_name_core =
        options.out_sam_read_id_number ? std::to_string(view.read_ordinal) : span_string(view.read_name);
    const std::string star_read_name = std::string(1, header_prefix) + read_name_core;

    if (!copy_string_to_cstr(star_read_name,
                             buffers->read_name_mates[0],
                             DEF_readNameLengthMax,
                             "read name",
                             error)) {
        return false;
    }
    for (uint32 imate = 1; imate < options.read_nends; ++imate) {
        if (!copy_string_to_cstr(star_read_name,
                                 buffers->read_name_mates[imate],
                                 DEF_readNameLengthMax,
                                 "mate read name",
                                 error)) {
            return false;
        }
    }

    *buffers->i_read_all = view.read_ordinal;
    *buffers->read_files_index = view.lane_index;
    *buffers->read_filter = view.read_filter;
    *buffers->read_file_type = has_quality ? 2 : 1;

    for (uint32 imate = 0; imate < options.read_nends; ++imate) {
        const CbqSegmentView& segment = view.segments[imate];
        if (segment.sequence.size > DEF_readSeqLengthMax) {
            std::ostringstream msg;
            msg << "CBQ STAR adapter mate " << (imate + 1)
                << " sequence length " << segment.sequence.size
                << " exceeds DEF_readSeqLengthMax=" << DEF_readSeqLengthMax;
            return set_error(error, msg.str());
        }
        if (segment.has_quality && segment.quality.size != segment.sequence.size) {
            std::ostringstream msg;
            msg << "CBQ STAR adapter mate " << (imate + 1)
                << " quality length " << segment.quality.size
                << " does not match sequence length " << segment.sequence.size;
            return set_error(error, msg.str());
        }

        if (!copy_span_to_cstr(segment.sequence,
                               buffers->read0[imate],
                               DEF_readSeqLengthMax + 1,
                               "sequence",
                               error)) {
            return false;
        }

        const uint length = static_cast<uint>(segment.sequence.size);
        if (segment.has_quality) {
            if (!copy_span_to_cstr(segment.quality,
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
            (*buffers->read_name_extra)[imate] = span_string(view.read_name_extra);
        } else {
            (*buffers->read_name_extra)[imate].clear();
        }

        if (clip_mates != nullptr && segment.has_quality) {
            (*clip_mates)[imate][0].clippedInfo = '+';
        }
        convertNucleotidesToNumbers(buffers->read0[imate], buffers->read1[imate], buffers->read_length[imate]);
        if (!options.trim_cutadapt_enabled && clip_mates != nullptr) {
            (*clip_mates)[imate][0].clip(buffers->read_length[imate], buffers->read1[imate]);
            (*clip_mates)[imate][1].clip(buffers->read_length[imate], buffers->read1[imate]);
        }
    }

    return true;
}

} // namespace input
} // namespace star

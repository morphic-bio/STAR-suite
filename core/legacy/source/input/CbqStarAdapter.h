#ifndef CODE_input_CbqStarAdapter
#define CODE_input_CbqStarAdapter

#include "IncludeDefine.h"
#include "ClipMate.h"
#include "input/CbqInputModule.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace star {
namespace input {

struct CbqStarAdapterOptions {
    uint32 read_nends = 0;
    bool out_sam_read_id_number = false;
    int out_qs_conversion_add = 0;
    bool trim_cutadapt_enabled = false;
    bool preserve_read_name_extra = false;
};

struct CbqStarReadBuffers {
    char** read_name_mates = nullptr;
    char** read0 = nullptr;
    char** read1 = nullptr;
    char** qual0 = nullptr;
    std::vector<std::string>* read_name_extra = nullptr;
    uint* read_length = nullptr;
    uint* read_length_original = nullptr;
    uint64* i_read_all = nullptr;
    uint32* read_files_index = nullptr;
    char* read_filter = nullptr;
    int* read_file_type = nullptr;
};

struct CbqStarChunkMate {
    CbqByteSpan sequence;
    CbqByteSpan quality;
    CbqPackedSequenceView packed_sequence;
    uint32_t length = 0;
    uint32_t original_length = 0;
    char clip_info_5p = 0;
    bool clip_info_5p_prepared = false;
    bool has_quality = true;
};

struct CbqStarChunkRead {
    CbqByteSpan read_name;
    CbqByteSpan read_name_extra;
    uint32_t lane_index = 0;
    uint64_t read_ordinal = 0;
    char read_filter = 'Y';
    uint32_t mate_offset = 0;
    uint32_t mate_count = 0;
    bool has_quality = true;
};

struct CbqStarChunk {
    std::vector<std::shared_ptr<const CbqBatchBacking>> backings;
    std::vector<CbqStarChunkRead> reads;
    std::vector<CbqStarChunkMate> mates;

    void clear();
    size_t read_count() const;
};

bool append_cbq_batch_to_star_chunk(const CbqReadBatchView& batch,
                                    uint32_t record_count,
                                    uint32_t read_nends,
                                    CbqStarChunk* chunk,
                                    std::string* error);

bool prepare_cbq_star_chunk_clip_info(CbqStarChunk* chunk,
                                      uint32_t mate_index,
                                      ClipMate* clip_mate,
                                      std::string* error);

bool load_cbq_read_view_into_star_mates(const CbqReadView& view,
                                        const CbqStarAdapterOptions& options,
                                        CbqStarReadBuffers* buffers,
                                        std::vector<std::vector<ClipMate>>* clip_mates,
                                        std::string* error);

bool load_cbq_star_chunk_read_into_star_mates(const CbqStarChunk& chunk,
                                              size_t read_index,
                                              const CbqStarAdapterOptions& options,
                                              CbqStarReadBuffers* buffers,
                                              std::vector<std::vector<ClipMate>>* clip_mates,
                                              std::string* error);

} // namespace input
} // namespace star

#endif

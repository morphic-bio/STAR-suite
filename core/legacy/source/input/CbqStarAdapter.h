#ifndef CODE_input_CbqStarAdapter
#define CODE_input_CbqStarAdapter

#include "IncludeDefine.h"
#include "ClipMate.h"
#include "input/CbqInputModule.h"

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

bool load_cbq_read_view_into_star_mates(const CbqReadView& view,
                                        const CbqStarAdapterOptions& options,
                                        CbqStarReadBuffers* buffers,
                                        std::vector<std::vector<ClipMate>>* clip_mates,
                                        std::string* error);

} // namespace input
} // namespace star

#endif

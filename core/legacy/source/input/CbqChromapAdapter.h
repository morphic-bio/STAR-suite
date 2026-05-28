#ifndef CODE_input_CbqChromapAdapter
#define CODE_input_CbqChromapAdapter

#include <string>

namespace star {
namespace input {

struct CbqChromapAdapterOptions {
    std::string read_pair_cbq;
    std::string barcode_cbq;
    std::string output_dir;
    std::string sample_name;
    bool require_name_match = true;
    bool allow_missing_qualities = false;
};

struct CbqChromapFastqPaths {
    std::string read1_fastq;
    std::string read2_fastq;
    std::string barcode_fastq;
};

bool materialize_cbq_chromap_fastqs(const CbqChromapAdapterOptions& options,
                                    CbqChromapFastqPaths* paths,
                                    std::string* error);

} // namespace input
} // namespace star

#endif

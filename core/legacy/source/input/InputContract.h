#ifndef CODE_input_InputContract
#define CODE_input_InputContract

#include <cstdint>
#include <string>
#include <vector>

namespace star {
namespace input {

enum class SourceFormat {
    Unknown,
    Fastx,
    Sam,
    Binseq
};

enum class InputStatus {
    Record,
    End,
    Error
};

struct InputMateRecord {
    std::string sequence;
    std::string quality;
    uint32_t original_length = 0;
    bool has_quality = true;
};

// Input modules emit records in module-defined order. The contract preserves
// record content and mate synchronization, but does not require this order to
// match the original source file order unless the source plan advertises that
// stronger guarantee.
struct InputRecord {
    std::string read_name;
    std::string read_name_extra;
    uint32_t lane_index = 0;
    // Ordinal in the module-emitted stream. This is not a
    // source-file ordinal unless preserves_source_order is true for the plan.
    uint64_t read_ordinal = 0;
    char read_filter = 'Y';
    uint32_t mate_count = 0;
    std::vector<InputMateRecord> mates;
};

struct InputSourcePlan {
    SourceFormat format = SourceFormat::Unknown;
    std::string module_name;
    std::vector<std::vector<std::string>> mate_files;
    std::vector<std::string> read_groups;
    std::vector<char> read_name_separator_chars;
    uint32_t read_files_n = 0;
    uint32_t mate_count = 0;
    bool preserves_source_order = false;
    bool uses_helper_stream = false;
    bool uses_internal_gzip = false;
    int out_qs_conversion_add = 0;
    std::string command_string;
    std::string prefix;
};

class InputModule {
public:
    virtual ~InputModule() {}

    virtual const char* name() const = 0;
    virtual bool configure(const InputSourcePlan& plan, std::string* error) = 0;
    virtual const InputSourcePlan& plan() const = 0;

    virtual bool supports_record_iteration() const {
        return false;
    }

    virtual bool open(std::string* error) {
        if (error != nullptr) {
            *error = "record iteration is not implemented by this input module";
        }
        return false;
    }

    virtual InputStatus next_record(InputRecord*, std::string* error) {
        if (error != nullptr) {
            *error = "record iteration is not implemented by this input module";
        }
        return InputStatus::Error;
    }

    virtual void close() {}
};

} // namespace input
} // namespace star

#endif

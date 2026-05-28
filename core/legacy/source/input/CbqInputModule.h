#ifndef CODE_input_CbqInputModule
#define CODE_input_CbqInputModule

#include "input/InputContract.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace star {
namespace input {

InputSourcePlan make_cbq_input_source_plan(
    const std::vector<std::vector<std::string>>& read_files_names,
    const std::vector<std::string>& read_groups,
    uint32_t mate_count);

struct CbqByteSpan {
    const char* data = nullptr;
    size_t size = 0;
};

enum class CbqSegmentRole {
    Unknown,
    Read1,
    Read2,
    Barcode,
    Feature,
    Reverse
};

struct CbqSegmentView {
    CbqSegmentRole role = CbqSegmentRole::Unknown;
    uint32_t source_index = 0;
    CbqByteSpan sequence;
    CbqByteSpan quality;
    uint32_t original_length = 0;
    bool has_quality = true;
};

struct CbqReadView {
    CbqByteSpan read_name;
    CbqByteSpan read_name_extra;
    uint32_t lane_index = 0;
    uint64_t read_ordinal = 0;
    uint64_t lane_read_ordinal = 0;
    char read_filter = 'Y';
    uint32_t segment_count = 0;
    const CbqSegmentView* segments = nullptr;
};

// Borrowed view over one decoded CBQ block. Pointers remain valid until the
// next CbqInputModule call that can load another block, or close().
struct CbqReadBatchView {
    const CbqReadView* records = nullptr;
    uint32_t record_count = 0;
    uint32_t lane_index = 0;
    bool preserves_source_order = true;
    bool backing_storage_owned_by_reader = true;
};

class CbqInputModule : public InputModule {
public:
    CbqInputModule();
    CbqInputModule(const CbqInputModule&) = delete;
    CbqInputModule& operator=(const CbqInputModule&) = delete;
    ~CbqInputModule() override;

    const char* name() const override;
    bool configure(const InputSourcePlan& plan, std::string* error) override;
    const InputSourcePlan& plan() const override;
    bool supports_record_iteration() const override;
    bool open(std::string* error) override;
    InputStatus next_batch(CbqReadBatchView* batch, std::string* error);
    InputStatus next_record(InputRecord* record, std::string* error) override;
    void close() override;

private:
    struct Impl;

    InputSourcePlan plan_;
    bool configured_;
    std::unique_ptr<Impl> impl_;
};

} // namespace input
} // namespace star

#endif

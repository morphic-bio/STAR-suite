#ifndef CODE_input_FastxInputModule
#define CODE_input_FastxInputModule

#include "input/InputContract.h"

#include <memory>
#include <string>
#include <vector>

namespace star {
namespace input {

InputSourcePlan make_fastx_input_source_plan(
    const std::vector<std::vector<std::string>>& read_files_names,
    const std::vector<std::string>& read_groups,
    const std::string& command_string,
    const std::string& prefix,
    bool uses_internal_gzip);

class FastxInputModule : public InputModule {
public:
    FastxInputModule();
    ~FastxInputModule() override;

    const char* name() const override;
    bool configure(const InputSourcePlan& plan, std::string* error) override;
    const InputSourcePlan& plan() const override;
    bool supports_record_iteration() const override;
    bool open(std::string* error) override;
    InputStatus next_record(InputRecord* record, std::string* error) override;
    void close() override;

private:
    struct Impl;

    bool open_lane(uint32_t lane_index, std::string* error);
    void close_lane();

    InputSourcePlan plan_;
    bool configured_;
    std::unique_ptr<Impl> impl_;
};

} // namespace input
} // namespace star

#endif

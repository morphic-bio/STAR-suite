#ifndef CODE_input_BinseqInputProbeModule
#define CODE_input_BinseqInputProbeModule

#include "input/InputContract.h"
#include "input/FastxInputModule.h"

#include <memory>
#include <string>

namespace star {
namespace input {

class BinseqInputProbeModule : public InputModule {
public:
    BinseqInputProbeModule();
    BinseqInputProbeModule(const BinseqInputProbeModule&) = delete;
    BinseqInputProbeModule& operator=(const BinseqInputProbeModule&) = delete;
    ~BinseqInputProbeModule() override;

    const char* name() const override;
    bool configure(const InputSourcePlan& plan, std::string* error) override;
    const InputSourcePlan& plan() const override;
    bool supports_record_iteration() const override;
    bool open(std::string* error) override;
    InputStatus next_record(InputRecord* record, std::string* error) override;
    void close() override;

private:
    bool decode_lane(uint32_t lane_index,
                     std::vector<std::vector<std::string>>& decoded_mate_files,
                     std::string* error);

    InputSourcePlan plan_;
    bool configured_;
    std::unique_ptr<FastxInputModule> fastx_;
};

} // namespace input
} // namespace star

#endif

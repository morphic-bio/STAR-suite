#ifndef CODE_input_CbqInputModule
#define CODE_input_CbqInputModule

#include "input/InputContract.h"

#include <memory>
#include <string>

namespace star {
namespace input {

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

#ifndef CODE_input_FastxProducerPool
#define CODE_input_FastxProducerPool

#include "input/InputContract.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace star {
namespace input {

// Concurrent FASTX source with deterministic lane-major delivery.
//
// Producers independently decompress and parse lanes into bounded batches.
// Consumers still observe exactly the InputModule record stream: every record
// from lane 0, then lane 1, and so on.  The separation keeps decompression and
// FASTQ parsing out of STAR's global input-chunk mutex without coupling the
// input layer to a particular mapping or feature pipeline.
class FastxProducerPool : public InputModule {
public:
    explicit FastxProducerPool(uint32_t producer_count = 0,
                               size_t records_per_batch = 2048,
                               size_t queued_batches_per_lane = 4);
    ~FastxProducerPool() override;

    FastxProducerPool(const FastxProducerPool&) = delete;
    FastxProducerPool& operator=(const FastxProducerPool&) = delete;

    const char* name() const override;
    bool configure(const InputSourcePlan& plan, std::string* error) override;
    const InputSourcePlan& plan() const override;
    bool supports_record_iteration() const override;
    bool open(std::string* error) override;
    InputStatus next_record(InputRecord* record, std::string* error) override;
    void close() override;

    uint32_t producer_count() const;
    size_t records_per_batch() const;
    size_t queued_batches_per_lane() const;

private:
    struct Impl;

    InputSourcePlan plan_;
    uint32_t requested_producer_count_;
    size_t records_per_batch_;
    size_t queued_batches_per_lane_;
    bool configured_;
    std::unique_ptr<Impl> impl_;
};

} // namespace input
} // namespace star

#endif

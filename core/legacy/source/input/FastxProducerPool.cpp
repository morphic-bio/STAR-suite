#include "input/FastxProducerPool.h"

#include "input/FastxInputModule.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <exception>
#include <mutex>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>

namespace star {
namespace input {
namespace {

bool set_error(std::string* error, const std::string& message) {
    if (error != nullptr) {
        *error = message;
    }
    return false;
}

struct RecordBatch {
    std::vector<InputRecord> records;
    size_t next = 0;
};

struct LaneState {
    std::deque<std::unique_ptr<RecordBatch>> ready;
    bool done = false;
};

InputSourcePlan one_lane_plan(const InputSourcePlan& plan, uint32_t lane_index) {
    InputSourcePlan lane_plan = plan;
    lane_plan.mate_files.clear();
    lane_plan.mate_files.resize(plan.mate_files.size());
    for (size_t imate = 0; imate < plan.mate_files.size(); ++imate) {
        lane_plan.mate_files[imate].push_back(plan.mate_files[imate][lane_index]);
    }
    lane_plan.read_groups.clear();
    if (!plan.read_groups.empty()) {
        lane_plan.read_groups.push_back(plan.read_groups[lane_index]);
    }
    lane_plan.read_files_n = 1;
    return lane_plan;
}

} // namespace

struct FastxProducerPool::Impl {
    std::mutex mutex;
    std::condition_variable data_ready;
    std::condition_variable queue_space;
    std::vector<std::unique_ptr<LaneState>> lanes;
    std::vector<std::thread> producers;
    std::atomic<uint32_t> next_lane{0};
    uint32_t delivery_lane = 0;
    uint64_t delivered_records = 0;
    uint32_t active_producer_count = 0;
    bool opened = false;
    bool stopping = false;
    bool failed = false;
    std::string failure;
};

FastxProducerPool::FastxProducerPool(uint32_t producer_count,
                                     size_t records_per_batch,
                                     size_t queued_batches_per_lane)
    : requested_producer_count_(producer_count),
      records_per_batch_(records_per_batch),
      queued_batches_per_lane_(queued_batches_per_lane),
      configured_(false),
      impl_(new Impl()) {
    plan_.format = SourceFormat::Fastx;
    plan_.module_name = "FastxProducerPool";
}

FastxProducerPool::~FastxProducerPool() {
    close();
}

const char* FastxProducerPool::name() const {
    return "FastxProducerPool";
}

bool FastxProducerPool::configure(const InputSourcePlan& plan, std::string* error) {
    close();
    if (records_per_batch_ == 0) {
        return set_error(error, "Fastx producer pool records_per_batch must be greater than zero");
    }
    if (queued_batches_per_lane_ == 0) {
        return set_error(error, "Fastx producer pool queued_batches_per_lane must be greater than zero");
    }

    // Use the sequential reader as the single source of validation rules so
    // both modes accept and reject identical source plans.
    FastxInputModule validator;
    if (!validator.configure(plan, error)) {
        return false;
    }

    plan_ = validator.plan();
    plan_.module_name = "FastxProducerPool";
    plan_.preserves_source_order = true;
    configured_ = true;
    return true;
}

const InputSourcePlan& FastxProducerPool::plan() const {
    return plan_;
}

bool FastxProducerPool::supports_record_iteration() const {
    return true;
}

bool FastxProducerPool::open(std::string* error) {
    if (!configured_) {
        return set_error(error, "FastxProducerPool must be configured before open()");
    }
    close();

    {
        std::lock_guard<std::mutex> lock(impl_->mutex);
        impl_->lanes.clear();
        impl_->lanes.reserve(plan_.read_files_n);
        for (uint32_t lane = 0; lane < plan_.read_files_n; ++lane) {
            impl_->lanes.emplace_back(new LaneState());
        }
        impl_->next_lane.store(0);
        impl_->delivery_lane = 0;
        impl_->delivered_records = 0;
        impl_->stopping = false;
        impl_->failed = false;
        impl_->failure.clear();
        impl_->opened = true;
        const uint32_t automatic_count = std::min<uint32_t>(4, plan_.read_files_n);
        impl_->active_producer_count = std::max<uint32_t>(
            1, std::min<uint32_t>(
                   requested_producer_count_ == 0 ? automatic_count : requested_producer_count_,
                   plan_.read_files_n));
    }

    const uint32_t producer_count_local = impl_->active_producer_count;
    try {
        for (uint32_t producer = 0; producer < producer_count_local; ++producer) {
            impl_->producers.emplace_back([this]() {
                try {
                for (;;) {
                    const uint32_t lane_index = impl_->next_lane.fetch_add(1);
                    if (lane_index >= plan_.read_files_n) {
                        return;
                    }

                    bool stopping = false;
                    {
                        std::lock_guard<std::mutex> lock(impl_->mutex);
                        stopping = impl_->stopping;
                    }
                    if (stopping) {
                        return;
                    }

                    FastxInputModule reader;
                    std::string reader_error;
                    const InputSourcePlan lane_plan = one_lane_plan(plan_, lane_index);
                    if (!reader.configure(lane_plan, &reader_error) ||
                        !reader.open(&reader_error)) {
                        std::lock_guard<std::mutex> lock(impl_->mutex);
                        if (!impl_->failed) {
                            std::ostringstream message;
                            message << "FASTX lane " << lane_index << ": " << reader_error;
                            impl_->failure = message.str();
                            impl_->failed = true;
                        }
                        impl_->stopping = true;
                        impl_->data_ready.notify_all();
                        impl_->queue_space.notify_all();
                        return;
                    }

                    bool lane_complete = false;
                    while (!lane_complete) {
                        std::unique_ptr<RecordBatch> batch(new RecordBatch());
                        batch->records.reserve(records_per_batch_);
                        for (size_t ii = 0; ii < records_per_batch_; ++ii) {
                            stopping = false;
                            {
                                std::lock_guard<std::mutex> lock(impl_->mutex);
                                stopping = impl_->stopping;
                            }
                            if (stopping) {
                                reader.close();
                                return;
                            }

                            InputRecord record;
                            reader_error.clear();
                            const InputStatus status = reader.next_record(&record, &reader_error);
                            if (status == InputStatus::Error) {
                                {
                                    std::lock_guard<std::mutex> lock(impl_->mutex);
                                    if (!impl_->failed) {
                                        std::ostringstream message;
                                        message << "FASTX lane " << lane_index << ": " << reader_error;
                                        impl_->failure = message.str();
                                        impl_->failed = true;
                                    }
                                    impl_->stopping = true;
                                    impl_->data_ready.notify_all();
                                    impl_->queue_space.notify_all();
                                }
                                reader.close();
                                return;
                            }
                            if (status == InputStatus::End) {
                                lane_complete = true;
                                break;
                            }
                            record.lane_index = lane_index;
                            batch->records.push_back(std::move(record));
                        }

                        if (!batch->records.empty()) {
                            std::unique_lock<std::mutex> lock(impl_->mutex);
                            LaneState& lane = *impl_->lanes[lane_index];
                            impl_->queue_space.wait(lock, [this, &lane]() {
                                return impl_->stopping ||
                                       lane.ready.size() < queued_batches_per_lane_;
                            });
                            if (impl_->stopping) {
                                lock.unlock();
                                reader.close();
                                return;
                            }
                            lane.ready.push_back(std::move(batch));
                            impl_->data_ready.notify_all();
                        }
                    }

                    reader.close();
                    {
                        std::lock_guard<std::mutex> lock(impl_->mutex);
                        impl_->lanes[lane_index]->done = true;
                        impl_->data_ready.notify_all();
                    }
                }
                } catch (const std::exception& ex) {
                    std::lock_guard<std::mutex> lock(impl_->mutex);
                    if (!impl_->failed) {
                        impl_->failure = std::string("FASTX producer exception: ") + ex.what();
                        impl_->failed = true;
                    }
                    impl_->stopping = true;
                    impl_->data_ready.notify_all();
                    impl_->queue_space.notify_all();
                } catch (...) {
                    std::lock_guard<std::mutex> lock(impl_->mutex);
                    if (!impl_->failed) {
                        impl_->failure = "FASTX producer failed with an unknown exception";
                        impl_->failed = true;
                    }
                    impl_->stopping = true;
                    impl_->data_ready.notify_all();
                    impl_->queue_space.notify_all();
                }
            });
        }
    } catch (const std::exception& ex) {
        close();
        return set_error(error, std::string("could not start FASTX producer thread: ") + ex.what());
    }

    return true;
}

InputStatus FastxProducerPool::next_record(InputRecord* record, std::string* error) {
    if (error != nullptr) {
        error->clear();
    }
    if (record == nullptr) {
        set_error(error, "next_record() requires a non-null InputRecord");
        return InputStatus::Error;
    }

    std::unique_lock<std::mutex> lock(impl_->mutex);
    if (!impl_->opened) {
        set_error(error, "FastxProducerPool is not open");
        return InputStatus::Error;
    }

    for (;;) {
        if (impl_->failed) {
            set_error(error, impl_->failure);
            return InputStatus::Error;
        }
        if (impl_->delivery_lane >= impl_->lanes.size()) {
            return InputStatus::End;
        }

        LaneState& lane = *impl_->lanes[impl_->delivery_lane];
        if (!lane.ready.empty()) {
            RecordBatch& batch = *lane.ready.front();
            *record = std::move(batch.records[batch.next++]);
            record->lane_index = impl_->delivery_lane;
            record->read_ordinal = ++impl_->delivered_records;
            if (batch.next >= batch.records.size()) {
                lane.ready.pop_front();
                impl_->queue_space.notify_all();
            }
            return InputStatus::Record;
        }
        if (lane.done) {
            ++impl_->delivery_lane;
            continue;
        }
        if (impl_->stopping) {
            return InputStatus::End;
        }
        impl_->data_ready.wait(lock);
    }
}

void FastxProducerPool::close() {
    if (!impl_) {
        return;
    }
    {
        std::lock_guard<std::mutex> lock(impl_->mutex);
        impl_->stopping = true;
        impl_->data_ready.notify_all();
        impl_->queue_space.notify_all();
    }
    for (std::thread& producer : impl_->producers) {
        if (producer.joinable()) {
            producer.join();
        }
    }
    impl_->producers.clear();
    {
        std::lock_guard<std::mutex> lock(impl_->mutex);
        impl_->lanes.clear();
        impl_->opened = false;
        impl_->active_producer_count = 0;
    }
}

uint32_t FastxProducerPool::producer_count() const {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    return impl_->active_producer_count;
}

size_t FastxProducerPool::records_per_batch() const {
    return records_per_batch_;
}

size_t FastxProducerPool::queued_batches_per_lane() const {
    return queued_batches_per_lane_;
}

} // namespace input
} // namespace star

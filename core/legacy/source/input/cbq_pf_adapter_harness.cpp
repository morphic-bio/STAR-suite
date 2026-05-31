#include "input/CbqInputModule.h"

#include "pf_api.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr size_t kPfLineLength = 1024;
constexpr size_t kPfSequenceCapacity = kPfLineLength - 1;

struct Args {
    std::string mode;
    std::string cbq;
    std::string barcode_fastq;
    std::string feature_fastq;
    std::string whitelist;
    std::string feature_ref;
    std::string output_dir;
    std::string sample_name = "sample";
    int barcode_length = 16;
    int umi_length = 12;
    int threads = 1;
    int search_threads = 1;
    int consumer_threads = 1;
    int read_buffer_lines = -1;
    int feature_offset = -1;
    std::string materialize_mode = "none";
    bool skip_emptydrops = true;
    bool skip_heatmaps = true;
};

void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " --mode fastq|cbq|cbq-direct-decode|cbq-direct-pf --whitelist WL --featureRef CSV "
        << "--outputDir OUT [--sampleName NAME]\n"
        << "  FASTQ mode: --barcodeFastq R1 --featureFastq R2\n"
        << "  CBQ mode:   --readFilesIn paired.cbq\n"
        << "  Direct CBQ decode mode: --readFilesIn paired.cbq [--consumerThreads N] "
        << "[--materializeMode none|sequence]\n"
        << "  Optional:   --readBufferLines N [--featureOffset N]\n";
}

bool parse_args(int argc, char** argv, Args* args) {
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto need_value = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                return nullptr;
            }
            return argv[++i];
        };

        if (key == "--mode") {
            const char* value = need_value("--mode");
            if (!value) return false;
            args->mode = value;
        } else if (key == "--readFilesIn") {
            const char* value = need_value("--readFilesIn");
            if (!value) return false;
            args->cbq = value;
        } else if (key == "--barcodeFastq") {
            const char* value = need_value("--barcodeFastq");
            if (!value) return false;
            args->barcode_fastq = value;
        } else if (key == "--featureFastq") {
            const char* value = need_value("--featureFastq");
            if (!value) return false;
            args->feature_fastq = value;
        } else if (key == "--whitelist") {
            const char* value = need_value("--whitelist");
            if (!value) return false;
            args->whitelist = value;
        } else if (key == "--featureRef") {
            const char* value = need_value("--featureRef");
            if (!value) return false;
            args->feature_ref = value;
        } else if (key == "--outputDir") {
            const char* value = need_value("--outputDir");
            if (!value) return false;
            args->output_dir = value;
        } else if (key == "--sampleName") {
            const char* value = need_value("--sampleName");
            if (!value) return false;
            args->sample_name = value;
        } else if (key == "--barcodeLength") {
            const char* value = need_value("--barcodeLength");
            if (!value) return false;
            args->barcode_length = std::atoi(value);
        } else if (key == "--umiLength") {
            const char* value = need_value("--umiLength");
            if (!value) return false;
            args->umi_length = std::atoi(value);
        } else if (key == "--threads") {
            const char* value = need_value("--threads");
            if (!value) return false;
            args->threads = std::atoi(value);
        } else if (key == "--searchThreads") {
            const char* value = need_value("--searchThreads");
            if (!value) return false;
            args->search_threads = std::atoi(value);
        } else if (key == "--consumerThreads") {
            const char* value = need_value("--consumerThreads");
            if (!value) return false;
            args->consumer_threads = std::atoi(value);
        } else if (key == "--readBufferLines") {
            const char* value = need_value("--readBufferLines");
            if (!value) return false;
            args->read_buffer_lines = std::atoi(value);
        } else if (key == "--featureOffset") {
            const char* value = need_value("--featureOffset");
            if (!value) return false;
            args->feature_offset = std::atoi(value);
        } else if (key == "--materializeMode") {
            const char* value = need_value("--materializeMode");
            if (!value) return false;
            args->materialize_mode = value;
        } else if (key == "--runEmptyDrops") {
            args->skip_emptydrops = false;
        } else if (key == "--runHeatmaps") {
            args->skip_heatmaps = false;
        } else if (key == "--help" || key == "-h") {
            return false;
        } else {
            std::cerr << "Unknown argument: " << key << "\n";
            return false;
        }
    }

    if (args->mode != "fastq" && args->mode != "cbq" &&
        args->mode != "cbq-direct-decode" && args->mode != "cbq-direct-pf") {
        std::cerr << "--mode must be fastq, cbq, cbq-direct-decode, or cbq-direct-pf\n";
        return false;
    }
    if (args->consumer_threads <= 0) {
        std::cerr << "--consumerThreads must be positive\n";
        return false;
    }
    if (args->mode == "cbq-direct-decode") {
        if (args->cbq.empty()) {
            std::cerr << "CBQ direct decode mode requires --readFilesIn\n";
            return false;
        }
        if (args->materialize_mode != "none" && args->materialize_mode != "sequence") {
            std::cerr << "--materializeMode must be none or sequence\n";
            return false;
        }
        return true;
    }
    if (args->whitelist.empty() || args->feature_ref.empty() || args->output_dir.empty()) {
        std::cerr << "--whitelist, --featureRef, and --outputDir are required\n";
        return false;
    }
    if (args->mode == "fastq" &&
        (args->barcode_fastq.empty() || args->feature_fastq.empty())) {
        std::cerr << "FASTQ mode requires --barcodeFastq and --featureFastq\n";
        return false;
    }
    if (args->mode == "cbq" && args->cbq.empty()) {
        std::cerr << "CBQ mode requires --readFilesIn\n";
        return false;
    }
    if (args->mode == "cbq-direct-pf" && args->cbq.empty()) {
        std::cerr << "CBQ direct PF mode requires --readFilesIn\n";
        return false;
    }
    return true;
}

pf_sequence_view to_pf_view(star::input::CbqByteSpan span) {
    pf_sequence_view view;
    view.data = span.data;
    view.length = span.size;
    return view;
}

pf_sequence_view to_pf_view_buffer(const char* data, size_t length) {
    pf_sequence_view view;
    view.data = data;
    view.length = length;
    return view;
}

star::input::InputSourcePlan make_cbq_pair_plan(const Args& args) {
    std::vector<std::vector<std::string>> read_files(1);
    read_files[0].push_back(args.cbq);
    star::input::InputSourcePlan plan =
        star::input::make_cbq_input_source_plan(read_files, std::vector<std::string>(), 2);
    plan.read_name_separator_chars.clear();
    plan.read_name_separator_chars.push_back(' ');
    return plan;
}

struct DirectDecodeStats {
    uint64_t records = 0;
    uint64_t batches = 0;
    uint64_t bases = 0;
    uint64_t read_ordinal_sum = 0;
    uint64_t lane_ordinal_sum = 0;
};

bool process_cbq_direct_range(const Args& args,
                              const star::input::InputSourcePlan& plan,
                              uint64_t first_record,
                              uint64_t record_count,
                              DirectDecodeStats* stats,
                              std::string* error) {
    star::input::CbqInputModule module;
    if (!module.configure(plan, error) ||
        !module.open_range(0, first_record, record_count, error)) {
        return false;
    }

    std::vector<char> sequence_storage;
    for (;;) {
        star::input::CbqReadBatchView batch;
        const star::input::InputStatus status = module.next_batch(&batch, error);
        if (status == star::input::InputStatus::Error) {
            module.close();
            return false;
        }
        if (status == star::input::InputStatus::End) {
            break;
        }

        ++stats->batches;
        stats->records += batch.record_count;
        if (args.materialize_mode == "sequence") {
            const size_t required_sequence_storage =
                static_cast<size_t>(batch.record_count) * 2U * kPfLineLength;
            if (sequence_storage.size() < required_sequence_storage) {
                sequence_storage.resize(required_sequence_storage);
            }
            for (uint32_t i = 0; i < batch.record_count; ++i) {
                const star::input::CbqReadView& view = batch.records[i];
                if (view.segment_count < 2 || view.segments == nullptr) {
                    if (error) *error = "CBQ direct decode expects paired CBQ records";
                    module.close();
                    return false;
                }
                stats->read_ordinal_sum += view.read_ordinal;
                stats->lane_ordinal_sum += view.lane_read_ordinal;

                char* barcode_sequence = sequence_storage.data() +
                    (static_cast<size_t>(i) * 2U * kPfLineLength);
                char* feature_sequence = barcode_sequence + kPfLineLength;
                size_t barcode_length = 0;
                size_t feature_length = 0;
                if (!star::input::materialize_cbq_segment_sequence_to_buffer(view.segments[0],
                                                                             barcode_sequence,
                                                                             kPfSequenceCapacity,
                                                                             &barcode_length,
                                                                             error) ||
                    !star::input::materialize_cbq_segment_sequence_to_buffer(view.segments[1],
                                                                             feature_sequence,
                                                                             kPfSequenceCapacity,
                                                                             &feature_length,
                                                                             error)) {
                    module.close();
                    return false;
                }
                stats->bases += static_cast<uint64_t>(barcode_length + feature_length);
            }
        }
    }

    module.close();
    return true;
}

int run_cbq_direct_decode(const Args& args) {
    const star::input::InputSourcePlan plan = make_cbq_pair_plan(args);
    std::string error;
    star::input::CbqInputModule metadata_module;
    if (!metadata_module.configure(plan, &error) ||
        !metadata_module.open_range(0, 0, std::numeric_limits<uint64_t>::max(), &error)) {
        std::cerr << "Failed to inspect CBQ index: " << error << "\n";
        return 1;
    }
    const uint64_t total_records = metadata_module.current_lane_record_count();
    metadata_module.close();

    const int requested_threads = args.consumer_threads > 0 ? args.consumer_threads : 1;
    const uint64_t chunk_size = total_records == 0
        ? 0
        : (total_records + static_cast<uint64_t>(requested_threads) - 1U) /
              static_cast<uint64_t>(requested_threads);

    struct Range {
        uint64_t first = 0;
        uint64_t count = 0;
    };
    std::vector<Range> ranges;
    for (int ithread = 0; ithread < requested_threads && chunk_size != 0; ++ithread) {
        const uint64_t first = static_cast<uint64_t>(ithread) * chunk_size;
        if (first >= total_records) {
            break;
        }
        Range range;
        range.first = first;
        range.count = std::min(chunk_size, total_records - first);
        ranges.push_back(range);
    }

    std::vector<DirectDecodeStats> per_thread(ranges.size());
    std::vector<std::string> errors(ranges.size());
    std::vector<int> ok(ranges.size(), 0);
    std::vector<std::thread> workers;
    workers.reserve(ranges.size());

    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for (size_t irange = 0; irange < ranges.size(); ++irange) {
        workers.push_back(std::thread([&, irange]() {
            ok[irange] = process_cbq_direct_range(args,
                                                  plan,
                                                  ranges[irange].first,
                                                  ranges[irange].count,
                                                  &per_thread[irange],
                                                  &errors[irange]) ? 1 : 0;
        }));
    }
    for (size_t iworker = 0; iworker < workers.size(); ++iworker) {
        workers[iworker].join();
    }
    const std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();

    for (size_t irange = 0; irange < ok.size(); ++irange) {
        if (!ok[irange]) {
            std::cerr << "CBQ direct decode worker failed";
            if (!errors[irange].empty()) {
                std::cerr << ": " << errors[irange];
            }
            std::cerr << "\n";
            return 1;
        }
    }

    DirectDecodeStats total;
    for (size_t irange = 0; irange < per_thread.size(); ++irange) {
        total.records += per_thread[irange].records;
        total.batches += per_thread[irange].batches;
        total.bases += per_thread[irange].bases;
        total.read_ordinal_sum += per_thread[irange].read_ordinal_sum;
        total.lane_ordinal_sum += per_thread[irange].lane_ordinal_sum;
    }
    const std::chrono::duration<double> elapsed = finish - start;
    std::cout << "processed_reads\t" << total.records << "\n";
    std::cout << "indexed_reads\t" << total_records << "\n";
    std::cout << "batches\t" << total.batches << "\n";
    std::cout << "bases_materialized\t" << total.bases << "\n";
    std::cout << "read_ordinal_sum\t" << total.read_ordinal_sum << "\n";
    std::cout << "lane_ordinal_sum\t" << total.lane_ordinal_sum << "\n";
    std::cout << "reader_threads_requested\t" << requested_threads << "\n";
    std::cout << "reader_threads_used\t" << ranges.size() << "\n";
    std::cout << "materialize_mode\t" << args.materialize_mode << "\n";
    std::cout << "elapsed_seconds\t" << elapsed.count() << "\n";
    return total.records == total_records ? 0 : 1;
}

struct CbqRange {
    uint64_t first = 0;
    uint64_t count = 0;
};

std::vector<CbqRange> make_ranges(uint64_t total_records, int requested_threads) {
    std::vector<CbqRange> ranges;
    const int threads = requested_threads > 0 ? requested_threads : 1;
    const uint64_t chunk_size = total_records == 0
        ? 0
        : (total_records + static_cast<uint64_t>(threads) - 1U) /
              static_cast<uint64_t>(threads);
    for (int ithread = 0; ithread < threads && chunk_size != 0; ++ithread) {
        const uint64_t first = static_cast<uint64_t>(ithread) * chunk_size;
        if (first >= total_records) {
            break;
        }
        CbqRange range;
        range.first = first;
        range.count = std::min(chunk_size, total_records - first);
        ranges.push_back(range);
    }
    return ranges;
}

bool inspect_cbq_records(const star::input::InputSourcePlan& plan,
                         uint64_t* total_records,
                         std::string* error) {
    star::input::CbqInputModule metadata_module;
    if (!metadata_module.configure(plan, error) ||
        !metadata_module.open_range(0, 0, std::numeric_limits<uint64_t>::max(), error)) {
        return false;
    }
    *total_records = metadata_module.current_lane_record_count();
    metadata_module.close();
    return true;
}

bool process_cbq_direct_pf_range(pf_context* ctx,
                                 pf_direct_range_job* job,
                                 const star::input::InputSourcePlan& plan,
                                 int worker_id,
                                 uint64_t first_record,
                                 uint64_t record_count,
                                 std::string* error) {
    star::input::CbqInputModule module;
    if (!module.configure(plan, error) ||
        !module.open_range(0, first_record, record_count, error)) {
        return false;
    }

    std::vector<pf_read_record_view> records;
    std::vector<char> sequence_storage;
    for (;;) {
        star::input::CbqReadBatchView batch;
        const star::input::InputStatus status = module.next_batch(&batch, error);
        if (status == star::input::InputStatus::Error) {
            module.close();
            return false;
        }
        if (status == star::input::InputStatus::End) {
            break;
        }

        records.clear();
        records.resize(batch.record_count);
        const size_t required_sequence_storage =
            static_cast<size_t>(batch.record_count) * 2U * kPfLineLength;
        if (sequence_storage.size() < required_sequence_storage) {
            sequence_storage.resize(required_sequence_storage);
        }
        for (uint32_t i = 0; i < batch.record_count; ++i) {
            const star::input::CbqReadView& view = batch.records[i];
            if (view.segment_count < 2 || view.segments == nullptr) {
                if (error) *error = "CBQ direct PF expects paired CBQ records";
                module.close();
                return false;
            }

            char* barcode_sequence = sequence_storage.data() +
                (static_cast<size_t>(i) * 2U * kPfLineLength);
            char* feature_sequence = barcode_sequence + kPfLineLength;
            size_t barcode_length = 0;
            size_t feature_length = 0;
            if (!star::input::materialize_cbq_segment_sequence_to_buffer(view.segments[0],
                                                                         barcode_sequence,
                                                                         kPfSequenceCapacity,
                                                                         &barcode_length,
                                                                         error) ||
                !star::input::materialize_cbq_segment_sequence_to_buffer(view.segments[1],
                                                                         feature_sequence,
                                                                         kPfSequenceCapacity,
                                                                         &feature_length,
                                                                         error)) {
                module.close();
                return false;
            }

            pf_read_record_view rec = {};
            rec.barcode_sequence = to_pf_view_buffer(barcode_sequence, barcode_length);
            rec.barcode_quality = to_pf_view(view.segments[0].quality);
            rec.feature_sequence = to_pf_view_buffer(feature_sequence, feature_length);
            rec.feature_quality = to_pf_view(view.segments[1].quality);
            records[i] = rec;
        }

        if (!records.empty()) {
            const pf_error err = pf_direct_range_process_record_views(
                job, worker_id, records.data(), batch.record_count);
            if (err != PF_OK) {
                if (error) {
                    const char* pf_error_message = pf_get_error(ctx);
                    *error = pf_error_message ? pf_error_message : "process_features direct CBQ batch failed";
                }
                module.close();
                return false;
            }
        }
    }

    module.close();
    return true;
}

bool process_cbq_direct_pf_records(pf_context* ctx,
                                   const Args& args,
                                   pf_stats* stats,
                                   std::string* error) {
    const star::input::InputSourcePlan plan = make_cbq_pair_plan(args);
    uint64_t total_records = 0;
    if (!inspect_cbq_records(plan, &total_records, error)) {
        return false;
    }

    const int requested_threads = args.consumer_threads > 0 ? args.consumer_threads : 1;
    const std::vector<CbqRange> ranges = make_ranges(total_records, requested_threads);

    pf_direct_range_job* job = nullptr;
    pf_error err = pf_direct_range_begin(ctx,
                                         args.output_dir.c_str(),
                                         args.sample_name.c_str(),
                                         requested_threads,
                                         2,
                                         &job);
    if (err != PF_OK) {
        if (error) {
            const char* pf_error_message = pf_get_error(ctx);
            *error = pf_error_message ? pf_error_message : "process_features direct range begin failed";
        }
        return false;
    }

    std::vector<std::thread> workers;
    std::vector<int> ok(ranges.size(), 0);
    std::vector<std::string> errors(ranges.size());
    workers.reserve(ranges.size());
    for (size_t irange = 0; irange < ranges.size(); ++irange) {
        workers.push_back(std::thread([&, irange]() {
            ok[irange] = process_cbq_direct_pf_range(ctx,
                                                     job,
                                                     plan,
                                                     static_cast<int>(irange),
                                                     ranges[irange].first,
                                                     ranges[irange].count,
                                                     &errors[irange]) ? 1 : 0;
        }));
    }
    for (size_t iworker = 0; iworker < workers.size(); ++iworker) {
        workers[iworker].join();
    }

    for (size_t irange = 0; irange < ok.size(); ++irange) {
        if (!ok[irange]) {
            if (error) {
                *error = errors[irange].empty()
                    ? "process_features direct CBQ worker failed"
                    : errors[irange];
            }
            pf_direct_range_abort(job);
            return false;
        }
    }

    err = pf_direct_range_end(job, stats);
    if (err != PF_OK) {
        if (error) {
            const char* pf_error_message = pf_get_error(ctx);
            *error = pf_error_message ? pf_error_message : "process_features direct range finish failed";
        }
        return false;
    }
    return true;
}

bool process_cbq_records(pf_context* ctx,
                         const Args& args,
                         pf_stats* stats,
                         std::string* error) {
    star::input::InputSourcePlan plan = make_cbq_pair_plan(args);

    star::input::CbqInputModule module;
    if (!module.configure(plan, error) || !module.open(error)) {
        return false;
    }

    pf_record_stream* stream = nullptr;
    pf_error err = pf_process_records_begin(ctx,
                                            args.output_dir.c_str(),
                                            args.sample_name.c_str(),
                                            &stream);
    if (err != PF_OK) {
        if (error) {
            const char* pf_error_message = pf_get_error(ctx);
            *error = pf_error_message ? pf_error_message : "process_features stream begin failed";
        }
        module.close();
        return false;
    }

    std::vector<pf_read_record_view> records;
    std::vector<char> sequence_storage;
    for (;;) {
        star::input::CbqReadBatchView batch;
        const star::input::InputStatus status = module.next_batch(&batch, error);
        if (status == star::input::InputStatus::Error) {
            pf_process_records_abort(stream);
            module.close();
            return false;
        }
        if (status == star::input::InputStatus::End) {
            break;
        }
        records.clear();
        records.resize(batch.record_count);
        const size_t required_sequence_storage =
            static_cast<size_t>(batch.record_count) * 2U * kPfLineLength;
        if (sequence_storage.size() < required_sequence_storage) {
            sequence_storage.resize(required_sequence_storage);
        }
        for (uint32_t i = 0; i < batch.record_count; ++i) {
            const star::input::CbqReadView& view = batch.records[i];
            if (view.segment_count < 2 || view.segments == nullptr) {
                if (error) *error = "CBQ PF adapter expects paired CBQ records";
                pf_process_records_abort(stream);
                module.close();
                return false;
            }
            pf_read_record_view rec = {};
            char* barcode_sequence = sequence_storage.data() +
                (static_cast<size_t>(i) * 2U * kPfLineLength);
            char* feature_sequence = barcode_sequence + kPfLineLength;
            size_t barcode_length = 0;
            size_t feature_length = 0;
            if (!star::input::materialize_cbq_segment_sequence_to_buffer(view.segments[0],
                                                                         barcode_sequence,
                                                                         kPfSequenceCapacity,
                                                                         &barcode_length,
                                                                         error) ||
                !star::input::materialize_cbq_segment_sequence_to_buffer(view.segments[1],
                                                                         feature_sequence,
                                                                         kPfSequenceCapacity,
                                                                         &feature_length,
                                                                         error)) {
                pf_process_records_abort(stream);
                module.close();
                return false;
            }
            rec.barcode_sequence = to_pf_view_buffer(barcode_sequence, barcode_length);
            rec.barcode_quality = to_pf_view(view.segments[0].quality);
            rec.feature_sequence = to_pf_view_buffer(feature_sequence, feature_length);
            rec.feature_quality = to_pf_view(view.segments[1].quality);
            records[i] = rec;
        }
        if (!records.empty()) {
            err = pf_process_record_views(stream, records.data(), batch.record_count);
            if (err != PF_OK) {
                if (error) {
                    const char* pf_error_message = pf_get_error(ctx);
                    *error = pf_error_message ? pf_error_message : "process_features CBQ batch failed";
                }
                pf_process_records_abort(stream);
                module.close();
                return false;
            }
        }
    }

    module.close();
    err = pf_process_records_end(stream, stats);
    if (err != PF_OK) {
        if (error) {
            const char* pf_error_message = pf_get_error(ctx);
            *error = pf_error_message ? pf_error_message : "process_features stream finish failed";
        }
        return false;
    }
    return true;
}

struct PfConfigDeleter {
    void operator()(pf_config* config) const { pf_config_destroy(config); }
};

struct PfContextDeleter {
    void operator()(pf_context* ctx) const { pf_destroy(ctx); }
};

int run_pf(const Args& args) {
    std::unique_ptr<pf_config, PfConfigDeleter> config(pf_config_create());
    if (!config) {
        std::cerr << "Failed to create process_features config\n";
        return 1;
    }
    pf_config_set_barcode_length(config.get(), args.barcode_length);
    pf_config_set_umi_length(config.get(), args.umi_length);
    pf_config_set_threads(config.get(), args.threads);
    pf_config_set_search_threads(config.get(), args.search_threads);
    pf_config_set_consumer_threads(config.get(), args.consumer_threads);
    if (args.read_buffer_lines > 0) {
        pf_config_set_read_buffer_lines(config.get(), args.read_buffer_lines);
    }
    if (args.feature_offset >= 0) {
        pf_config_set_feature_offset(config.get(), args.feature_offset);
    }
    pf_config_set_skip_emptydrops(config.get(), args.skip_emptydrops ? 1 : 0);
    pf_config_set_skip_heatmaps(config.get(), args.skip_heatmaps ? 1 : 0);

    std::unique_ptr<pf_context, PfContextDeleter> ctx(pf_init(config.get()));
    if (!ctx) {
        std::cerr << "Failed to initialize process_features context\n";
        return 1;
    }
    pf_error err = pf_load_feature_ref(ctx.get(), args.feature_ref.c_str());
    if (err != PF_OK) {
        std::cerr << "Failed to load feature reference: " << pf_get_error(ctx.get()) << "\n";
        return 1;
    }
    err = pf_load_whitelist(ctx.get(), args.whitelist.c_str());
    if (err != PF_OK) {
        std::cerr << "Failed to load whitelist: " << pf_get_error(ctx.get()) << "\n";
        return 1;
    }

    pf_stats stats;
    if (args.mode == "fastq") {
        const char* barcode_fastqs[] = {args.barcode_fastq.c_str()};
        const char* feature_fastqs[] = {args.feature_fastq.c_str()};
        err = pf_process_fastqs(ctx.get(),
                                barcode_fastqs,
                                feature_fastqs,
                                1,
                                args.output_dir.c_str(),
                                args.sample_name.c_str(),
                                &stats);
    } else if (args.mode == "cbq-direct-pf") {
        std::string input_error;
        if (!process_cbq_direct_pf_records(ctx.get(), args, &stats, &input_error)) {
            std::cerr << "Failed to process direct CBQ records: " << input_error << "\n";
            return 1;
        }
        err = PF_OK;
    } else {
        std::string input_error;
        if (!process_cbq_records(ctx.get(), args, &stats, &input_error)) {
            std::cerr << "Failed to process CBQ records: " << input_error << "\n";
            return 1;
        }
        err = PF_OK;
    }

    if (err != PF_OK) {
        const char* error = pf_get_error(ctx.get());
        std::cerr << "process_features processing failed";
        if (error) std::cerr << ": " << error;
        std::cerr << "\n";
        return 1;
    }

    std::cout << "processed_reads\t" << stats.total_reads << "\n";
    std::cout << "matched_reads\t" << stats.matched_reads << "\n";
    std::cout << "features\t" << stats.total_features << "\n";
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, &args)) {
        usage(argv[0]);
        return 2;
    }
    if (args.mode == "cbq-direct-decode") {
        return run_cbq_direct_decode(args);
    }
    return run_pf(args);
}

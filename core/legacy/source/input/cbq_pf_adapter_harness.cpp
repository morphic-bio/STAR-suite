#include "input/CbqInputModule.h"

#include "pf_api.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

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
    bool skip_emptydrops = true;
    bool skip_heatmaps = true;
};

void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " --mode fastq|cbq --whitelist WL --featureRef CSV "
        << "--outputDir OUT [--sampleName NAME]\n"
        << "  FASTQ mode: --barcodeFastq R1 --featureFastq R2\n"
        << "  CBQ mode:   --readFilesIn paired.cbq\n";
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

    if (args->mode != "fastq" && args->mode != "cbq") {
        std::cerr << "--mode must be fastq or cbq\n";
        return false;
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
    return true;
}

std::string span_to_string(star::input::CbqByteSpan span) {
    if (span.data == nullptr || span.size == 0) {
        return std::string();
    }
    return std::string(span.data, span.size);
}

struct OwnedPfRecord {
    std::string barcode_sequence;
    std::string barcode_quality;
    std::string feature_sequence;
    std::string feature_quality;
    pf_read_record record;
};

bool load_cbq_records(const std::string& cbq,
                      std::vector<OwnedPfRecord>* owned_records,
                      std::string* error) {
    std::vector<std::vector<std::string>> read_files(1);
    read_files[0].push_back(cbq);
    star::input::InputSourcePlan plan =
        star::input::make_cbq_input_source_plan(read_files, std::vector<std::string>(), 2);
    plan.read_name_separator_chars.clear();
    plan.read_name_separator_chars.push_back(' ');

    star::input::CbqInputModule module;
    if (!module.configure(plan, error) || !module.open(error)) {
        return false;
    }

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
        for (uint32_t i = 0; i < batch.record_count; ++i) {
            const star::input::CbqReadView& view = batch.records[i];
            if (view.segment_count < 2 || view.segments == nullptr) {
                if (error) *error = "CBQ PF adapter expects paired CBQ records";
                module.close();
                return false;
            }
            OwnedPfRecord rec;
            rec.barcode_sequence = span_to_string(view.segments[0].sequence);
            rec.barcode_quality = span_to_string(view.segments[0].quality);
            rec.feature_sequence = span_to_string(view.segments[1].sequence);
            rec.feature_quality = span_to_string(view.segments[1].quality);
            rec.record.barcode_sequence = rec.barcode_sequence.c_str();
            rec.record.barcode_quality = rec.barcode_quality.c_str();
            rec.record.feature_sequence = rec.feature_sequence.c_str();
            rec.record.feature_quality = rec.feature_quality.c_str();
            rec.record.feature_sequence2 = nullptr;
            rec.record.feature_quality2 = nullptr;
            owned_records->push_back(std::move(rec));
        }
    }

    module.close();

    for (OwnedPfRecord& rec : *owned_records) {
        rec.record.barcode_sequence = rec.barcode_sequence.c_str();
        rec.record.barcode_quality = rec.barcode_quality.c_str();
        rec.record.feature_sequence = rec.feature_sequence.c_str();
        rec.record.feature_quality = rec.feature_quality.c_str();
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
    } else {
        std::vector<OwnedPfRecord> owned_records;
        std::string input_error;
        if (!load_cbq_records(args.cbq, &owned_records, &input_error)) {
            std::cerr << "Failed to load CBQ records: " << input_error << "\n";
            return 1;
        }
        std::vector<pf_read_record> records;
        records.reserve(owned_records.size());
        for (const OwnedPfRecord& rec : owned_records) {
            records.push_back(rec.record);
        }
        err = pf_process_records(ctx.get(),
                                 records.data(),
                                 records.size(),
                                 args.output_dir.c_str(),
                                 args.sample_name.c_str(),
                                 &stats);
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
    return run_pf(args);
}

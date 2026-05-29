#include "input/CbqInputModule.h"
#include "input/CbqStarAdapter.h"
#include "SequenceFuns.h"

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using star::input::CbqByteSpan;
using star::input::CbqInputModule;
using star::input::CbqReadBatchView;
using star::input::CbqReadView;
using star::input::CbqSegmentView;
using star::input::CbqStarAdapterOptions;
using star::input::CbqStarReadBuffers;
using star::input::InputSourcePlan;
using star::input::InputStatus;
using star::input::SourceFormat;

namespace {

bool starts_with_option(const std::string& value) {
    return value.size() >= 2 && value[0] == '-' && value[1] == '-';
}

void split_string(const std::string& value, char delimiter, std::vector<std::string>& out) {
    out.clear();
    std::string item;
    std::istringstream in(value);
    while (std::getline(in, item, delimiter)) {
        if (!item.empty()) {
            out.push_back(item);
        }
    }
}

std::vector<std::string> split_tab_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream in(line);
    while (std::getline(in, field, '\t')) {
        fields.push_back(field);
    }
    return fields;
}

std::string read_group_id(const std::string& rg_line) {
    std::string line = rg_line;
    if (line.substr(0, 3) != "ID:") {
        line = "ID:" + line;
    }
    const size_t tab = line.find('\t');
    return line.substr(3, tab == std::string::npos ? std::string::npos : tab - 3);
}

std::string span_string(CbqByteSpan span) {
    if (span.data == nullptr || span.size == 0) {
        return std::string();
    }
    return std::string(span.data, span.size);
}

struct HarnessOptions {
    std::vector<std::string> read_files_in;
    std::string read_files_manifest = "-";
    std::string read_files_prefix;
    std::vector<char> read_name_separator_chars{' '};
    uint32_t mate_count = 2;
    std::string mode = "direct";
    uint64_t benchmark_iterations = 5;
    uint64_t benchmark_max_records = 1000000;
};

void usage(std::ostream& out) {
    out << "Usage: cbq_star_adapter_harness [options]\n"
        << "  --readFilesIn sample.cbq[,lane2.cbq]\n"
        << "  --readFilesManifest manifest.tsv  # CBQ<TAB>-<TAB>ReadGroup\n"
        << "  --mateCount 1|2\n"
        << "  --readFilesPrefix PREFIX\n"
        << "  --readNameSeparator space|none|CHAR[,CHAR...]\n"
        << "  --mode direct|reference|benchmark\n"
        << "  --iterations N          # benchmark mode only, default 5\n"
        << "  --maxRecords N          # benchmark mode only, 0 means all, default 1000000\n";
}

std::vector<char> parse_separators(const std::string& value) {
    std::vector<char> separators;
    std::vector<std::string> tokens;
    split_string(value, ',', tokens);
    for (const std::string& token : tokens) {
        if (token == "space") {
            separators.push_back(' ');
        } else if (token == "none") {
            continue;
        } else if (token.size() == 1) {
            separators.push_back(token[0]);
        } else {
            throw std::runtime_error("unsupported --readNameSeparator value: " + token);
        }
    }
    return separators;
}

HarnessOptions parse_args(int argc, char* argv[]) {
    HarnessOptions opts;

    for (int ii = 1; ii < argc; ++ii) {
        const std::string arg = argv[ii];
        if (arg == "--help" || arg == "-h") {
            usage(std::cout);
            std::exit(0);
        } else if (arg == "--readFilesIn") {
            if (++ii >= argc || starts_with_option(argv[ii])) {
                throw std::runtime_error("--readFilesIn requires one comma-separated CBQ lane list");
            }
            opts.read_files_in.push_back(argv[ii]);
        } else if (arg == "--readFilesManifest") {
            if (++ii >= argc) {
                throw std::runtime_error("--readFilesManifest requires a value");
            }
            opts.read_files_manifest = argv[ii];
        } else if (arg == "--readFilesPrefix") {
            if (++ii >= argc) {
                throw std::runtime_error("--readFilesPrefix requires a value");
            }
            opts.read_files_prefix = argv[ii] == std::string("-") ? "" : argv[ii];
        } else if (arg == "--mateCount") {
            if (++ii >= argc) {
                throw std::runtime_error("--mateCount requires a value");
            }
            opts.mate_count = static_cast<uint32_t>(std::stoul(argv[ii]));
            if (opts.mate_count != 1 && opts.mate_count != 2) {
                throw std::runtime_error("--mateCount must be 1 or 2");
            }
        } else if (arg == "--readNameSeparator") {
            if (++ii >= argc) {
                throw std::runtime_error("--readNameSeparator requires a value");
            }
            opts.read_name_separator_chars = parse_separators(argv[ii]);
        } else if (arg == "--mode") {
            if (++ii >= argc) {
                throw std::runtime_error("--mode requires a value");
            }
            opts.mode = argv[ii];
            if (opts.mode != "direct" && opts.mode != "reference" && opts.mode != "benchmark") {
                throw std::runtime_error("--mode must be direct, reference, or benchmark");
            }
        } else if (arg == "--iterations") {
            if (++ii >= argc) {
                throw std::runtime_error("--iterations requires a value");
            }
            opts.benchmark_iterations = std::stoull(argv[ii]);
            if (opts.benchmark_iterations == 0) {
                throw std::runtime_error("--iterations must be > 0");
            }
        } else if (arg == "--maxRecords") {
            if (++ii >= argc) {
                throw std::runtime_error("--maxRecords requires a value");
            }
            opts.benchmark_max_records = std::stoull(argv[ii]);
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    return opts;
}

InputSourcePlan build_plan_from_manifest(const HarnessOptions& opts) {
    std::ifstream manifest(opts.read_files_manifest.c_str());
    if (!manifest.good()) {
        throw std::runtime_error("could not open manifest: " + opts.read_files_manifest);
    }

    std::vector<std::string> cbq_files;
    std::vector<std::string> read_groups;
    std::string line;
    while (std::getline(manifest, line)) {
        if (line.find_first_not_of(" \t") == std::string::npos) {
            continue;
        }
        const std::vector<std::string> fields = split_tab_line(line);
        if (fields.size() < 3) {
            throw std::runtime_error("CBQ manifest rows require CBQ<TAB>-<TAB>ReadGroup");
        }
        if (fields[1] != "-") {
            throw std::runtime_error("CBQ manifest second column must be '-' for records packed in one CBQ");
        }
        cbq_files.push_back(opts.read_files_prefix + fields[0]);
        read_groups.push_back(read_group_id(fields[2]));
    }
    if (cbq_files.empty()) {
        throw std::runtime_error("manifest did not contain any CBQ files");
    }

    InputSourcePlan plan;
    plan.format = SourceFormat::Binseq;
    plan.module_name = "Cbq";
    plan.mate_files.push_back(cbq_files);
    plan.read_groups = read_groups;
    plan.read_files_n = static_cast<uint32_t>(cbq_files.size());
    plan.mate_count = opts.mate_count;
    plan.read_name_separator_chars = opts.read_name_separator_chars;
    return plan;
}

InputSourcePlan build_plan(const HarnessOptions& opts) {
    if (opts.read_files_manifest != "-") {
        return build_plan_from_manifest(opts);
    }
    if (opts.read_files_in.size() != 1) {
        throw std::runtime_error("--readFilesIn requires one comma-separated CBQ lane list");
    }

    std::vector<std::string> cbq_files;
    split_string(opts.read_files_in.front(), ',', cbq_files);
    if (cbq_files.empty()) {
        throw std::runtime_error("--readFilesIn did not contain any CBQ files");
    }
    for (std::string& path : cbq_files) {
        path = opts.read_files_prefix + path;
    }

    InputSourcePlan plan;
    plan.format = SourceFormat::Binseq;
    plan.module_name = "Cbq";
    plan.mate_files.push_back(cbq_files);
    plan.read_files_n = static_cast<uint32_t>(cbq_files.size());
    plan.mate_count = opts.mate_count;
    plan.read_name_separator_chars = opts.read_name_separator_chars;
    return plan;
}

struct StarReadState {
    explicit StarReadState(uint32_t nends_in) : nends(nends_in) {
        read_name_storage.resize(nends);
        read0_storage.resize(nends);
        read1_storage.resize(nends);
        qual0_storage.resize(nends);
        read_name_mates.resize(nends);
        read0.resize(nends);
        read1.resize(nends);
        qual0.resize(nends);
        read_name_extra.resize(nends);
        read_length.assign(nends, 0);
        read_length_original.assign(nends, 0);
        clip_mates.resize(nends);

        for (uint32_t imate = 0; imate < nends; ++imate) {
            read_name_storage[imate].assign(DEF_readNameLengthMax, '\0');
            read0_storage[imate].assign(DEF_readSeqLengthMax + 1, '\0');
            read1_storage[imate].assign(DEF_readSeqLengthMax + 1, '\0');
            qual0_storage[imate].assign(DEF_readSeqLengthMax + 1, '\0');
            read_name_mates[imate] = read_name_storage[imate].data();
            read0[imate] = read0_storage[imate].data();
            read1[imate] = read1_storage[imate].data();
            qual0[imate] = qual0_storage[imate].data();
            clip_mates[imate].resize(2);
            for (uint32_t iclip = 0; iclip < 2; ++iclip) {
                clip_mates[imate][iclip].type = -1;
                clip_mates[imate][iclip].N = 0;
                clip_mates[imate][iclip].NafterAd = 0;
                clip_mates[imate][iclip].adMMp = 0.0;
                clip_mates[imate][iclip].clippedInfo = 0;
                clip_mates[imate][iclip].clippedAdN = 0;
                clip_mates[imate][iclip].clippedAdMM = 0;
                clip_mates[imate][iclip].clippedN = 0;
            }
        }
    }

    uint32_t nends = 0;
    std::vector<std::vector<char>> read_name_storage;
    std::vector<std::vector<char>> read0_storage;
    std::vector<std::vector<char>> read1_storage;
    std::vector<std::vector<char>> qual0_storage;
    std::vector<char*> read_name_mates;
    std::vector<char*> read0;
    std::vector<char*> read1;
    std::vector<char*> qual0;
    std::vector<std::string> read_name_extra;
    std::vector<uint> read_length;
    std::vector<uint> read_length_original;
    std::vector<std::vector<ClipMate>> clip_mates;
    uint64 i_read_all = 0;
    uint32 read_files_index = 0;
    char read_filter = 'N';
    int read_file_type = 0;
};

void apply_quality_conversion(char* quality, uint length, int add) {
    if (add == 0) {
        return;
    }
    for (uint ii = 0; ii < length; ++ii) {
        int qs = static_cast<int>(quality[ii]) + add;
        if (qs < 33) {
            qs = 33;
        } else if (qs > 126) {
            qs = 126;
        }
        quality[ii] = static_cast<char>(qs);
    }
}

void copy_string_checked(const std::string& value, char* dest, size_t capacity, const std::string& label) {
    if (value.size() + 1 > capacity) {
        std::ostringstream msg;
        msg << label << " length " << value.size() << " exceeds capacity " << capacity;
        throw std::runtime_error(msg.str());
    }
    if (!value.empty()) {
        std::memcpy(dest, value.data(), value.size());
    }
    dest[value.size()] = '\0';
}

void copy_span_checked(CbqByteSpan span, char* dest, size_t capacity, const std::string& label) {
    if (span.size + 1 > capacity) {
        std::ostringstream msg;
        msg << label << " length " << span.size << " exceeds capacity " << capacity;
        throw std::runtime_error(msg.str());
    }
    if (span.size != 0) {
        std::memcpy(dest, span.data, span.size);
    }
    dest[span.size] = '\0';
}

void load_reference_state(const CbqReadView& view, StarReadState* state) {
    const bool has_quality = view.segments[0].has_quality;
    const char prefix = has_quality ? '@' : '>';

    std::ostringstream header;
    header << prefix << span_string(view.read_name)
           << ' ' << view.read_ordinal
           << ' ' << view.read_filter
           << ' ' << view.lane_index;

    std::istringstream header_stream(header.str());
    std::string parsed_name;
    header_stream >> parsed_name;
    header_stream >> state->i_read_all >> state->read_filter >> state->read_files_index;
    std::string ignored_extra;
    header_stream >> std::ws;
    std::getline(header_stream, ignored_extra);
    state->read_file_type = has_quality ? 2 : 1;

    for (uint32_t imate = 0; imate < state->nends; ++imate) {
        const CbqSegmentView& segment = view.segments[imate];
        std::string sequence;
        materialize_cbq_segment_sequence(segment, &sequence);
        copy_string_checked(parsed_name, state->read_name_mates[imate], DEF_readNameLengthMax, "read name");
        copy_string_checked(sequence, state->read0[imate], DEF_readSeqLengthMax + 1, "sequence");

        const uint length = static_cast<uint>(sequence.size());
        if (segment.has_quality) {
            copy_span_checked(segment.quality, state->qual0[imate], DEF_readSeqLengthMax + 1, "quality");
            apply_quality_conversion(state->qual0[imate], length, 0);
        } else {
            std::fill(state->qual0[imate], state->qual0[imate] + length, 'A');
            state->qual0[imate][length] = '\0';
        }

        state->read_length[imate] = length;
        state->read_length_original[imate] = length;
        state->read_name_extra[imate].clear();
        if (segment.has_quality) {
            state->clip_mates[imate][0].clippedInfo = '+';
        }
        convertNucleotidesToNumbers(state->read0[imate], state->read1[imate], state->read_length[imate]);
        state->clip_mates[imate][0].clip(state->read_length[imate], state->read1[imate]);
        state->clip_mates[imate][1].clip(state->read_length[imate], state->read1[imate]);
    }
}

void load_direct_state(const CbqReadView& view, StarReadState* state) {
    CbqStarAdapterOptions options;
    options.read_nends = state->nends;
    options.out_sam_read_id_number = false;
    options.out_qs_conversion_add = 0;
    options.trim_cutadapt_enabled = false;
    options.preserve_read_name_extra = false;

    CbqStarReadBuffers buffers;
    buffers.read_name_mates = state->read_name_mates.data();
    buffers.read0 = state->read0.data();
    buffers.read1 = state->read1.data();
    buffers.qual0 = state->qual0.data();
    buffers.read_name_extra = &state->read_name_extra;
    buffers.read_length = state->read_length.data();
    buffers.read_length_original = state->read_length_original.data();
    buffers.i_read_all = &state->i_read_all;
    buffers.read_files_index = &state->read_files_index;
    buffers.read_filter = &state->read_filter;
    buffers.read_file_type = &state->read_file_type;

    std::string error;
    if (!star::input::load_cbq_read_view_into_star_mates(
            view, options, &buffers, &state->clip_mates, &error)) {
        throw std::runtime_error("direct adapter failed: " + error);
    }
}

void write_u8(std::ostream& out, uint8_t value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void write_u32(std::ostream& out, uint32_t value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void write_u64(std::ostream& out, uint64_t value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void write_string(std::ostream& out, const std::string& value) {
    write_u64(out, static_cast<uint64_t>(value.size()));
    if (!value.empty()) {
        out.write(value.data(), static_cast<std::streamsize>(value.size()));
    }
}

void write_c_string(std::ostream& out, const char* value) {
    write_string(out, value == nullptr ? std::string() : std::string(value));
}

void emit_state(const StarReadState& state, std::ostream& out) {
    write_u32(out, state.nends);
    write_u64(out, state.i_read_all);
    write_u32(out, state.read_files_index);
    write_u8(out, static_cast<uint8_t>(state.read_filter));
    write_u32(out, static_cast<uint32_t>(state.read_file_type));

    for (uint32_t imate = 0; imate < state.nends; ++imate) {
        write_c_string(out, state.read_name_mates[imate]);
        write_string(out, state.read_name_extra[imate]);
        write_u32(out, static_cast<uint32_t>(state.read_length[imate]));
        write_u32(out, static_cast<uint32_t>(state.read_length_original[imate]));

        const uint32_t original_length = static_cast<uint32_t>(state.read_length_original[imate]);
        out.write(state.read0[imate], static_cast<std::streamsize>(original_length + 1));
        out.write(state.qual0[imate], static_cast<std::streamsize>(original_length + 1));
        out.write(state.read1[imate], static_cast<std::streamsize>(state.read_length[imate]));

        write_u8(out, static_cast<uint8_t>(state.clip_mates[imate][0].clippedInfo));
        write_u8(out, static_cast<uint8_t>(state.clip_mates[imate][1].clippedInfo));
        write_u32(out, state.clip_mates[imate][0].clippedN);
        write_u32(out, state.clip_mates[imate][1].clippedN);
    }
}

struct PreloadedBatch {
    CbqReadBatchView batch;
};

struct PreloadedData {
    std::vector<PreloadedBatch> batches;
    uint64_t record_count = 0;
    uint64_t mate_count = 0;
    uint64_t base_count = 0;
};

double elapsed_seconds(std::chrono::steady_clock::time_point begin,
                       std::chrono::steady_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count();
}

void hash_mix(uint64_t* hash, uint64_t value) {
    *hash ^= value;
    *hash *= 1099511628211ULL;
}

void hash_bytes(uint64_t* hash, const char* data, size_t size) {
    for (size_t ii = 0; ii < size; ++ii) {
        hash_mix(hash, static_cast<unsigned char>(data[ii]));
    }
}

void hash_sampled_bytes(uint64_t* hash, const char* data, size_t size) {
    if (data == nullptr || size == 0) {
        hash_mix(hash, 0);
        return;
    }
    hash_mix(hash, size);
    if (size <= 16) {
        hash_bytes(hash, data, size);
        return;
    }
    hash_bytes(hash, data, 8);
    hash_bytes(hash, data + size - 8, 8);
}

uint64_t hash_state_summary(uint64_t hash, const StarReadState& state) {
    hash_mix(&hash, state.i_read_all);
    hash_mix(&hash, state.read_files_index);
    hash_mix(&hash, static_cast<unsigned char>(state.read_filter));
    hash_mix(&hash, static_cast<uint64_t>(state.read_file_type));
    for (uint32_t imate = 0; imate < state.nends; ++imate) {
        const uint length = state.read_length[imate];
        hash_mix(&hash, length);
        hash_mix(&hash, state.read_length_original[imate]);
        hash_sampled_bytes(&hash, state.read_name_mates[imate], std::strlen(state.read_name_mates[imate]));
        hash_sampled_bytes(&hash, state.read0[imate], length);
        hash_sampled_bytes(&hash, state.read1[imate], length);
        hash_sampled_bytes(&hash, state.qual0[imate], length);
        hash_mix(&hash, static_cast<unsigned char>(state.clip_mates[imate][0].clippedInfo));
        hash_mix(&hash, static_cast<unsigned char>(state.clip_mates[imate][1].clippedInfo));
    }
    return hash;
}

PreloadedData preload_decoded_batches(CbqInputModule* module,
                                      uint64_t max_records,
                                      uint32_t mate_count) {
    PreloadedData data;
    CbqReadBatchView batch;
    std::string error;
    while (true) {
        error.clear();
        const InputStatus status = module->next_batch(&batch, &error);
        if (status == InputStatus::End) {
            break;
        }
        if (status == InputStatus::Error) {
            throw std::runtime_error("read failed: " + error);
        }

        CbqReadBatchView stored = batch;
        if (max_records != 0 && data.record_count + stored.record_count > max_records) {
            stored.record_count = static_cast<uint32_t>(max_records - data.record_count);
        }
        if (stored.record_count == 0) {
            break;
        }

        for (uint32_t irecord = 0; irecord < stored.record_count; ++irecord) {
            const CbqReadView& record = stored.records[irecord];
            const uint32_t nends = std::min<uint32_t>(record.segment_count, mate_count);
            data.mate_count += nends;
            for (uint32_t imate = 0; imate < nends; ++imate) {
                data.base_count += cbq_segment_sequence_length(record.segments[imate]);
            }
        }

        data.record_count += stored.record_count;
        data.batches.push_back(PreloadedBatch{stored});
        if (max_records != 0 && data.record_count >= max_records) {
            break;
        }
    }
    return data;
}

int run_benchmark(const HarnessOptions& opts, const InputSourcePlan& plan) {
    CbqInputModule module;
    std::string error;
    if (!module.configure(plan, &error)) {
        std::cerr << "configure failed: " << error << "\n";
        return 2;
    }
    if (!module.open(&error)) {
        std::cerr << "open failed: " << error << "\n";
        return 2;
    }

    const auto preload_begin = std::chrono::steady_clock::now();
    PreloadedData data = preload_decoded_batches(&module, opts.benchmark_max_records, opts.mate_count);
    const auto preload_end = std::chrono::steady_clock::now();
    module.close();

    if (data.record_count == 0) {
        std::cerr << "benchmark loaded zero records\n";
        return 4;
    }

    StarReadState state(opts.mate_count);
    uint64_t checksum = 1469598103934665603ULL;
    const auto adapter_begin = std::chrono::steady_clock::now();
    for (uint64_t iter = 0; iter < opts.benchmark_iterations; ++iter) {
        for (const PreloadedBatch& loaded : data.batches) {
            for (uint32_t irecord = 0; irecord < loaded.batch.record_count; ++irecord) {
                load_direct_state(loaded.batch.records[irecord], &state);
                checksum = hash_state_summary(checksum, state);
            }
        }
    }
    const auto adapter_end = std::chrono::steady_clock::now();

    const double preload_seconds = elapsed_seconds(preload_begin, preload_end);
    const double adapter_seconds = elapsed_seconds(adapter_begin, adapter_end);
    const uint64_t adapter_records = data.record_count * opts.benchmark_iterations;
    const uint64_t adapter_bases = data.base_count * opts.benchmark_iterations;

    std::cout << std::fixed << std::setprecision(6)
              << "mode=benchmark\n"
              << "records_loaded=" << data.record_count << "\n"
              << "mates_loaded=" << data.mate_count << "\n"
              << "bases_loaded=" << data.base_count << "\n"
              << "batches_loaded=" << data.batches.size() << "\n"
              << "iterations=" << opts.benchmark_iterations << "\n"
              << "preload_seconds=" << preload_seconds << "\n"
              << "adapter_seconds=" << adapter_seconds << "\n"
              << "adapter_records=" << adapter_records << "\n"
              << "adapter_bases=" << adapter_bases << "\n"
              << "adapter_records_per_second=" << (static_cast<double>(adapter_records) / adapter_seconds) << "\n"
              << "adapter_bases_per_second=" << (static_cast<double>(adapter_bases) / adapter_seconds) << "\n"
              << "checksum=" << checksum << "\n";
    return 0;
}

int run(const HarnessOptions& opts) {
    InputSourcePlan plan = build_plan(opts);
    if (opts.mode == "benchmark") {
        return run_benchmark(opts, plan);
    }

    CbqInputModule module;
    std::string error;
    if (!module.configure(plan, &error)) {
        std::cerr << "configure failed: " << error << "\n";
        return 2;
    }
    if (!module.open(&error)) {
        std::cerr << "open failed: " << error << "\n";
        return 2;
    }

    CbqReadBatchView batch;
    while (true) {
        error.clear();
        const InputStatus status = module.next_batch(&batch, &error);
        if (status == InputStatus::End) {
            break;
        }
        if (status == InputStatus::Error) {
            std::cerr << "read failed: " << error << "\n";
            return 3;
        }
        for (uint32_t irecord = 0; irecord < batch.record_count; ++irecord) {
            StarReadState state(opts.mate_count);
            const CbqReadView& record = batch.records[irecord];
            if (opts.mode == "direct") {
                load_direct_state(record, &state);
            } else {
                load_reference_state(record, &state);
            }
            emit_state(state, std::cout);
        }
    }
    module.close();
    return 0;
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        return run(parse_args(argc, argv));
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << "\n\n";
        usage(std::cerr);
        return 1;
    }
}

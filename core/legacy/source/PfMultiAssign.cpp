#include "PfMultiAssign.h"
#include "GlobalVariables.h"
#include "input/CbqInputModule.h"
#include "pf_api.h"
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <limits>
#include <chrono>
#include <thread>
#include <sys/stat.h>
#include <stdexcept>
using std::cerr;
using std::endl;

namespace {
ThreadControl::PermitHookContext kFeaturePermitHookContext{ThreadControl::PermitDomain::FEATURE};
constexpr size_t kPfLineLength = 1024;
constexpr size_t kPfSequenceCapacity = kPfLineLength - 1;
}

extern "C" uint64_t pfStarDynamicPermitAcquire(void *hookCtx) {
    const ThreadControl::PermitHookContext *permitCtx =
        static_cast<const ThreadControl::PermitHookContext *>(hookCtx);
    const ThreadControl::PermitDomain domain =
        (permitCtx == nullptr) ? ThreadControl::PermitDomain::MAP : permitCtx->domain;
    return g_threadChunks.mapPermitAcquireForDomain(domain);
}

extern "C" void pfStarDynamicPermitRelease(
    void *hookCtx,
    uint64_t waitNs,
    uint64_t workUnits,
    uint64_t workBytes,
    uint64_t workNs
) {
    const ThreadControl::PermitHookContext *permitCtx =
        static_cast<const ThreadControl::PermitHookContext *>(hookCtx);
    const ThreadControl::PermitDomain domain =
        (permitCtx == nullptr) ? ThreadControl::PermitDomain::MAP : permitCtx->domain;
    g_threadChunks.mapPermitReleaseForDomain(domain, waitNs, workUnits, workBytes, workNs);
}

namespace PfMultiAssign {

namespace {

static bool fileExists(const string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

static bool dirExists(const string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

static string lowerCopyLocal(string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

static bool hasSuffix(const string& value, const string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static string stripOuterQuotes(const string& value) {
    size_t first = value.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        return "";
    }
    size_t last = value.find_last_not_of(" \t\r\n");
    string stripped = value.substr(first, last - first + 1);
    if (stripped.size() >= 2 &&
        ((stripped.front() == '"' && stripped.back() == '"') ||
         (stripped.front() == '\'' && stripped.back() == '\''))) {
        stripped = stripped.substr(1, stripped.size() - 2);
    }
    return stripped;
}

static bool isCbqPath(const string& path) {
    return hasSuffix(lowerCopyLocal(path), ".cbq");
}

static vector<string> splitCommaList(const string& value) {
    vector<string> out;
    std::stringstream ss(value);
    string item;
    while (std::getline(ss, item, ',')) {
        size_t first = item.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t last = item.find_last_not_of(" \t\r\n");
        const string stripped = stripOuterQuotes(item.substr(first, last - first + 1));
        if (!stripped.empty()) {
            out.push_back(stripped);
        }
    }
    return out;
}

static vector<string> listCbqFilesInDirectory(const string& dirPath) {
    vector<string> files;
    DIR* dir = opendir(dirPath.c_str());
    if (dir == nullptr) {
        return files;
    }
    struct dirent* entry = nullptr;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_name[0] == '.') {
            continue;
        }
        const string full = dirPath + "/" + entry->d_name;
        if (!isCbqPath(full)) {
            continue;
        }
        struct stat st;
        if (stat(full.c_str(), &st) == 0 && S_ISREG(st.st_mode)) {
            files.push_back(full);
        }
    }
    closedir(dir);
    std::sort(files.begin(), files.end());
    return files;
}

static bool resolveCbqSources(const string& source, vector<string>* cbqPaths) {
    if (cbqPaths == nullptr) {
        return false;
    }
    cbqPaths->clear();
    const string normalizedSource = stripOuterQuotes(source);
    if (dirExists(normalizedSource)) {
        *cbqPaths = listCbqFilesInDirectory(normalizedSource);
        return !cbqPaths->empty();
    }

    vector<string> candidates = splitCommaList(normalizedSource);
    if (candidates.empty()) {
        return false;
    }
    for (const string& path : candidates) {
        if (!isCbqPath(path) || !fileExists(path)) {
            cbqPaths->clear();
            return false;
        }
        cbqPaths->push_back(path);
    }
    return !cbqPaths->empty();
}

static bool looksLikeMultiColumnWhitelist(const string& whitelistPath) {
    std::ifstream in(whitelistPath.c_str());
    if (!in.is_open()) {
        return false;
    }
    string line;
    while (std::getline(in, line)) {
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        return line.find('\t') != string::npos || line.find(',') != string::npos;
    }
    return false;
}

static string upperCopy(string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return value;
}

static bool pathContainsToken(const string& path, const string& token) {
    const string up = upperCopy(path);
    return up.find(token) != string::npos;
}

static char complementBase(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    }
}

static string translateNxtMiddleTwoBases(string barcode) {
    barcode = upperCopy(barcode);
    if (barcode.size() >= 9) {
        barcode[7] = complementBase(barcode[7]);
        barcode[8] = complementBase(barcode[8]);
    }
    return barcode;
}

static bool isValidBarcodeSeq(const string& seq) {
    if (seq.empty()) {
        return false;
    }
    for (unsigned char c : seq) {
        unsigned char u = static_cast<unsigned char>(std::toupper(c));
        if (!(u == 'A' || u == 'C' || u == 'G' || u == 'T' || u == 'N')) {
            return false;
        }
    }
    return true;
}

static uint64 countValidWhitelistRows(const string& whitelistPath) {
    std::ifstream in(whitelistPath.c_str());
    if (!in.is_open()) {
        return 0;
    }
    string line;
    uint64 count = 0;
    while (std::getline(in, line)) {
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t end = line.find_first_of("\t, \r\n", first);
        string token = (end == string::npos) ? line.substr(first) : line.substr(first, end - first);
        if (isValidBarcodeSeq(token)) {
            ++count;
        }
    }
    return count;
}

static string inferOneColumnNamespace(const string& whitelistPath, bool& confident) {
    confident = true;
    if (pathContainsToken(whitelistPath, "NXT")) {
        return "NXT";
    }
    if (pathContainsToken(whitelistPath, "TRU")) {
        return "TRU";
    }
    confident = false;
    return "UNKNOWN";
}

static string inferTwoColumnNamespace(const string& whitelistPath, bool& confident) {
    confident = false;
    if (pathContainsToken(whitelistPath, "NXT")) {
        confident = true;
        return "NXT";
    }
    if (pathContainsToken(whitelistPath, "TRU")) {
        confident = true;
        return "TRU";
    }

    std::ifstream in(whitelistPath.c_str());
    if (!in.is_open()) {
        return "UNKNOWN";
    }

    string line;
    uint64 sampled = 0;
    uint64 twoCol = 0;
    uint64 matchedRule = 0;
    const uint64 kMaxSample = 10000;
    while (sampled < kMaxSample && std::getline(in, line)) {
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        sampled++;

        size_t delim1 = line.find_first_of("\t, \r\n", first);
        if (delim1 == string::npos) {
            continue;
        }
        string col1 = line.substr(first, delim1 - first);
        size_t second = line.find_first_not_of("\t, \r\n", delim1);
        if (second == string::npos) {
            continue;
        }
        size_t delim2 = line.find_first_of("\t, \r\n", second);
        string col2 = (delim2 == string::npos) ? line.substr(second)
                                               : line.substr(second, delim2 - second);
        if (!isValidBarcodeSeq(col1) || !isValidBarcodeSeq(col2)) {
            continue;
        }
        twoCol++;
        if (translateNxtMiddleTwoBases(col1) == upperCopy(col2)) {
            matchedRule++;
        }
    }

    if (twoCol > 0) {
        const double frac = static_cast<double>(matchedRule) / static_cast<double>(twoCol);
        if (frac >= 0.80) {
            // Complement rule is symmetric (translate(A)==B iff translate(B)==A),
            // so content alone cannot distinguish COL1=NXT from COL1=TRU.
            // Return UNKNOWN to force explicit --crChemistry override.
            return "UNKNOWN";
        }
    }
    return "UNKNOWN";
}

static WhitelistNormalizationResult normalizeWhitelistInternal(const string& whitelistPath,
                                                               const string& assignOut) {
    WhitelistNormalizationResult result;
    result.sourcePath = whitelistPath;
    result.normalizedPath = whitelistPath;

    if (!looksLikeMultiColumnWhitelist(whitelistPath)) {
        result.assignmentNamespace = inferOneColumnNamespace(whitelistPath, result.namespaceConfidence);
        result.normalizedRowCount = countValidWhitelistRows(whitelistPath);
        return result;
    }

    result.hasTwoColumnSource = true;
    result.assignmentNamespace = inferTwoColumnNamespace(whitelistPath, result.namespaceConfidence);

    std::ifstream in(whitelistPath.c_str());
    if (!in.is_open()) {
        return result;
    }

    const string normalizedPath = assignOut + "/whitelist.normalized.txt";
    std::ofstream out(normalizedPath.c_str());
    if (!out.is_open()) {
        return result;
    }

    string line;
    uint64 emitted = 0;
    while (std::getline(in, line)) {
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t end = line.find_first_of("\t, \r\n", first);
        string token = (end == string::npos) ? line.substr(first) : line.substr(first, end - first);
        if (!isValidBarcodeSeq(token)) {
            continue;
        }
        out << token << "\n";
        ++emitted;
    }

    if (emitted == 0) {
        return result;
    }
    result.normalizedPath = normalizedPath;
    result.normalizedRowCount = emitted;
    return result;
}

static string pfErrorCodeString(pf_error err) {
    switch (err) {
        case PF_OK: return "PF_OK";
        case PF_ERR_INVALID_ARG: return "PF_ERR_INVALID_ARG";
        case PF_ERR_FILE_NOT_FOUND: return "PF_ERR_FILE_NOT_FOUND";
        case PF_ERR_PARSE_ERROR: return "PF_ERR_PARSE_ERROR";
        case PF_ERR_OUT_OF_MEMORY: return "PF_ERR_OUT_OF_MEMORY";
        case PF_ERR_IO_ERROR: return "PF_ERR_IO_ERROR";
        case PF_ERR_NOT_INITIALIZED: return "PF_ERR_NOT_INITIALIZED";
        case PF_ERR_ALREADY_INITIALIZED: return "PF_ERR_ALREADY_INITIALIZED";
        case PF_ERR_OFFSET_CONFLICT: return "PF_ERR_OFFSET_CONFLICT";
        case PF_ERR_MULTI_OFFSET_DETECTED: return "PF_ERR_MULTI_OFFSET_DETECTED";
        case PF_ERR_ALLOC: return "PF_ERR_ALLOC";
        case PF_ERR_NAMESPACE: return "PF_ERR_NAMESPACE";
        default: return "PF_ERR_UNKNOWN";
    }
}

static void applyAssignOptions(pf_config* cfg, const AssignOptions& options) {
    // Match assignBarcodes CLI default for min_counts.
    pf_config_set_min_counts(cfg, 0);

    if (options.maxHammingDistance >= 0) {
        pf_config_set_max_hamming_distance(cfg, options.maxHammingDistance);
    }
    if (options.featureConstantOffset >= 0) {
        pf_config_set_feature_offset(cfg, options.featureConstantOffset);
    }
    if (options.limitSearch != -2) {
        pf_config_set_limit_search(cfg, options.limitSearch);
    }
    if (options.minCounts >= 0) {
        pf_config_set_min_counts(cfg, options.minCounts);
    }
    if (options.maxBarcodeMismatches >= 0) {
        pf_config_set_max_barcode_mismatches(cfg, options.maxBarcodeMismatches);
    }
    if (options.featureN >= 0) {
        pf_config_set_max_feature_n(cfg, options.featureN);
    }
    if (options.barcodeN >= 0) {
        pf_config_set_max_barcode_n(cfg, options.barcodeN);
    }
    if (options.consumerThreadsPerSet > 0) {
        pf_config_set_consumer_threads(cfg, options.consumerThreadsPerSet);
    }
    if (options.searchThreads > 0) {
        pf_config_set_search_threads(cfg, options.searchThreads);
    }
    if (options.readBufferLines > 0) {
        pf_config_set_read_buffer_lines(cfg, options.readBufferLines);
    }
    if (options.minPosterior >= 0.0) {
        pf_config_set_min_posterior(cfg, options.minPosterior);
    }
    if (options.maxReads > 0) {
        pf_config_set_max_reads(cfg, options.maxReads);
    }
    if (options.legacyCbRescue) {
        pf_config_set_legacy_cb_rescue(cfg, 1);
    }
    if (options.translateNxt) {
        pf_config_set_translate_nxt(cfg, 1);
    }
    if (options.autodetectChemistry) {
        pf_config_set_autodetect_chemistry(cfg, 1);
        pf_config_set_autodetect_chemistry_reads(cfg, options.autodetectChemistryReads);
        pf_config_set_autodetect_chemistry_min_hits(cfg, options.autodetectChemistryMinHits);
    }
    if (options.probeOnly) {
        pf_config_set_probe_only(cfg, 1);
    }
    if (options.skipQcOutputs) {
        pf_config_set_skip_qc_outputs(cfg, 1);
    }
    if (options.enableStarDynamicPermitHooks) {
        pf_config_set_permit_hooks(
            cfg,
            pfStarDynamicPermitAcquire,
            pfStarDynamicPermitRelease,
            &kFeaturePermitHookContext
        );
    }
    if (options.allowUnionWhitelist) {
        pf_config_set_allow_union_whitelist(cfg, 1);
    }
    if (options.sourceNamespace != "UNKNOWN") {
        pf_config_set_source_namespace(cfg, pf_namespace_from_string(options.sourceNamespace.c_str()));
    }
    if (options.targetNamespace != "UNKNOWN") {
        pf_config_set_target_namespace(cfg, pf_namespace_from_string(options.targetNamespace.c_str()));
    }
    if (options.useFeatureAnchorSearch) {
        pf_config_set_use_feature_anchor_search(cfg, 1);
    }
    if (options.requireFeatureAnchorMatch) {
        pf_config_set_require_feature_anchor_match(cfg, 1);
    }
    if (options.featureModeBootstrapReads > 0) {
        pf_config_set_feature_mode_bootstrap_reads(cfg, options.featureModeBootstrapReads);
    }
    if (options.useHotHash) {
        pf_config_set_use_hot_hash(cfg, 1);
    }
    if (options.skipHeatmaps) {
        pf_config_set_skip_heatmaps(cfg, 1);
    }
    if (options.adtMexOutput) {
        pf_config_set_adt_mex_output(cfg, 1);
        pf_config_set_skip_emptydrops(cfg, 1);
    }
    if (options.useSplitReadLayout) {
        pf_config_set_split_read_layout(cfg, &options.splitReadLayout);
    } else {
        pf_config_clear_split_read_layout(cfg);
    }
}

static string pfErrorMessage(pf_context* ctx, pf_error err, const string& stage) {
    std::ostringstream oss;
    oss << "pf_api failed at " << stage << " (" << pfErrorCodeString(err) << ")";
    const char* errMsg = pf_get_error(ctx);
    if (errMsg != nullptr && errMsg[0] != '\0') {
        oss << ": " << errMsg;
    }
    return oss.str();
}

static pf_sequence_view pfViewFromCbqSpan(star::input::CbqByteSpan span) {
    pf_sequence_view view;
    view.data = span.data;
    view.length = span.size;
    return view;
}

static pf_sequence_view pfViewFromBuffer(const char* data, size_t length) {
    pf_sequence_view view;
    view.data = data;
    view.length = length;
    return view;
}

static string normalizeCbqMode(string mode) {
    mode = lowerCopyLocal(stripOuterQuotes(mode));
    if (mode.empty() || mode == "-") {
        return "auto";
    }
    if (mode == "direct") {
        return "range";
    }
    return mode;
}

static star::input::InputSourcePlan makeCbqPairPlan(const vector<string>& cbqPaths) {
    std::vector<std::vector<std::string>> readFiles(1);
    readFiles[0] = cbqPaths;
    star::input::InputSourcePlan plan =
        star::input::make_cbq_input_source_plan(readFiles, std::vector<std::string>(), 2);
    plan.read_name_separator_chars.clear();
    plan.read_name_separator_chars.push_back(' ');
    return plan;
}

struct CbqLaneRange {
    uint32_t lane = 0;
    uint64_t first = 0;
    uint64_t count = 0;
};

struct CbqLaneCount {
    uint32_t lane = 0;
    uint64_t count = 0;
};

static bool directCbqSettingsSupported(const AssignOptions& options,
                                       string* fallbackReason) {
    if (options.legacyCbRescue) {
        if (fallbackReason != nullptr) {
            *fallbackReason = "legacy_cb_rescue";
        }
        return false;
    }
    if (options.featureModeBootstrapReads > 0) {
        if (fallbackReason != nullptr) {
            *fallbackReason = "feature_mode_bootstrap";
        }
        return false;
    }
    return true;
}

static bool inspectCbqLaneCounts(const star::input::InputSourcePlan& plan,
                                 long long maxReads,
                                 vector<CbqLaneCount>* lanes,
                                 uint64_t* totalRecords,
                                 string* errorMessage) {
    lanes->clear();
    *totalRecords = 0;
    uint64_t remaining = maxReads > 0
        ? static_cast<uint64_t>(maxReads)
        : std::numeric_limits<uint64_t>::max();

    for (uint32_t lane = 0; lane < plan.read_files_n; ++lane) {
        std::string inputError;
        star::input::CbqInputModule module;
        if (!module.configure(plan, &inputError) ||
            !module.open_range(lane, 0, std::numeric_limits<uint64_t>::max(), &inputError)) {
            if (errorMessage != nullptr) {
                *errorMessage = inputError;
            }
            return false;
        }
        const uint64_t laneRecords = module.current_lane_record_count();
        module.close();

        const uint64_t allowed = remaining == std::numeric_limits<uint64_t>::max()
            ? laneRecords
            : std::min(laneRecords, remaining);
        CbqLaneCount laneCount;
        laneCount.lane = lane;
        laneCount.count = allowed;
        lanes->push_back(laneCount);
        *totalRecords += allowed;
        if (remaining != std::numeric_limits<uint64_t>::max()) {
            remaining -= allowed;
            if (remaining == 0) {
                for (uint32_t rest = lane + 1; rest < plan.read_files_n; ++rest) {
                    CbqLaneCount zero;
                    zero.lane = rest;
                    zero.count = 0;
                    lanes->push_back(zero);
                }
                break;
            }
        }
    }
    return true;
}

static vector<vector<CbqLaneRange>> makeDirectWorkerRanges(const vector<CbqLaneCount>& lanes,
                                                           uint64_t totalRecords,
                                                           int nworkers) {
    vector<vector<CbqLaneRange>> workerRanges(static_cast<size_t>(nworkers));
    if (nworkers <= 0 || totalRecords == 0) {
        return workerRanges;
    }

    const uint64_t chunkSize =
        (totalRecords + static_cast<uint64_t>(nworkers) - 1U) /
        static_cast<uint64_t>(nworkers);
    for (int worker = 0; worker < nworkers; ++worker) {
        const uint64_t globalFirst = static_cast<uint64_t>(worker) * chunkSize;
        if (globalFirst >= totalRecords) {
            break;
        }
        const uint64_t globalEnd = std::min(totalRecords, globalFirst + chunkSize);
        uint64_t laneGlobalFirst = 0;
        for (const CbqLaneCount& lane : lanes) {
            const uint64_t laneGlobalEnd = laneGlobalFirst + lane.count;
            if (lane.count > 0 && globalFirst < laneGlobalEnd && globalEnd > laneGlobalFirst) {
                const uint64_t overlapFirst = std::max(globalFirst, laneGlobalFirst);
                const uint64_t overlapEnd = std::min(globalEnd, laneGlobalEnd);
                CbqLaneRange range;
                range.lane = lane.lane;
                range.first = overlapFirst - laneGlobalFirst;
                range.count = overlapEnd - overlapFirst;
                workerRanges[static_cast<size_t>(worker)].push_back(range);
            }
            laneGlobalFirst = laneGlobalEnd;
            if (laneGlobalFirst >= globalEnd) {
                break;
            }
        }
    }
    return workerRanges;
}

static bool processCbqDirectRangeWorker(pf_context* ctx,
                                        pf_direct_range_job* job,
                                        const star::input::InputSourcePlan& plan,
                                        int workerId,
                                        const vector<CbqLaneRange>& ranges,
                                        string* errorMessage) {
    std::vector<pf_read_record_view> records;
    std::vector<char> sequenceStorage;
    for (const CbqLaneRange& range : ranges) {
        if (range.count == 0) {
            continue;
        }
        std::string inputError;
        star::input::CbqInputModule module;
        if (!module.configure(plan, &inputError) ||
            !module.open_range(range.lane, range.first, range.count, &inputError)) {
            if (errorMessage != nullptr) {
                *errorMessage = "CBQ direct range open error: " + inputError;
            }
            return false;
        }

        for (;;) {
            star::input::CbqReadBatchView batch;
            const star::input::InputStatus status = module.next_batch(&batch, &inputError);
            if (status == star::input::InputStatus::Error) {
                module.close();
                if (errorMessage != nullptr) {
                    *errorMessage = "CBQ direct range read error: " + inputError;
                }
                return false;
            }
            if (status == star::input::InputStatus::End) {
                break;
            }
            if (batch.record_count == 0) {
                continue;
            }

            records.clear();
            records.resize(batch.record_count);
            const size_t requiredSequenceStorage =
                static_cast<size_t>(batch.record_count) * 2U * kPfLineLength;
            if (sequenceStorage.size() < requiredSequenceStorage) {
                sequenceStorage.resize(requiredSequenceStorage);
            }

            for (uint32_t i = 0; i < batch.record_count; ++i) {
                const star::input::CbqReadView& view = batch.records[i];
                if (view.segment_count < 2 || view.segments == nullptr) {
                    module.close();
                    if (errorMessage != nullptr) {
                        *errorMessage = "CBQ direct range expects paired CBQ records";
                    }
                    return false;
                }

                pf_read_record_view rec = {};
                char* barcodeSequence = sequenceStorage.data() +
                    (static_cast<size_t>(i) * 2U * kPfLineLength);
                char* featureSequence = barcodeSequence + kPfLineLength;
                size_t barcodeLength = 0;
                size_t featureLength = 0;
                if (!star::input::materialize_cbq_segment_sequence_to_buffer(
                        view.segments[0], barcodeSequence, kPfSequenceCapacity,
                        &barcodeLength, &inputError) ||
                    !star::input::materialize_cbq_segment_sequence_to_buffer(
                        view.segments[1], featureSequence, kPfSequenceCapacity,
                        &featureLength, &inputError)) {
                    module.close();
                    if (errorMessage != nullptr) {
                        *errorMessage = "CBQ direct range sequence decode error: " + inputError;
                    }
                    return false;
                }

                rec.barcode_sequence = pfViewFromBuffer(barcodeSequence, barcodeLength);
                rec.barcode_quality = pfViewFromCbqSpan(view.segments[0].quality);
                rec.feature_sequence = pfViewFromBuffer(featureSequence, featureLength);
                rec.feature_quality = pfViewFromCbqSpan(view.segments[1].quality);
                records[i] = rec;
            }

            pf_error err = pf_direct_range_process_record_views(
                job, workerId, records.data(), batch.record_count);
            if (err != PF_OK) {
                module.close();
                if (errorMessage != nullptr) {
                    const char* pfMsg = pf_get_error(ctx);
                    *errorMessage = "CBQ direct range PF worker error";
                    if (pfMsg != nullptr && pfMsg[0] != '\0') {
                        *errorMessage += ": ";
                        *errorMessage += pfMsg;
                    }
                }
                return false;
            }
        }
        module.close();
    }
    return true;
}

static pf_error processCbqSourcesDirect(pf_context* ctx,
                                        const vector<string>& cbqPaths,
                                        const string& outputDir,
                                        const string& sampleName,
                                        long long maxReads,
                                        int requestedWorkers,
                                        pf_stats* stats,
                                        string* errorMessage) {
    const star::input::InputSourcePlan plan = makeCbqPairPlan(cbqPaths);
    vector<CbqLaneCount> laneCounts;
    uint64_t totalRecords = 0;
    if (!inspectCbqLaneCounts(plan, maxReads, &laneCounts, &totalRecords, errorMessage)) {
        return PF_ERR_IO_ERROR;
    }

    const int nworkers = std::max(1, requestedWorkers);
    const vector<vector<CbqLaneRange>> workerRanges =
        makeDirectWorkerRanges(laneCounts, totalRecords, nworkers);

    pf_direct_range_job* job = nullptr;
    pf_error err = pf_direct_range_begin(ctx, outputDir.c_str(),
                                         sampleName.c_str(), nworkers, 2, &job);
    if (err != PF_OK) {
        return err;
    }

    vector<std::thread> workers;
    vector<int> ok(static_cast<size_t>(nworkers), 1);
    vector<string> workerErrors(static_cast<size_t>(nworkers));
    workers.reserve(static_cast<size_t>(nworkers));
    for (int worker = 0; worker < nworkers; ++worker) {
        workers.push_back(std::thread([&, worker]() {
            ok[static_cast<size_t>(worker)] =
                processCbqDirectRangeWorker(ctx, job, plan, worker,
                                            workerRanges[static_cast<size_t>(worker)],
                                            &workerErrors[static_cast<size_t>(worker)])
                ? 1 : 0;
        }));
    }
    for (std::thread& worker : workers) {
        worker.join();
    }

    for (int worker = 0; worker < nworkers; ++worker) {
        if (!ok[static_cast<size_t>(worker)]) {
            if (errorMessage != nullptr) {
                *errorMessage = workerErrors[static_cast<size_t>(worker)].empty()
                    ? "CBQ direct range worker failed"
                    : workerErrors[static_cast<size_t>(worker)];
            }
            pf_direct_range_abort(job);
            return PF_ERR_IO_ERROR;
        }
    }

    return pf_direct_range_end(job, stats);
}

static pf_error processCbqSources(pf_context* ctx,
                                  const vector<string>& cbqPaths,
                                  const string& outputDir,
                                  const string& sampleName,
                                  long long maxReads,
                                  pf_stats* stats,
                                  string* errorMessage) {
    star::input::InputSourcePlan plan = makeCbqPairPlan(cbqPaths);

    std::string inputError;
    star::input::CbqInputModule module;
    if (!module.configure(plan, &inputError) || !module.open(&inputError)) {
        if (errorMessage != nullptr) {
            *errorMessage = "CBQ feature assignment input error: " + inputError;
        }
        return PF_ERR_IO_ERROR;
    }

    pf_record_stream* stream = nullptr;
    pf_error err = pf_process_records_begin(ctx, outputDir.c_str(),
                                            sampleName.c_str(), &stream);
    if (err != PF_OK) {
        module.close();
        return err;
    }

    long long emitted = 0;
    std::vector<pf_read_record_view> records;
    std::vector<char> sequenceStorage;
    for (;;) {
        star::input::CbqReadBatchView batch;
        const star::input::InputStatus status = module.next_batch(&batch, &inputError);
        if (status == star::input::InputStatus::Error) {
            pf_process_records_abort(stream);
            module.close();
            if (errorMessage != nullptr) {
                *errorMessage = "CBQ feature assignment read error: " + inputError;
            }
            return PF_ERR_IO_ERROR;
        }
        if (status == star::input::InputStatus::End) {
            break;
        }

        uint32_t recordsToProcess = batch.record_count;
        if (maxReads > 0) {
            const long long remaining = maxReads - emitted;
            if (remaining <= 0) {
                break;
            }
            if (static_cast<long long>(recordsToProcess) > remaining) {
                recordsToProcess = static_cast<uint32_t>(remaining);
            }
        }
        if (recordsToProcess == 0) {
            continue;
        }

        records.clear();
        records.resize(recordsToProcess);
        const size_t requiredSequenceStorage =
            static_cast<size_t>(recordsToProcess) * 2U * kPfLineLength;
        if (sequenceStorage.size() < requiredSequenceStorage) {
            sequenceStorage.resize(requiredSequenceStorage);
        }
        for (uint32_t i = 0; i < recordsToProcess; ++i) {
            const star::input::CbqReadView& view = batch.records[i];
            if (view.segment_count < 2 || view.segments == nullptr) {
                pf_process_records_abort(stream);
                module.close();
                if (errorMessage != nullptr) {
                    *errorMessage = "CBQ feature assignment expects paired CBQ records";
                }
                return PF_ERR_PARSE_ERROR;
            }

            pf_read_record_view rec = {};
            char* barcodeSequence = sequenceStorage.data() +
                (static_cast<size_t>(i) * 2U * kPfLineLength);
            char* featureSequence = barcodeSequence + kPfLineLength;
            size_t barcodeLength = 0;
            size_t featureLength = 0;
            if (!star::input::materialize_cbq_segment_sequence_to_buffer(
                    view.segments[0], barcodeSequence, kPfSequenceCapacity,
                    &barcodeLength, &inputError) ||
                !star::input::materialize_cbq_segment_sequence_to_buffer(
                    view.segments[1], featureSequence, kPfSequenceCapacity,
                    &featureLength, &inputError)) {
                pf_process_records_abort(stream);
                module.close();
                if (errorMessage != nullptr) {
                    *errorMessage = "CBQ feature assignment sequence decode error: " + inputError;
                }
                return PF_ERR_PARSE_ERROR;
            }

            rec.barcode_sequence = pfViewFromBuffer(barcodeSequence, barcodeLength);
            rec.barcode_quality = pfViewFromCbqSpan(view.segments[0].quality);
            rec.feature_sequence = pfViewFromBuffer(featureSequence, featureLength);
            rec.feature_quality = pfViewFromCbqSpan(view.segments[1].quality);
            records[i] = rec;
        }

        err = pf_process_record_views(stream, records.data(), recordsToProcess);
        if (err != PF_OK) {
            pf_process_records_abort(stream);
            module.close();
            return err;
        }
        emitted += static_cast<long long>(recordsToProcess);
    }

    module.close();
    return pf_process_records_end(stream, stats);
}

static void writeApiRunSummary(const string& assignOut,
                               const WhitelistNormalizationResult& whitelistInfo,
                               const string& featureRef,
                               const string& fastqDir,
                               const string& inputFormat,
                               const string& cbqModeRequested,
                               const string& cbqModeEffective,
                               const string& cbqModeFallbackReason,
                               const AssignOptions& options,
                               const pf_stats& stats,
                               const ThreadControl::MapPermitSnapshot* permitBefore,
                               const ThreadControl::MapPermitSnapshot* permitAfter) {
    std::ofstream out((assignOut + "/assignBarcodes.api_run.txt").c_str());
    if (!out.is_open()) {
        return;
    }

    out << "mode=in_process_pf_api\n";
    out << "whitelist=" << whitelistInfo.normalizedPath << "\n";
    out << "whitelist_source=" << whitelistInfo.sourcePath << "\n";
    out << "whitelist_has_two_columns=" << (whitelistInfo.hasTwoColumnSource ? 1 : 0) << "\n";
    out << "whitelist_assignment_namespace=" << whitelistInfo.assignmentNamespace << "\n";
    out << "whitelist_namespace_confidence=" << (whitelistInfo.namespaceConfidence ? 1 : 0) << "\n";
    out << "whitelist_normalized_rows=" << whitelistInfo.normalizedRowCount << "\n";
    out << "feature_ref=" << featureRef << "\n";
    out << "input_format=" << inputFormat << "\n";
    out << "cbq_mode_requested=" << cbqModeRequested << "\n";
    out << "cbq_mode_effective=" << cbqModeEffective << "\n";
    out << "cbq_mode_fallback_reason=" << cbqModeFallbackReason << "\n";
    out << "fastq_dir=" << fastqDir << "\n";
    out << "maxHammingDistance=" << options.maxHammingDistance << "\n";
    out << "featureConstantOffset=" << options.featureConstantOffset << "\n";
    out << "limitSearch=" << options.limitSearch << "\n";
    out << "minCounts=" << options.minCounts << "\n";
    out << "maxBarcodeMismatches=" << options.maxBarcodeMismatches << "\n";
    out << "featureN=" << options.featureN << "\n";
    out << "barcodeN=" << options.barcodeN << "\n";
    out << "consumerThreadsPerSet=" << options.consumerThreadsPerSet << "\n";
    out << "searchThreads=" << options.searchThreads << "\n";
    out << "readBufferLines=" << options.readBufferLines << "\n";
    out << "minPosterior=" << options.minPosterior << "\n";
    out << "maxReads=" << options.maxReads << "\n";
    out << "legacyCbRescue=" << (options.legacyCbRescue ? 1 : 0) << "\n";
    out << "translateNxt=" << (options.translateNxt ? 1 : 0) << "\n";
    out << "probeOnly=" << (options.probeOnly ? 1 : 0) << "\n";
    out << "skipQcOutputs=" << (options.skipQcOutputs ? 1 : 0) << "\n";
    out << "adtMexOutput=" << (options.adtMexOutput ? 1 : 0) << "\n";
    out << "output_mode=" << (options.adtMexOutput ? "adt_mex" : "default") << "\n";
    out << "enableStarDynamicPermitHooks=" << (options.enableStarDynamicPermitHooks ? 1 : 0) << "\n";
    out << "filteredBarcodesPath=" << options.filteredBarcodesPath << "\n";
    out << "stats.total_reads=" << stats.total_reads << "\n";
    out << "stats.matched_reads=" << stats.matched_reads << "\n";
    out << "stats.unmatched_reads=" << stats.unmatched_reads << "\n";
    out << "stats.total_deduped_counts=" << stats.total_deduped_counts << "\n";
    out << "stats.processing_time_sec=" << stats.processing_time_sec << "\n";

    if (permitBefore != nullptr && permitAfter != nullptr) {
        auto delta = [](uint64_t after, uint64_t before) -> uint64_t {
            return after >= before ? (after - before) : 0;
        };
        out << "dynamicPermitDelta.acquires="
            << delta(permitAfter->acquireCalls, permitBefore->acquireCalls) << "\n";
        out << "dynamicPermitDelta.retunes="
            << delta(permitAfter->retuneCalls, permitBefore->retuneCalls) << "\n";
        out << "dynamicPermitDelta.blockedAcquires="
            << delta(permitAfter->blockedAcquireCalls, permitBefore->blockedAcquireCalls) << "\n";
        out << "dynamicPermitDelta.waitTimeoutEvents="
            << delta(permitAfter->waitTimeoutEvents, permitBefore->waitTimeoutEvents) << "\n";
        out << "dynamicPermitDelta.stallWarnEvents="
            << delta(permitAfter->stallWarnEvents, permitBefore->stallWarnEvents) << "\n";
        out << "dynamicPermitDelta.waitNs="
            << delta(permitAfter->waitNsTotal, permitBefore->waitNsTotal) << "\n";
        out << "dynamicPermitDelta.workUnits="
            << delta(permitAfter->workUnitsTotal, permitBefore->workUnitsTotal) << "\n";
        out << "dynamicPermitDelta.workBytes="
            << delta(permitAfter->workBytesTotal, permitBefore->workBytesTotal) << "\n";
        out << "dynamicPermitDelta.workNs="
            << delta(permitAfter->workNsTotal, permitBefore->workNsTotal) << "\n";
        out << "dynamicPermitDelta.map.acquires="
            << delta(permitAfter->mapDomain.acquireCalls, permitBefore->mapDomain.acquireCalls) << "\n";
        out << "dynamicPermitDelta.map.waitNs="
            << delta(permitAfter->mapDomain.waitNsTotal, permitBefore->mapDomain.waitNsTotal) << "\n";
        out << "dynamicPermitDelta.map.workUnits="
            << delta(permitAfter->mapDomain.workUnitsTotal, permitBefore->mapDomain.workUnitsTotal) << "\n";
        out << "dynamicPermitDelta.map.workBytes="
            << delta(permitAfter->mapDomain.workBytesTotal, permitBefore->mapDomain.workBytesTotal) << "\n";
        out << "dynamicPermitDelta.map.workNs="
            << delta(permitAfter->mapDomain.workNsTotal, permitBefore->mapDomain.workNsTotal) << "\n";
        out << "dynamicPermitDelta.feature.acquires="
            << delta(permitAfter->featureDomain.acquireCalls, permitBefore->featureDomain.acquireCalls) << "\n";
        out << "dynamicPermitDelta.feature.waitNs="
            << delta(permitAfter->featureDomain.waitNsTotal, permitBefore->featureDomain.waitNsTotal) << "\n";
        out << "dynamicPermitDelta.feature.workUnits="
            << delta(permitAfter->featureDomain.workUnitsTotal, permitBefore->featureDomain.workUnitsTotal) << "\n";
        out << "dynamicPermitDelta.feature.workBytes="
            << delta(permitAfter->featureDomain.workBytesTotal, permitBefore->featureDomain.workBytesTotal) << "\n";
        out << "dynamicPermitDelta.feature.workNs="
            << delta(permitAfter->featureDomain.workNsTotal, permitBefore->featureDomain.workNsTotal) << "\n";
        out << "dynamicPermitAfter.configuredPermits=" << permitAfter->configuredPermits << "\n";
        out << "dynamicPermitAfter.targetPermits=" << permitAfter->targetPermits << "\n";
        out << "dynamicPermitAfter.availablePermits=" << permitAfter->availablePermits << "\n";
        out << "dynamicPermitAfter.inUsePermits=" << permitAfter->inUsePermits << "\n";
        out << "dynamicPermitAfter.waiters.current=" << permitAfter->currentWaiters << "\n";
        out << "dynamicPermitAfter.waiters.max=" << permitAfter->maxWaiters << "\n";
        out << "dynamicPermitAfter.blockedAcquires=" << permitAfter->blockedAcquireCalls << "\n";
        out << "dynamicPermitAfter.waitTimeoutEvents=" << permitAfter->waitTimeoutEvents << "\n";
        out << "dynamicPermitAfter.stallWarnEvents=" << permitAfter->stallWarnEvents << "\n";
        out << "dynamicPermitAfter.lastReleaseAgoNs=" << permitAfter->lastReleaseAgoNs << "\n";
    }
}

} // namespace

WhitelistNormalizationResult normalizeWhitelistForAssign(const string& whitelistPath,
                                                         const string& assignOut) {
    return normalizeWhitelistInternal(whitelistPath, assignOut);
}

WhitelistNormalizationResult normalizeWhitelistToNamespace(
    const string& whitelistPath,
    const string& assignOut,
    const string& desiredNamespace)
{
    WhitelistNormalizationResult base = normalizeWhitelistInternal(whitelistPath, assignOut);

    if (desiredNamespace != "NXT" && desiredNamespace != "TRU") {
        return base;
    }
    if (base.assignmentNamespace == desiredNamespace) {
        return base;
    }
    if (base.assignmentNamespace != "NXT" && base.assignmentNamespace != "TRU") {
        return base;
    }

    const string srcPath = base.normalizedPath;
    const string translatedPath = assignOut + "/whitelist.ns_" + desiredNamespace + ".txt";

    std::ifstream in(srcPath.c_str());
    if (!in.is_open()) {
        return base;
    }
    std::ofstream out(translatedPath.c_str());
    if (!out.is_open()) {
        return base;
    }

    string line;
    uint64 emitted = 0;
    while (std::getline(in, line)) {
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t end = line.find_first_of("\t, \r\n", first);
        string token = (end == string::npos) ? line.substr(first)
                                             : line.substr(first, end - first);
        if (!isValidBarcodeSeq(token)) {
            continue;
        }
        out << translateNxtMiddleTwoBases(token) << "\n";
        ++emitted;
    }

    if (emitted == 0) {
        return base;
    }

    WhitelistNormalizationResult result = base;
    result.normalizedPath = translatedPath;
    result.assignmentNamespace = desiredNamespace;
    result.normalizedRowCount = emitted;
    return result;
}

AssignResult runAssignBarcodes(const string& whitelist,
                     const string& featureRef, const string& fastqDir,
                     const string& assignOut,
                     const AssignOptions& options) {
    vector<string> cbqPaths;
    const bool useCbqInput = resolveCbqSources(fastqDir, &cbqPaths);
    const string cbqModeRequested = normalizeCbqMode(options.cbqMode);
    if (cbqModeRequested != "auto" &&
        cbqModeRequested != "stream" &&
        cbqModeRequested != "range") {
        throw runtime_error("Invalid CBQ assignment mode '" + options.cbqMode +
                            "'; expected auto, stream, or range");
    }

    // Check inputs exist
    if (!fileExists(whitelist)) {
        ostringstream err;
        err << "Whitelist file not found: " << whitelist;
        throw runtime_error(err.str());
    }
    if (!fileExists(featureRef)) {
        ostringstream err;
        err << "Feature reference file not found: " << featureRef;
        throw runtime_error(err.str());
    }
    if (!useCbqInput && !dirExists(fastqDir)) {
        ostringstream err;
        err << "FASTQ directory or CBQ source not found: " << fastqDir;
        throw runtime_error(err.str());
    }
    
    // Create output directory
    string cmd = "mkdir -p \"" + assignOut + "\"";
    int ret = system(cmd.c_str());
    if (ret != 0) {
        ostringstream err;
        err << "Failed to create output directory: " << assignOut;
        throw runtime_error(err.str());
    }

    WhitelistNormalizationResult whitelistInfo = normalizeWhitelistInternal(whitelist, assignOut);
    const string whitelistForAssign = whitelistInfo.normalizedPath;

    pf_config* cfg = pf_config_create();
    if (cfg == nullptr) {
        throw runtime_error("Failed to create pf_api config");
    }
    applyAssignOptions(cfg, options);

    pf_context* ctx = pf_init(cfg);
    pf_config_destroy(cfg);
    if (ctx == nullptr) {
        throw runtime_error("Failed to initialize pf_api context");
    }

    pf_error err = pf_load_feature_ref(ctx, featureRef.c_str());
    if (err != PF_OK) {
        string msg = pfErrorMessage(ctx, err, "pf_load_feature_ref");
        pf_destroy(ctx);
        throw runtime_error(msg);
    }

    err = pf_load_whitelist(ctx, whitelistForAssign.c_str());
    if (err != PF_OK) {
        string msg = pfErrorMessage(ctx, err, "pf_load_whitelist");
        pf_destroy(ctx);
        throw runtime_error(msg);
    }

    if (!options.filteredBarcodesPath.empty()) {
        if (!fileExists(options.filteredBarcodesPath)) {
            pf_destroy(ctx);
            throw runtime_error("Filtered barcodes file not found: " + options.filteredBarcodesPath);
        }
        err = pf_load_filtered_barcodes(ctx, options.filteredBarcodesPath.c_str());
        if (err != PF_OK) {
            string msg = pfErrorMessage(ctx, err, "pf_load_filtered_barcodes");
            pf_destroy(ctx);
            throw runtime_error(msg);
        }
    }

    ThreadControl::MapPermitSnapshot permitBefore{};
    ThreadControl::MapPermitSnapshot permitAfter{};
    const bool capturePermitDelta = options.enableStarDynamicPermitHooks && g_threadChunks.mapPermitEnabled();
    if (capturePermitDelta) {
        permitBefore = g_threadChunks.mapPermitSnapshot();
    }

    pf_stats stats = {};
    string cbqErrorMessage;
    string inputFormat = useCbqInput ? "cbq_stream" : "fastq";
    string cbqModeEffective = useCbqInput ? "stream" : "";
    string cbqModeFallbackReason;
    if (useCbqInput) {
        const string sampleName = options.sampleName.empty() ? "sample" : options.sampleName;
        bool runStreamCbq = (cbqModeRequested == "stream");
        if (!runStreamCbq) {
            string unsupportedReason;
            if (!directCbqSettingsSupported(options, &unsupportedReason)) {
                cbqModeFallbackReason = unsupportedReason;
                if (cbqModeRequested == "range") {
                    pf_destroy(ctx);
                    throw runtime_error("CBQ range mode is unsupported for this PF configuration: " +
                                        unsupportedReason);
                }
                runStreamCbq = true;
            } else {
                const star::input::InputSourcePlan preflightPlan = makeCbqPairPlan(cbqPaths);
                vector<CbqLaneCount> preflightLaneCounts;
                uint64_t preflightRecords = 0;
                string preflightError;
                if (!inspectCbqLaneCounts(preflightPlan, options.maxReads,
                                          &preflightLaneCounts, &preflightRecords,
                                          &preflightError)) {
                    cbqModeFallbackReason = "cbq_index_unavailable: " + preflightError;
                    if (cbqModeRequested == "range") {
                        pf_destroy(ctx);
                        throw runtime_error("CBQ range mode requires indexed CBQ input: " +
                                            preflightError);
                    }
                    runStreamCbq = true;
                } else {
                    err = processCbqSourcesDirect(ctx, cbqPaths, assignOut, sampleName,
                                                  options.maxReads,
                                                  options.consumerThreadsPerSet > 0
                                                      ? options.consumerThreadsPerSet
                                                      : 1,
                                                  &stats, &cbqErrorMessage);
                    inputFormat = "cbq_range";
                    cbqModeEffective = "range";
                }
            }
        }
        if (runStreamCbq) {
            err = processCbqSources(ctx, cbqPaths, assignOut, sampleName,
                                    options.maxReads, &stats, &cbqErrorMessage);
            inputFormat = "cbq_stream";
            cbqModeEffective = "stream";
        }
    } else if (options.useSplitReadLayout) {
        pf_split_read_metrics splitMetrics = {};
        err = pf_process_split_fastq_dir(ctx, fastqDir.c_str(), assignOut.c_str(),
                                         &stats, &splitMetrics);
        inputFormat = "split_read";
    } else {
        err = pf_process_fastq_dir(ctx, fastqDir.c_str(), assignOut.c_str(), &stats);
    }
    if (err != PF_OK) {
        string msg;
        if (useCbqInput && !cbqErrorMessage.empty()) {
            msg = "pf_api failed at process_cbq_sources (" + pfErrorCodeString(err) + "): " +
                  cbqErrorMessage;
        } else {
            msg = pfErrorMessage(ctx, err,
                                 useCbqInput ? "process_cbq_sources" : "pf_process_fastq_dir");
        }
        pf_destroy(ctx);
        throw runtime_error(msg);
    }

    if (capturePermitDelta) {
        permitAfter = g_threadChunks.mapPermitSnapshot();
    }

    writeApiRunSummary(
        assignOut,
        whitelistInfo,
        featureRef,
        fastqDir,
        inputFormat,
        cbqModeRequested,
        cbqModeEffective,
        cbqModeFallbackReason,
        options,
        stats,
        capturePermitDelta ? &permitBefore : nullptr,
        capturePermitDelta ? &permitAfter : nullptr
    );

    AssignResult result;
    result.returnCode = 0;
    result.whitelistNormalization = whitelistInfo;
    result.inputFormat = inputFormat;
    result.cbqModeRequested = cbqModeRequested;
    result.cbqModeEffective = cbqModeEffective;
    result.cbqModeFallbackReason = cbqModeFallbackReason;
    if (options.autodetectChemistry) {
        const char* mode = pf_get_detected_match_mode(ctx);
        result.detectedMatchMode = (mode != nullptr) ? mode : "UNKNOWN";
    }

    pf_destroy(ctx);

    return result;
}

int processFeatureLibraries(const PfMultiConfig::Config& config,
                            const string& whitelist,
                            const string& featureRef, const map<string, string>& fastqMap,
                            const string& fastqRoot, const string& outPrefix,
                            const vector<string>& featureTypes,
                            const AssignOptions& options) {
    for (const auto& featureType : featureTypes) {
        vector<PfMultiConfig::LibraryEntry> libs = config.getFeatureLibraries(featureType);
        if (libs.empty()) {
            continue; // Skip if no libraries of this type
        }
        
        // Process first library of this type (could extend to handle multiple)
        const auto& lib = libs[0];
        string resolvedFastq = PfMultiConfig::resolveFastqDir(lib.fastqs, fastqRoot, fastqMap);
        
        // Create output directory name from feature type
        string featureTypeDir = featureType;
        // Replace spaces and special chars with underscores
        for (char& c : featureTypeDir) {
            if (c == ' ' || c == '/' || c == '\\') {
                c = '_';
            }
        }
        string assignOut = outPrefix + "/cr_assign/" + featureTypeDir;
        
        try {
            AssignResult result = runAssignBarcodes(whitelist, featureRef, resolvedFastq, assignOut, options);
            if (result.returnCode != 0) {
                cerr << "ERROR processing " << featureType << ": returnCode=" << result.returnCode << endl;
                return 1;
            }
        } catch (const exception& e) {
            cerr << "ERROR processing " << featureType << ": " << e.what() << endl;
            return 1;
        }
    }
    
    return 0;
}

void waitForFeaturePermitInterface(bool hooksEnabled) {
    if (!hooksEnabled) {
        return;
    }
    if (g_threadChunks.mapPermitEnabled()) {
        return;
    }
    constexpr int kPollIntervalMs = 10;
    constexpr int kTimeoutMs = 30 * 60 * 1000;
    const auto deadline = std::chrono::steady_clock::now() +
                          std::chrono::milliseconds(kTimeoutMs);
    while (!g_threadChunks.mapPermitEnabled()) {
        if (std::chrono::steady_clock::now() >= deadline) {
            throw std::runtime_error(
                "table import: timed out waiting for dynamic FEATURE permit interface");
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(kPollIntervalMs));
    }
}

bool acquireFeaturePermitChunk(bool enabled, uint64_t& waitNsOut) {
    if (!enabled || !g_threadChunks.mapPermitEnabled()) {
        waitNsOut = 0;
        return false;
    }
    const auto t0 = std::chrono::steady_clock::now();
    g_threadChunks.mapPermitAcquireForDomain(ThreadControl::PermitDomain::FEATURE);
    const auto t1 = std::chrono::steady_clock::now();
    waitNsOut = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    return true;
}

void releaseFeaturePermitChunk(bool enabled,
                             uint64_t waitNs,
                             uint64_t workUnits,
                             uint64_t workBytes) {
    if (!enabled || !g_threadChunks.mapPermitEnabled()) {
        return;
    }
    g_threadChunks.mapPermitReleaseForDomain(
        ThreadControl::PermitDomain::FEATURE, waitNs, workUnits, workBytes, 0);
}

bool featurePermitTelemetryEnabled(bool hooksEnabled) {
    return hooksEnabled && g_threadChunks.mapPermitEnabled();
}

ThreadControl::MapPermitSnapshot featurePermitSnapshot() {
    return g_threadChunks.mapPermitSnapshot();
}

} // namespace PfMultiAssign

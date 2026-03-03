#include "PfMultiAssign.h"
#include "GlobalVariables.h"
#include "pf_api.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <stdexcept>
using std::cerr;
using std::endl;

namespace {
ThreadControl::PermitHookContext kFeaturePermitHookContext{ThreadControl::PermitDomain::FEATURE};
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
    if (options.enableStarDynamicPermitHooks) {
        pf_config_set_permit_hooks(
            cfg,
            pfStarDynamicPermitAcquire,
            pfStarDynamicPermitRelease,
            &kFeaturePermitHookContext
        );
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

static void writeApiRunSummary(const string& assignOut,
                               const WhitelistNormalizationResult& whitelistInfo,
                               const string& featureRef,
                               const string& fastqDir,
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
    out << "minPosterior=" << options.minPosterior << "\n";
    out << "maxReads=" << options.maxReads << "\n";
    out << "legacyCbRescue=" << (options.legacyCbRescue ? 1 : 0) << "\n";
    out << "translateNxt=" << (options.translateNxt ? 1 : 0) << "\n";
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

AssignResult runAssignBarcodes(const string& whitelist,
                     const string& featureRef, const string& fastqDir,
                     const string& assignOut,
                     const AssignOptions& options) {
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
    if (!dirExists(fastqDir)) {
        ostringstream err;
        err << "FASTQ directory not found: " << fastqDir;
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
    err = pf_process_fastq_dir(ctx, fastqDir.c_str(), assignOut.c_str(), &stats);
    if (err != PF_OK) {
        string msg = pfErrorMessage(ctx, err, "pf_process_fastq_dir");
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
        options,
        stats,
        capturePermitDelta ? &permitBefore : nullptr,
        capturePermitDelta ? &permitAfter : nullptr
    );

    AssignResult result;
    result.returnCode = 0;
    result.whitelistNormalization = whitelistInfo;
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

} // namespace PfMultiAssign

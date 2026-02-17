#include "PfMultiAssign.h"
#include "pf_api.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cctype>
#include <sys/stat.h>
#include <stdexcept>
using std::cerr;
using std::endl;

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

static string normalizeWhitelistIfNeeded(const string& whitelistPath, const string& assignOut) {
    if (!looksLikeMultiColumnWhitelist(whitelistPath)) {
        return whitelistPath;
    }

    std::ifstream in(whitelistPath.c_str());
    if (!in.is_open()) {
        return whitelistPath;
    }

    string normalizedPath = assignOut + "/whitelist.normalized.txt";
    std::ofstream out(normalizedPath.c_str());
    if (!out.is_open()) {
        return whitelistPath;
    }

    string line;
    size_t emitted = 0;
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
        return whitelistPath;
    }
    return normalizedPath;
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
                               const string& whitelistPath,
                               const string& featureRef,
                               const string& fastqDir,
                               const AssignOptions& options,
                               const pf_stats& stats) {
    std::ofstream out((assignOut + "/assignBarcodes.api_run.txt").c_str());
    if (!out.is_open()) {
        return;
    }

    out << "mode=in_process_pf_api\n";
    out << "whitelist=" << whitelistPath << "\n";
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
    out << "filteredBarcodesPath=" << options.filteredBarcodesPath << "\n";
    out << "stats.total_reads=" << stats.total_reads << "\n";
    out << "stats.matched_reads=" << stats.matched_reads << "\n";
    out << "stats.unmatched_reads=" << stats.unmatched_reads << "\n";
    out << "stats.total_deduped_counts=" << stats.total_deduped_counts << "\n";
    out << "stats.processing_time_sec=" << stats.processing_time_sec << "\n";
}

} // namespace

int runAssignBarcodes(const string& whitelist,
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

    string whitelistForAssign = normalizeWhitelistIfNeeded(whitelist, assignOut);

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

    pf_stats stats = {};
    err = pf_process_fastq_dir(ctx, fastqDir.c_str(), assignOut.c_str(), &stats);
    if (err != PF_OK) {
        string msg = pfErrorMessage(ctx, err, "pf_process_fastq_dir");
        pf_destroy(ctx);
        throw runtime_error(msg);
    }

    writeApiRunSummary(assignOut, whitelistForAssign, featureRef, fastqDir, options, stats);
    pf_destroy(ctx);

    return 0;
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
            runAssignBarcodes(whitelist, featureRef, resolvedFastq, assignOut, options);
        } catch (const exception& e) {
            cerr << "ERROR processing " << featureType << ": " << e.what() << endl;
            return 1;
        }
    }
    
    return 0;
}

} // namespace PfMultiAssign

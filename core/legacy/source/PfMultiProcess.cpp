#include "PfMultiProcess.h"
#include "Solo.h"
#include "PfMultiConfig.h"
#include "PfMultiAssign.h"
#include "PfMultiMexStub.h"
#include "PfMultiMerge.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"
#include "serviceFuns.cpp"
#include "TimeFunctions.h"
#include "call_features.h"
#include "MexWriter.h"
#include "PfMultiFeatureSpecs.h"
#include <sys/stat.h>
#include <dirent.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <cctype>
#include <cmath>
#include <chrono>
#include <thread>
#include <atomic>
#include <limits>
#include <set>
#include <unordered_set>
#include <zlib.h>
using std::cerr;
using std::endl;

struct AnchoredReadEstimate {
    bool valid = false;
    uint64_t estimatedReads = 0;
    uint64_t anchorReads = 0;
    uint64_t anchorBytes = 0;
    uint64_t totalBytes = 0;
    string anchorFile;
    size_t fileCount = 0;
};

struct PfPreparedFeatureLibrary {
    string libraryType;
    string featureRefType;
    string resolvedFastq;
    string sampleName;
    string libraryId;             // from star_library_id (auto-generated if absent)
    string assignOut;
    string featureRefPath;
    string whitelistPath;
    string effectiveChem;
    string resolvedChemRequest;   // "auto", "nxt", or "tru" — after column > global resolution
    bool explicitChem = false;    // true when star_chemistry column provided NXT/TRU
    bool usedFilteredRef = false;
    int starMaxHamming = -1;    // Per-library max Hamming override (-1 = use global)
    bool useSplitReadLayout = false;
    bool adtMexOutput = false;
    pf_split_read_layout splitReadLayout = {};
    string starBarcodeOutputMap;
    string starFeatureSearchMode;
    AnchoredReadEstimate featureEstimate;
};

struct FilteredBarcodeNormalizationStats {
    uint64_t inSet = 0;
    uint64_t translatedToSet = 0;
    uint64_t unmatched = 0;
    uint64_t dedupDropped = 0;
    uint64_t outputCount = 0;
};

struct PfLibraryNamespaceContext {
    string libraryId;
    string effectiveReadNamespace;       // NXT | TRU
    string assignmentWhitelistNamespace; // NXT | TRU | UNKNOWN
    bool translateNxtForAssign = false;
    string outputNamespace;              // NXT | TRU (global output namespace)
    string autoDetectMatchMode;          // RAW_MATCH | TRANSLATED_MATCH | AMBIGUOUS | UNKNOWN
    bool isChemistryExplicit = false;
    bool isNamespaceConfident = false;
    bool allowUnionWhitelist = false;
    FilteredBarcodeNormalizationStats normalizationStats;
};

struct PfMultiPreparedContext {
    PfMultiConfig::Config config;
    string featureRef;
    string whitelist;
    string soloWhitelist;
    string requestedChem;
    string requestedOutputChem;
    string inferredChem;
    string inferredReason;
    bool inferredChemConfident = false;
    string effectiveChem;
    bool hasTwoColumnWhitelist = false;
    string soloInferredChem;
    string soloInferredReason;
    bool soloInferredChemConfident = false;
    string soloEffectiveChem;
    bool soloHasTwoColumnWhitelist = false;
    string soloOutputNamespace;  // actual namespace of Solo's raw MEX barcodes
    string outputChem;
    string outPrefix;
    string crAssignRoot;
    AnchoredReadEstimate mapEstimate;
    vector<PfPreparedFeatureLibrary> featureLibraries;
    string prepLog;
};

class PfMultiPreloadHandle;

static PfMultiPreparedContext resolvePfMultiPreparedContext(
    Parameters& P,
    const std::shared_ptr<PfMultiPreloadHandle>& preload);

class PfMultiPreloadHandle {
public:
    ~PfMultiPreloadHandle() {
        if (worker.joinable()) {
            worker.join();
        }
    }

private:
    friend std::shared_ptr<PfMultiPreloadHandle> startPfMultiConfigPreload(const Parameters& P);
    friend PfMultiPreparedContext resolvePfMultiPreparedContext(
        Parameters& P,
        const std::shared_ptr<PfMultiPreloadHandle>& preload);
    friend std::shared_ptr<PfMultiAssignPhaseResult> runPfMultiAssignPhase(
        Parameters& P,
        const std::shared_ptr<PfMultiPreloadHandle>& preload);

    std::thread worker;
    bool started = false;
    bool finished = false;
    bool success = false;
    double preloadSec = 0.0;
    string error;
    PfMultiPreparedContext context;
};

struct PfMultiFeatureRun {
    string featureType;
    string assignOut;
    string featureRefPath;
    string whitelistPath;
    string barcodeOutputMapPath;
    string effectiveChem;
    string effectiveReadNamespace;
    string assignmentWhitelistNamespace;
    bool translateNxtForAssign = false;
    string featureMexOutputNamespace;
    string outputNamespace;
    bool namespaceConfident = false;
    string detectedMatchMode;
    PfMultiAssign::WhitelistNormalizationResult whitelistNormalization;
    FilteredBarcodeNormalizationStats normalizationStats;
    string libraryId;
    string sampleName;
    string resolvedFastq;
    string resolvedChemRequest;
    bool explicitChem = false;
    bool adtMexOutput = false;
    int returnCode = 0;
};

struct PfMultiFeatureMexEntry {
    PfMultiMerge::MexData data;
    string effectiveChem;
    string featureMexOutputNamespace;
    string featureType;
    string libraryId;
};

struct PfMultiAssignPhaseResult {
    PfMultiPreparedContext prepared;
    vector<PfMultiFeatureRun> featureRuns;
    vector<PfMultiFeatureMexEntry> featureMexEntries;
    bool usedExplicitAssignFilteredBarcodes = false;
};

class PfFeatureAssignHandle {
public:
    ~PfFeatureAssignHandle() {
        if (worker.joinable()) {
            worker.join();
        }
    }

    void join() {
        if (worker.joinable()) {
            worker.join();
        }
    }

    bool isFinished() const { return finished; }
    bool isSuccess() const { return success; }
    const string& getError() const { return error; }
    const std::shared_ptr<PfMultiAssignPhaseResult>& getResult() const { return result; }

private:
    friend std::shared_ptr<PfFeatureAssignHandle> startPfFeatureAssignment(
        Parameters& P,
        const std::shared_ptr<PfMultiPreloadHandle>& preload);
    friend int processPfMultiConfig(Parameters& P,
                                    const Solo* solo,
                                    const std::shared_ptr<PfMultiPreloadHandle>& preload,
                                    const std::shared_ptr<PfFeatureAssignHandle>& assignHandle);

    std::thread worker;
    bool started = false;
    bool finished = false;
    bool success = false;
    string error;
    std::shared_ptr<PfMultiAssignPhaseResult> result;
};

using PfMultiFeatureSpecs::FeatureSpec;
using PfMultiFeatureSpecs::buildFeatureSpecsFromConfig;

namespace {

static string normalizeType(const string& input) {
    string out;
    out.reserve(input.size());
    for (unsigned char c : input) {
        if (std::isalnum(c)) {
            out.push_back(static_cast<char>(std::tolower(c)));
        }
    }
    return out;
}

static string trimCopy(const string& input) {
    size_t first = input.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        return "";
    }
    size_t last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

static string lowerCopy(string input) {
    std::transform(input.begin(), input.end(), input.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return input;
}

static bool isUnsetToken(const string& input) {
    string token = lowerCopy(trimCopy(input));
    return token.empty() || token == "-" || token == "none";
}

static string basenameOf(const string& path);

static string upperCopy(string input) {
    std::transform(input.begin(), input.end(), input.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return input;
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

static string normalizeBarcodeToTru(string barcode, const string& fromChem) {
    string norm = upperCopy(fromChem);
    if (norm == "NXT" && barcode.size() >= 9) {
        barcode[7] = complementBase(barcode[7]);
        barcode[8] = complementBase(barcode[8]);
    }
    return barcode;
}

static void normalizeMexBarcodesToTru(PfMultiMerge::MexData& data, const string& fromChem) {
    for (auto& bc : data.barcodes) {
        bc = normalizeBarcodeToTru(bc, fromChem);
    }
}

static void normalizeBarcodeVecToTru(vector<string>& barcodes, const string& fromChem) {
    for (auto& bc : barcodes) {
        bc = normalizeBarcodeToTru(bc, fromChem);
    }
}

static vector<string> splitWhitelistColumns(const string& line) {
    vector<string> fields;
    string current;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '\t' || c == ',' || std::isspace(static_cast<unsigned char>(c))) {
            if (!current.empty()) {
                fields.push_back(current);
                current.clear();
                if (fields.size() >= 2) {
                    break;
                }
            }
        } else {
            current.push_back(c);
        }
    }
    if (!current.empty() && fields.size() < 2) {
        fields.push_back(current);
    }
    return fields;
}

static string detectChemistryFromWhitelistPath(const string& whitelistPath, string& reason,
                                                bool& confident) {
    confident = false;
    if (isUnsetToken(whitelistPath)) {
        reason = "whitelist path is unset";
        return "UNKNOWN";
    }

    bool complementRuleMatched = false;
    string complementRuleReason;

    // Prefer content-based detection first for 2-column NXT/TRU translation lists.
    std::ifstream in(whitelistPath.c_str());
    if (in.is_open()) {
        const size_t kMaxScanRows = 10000;
        size_t scanned = 0;
        size_t twoColumnRows = 0;
        size_t nxtRuleRows = 0;
        string line;

        while (scanned < kMaxScanRows && std::getline(in, line)) {
            string trimmed = trimCopy(line);
            if (trimmed.empty() || trimmed[0] == '#') {
                continue;
            }

            vector<string> fields = splitWhitelistColumns(trimmed);
            if (fields.empty()) {
                continue;
            }
            scanned++;

            if (fields.size() >= 2 && isValidBarcodeSeq(fields[0]) && isValidBarcodeSeq(fields[1])) {
                twoColumnRows++;
                const string col1 = upperCopy(fields[0]);
                const string col2 = upperCopy(fields[1]);
                // Symmetric: complement is an involution so translate(A)==B iff translate(B)==A,
                // meaning column order does not matter.
                if (translateNxtMiddleTwoBases(col1) == col2) {
                    nxtRuleRows++;
                }
            }
        }

        if (twoColumnRows > 0) {
            const double nxtFraction = static_cast<double>(nxtRuleRows) / static_cast<double>(twoColumnRows);
            std::ostringstream msg;
            msg << "whitelist has two barcode columns (" << twoColumnRows
                << " sampled rows); center-2bp complement rule matched "
                << nxtRuleRows << "/" << twoColumnRows
                << " (" << static_cast<int>(nxtFraction * 100.0 + 0.5) << "%)";
            if (nxtFraction >= 0.8) {
                complementRuleMatched = true;
                complementRuleReason = msg.str();
            } else {
                msg << " — below 80% threshold, falling through to filename/default";
                reason = msg.str();
            }
        }
    }

    // Filename hints provide orientation confirmation for 2-column whitelists.
    string base = lowerCopy(basenameOf(whitelistPath));
    if (base.find("nxt") != string::npos) {
        reason = complementRuleMatched
            ? (complementRuleReason + "; confirmed by filename 'nxt'")
            : "whitelist filename contains 'nxt'";
        confident = true;
        return "NXT";
    }
    if (base.find("tru") != string::npos) {
        reason = complementRuleMatched
            ? (complementRuleReason + "; confirmed by filename 'tru'")
            : "whitelist filename contains 'tru'";
        confident = true;
        return "TRU";
    }

    // Complement rule matched but no filename confirmation. The rule is
    // symmetric (translate(A)==B iff translate(B)==A), so content alone
    // cannot distinguish COL1=NXT from COL1=TRU orientation.
    if (complementRuleMatched) {
        reason = complementRuleReason
            + "; ERROR: no filename hint to confirm NXT/TRU orientation"
              " — cannot determine column order from content alone";
        confident = false;
        return "UNKNOWN";
    }

    reason = "no chemistry marker found in whitelist content or filename";
    confident = false;
    return "UNKNOWN";
}

static string oppositeNamespace(const string& ns) {
    if (ns == "NXT") return "TRU";
    if (ns == "TRU") return "NXT";
    return ns;
}

static bool isKnownNamespace(const string& ns) {
    return ns == "NXT" || ns == "TRU";
}

static bool whitelistHasTwoColumns(const string& whitelistPath) {
    std::ifstream in(whitelistPath.c_str());
    if (!in.is_open()) return false;
    string line;
    const int maxProbe = 20; // scan up to 20 non-comment lines before giving up
    int probed = 0;
    while (std::getline(in, line) && probed < maxProbe) {
        string trimmed = trimCopy(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        probed++;
        vector<string> fields = splitWhitelistColumns(trimmed);
        if (fields.size() >= 2 && isValidBarcodeSeq(fields[0]) && isValidBarcodeSeq(fields[1])) {
            return true;
        }
    }
    return false;
}

static void enforceUniqueBarcodesOrThrow(const vector<string>& barcodes,
                                         const string& contextLabel) {
    std::unordered_set<string> seen;
    seen.reserve(barcodes.size() * 2 + 1);
    vector<string> dupPreview;
    dupPreview.reserve(10);
    uint64_t dupCount = 0;
    for (const auto& bc : barcodes) {
        const string token = upperCopy(trimCopy(bc));
        if (token.empty()) {
            continue;
        }
        if (!seen.insert(token).second) {
            ++dupCount;
            if (dupPreview.size() < 10) {
                dupPreview.push_back(token);
            }
        }
    }
    if (dupCount == 0) {
        return;
    }
    std::ostringstream oss;
    oss << "Duplicate barcodes detected in " << contextLabel
        << " (duplicates=" << dupCount << ").";
    if (!dupPreview.empty()) {
        oss << " Example duplicates:";
        for (const auto& d : dupPreview) {
            oss << " " << d;
        }
    }
    throw runtime_error(oss.str());
}

static void enforceUniqueBarcodeFileOrThrow(const string& barcodePath,
                                            const string& contextLabel) {
    std::ifstream in(barcodePath.c_str());
    if (!in.is_open()) {
        throw runtime_error("Failed to read barcode file for uniqueness check: " + barcodePath);
    }
    std::unordered_set<string> seen;
    vector<string> dupPreview;
    dupPreview.reserve(10);
    uint64_t dupCount = 0;
    string line;
    while (std::getline(in, line)) {
        string token = upperCopy(trimCopy(line));
        if (token.empty()) {
            continue;
        }
        if (!seen.insert(token).second) {
            ++dupCount;
            if (dupPreview.size() < 10) {
                dupPreview.push_back(token);
            }
        }
    }
    if (dupCount == 0) {
        return;
    }
    std::ostringstream oss;
    oss << "Duplicate filtered barcodes detected before pf_api load in "
        << contextLabel << " (" << barcodePath << ", duplicates=" << dupCount << ").";
    if (!dupPreview.empty()) {
        oss << " Example duplicates:";
        for (const auto& d : dupPreview) {
            oss << " " << d;
        }
    }
    throw runtime_error(oss.str());
}

static string parseChemistryToken(const string& rawValue, const string& flagName) {
    string value = lowerCopy(trimCopy(rawValue));
    if (value.empty() || value == "-" || value == "none") {
        value = "auto";
    }
    if (value != "auto" && value != "nxt" && value != "tru") {
        throw runtime_error("Invalid " + flagName + " value '" + rawValue +
                            "' (allowed: auto|NXT|TRU)");
    }
    return value;
}

static string sanitizeDirName(const string& input) {
    string out = input;
    for (char& c : out) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_')) {
            c = '_';
        }
    }
    return out;
}

static string basenameOf(const string& path) {
    size_t pos = path.find_last_of("/\\");
    return (pos == string::npos) ? path : path.substr(pos + 1);
}

enum class PfPermitControllerMode {
    Off = 0,
    Shadow = 1,
    Active = 2,
    Eta = 3,
    Chunked = 4
};

static PfPermitControllerMode parsePfPermitControllerMode(const string& rawMode) {
    const string mode = lowerCopy(trimCopy(rawMode));
    if (mode == "chunked" || mode == "chunks" || mode == "4") {
        return PfPermitControllerMode::Chunked;
    }
    if (mode == "eta" || mode == "3") {
        return PfPermitControllerMode::Eta;
    }
    if (mode == "active" || mode == "2") {
        return PfPermitControllerMode::Active;
    }
    if (mode == "shadow" || mode == "1") {
        return PfPermitControllerMode::Shadow;
    }
    return PfPermitControllerMode::Off;
}

static string pfPermitControllerModeName(PfPermitControllerMode mode) {
    switch (mode) {
        case PfPermitControllerMode::Shadow:
            return "shadow";
        case PfPermitControllerMode::Active:
            return "active";
        case PfPermitControllerMode::Eta:
            return "eta";
        case PfPermitControllerMode::Chunked:
            return "chunked";
        case PfPermitControllerMode::Off:
        default:
            return "off";
    }
}

static int normalizePfSearchThreadQuantum(int rawSearchThreads) {
    if (rawSearchThreads <= 1) {
        return 1;
    }
    if (rawSearchThreads <= 2) {
        return 2;
    }
    return 4;
}

static string joinIntVector(const vector<int>& values) {
    if (values.empty()) {
        return "";
    }
    std::ostringstream oss;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            oss << ",";
        }
        oss << values[i];
    }
    return oss.str();
}

static bool hasSuffixCaseInsensitive(const string& value, const string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    const size_t offset = value.size() - suffix.size();
    for (size_t i = 0; i < suffix.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(value[offset + i])) !=
            std::tolower(static_cast<unsigned char>(suffix[i]))) {
            return false;
        }
    }
    return true;
}

static bool isFastqPath(const string& path) {
    return hasSuffixCaseInsensitive(path, ".fastq") ||
           hasSuffixCaseInsensitive(path, ".fq") ||
           hasSuffixCaseInsensitive(path, ".fastq.gz") ||
           hasSuffixCaseInsensitive(path, ".fq.gz");
}

static bool isGzipPath(const string& path) {
    return hasSuffixCaseInsensitive(path, ".gz");
}

static uint64_t fileSizeBytes(const string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0 || st.st_size <= 0) {
        return 0;
    }
    return static_cast<uint64_t>(st.st_size);
}

static uint64_t countNewlinesInPlainFile(const string& path, bool* okOut) {
    ifstream in(path.c_str(), std::ios::binary);
    if (!in.is_open()) {
        if (okOut != nullptr) {
            *okOut = false;
        }
        return 0;
    }
    uint64_t lines = 0;
    char buffer[1 << 20];
    while (in.good()) {
        in.read(buffer, sizeof(buffer));
        const std::streamsize n = in.gcount();
        for (std::streamsize i = 0; i < n; ++i) {
            if (buffer[i] == '\n') {
                ++lines;
            }
        }
    }
    if (okOut != nullptr) {
        *okOut = true;
    }
    return lines;
}

static uint64_t countNewlinesInGzipFile(const string& path, bool* okOut) {
    gzFile gz = gzopen(path.c_str(), "rb");
    if (gz == nullptr) {
        if (okOut != nullptr) {
            *okOut = false;
        }
        return 0;
    }

    uint64_t lines = 0;
    char buffer[1 << 20];
    int nRead = 0;
    while ((nRead = gzread(gz, buffer, sizeof(buffer))) > 0) {
        for (int i = 0; i < nRead; ++i) {
            if (buffer[i] == '\n') {
                ++lines;
            }
        }
    }

    bool ok = (nRead >= 0);
    if (gzclose(gz) != Z_OK) {
        ok = false;
    }
    if (okOut != nullptr) {
        *okOut = ok;
    }
    return ok ? lines : 0;
}

static bool countFastqReadsExact(const string& path, uint64_t* readsOut) {
    if (readsOut == nullptr) {
        return false;
    }
    bool ok = false;
    const uint64_t lines = isGzipPath(path)
        ? countNewlinesInGzipFile(path, &ok)
        : countNewlinesInPlainFile(path, &ok);
    if (!ok || lines == 0) {
        return false;
    }
    *readsOut = lines / 4;
    return true;
}

static vector<string> listFastqFilesInDirectory(const string& dirPath) {
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
        struct stat st;
        if (stat(full.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) {
            continue;
        }
        if (isFastqPath(full)) {
            files.push_back(full);
        }
    }
    closedir(dir);
    std::sort(files.begin(), files.end());
    return files;
}

static vector<string> filterFastqFilesByMarker(const vector<string>& files, const string& marker) {
    if (marker.empty()) {
        return files;
    }
    const string markerLower = lowerCopy(marker);
    vector<string> out;
    for (const auto& path : files) {
        const string base = lowerCopy(basenameOf(path));
        if (base.find(markerLower) != string::npos) {
            out.push_back(path);
        }
    }
    return out;
}

struct PfMultiPreloadInput {
    string pfMultiConfig;
    string crFeatureRef;
    string crWhitelist;
    string crChemistry;
    string crOutputChemistry;
    string crFastqRoot;
    vector<string> crFastqMap;
    vector<vector<string>> readFilesNames;
    int soloType = 0;
    uint32 soloBarcodeRead = 0;
    vector<string> soloCbWhitelist;
    string outFileNamePrefix;
};

static PfMultiPreloadInput makePfMultiPreloadInput(const Parameters& P) {
    PfMultiPreloadInput input;
    input.pfMultiConfig = P.pfMulti.pfMultiConfig;
    input.crFeatureRef = P.pfMulti.crFeatureRef;
    input.crWhitelist = P.pfMulti.crWhitelist;
    input.crChemistry = P.pfMulti.crChemistry;
    input.crOutputChemistry = P.pfMulti.crOutputChemistry;
    input.crFastqRoot = P.pfMulti.crFastqRoot;
    input.crFastqMap = P.pfMulti.crFastqMap;
    input.readFilesNames = P.readFilesNames;
    input.soloType = P.pSolo.type;
    input.soloBarcodeRead = P.pSolo.barcodeRead;
    input.soloCbWhitelist = P.pSolo.soloCBwhitelist;
    input.outFileNamePrefix = P.outFileNamePrefix;
    return input;
}

static vector<string> selectMapPrimaryFastqs(const PfMultiPreloadInput& input) {
    vector<string> explicitR1;
    for (const auto& mateFiles : input.readFilesNames) {
        for (const auto& path : mateFiles) {
            const string base = lowerCopy(basenameOf(path));
            if (base.find("_r1_") != string::npos || base.find(".r1.") != string::npos) {
                explicitR1.push_back(path);
            }
        }
    }
    if (!explicitR1.empty()) {
        return explicitR1;
    }

    if (input.soloType != 0 && input.soloBarcodeRead < input.readFilesNames.size()) {
        return input.readFilesNames.at(input.soloBarcodeRead);
    }
    if (!input.readFilesNames.empty()) {
        return input.readFilesNames.at(0);
    }
    return {};
}

static AnchoredReadEstimate estimateReadsAnchored(const vector<string>& inputFiles,
                                                 const string& label,
                                                 ostream& logMain) {
    AnchoredReadEstimate estimate;
    vector<string> files;
    files.reserve(inputFiles.size());
    for (const auto& path : inputFiles) {
        const uint64_t sizeBytes = fileSizeBytes(path);
        if (sizeBytes > 0) {
            files.push_back(path);
            estimate.totalBytes += sizeBytes;
        }
    }
    estimate.fileCount = files.size();
    if (files.empty() || estimate.totalBytes == 0) {
        logMain << "pf-dynamic-controller: estimator[" << label
                << "] unavailable (no FASTQ files)\n";
        return estimate;
    }

    std::sort(files.begin(), files.end(),
              [](const string& a, const string& b) { return fileSizeBytes(a) < fileSizeBytes(b); });
    estimate.anchorFile = files.front();
    estimate.anchorBytes = fileSizeBytes(estimate.anchorFile);
    if (estimate.anchorBytes == 0) {
        logMain << "pf-dynamic-controller: estimator[" << label
                << "] anchor has zero size: " << estimate.anchorFile << "\n";
        return estimate;
    }

    const auto t0 = std::chrono::steady_clock::now();
    uint64_t anchorReads = 0;
    if (!countFastqReadsExact(estimate.anchorFile, &anchorReads)) {
        logMain << "pf-dynamic-controller: estimator[" << label
                << "] failed to count anchor reads: " << estimate.anchorFile << "\n";
        return estimate;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double anchorCountSec =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();

    estimate.anchorReads = anchorReads;
    const long double scaled = static_cast<long double>(estimate.anchorReads) *
                               static_cast<long double>(estimate.totalBytes) /
                               static_cast<long double>(estimate.anchorBytes);
    uint64_t estimatedReads = static_cast<uint64_t>(std::llround(scaled));
    if (estimatedReads == 0 && estimate.anchorReads > 0) {
        estimatedReads = estimate.anchorReads;
    }
    estimate.estimatedReads = estimatedReads;
    estimate.valid = (estimate.estimatedReads > 0);

    logMain << "pf-dynamic-controller: estimator[" << label
            << "] files=" << estimate.fileCount
            << ", anchor=" << estimate.anchorFile
            << ", anchorBytes=" << estimate.anchorBytes
            << ", anchorReads=" << estimate.anchorReads
            << ", totalBytes=" << estimate.totalBytes
            << ", estimatedReads=" << estimate.estimatedReads
            << ", anchorCountSec=" << anchorCountSec
            << "\n";

    return estimate;
}

static uint64_t inflateEstimateUntilNonNegativeRemaining(uint64_t estimateTotal, uint64_t done) {
    if (estimateTotal == 0) {
        estimateTotal = done;
    }
    while (estimateTotal < done) {
        const uint64_t bump = std::max<uint64_t>(1, estimateTotal / 10);
        const uint64_t next = estimateTotal + bump;
        if (next <= estimateTotal) {
            estimateTotal = done;
            break;
        }
        estimateTotal = next;
    }
    return estimateTotal;
}

static string findAssignOutputDir(const string& baseDir) {
    struct stat st;
    string matrixPath = baseDir + "/matrix.mtx";
    string matrixGz = baseDir + "/matrix.mtx.gz";
    if (stat(matrixPath.c_str(), &st) == 0 || stat(matrixGz.c_str(), &st) == 0) {
        return baseDir;
    }

    DIR* dir = opendir(baseDir.c_str());
    if (!dir) {
        return baseDir;
    }
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_name[0] == '.') {
            continue;
        }
        string sub = baseDir + "/" + entry->d_name;
        if (stat(sub.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            string subMatrix = sub + "/matrix.mtx";
            string subMatrixGz = sub + "/matrix.mtx.gz";
            if (stat(subMatrix.c_str(), &st) == 0 || stat(subMatrixGz.c_str(), &st) == 0) {
                closedir(dir);
                return sub;
            }
        }
    }
    closedir(dir);
    return baseDir;
}

static bool hasMexFiles(const string& dirPath) {
    struct stat st;
    string features = dirPath + "/features.tsv";
    string featuresGz = dirPath + "/features.tsv.gz";
    return (stat(features.c_str(), &st) == 0) || (stat(featuresGz.c_str(), &st) == 0);
}

static bool filterFeatureRefCsv(const string& inputPath, const string& featureType,
                                const string& outputPath) {
    ifstream in(inputPath);
    if (!in.is_open()) {
        return false;
    }

    string headerLine;
    if (!getline(in, headerLine)) {
        return false;
    }

    vector<string> headers;
    {
        string field;
        bool inQuotes = false;
        for (char c : headerLine) {
            if (c == '"') {
                inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                headers.push_back(field);
                field.clear();
            } else {
                field += c;
            }
        }
        headers.push_back(field);
    }

    int typeIdx = -1;
    for (size_t i = 0; i < headers.size(); ++i) {
        string h = headers[i];
        std::transform(h.begin(), h.end(), h.begin(), ::tolower);
        if (h == "feature_type" || h == "type") {
            typeIdx = static_cast<int>(i);
            break;
        }
    }

    if (typeIdx < 0) {
        return false;
    }

    string normalizedTarget = normalizeType(featureType);
    ofstream out(outputPath);
    if (!out.is_open()) {
        return false;
    }
    out << headerLine << "\n";

    string line;
    bool wroteAny = false;
    while (getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        vector<string> fields;
        string field;
        bool inQuotes = false;
        for (char c : line) {
            if (c == '"') {
                inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                fields.push_back(field);
                field.clear();
            } else {
                field += c;
            }
        }
        fields.push_back(field);
        if (static_cast<int>(fields.size()) <= typeIdx) {
            continue;
        }
        string value = fields[typeIdx];
        if (!value.empty() && value.front() == '"' && value.back() == '"') {
            value = value.substr(1, value.size() - 2);
        }
        string normalizedValue = normalizeType(value);
        if (normalizedValue == normalizedTarget) {
            out << line << "\n";
            wroteAny = true;
        }
    }
    return wroteAny;
}

static bool getFilteredBarcodesFromSolo(const Solo* solo,
                                        const Parameters& P,
                                        vector<string>& out,
                                        bool useOutputNamespace) {
    if (!solo || !solo->soloFeat) {
        return false;
    }
    int32 featIdx = -1;
    if (P.pSolo.crGexFeature == ParametersSolo::CrGexGene) {
        featIdx = P.pSolo.featureInd[SoloFeatureTypes::Gene];
    } else if (P.pSolo.crGexFeature == ParametersSolo::CrGexGeneFull) {
        featIdx = P.pSolo.featureInd[SoloFeatureTypes::GeneFull];
    } else {
        // Auto: prefer Gene, fallback to GeneFull
        featIdx = P.pSolo.featureInd[SoloFeatureTypes::Gene];
        if (featIdx < 0) {
            featIdx = P.pSolo.featureInd[SoloFeatureTypes::GeneFull];
        }
    }
    if (featIdx < 0 || solo->soloFeat[featIdx] == nullptr) {
        return false;
    }

    const SoloFeature* gex = solo->soloFeat[featIdx];
    if (gex->filteredCells.filtVecBool.empty()) {
        return false;
    }

    out.clear();
    out.reserve(gex->filteredCells.nCells);
    for (uint32 icb = 0; icb < gex->nCB; icb++) {
        if (gex->filteredCells.filtVecBool[icb]) {
            uint32 wlIdx = gex->indCB[icb];
            if (useOutputNamespace && wlIdx < gex->pSolo.cbWLstrOut.size()) {
                out.push_back(gex->pSolo.cbWLstrOut[wlIdx]);
            } else if (wlIdx < gex->pSolo.cbWLstr.size()) {
                out.push_back(gex->pSolo.cbWLstr[wlIdx]);
            }
        }
    }
    return !out.empty();
}

static bool loadWhitelistBarcodeSet(const string& whitelistPath,
                                    std::unordered_set<string>& outSet,
                                    uint64_t& rowsRead,
                                    uint64_t& invalidRows) {
    outSet.clear();
    rowsRead = 0;
    invalidRows = 0;

    std::ifstream in(whitelistPath.c_str());
    if (!in.is_open()) {
        return false;
    }

    string line;
    while (std::getline(in, line)) {
        string trimmed = trimCopy(line);
        if (trimmed.empty()) {
            continue;
        }
        size_t first = trimmed.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t end = trimmed.find_first_of("\t, \r\n", first);
        string token = (end == string::npos) ? trimmed.substr(first)
                                             : trimmed.substr(first, end - first);
        token = upperCopy(token);
        rowsRead++;
        if (!isValidBarcodeSeq(token)) {
            invalidRows++;
            continue;
        }
        outSet.insert(token);
    }

    return !outSet.empty();
}

static bool loadBarcodeListFromFile(const string& path,
                                    vector<string>& out,
                                    uint64_t& rowsRead,
                                    uint64_t& invalidRows) {
    out.clear();
    rowsRead = 0;
    invalidRows = 0;

    std::ifstream in(path.c_str());
    if (!in.is_open()) {
        return false;
    }

    string line;
    while (std::getline(in, line)) {
        string trimmed = trimCopy(line);
        if (trimmed.empty()) {
            continue;
        }
        size_t first = trimmed.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t end = trimmed.find_first_of("\t, \r\n", first);
        string token = (end == string::npos) ? trimmed.substr(first)
                                             : trimmed.substr(first, end - first);
        token = upperCopy(token);
        rowsRead++;
        if (!isValidBarcodeSeq(token)) {
            invalidRows++;
            continue;
        }
        out.push_back(token);
    }

    return !out.empty();
}

// Normalize filtered barcodes into the whitelist (assignment) namespace.
//
// sourceNamespace and whitelistNamespace must both be resolved ("NXT" or
// "TRU").  Uses deterministic mapping only.
static vector<string> normalizeFilteredBarcodesForAssignNamespace(
    const vector<string>& inputBarcodes,
    const std::unordered_set<string>& whitelistSet,
    FilteredBarcodeNormalizationStats& stats,
    const string& sourceNamespace = "",
    const string& whitelistNamespace = "") {
    stats = FilteredBarcodeNormalizationStats{};
    vector<string> output;
    output.reserve(inputBarcodes.size());
    std::unordered_set<string> seen;
    seen.reserve(inputBarcodes.size());

    if (!isKnownNamespace(sourceNamespace) || !isKnownNamespace(whitelistNamespace)) {
        throw runtime_error(
            "normalizeFilteredBarcodesForAssignNamespace: both sourceNamespace ('"
            + sourceNamespace + "') and whitelistNamespace ('" + whitelistNamespace
            + "') must be resolved to NXT or TRU.  Heuristic membership-based "
            "fallback has been removed per strict namespace policy.");
    }
    const bool needsTranslation = (sourceNamespace != whitelistNamespace);

    for (const auto& rawBc : inputBarcodes) {
        const string bc = upperCopy(rawBc);
        string candidate = needsTranslation ? translateNxtMiddleTwoBases(bc) : bc;
        string selected;

        if (whitelistSet.find(candidate) != whitelistSet.end()) {
            selected = candidate;
            if (needsTranslation) {
                stats.translatedToSet++;
            } else {
                stats.inSet++;
            }
        } else {
            stats.unmatched++;
        }

        if (!selected.empty()) {
            if (seen.insert(selected).second) {
                output.push_back(selected);
            } else {
                stats.dedupDropped++;
            }
        }
    }
    stats.outputCount = output.size();
    return output;
}

static void appendAssignNormalizationStats(const string& assignOut,
                                           const PfLibraryNamespaceContext& nsCtx,
                                           const PfMultiAssign::WhitelistNormalizationResult& wlInfo) {
    string runSummaryPath = assignOut + "/assignBarcodes.api_run.txt";
    std::ofstream out(runSummaryPath.c_str(), std::ios::app);
    if (!out.is_open()) {
        return;
    }
    out << "namespace_context.library_id=" << nsCtx.libraryId << "\n";
    out << "namespace_context.effective_read_namespace=" << nsCtx.effectiveReadNamespace << "\n";
    out << "namespace_context.assignment_whitelist_namespace=" << nsCtx.assignmentWhitelistNamespace << "\n";
    out << "namespace_context.translate_nxt_for_assign=" << (nsCtx.translateNxtForAssign ? 1 : 0) << "\n";
    out << "namespace_context.output_namespace=" << nsCtx.outputNamespace << "\n";
    out << "namespace_context.detected_match_mode=" << nsCtx.autoDetectMatchMode << "\n";
    out << "namespace_context.chemistry_explicit=" << (nsCtx.isChemistryExplicit ? 1 : 0) << "\n";
    out << "namespace_context.namespace_confident=" << (nsCtx.isNamespaceConfident ? 1 : 0) << "\n";
    out << "namespace_context.whitelist_source=" << wlInfo.sourcePath << "\n";
    out << "namespace_context.whitelist_normalized=" << wlInfo.normalizedPath << "\n";
    out << "namespace_context.whitelist_rows=" << wlInfo.normalizedRowCount << "\n";
    out << "barcode_normalization.in_set=" << nsCtx.normalizationStats.inSet << "\n";
    out << "barcode_normalization.translated_to_set=" << nsCtx.normalizationStats.translatedToSet << "\n";
    out << "barcode_normalization.unmatched=" << nsCtx.normalizationStats.unmatched << "\n";
    out << "barcode_normalization.dedup_dropped=" << nsCtx.normalizationStats.dedupDropped << "\n";
    out << "barcode_normalization.output_count=" << nsCtx.normalizationStats.outputCount << "\n";
    out << "namespace_policy.strict_exact_only=1\n";
    out << "namespace_policy.allow_union_whitelist="
        << (nsCtx.allowUnionWhitelist ? 1 : 0) << "\n";
}

// knownFeatureRefTypeMap() and buildFeatureSpecsFromConfig() are in PfMultiFeatureSpecs.h

static PfMultiPreparedContext buildPfMultiPreparedContext(const PfMultiPreloadInput& input) {
    PfMultiPreparedContext context;
    std::ostringstream prepLog;

    context.config = PfMultiConfig::parseConfig(input.pfMultiConfig);
    if (context.config.libraries.empty()) {
        throw runtime_error("No libraries found in multi config");
    }

    context.featureRef = input.crFeatureRef;
    if (isUnsetToken(context.featureRef)) {
        context.featureRef = context.config.featureRef;
    }
    // Global featureRef is only required if at least one non-GEX library
    // lacks a per-library star_feature_ref.
    if (isUnsetToken(context.featureRef)) {
        bool anyNeedsGlobal = false;
        for (const auto& lib : context.config.libraries) {
            string norm = lib.normalizedFeatureType();
            bool isGex = (norm == "geneexpression" || norm == "gex");
            if (!isGex && lib.starFeatureRef.empty()) {
                anyNeedsGlobal = true;
                break;
            }
        }
        if (anyNeedsGlobal) {
            throw runtime_error("Feature reference not provided: use --crFeatureRef, "
                "set [feature] ref in config, or provide star_feature_ref for every feature library");
        }
    }

    context.whitelist = input.crWhitelist;
    if (isUnsetToken(context.whitelist)) {
        if (!input.soloCbWhitelist.empty() && !isUnsetToken(input.soloCbWhitelist[0])) {
            context.whitelist = input.soloCbWhitelist[0];
        }
    }
    if (isUnsetToken(context.whitelist)) {
        throw runtime_error("Whitelist not provided (use --crWhitelist or --soloCBwhitelist)");
    }

    context.soloWhitelist = context.whitelist;
    if (!input.soloCbWhitelist.empty() && !isUnsetToken(input.soloCbWhitelist[0])) {
        context.soloWhitelist = input.soloCbWhitelist[0];
    }

    context.requestedChem = parseChemistryToken(input.crChemistry, "--crChemistry");
    context.requestedOutputChem = parseChemistryToken(input.crOutputChemistry, "--crOutputChemistry");
    context.inferredChem = detectChemistryFromWhitelistPath(context.whitelist, context.inferredReason,
                                                            context.inferredChemConfident);
    context.effectiveChem = context.inferredChem;
    if (context.requestedChem == "nxt") {
        context.effectiveChem = "NXT";
    } else if (context.requestedChem == "tru") {
        context.effectiveChem = "TRU";
    }
    if (context.requestedChem == "auto" && !isKnownNamespace(context.effectiveChem)) {
        prepLog
            << "NOTICE: whitelist chemistry namespace unresolved for "
            << context.whitelist
            << ". Continuing with per-library autodetect/identity handling; "
            << "mixed NXT/TRU translation remains disabled until a namespace is "
            << "explicitly known.\n";
    }
    if (!context.inferredChemConfident && isKnownNamespace(context.effectiveChem)) {
        prepLog
            << "NOTICE: whitelist chemistry '" << context.effectiveChem
            << "' set via explicit --crChemistry (" << context.inferredReason
            << "). Content-based auto-detection could not determine orientation.\n";
    }

    context.soloInferredChem = detectChemistryFromWhitelistPath(
        context.soloWhitelist, context.soloInferredReason, context.soloInferredChemConfident);
    context.soloEffectiveChem = context.soloInferredChem;

    // Solo writes barcodes.tsv using cbWLstrOut (COL2 for 2-column whitelists).
    // soloEffectiveChem represents COL1 (matching namespace), so the output
    // namespace is the opposite for 2-column files, regardless of NXT-first vs
    // TRU-first.
    context.hasTwoColumnWhitelist = whitelistHasTwoColumns(context.whitelist);
    context.soloHasTwoColumnWhitelist = whitelistHasTwoColumns(context.soloWhitelist);
    if (context.soloHasTwoColumnWhitelist && isKnownNamespace(context.soloEffectiveChem)) {
        context.soloOutputNamespace = oppositeNamespace(context.soloEffectiveChem);
    } else {
        context.soloOutputNamespace = context.soloEffectiveChem;
    }

    // CR-compat default output namespace is TRU.
    context.outputChem = "TRU";
    if (context.requestedOutputChem == "nxt") {
        context.outputChem = "NXT";
    } else if (context.requestedOutputChem == "tru") {
        context.outputChem = "TRU";
    }

    prepLog << "pf-multi chemistry: requested=" << context.requestedChem
            << " inferred=" << context.inferredChem
            << " (confident=" << (context.inferredChemConfident ? "yes" : "no")
            << "; " << context.inferredReason << ")"
            << " effective=" << context.effectiveChem
            << " hasTwoColumnWL=" << (context.hasTwoColumnWhitelist ? "yes" : "no")
            << " soloWhitelist=" << context.soloWhitelist
            << " soloInferred=" << context.soloInferredChem
            << " (confident=" << (context.soloInferredChemConfident ? "yes" : "no")
            << "; " << context.soloInferredReason << ")"
            << " soloEffective=" << context.soloEffectiveChem
            << " soloHasTwoColumnWL=" << (context.soloHasTwoColumnWhitelist ? "yes" : "no")
            << " soloOutputNamespace=" << context.soloOutputNamespace
            << " output_requested=" << context.requestedOutputChem
            << " output_effective=" << context.outputChem
            << "\n";

    // Per-library GEX chemistry: star_chemistry column on GEX row overrides the
    // Solo/GEX namespace interpretation, but not the feature whitelist namespace.
    {
        vector<PfMultiConfig::LibraryEntry> gexLibs = context.config.getGexLibraries();
        for (const auto& gex : gexLibs) {
            if (!gex.starChemistry.empty() && gex.starChemistry != "auto") {
                string prevSoloEffective = context.soloEffectiveChem;
                if (gex.starChemistry == "nxt") {
                    context.soloEffectiveChem = "NXT";
                } else if (gex.starChemistry == "tru") {
                    context.soloEffectiveChem = "TRU";
                }
                // Recompute soloOutputNamespace after override
                if (context.soloHasTwoColumnWhitelist && isKnownNamespace(context.soloEffectiveChem)) {
                    context.soloOutputNamespace = oppositeNamespace(context.soloEffectiveChem);
                } else {
                    context.soloOutputNamespace = context.soloEffectiveChem;
                }
                prepLog << "  GEX star_chemistry=" << gex.starChemistry
                        << " overrides Solo namespace: soloEffectiveChem " << prevSoloEffective
                        << " → " << context.soloEffectiveChem
                        << ", soloOutputNamespace=" << context.soloOutputNamespace
                        << " (user-specified)\n";
                break;
            }
        }
    }

    const map<string, string> fastqMap = PfMultiConfig::parseFastqMap(input.crFastqMap);

    context.outPrefix = input.outFileNamePrefix;
    while (!context.outPrefix.empty() && context.outPrefix.back() == '/') {
        context.outPrefix.pop_back();
    }
    context.crAssignRoot = context.outPrefix + "/cr_assign";

    const string mkCrAssignRootCmd = "mkdir -p \"" + context.crAssignRoot + "\"";
    if (system(mkCrAssignRootCmd.c_str()) != 0) {
        throw runtime_error("Failed to create CR assign root directory: " + context.crAssignRoot);
    }

    context.mapEstimate =
        estimateReadsAnchored(selectMapPrimaryFastqs(input), "map", prepLog);

    const vector<FeatureSpec> featureSpecs = buildFeatureSpecsFromConfig(context.config, prepLog);
    for (const auto& spec : featureSpecs) {
        const vector<PfMultiConfig::LibraryEntry> libs =
            context.config.getFeatureLibraries(spec.libraryType);
        if (libs.empty()) {
            continue;
        }

        const string featureDir = sanitizeDirName(spec.libraryType);
        const string assignBase = context.outPrefix + "/cr_assign/" + featureDir;

        for (const auto& lib : libs) {
            PfPreparedFeatureLibrary prepared;
            prepared.libraryType = spec.libraryType;
            prepared.featureRefType = spec.featureRefType;
            prepared.adtMexOutput =
                PfMultiFeatureSpecs::shouldEmitAdtMexOutput(
                    spec.featureRefType, spec.libraryType);
            prepared.resolvedFastq =
                PfMultiConfig::resolveFastqDir(lib.fastqs, input.crFastqRoot, fastqMap);
            prepared.sampleName = lib.sample.empty()
                ? basenameOf(prepared.resolvedFastq)
                : lib.sample;
            prepared.libraryId = lib.starLibraryId;
            prepared.starMaxHamming = lib.starMaxHamming;
            prepared.useSplitReadLayout = lib.hasSplitReadLayout();
            if (prepared.useSplitReadLayout) {
                prepared.splitReadLayout = PfMultiConfig::buildSplitReadLayout(lib);
                prepared.starBarcodeOutputMap = lib.starBarcodeOutputMap;
                prepared.starFeatureSearchMode = lib.starFeatureSearchMode;
            }
            prepared.assignOut =
                assignBase + "/" + sanitizeDirName(prepared.libraryId);

            const string mkAssignOutCmd = "mkdir -p \"" + prepared.assignOut + "\"";
            if (system(mkAssignOutCmd.c_str()) != 0) {
                throw runtime_error("Failed to create assign output directory: " + prepared.assignOut);
            }

            // Feature ref precedence: star_feature_ref > global [feature] ref
            if (!lib.starFeatureRef.empty()) {
                prepared.featureRefPath = lib.starFeatureRef;
                prepared.usedFilteredRef = false;
                prepLog << "  using per-library star_feature_ref=" << lib.starFeatureRef
                        << " (filter skipped)\n";
            } else {
                const string filteredRef = prepared.assignOut + "/feature_reference.filtered.csv";
                prepared.usedFilteredRef =
                    filterFeatureRefCsv(context.featureRef, spec.featureRefType, filteredRef);
                prepared.featureRefPath = prepared.usedFilteredRef ? filteredRef : context.featureRef;
            }
            prepared.whitelistPath = lib.starWhitelist.empty() ? context.whitelist : lib.starWhitelist;

            // Per-library chemistry: star_chemistry column > --crChemistry flag
            string libChemRequest = lib.starChemistry;
            if (libChemRequest.empty()) {
                libChemRequest = context.requestedChem;
            }
            prepared.resolvedChemRequest = libChemRequest;
            if (libChemRequest == "nxt") {
                prepared.effectiveChem = "NXT";
                prepared.explicitChem = true;
            } else if (libChemRequest == "tru") {
                prepared.effectiveChem = "TRU";
                prepared.explicitChem = true;
            } else {
                prepared.effectiveChem = context.effectiveChem;
                prepared.explicitChem = false;
            }
            prepLog << "  library " << prepared.libraryId
                    << " (" << prepared.sampleName << "/" << spec.libraryType << ")"
                    << ": star_chemistry=" << (lib.starChemistry.empty() ? "(empty)" : lib.starChemistry)
                    << ", global=" << context.requestedChem
                    << ", resolved=" << libChemRequest
                    << ", effectiveChem=" << prepared.effectiveChem
                    << (prepared.explicitChem ? " (explicit)" : " (auto-detect eligible)")
                    << ", star_max_hamming=" << (prepared.starMaxHamming >= 0 ? std::to_string(prepared.starMaxHamming) : "(global)")
                    << ", whitelist=" << prepared.whitelistPath
                    << ", featureRef=" << prepared.featureRefPath
                    << ", adt_mex_output=" << (prepared.adtMexOutput ? "yes" : "no")
                    << "\n";

            if (!prepared.usedFilteredRef && lib.starFeatureRef.empty()) {
                prepLog << "WARNING: feature reference not filtered for "
                        << spec.libraryType << "; using full reference\n";
            } else if (!prepared.usedFilteredRef && !lib.starFeatureRef.empty()) {
                prepLog << "NOTICE: per-library star_feature_ref provided for "
                        << spec.libraryType << " (library_id="
                        << prepared.libraryId << "); type-based filtering skipped as expected\n";
            }

            vector<string> allFastqs = listFastqFilesInDirectory(prepared.resolvedFastq);
            vector<string> featurePrimaryFastqs = filterFastqFilesByMarker(allFastqs, "_R1_");
            if (featurePrimaryFastqs.empty()) {
                featurePrimaryFastqs.swap(allFastqs);
            }
            prepared.featureEstimate =
                estimateReadsAnchored(featurePrimaryFastqs, "feature", prepLog);

            context.featureLibraries.push_back(std::move(prepared));
        }
    }

    if (context.featureLibraries.empty()) {
        throw runtime_error("No feature libraries found in multi config");
    }

    context.prepLog = prepLog.str();
    return context;
}

/**
 * Run CRISPR feature calling on filtered MEX with CRISPR Guide Capture features.
 * 
 * @param filteredMexDir Directory containing filtered_feature_bc_matrix
 * @param outputDir Output directory for crispr_analysis/
 * @param minUmi Minimum UMI threshold for GMM calling
 * @param logStream Log output stream
 * @return 0 on success, -1 on failure
 */
static int runCrisprFeatureCalling(const string& filteredMexDir, const string& outputDir,
                                    int minUmi, ostream& logStream) {
    logStream << timeMonthDayTime() << " ..... starting CRISPR feature calling\n";
    
    // Step 1: Read the filtered MEX
    PfMultiMerge::MexData mexData;
    try {
        mexData = PfMultiMerge::readMex(filteredMexDir);
    } catch (const exception& e) {
        logStream << "ERROR: Failed to read filtered MEX for CRISPR calling: " << e.what() << "\n";
        return -1;
    }
    
    // Step 2: Filter to CRISPR Guide Capture features only
    PfMultiMerge::MexData crisprData = PfMultiMerge::filterByFeatureType(mexData, "CRISPR Guide Capture");
    
    if (crisprData.features.empty()) {
        logStream << "NOTICE: No CRISPR Guide Capture features found, skipping feature calling\n";
        return 0;
    }
    
    logStream << "  CRISPR features found: " << crisprData.features.size() << "\n";
    logStream << "  Barcodes: " << crisprData.barcodes.size() << "\n";
    logStream << "  Non-zero entries: " << crisprData.triplets.size() << "\n";
    
    // Step 3: Write CRISPR-only MEX to temporary directory
    string tempMexDir = outputDir + "/.crispr_mex_tmp";
    string mkdirCmd = "mkdir -p \"" + tempMexDir + "\"";
    if (system(mkdirCmd.c_str()) != 0) {
        logStream << "ERROR: Failed to create temp directory: " << tempMexDir << "\n";
        return -1;
    }
    
    // Convert MexData to MexWriter format
    vector<MexWriter::Feature> features;
    for (size_t i = 0; i < crisprData.features.size(); ++i) {
        string name = (i < crisprData.featureNames.size()) ? crisprData.featureNames[i] : crisprData.features[i];
        string type = (i < crisprData.featureTypes.size()) ? crisprData.featureTypes[i] : "CRISPR Guide Capture";
        features.emplace_back(crisprData.features[i], name, type);
    }
    
    // Write MEX (uncompressed for call_features compatibility)
    string mexPrefix = tempMexDir + "/";
    int ret = MexWriter::writeMex(mexPrefix, crisprData.barcodes, features, crisprData.triplets, -1);
    if (ret != 0) {
        logStream << "ERROR: Failed to write temporary CRISPR MEX\n";
        return -1;
    }
    
    // Step 4: Run GMM feature calling with min_umi=10 (CR-compatible default)
    string crisprAnalysisDir = outputDir + "/crispr_analysis";
    mkdirCmd = "mkdir -p \"" + crisprAnalysisDir + "\"";
    if (system(mkdirCmd.c_str()) != 0) {
        logStream << "ERROR: Failed to create crispr_analysis directory\n";
        string rmCmd = "rm -rf \"" + tempMexDir + "\"";
        system(rmCmd.c_str());
        return -1;
    }
    
    cf_gmm_config *gmm_cfg = cf_gmm_config_create();
    if (!gmm_cfg) {
        logStream << "ERROR: Failed to create GMM config\n";
        return -1;
    }
    gmm_cfg->min_umi_threshold = minUmi;
    
    logStream << "  Calling mode: GMM (CR9-compatible)\n";
    logStream << "  min_umi: " << gmm_cfg->min_umi_threshold << "\n";
    logStream << "  Output: " << crisprAnalysisDir << "\n";
    
    ret = cf_process_mex_dir_gmm(tempMexDir.c_str(), crisprAnalysisDir.c_str(), gmm_cfg);
    cf_gmm_config_destroy(gmm_cfg);
    
    if (ret != 0) {
        logStream << "ERROR: CRISPR feature calling failed\n";
        // Cleanup temp dir
        string rmCmd = "rm -rf \"" + tempMexDir + "\"";
        system(rmCmd.c_str());
        return -1;
    }
    
    // Step 5: Cleanup temporary MEX directory
    string rmCmd = "rm -rf \"" + tempMexDir + "\"";
    system(rmCmd.c_str());
    
    logStream << timeMonthDayTime() << " ..... finished CRISPR feature calling\n";
    logStream << "  Output files:\n";
    logStream << "    " << crisprAnalysisDir << "/protospacer_calls_per_cell.csv\n";
    logStream << "    " << crisprAnalysisDir << "/protospacer_calls_summary.csv\n";
    logStream << "    " << crisprAnalysisDir << "/protospacer_umi_thresholds.csv\n";
    logStream << "    " << crisprAnalysisDir << "/protospacer_umi_thresholds.json\n";
    
    return 0;
}

} // namespace

static PfMultiPreparedContext resolvePfMultiPreparedContext(
    Parameters& P,
    const std::shared_ptr<PfMultiPreloadHandle>& preload) {
    PfMultiPreparedContext prepared;
    bool usedPreload = false;
    if (preload && preload->started) {
        if (preload->worker.joinable()) {
            preload->worker.join();
        }
        if (preload->success) {
            prepared = preload->context;
            usedPreload = true;
            P.inOut->logMain << "pf-preload: consumed async preparation (sec="
                             << preload->preloadSec
                             << ", featureLibraries=" << prepared.featureLibraries.size()
                             << ")\n";
        } else {
            P.inOut->logMain << "WARNING: pf-preload failed; using synchronous preparation. cause="
                             << (preload->error.empty() ? string("unknown") : preload->error)
                             << "\n";
        }
    }
    if (!usedPreload) {
        const auto prepStart = std::chrono::steady_clock::now();
        prepared = buildPfMultiPreparedContext(makePfMultiPreloadInput(P));
        const auto prepEnd = std::chrono::steady_clock::now();
        const double prepSec =
            std::chrono::duration_cast<std::chrono::duration<double>>(prepEnd - prepStart).count();
        P.inOut->logMain << "pf-preload: synchronous preparation complete (sec="
                         << prepSec
                         << ", featureLibraries=" << prepared.featureLibraries.size()
                         << ")\n";
    }
    if (!prepared.prepLog.empty()) {
        P.inOut->logMain << prepared.prepLog;
    }
    return prepared;
}

static string stripBarcodeSuffix(const string& barcode) {
    size_t dashPos = barcode.find_last_of('-');
    if (dashPos != string::npos && dashPos < barcode.size() - 1) {
        bool digits = true;
        for (size_t i = dashPos + 1; i < barcode.size(); ++i) {
            if (!std::isdigit(static_cast<unsigned char>(barcode[i]))) {
                digits = false;
                break;
            }
        }
        if (digits) {
            return barcode.substr(0, dashPos);
        }
    }
    return barcode;
}

static void writeDeferredFilteredAssignOutput(const string& assignOut,
                                              const vector<string>& filteredGexBarcodes,
                                              const string& featureBarcodeNamespace,
                                              const string& whitelistPath,
                                              ofstream& logStream) {
    if (filteredGexBarcodes.empty()) {
        return;
    }

    PfMultiMerge::MexData rawData = PfMultiMerge::readMex(assignOut);
    vector<string> nativeBarcodes;
    {
        const string nativeBarcodePath = assignOut + "/barcodes.txt";
        std::ifstream in(nativeBarcodePath.c_str());
        string line;
        while (std::getline(in, line)) {
            line = trimCopy(line);
            if (!line.empty()) {
                nativeBarcodes.push_back(line);
            }
        }
        if (!nativeBarcodes.empty() && nativeBarcodes.size() != rawData.barcodes.size()) {
            logStream << "WARNING: deferred filtered assign output for " << assignOut
                      << " found mismatched barcodes.txt/barcodes.tsv sizes (native="
                      << nativeBarcodes.size() << ", mex=" << rawData.barcodes.size()
                      << "); falling back to MEX barcode order only\n";
            nativeBarcodes.clear();
        }
        if (nativeBarcodes.empty()) {
            nativeBarcodes = rawData.barcodes;
        }
    }
    std::unordered_set<string> filteredSet;
    filteredSet.reserve(filteredGexBarcodes.size());
    for (const auto& bc : filteredGexBarcodes) {
        filteredSet.insert(stripBarcodeSuffix(upperCopy(bc)));
    }

    std::vector<string> subsetBarcodes;
    subsetBarcodes.reserve(filteredSet.size());
    std::vector<uint32_t> oldToCompact(rawData.barcodes.size(), std::numeric_limits<uint32_t>::max());
    for (size_t i = 0; i < rawData.barcodes.size(); ++i) {
        const string nativeBarcode = stripBarcodeSuffix(upperCopy(nativeBarcodes[i]));
        const string mexBarcode = stripBarcodeSuffix(upperCopy(rawData.barcodes[i]));
        const string normalizedMexBarcode = normalizeBarcodeToTru(mexBarcode, featureBarcodeNamespace);
        if (filteredSet.find(nativeBarcode) != filteredSet.end()
            || filteredSet.find(mexBarcode) != filteredSet.end()
            || filteredSet.find(normalizedMexBarcode) != filteredSet.end()) {
            oldToCompact[i] = static_cast<uint32_t>(subsetBarcodes.size());
            subsetBarcodes.push_back(nativeBarcodes[i]);
        }
    }

    const string filteredDir = assignOut + "/filtered";
    const string resetCmd = "rm -rf \"" + filteredDir + "\" && mkdir -p \"" + filteredDir + "\"";
    if (system(resetCmd.c_str()) != 0) {
        throw runtime_error("Failed to reset deferred filtered assign output dir: " + filteredDir);
    }
    if (subsetBarcodes.empty()) {
        logStream << "WARNING: deferred filtered assign output for " << assignOut
                  << " matched zero GEX barcodes\n";
        return;
    }

    std::vector<MexWriter::Triplet> subsetTriplets;
    subsetTriplets.reserve(rawData.triplets.size());
    uint64_t dedupedCount = 0;
    for (const auto& t : rawData.triplets) {
        if (t.cell_idx < oldToCompact.size() && oldToCompact[t.cell_idx] != std::numeric_limits<uint32_t>::max()) {
            MexWriter::Triplet newT;
            newT.cell_idx = oldToCompact[t.cell_idx];
            newT.gene_idx = t.gene_idx;
            newT.count = t.count;
            subsetTriplets.push_back(newT);
            dedupedCount += t.count;
        }
    }
    std::sort(subsetTriplets.begin(), subsetTriplets.end(),
              [](const MexWriter::Triplet& a, const MexWriter::Triplet& b) {
                  if (a.cell_idx != b.cell_idx) {
                      return a.cell_idx < b.cell_idx;
                  }
                  return a.gene_idx < b.gene_idx;
              });

    std::vector<MexWriter::Feature> features;
    features.reserve(rawData.features.size());
    for (size_t i = 0; i < rawData.features.size(); ++i) {
        const string name = (i < rawData.featureNames.size()) ? rawData.featureNames[i] : rawData.features[i];
        const string type = (i < rawData.featureTypes.size()) ? rawData.featureTypes[i] : "Gene Expression";
        features.emplace_back(rawData.features[i], name, type);
    }
    if (MexWriter::writeMex(filteredDir + "/", subsetBarcodes, features, subsetTriplets, -1) != 0) {
        throw runtime_error("Failed to write deferred filtered MEX for assign output: " + assignOut);
    }

    {
        std::ofstream out((filteredDir + "/barcodes.txt").c_str());
        if (!out.is_open()) {
            throw runtime_error("Failed to write filtered barcodes.txt: " + filteredDir);
        }
        for (const auto& bc : subsetBarcodes) {
            out << bc << "\n";
        }
    }
    if (!PfMultiMexStub::copyBarcodesTsv(filteredDir + "/barcodes.txt",
                                         filteredDir + "/barcodes.tsv",
                                         true,
                                         whitelistPath)) {
        throw runtime_error("Failed to regenerate filtered barcodes.tsv: " + filteredDir);
    }
    {
        std::ofstream out((filteredDir + "/features.txt").c_str());
        if (!out.is_open()) {
            throw runtime_error("Failed to write filtered features.txt: " + filteredDir);
        }
        for (const auto& featureName : rawData.featureNames) {
            out << featureName << "\n";
        }
    }
    {
        const string rawFeaturePerCell = assignOut + "/feature_per_cell.csv";
        std::ifstream in(rawFeaturePerCell.c_str());
        if (in.is_open()) {
            std::ofstream out((filteredDir + "/feature_per_cell.csv").c_str());
            if (!out.is_open()) {
                throw runtime_error("Failed to write filtered feature_per_cell.csv: " + filteredDir);
            }
            string line;
            if (std::getline(in, line)) {
                out << line << "\n";
            }
            std::unordered_set<string> subsetNativeSet;
            subsetNativeSet.reserve(subsetBarcodes.size());
            for (const auto& bc : subsetBarcodes) {
                subsetNativeSet.insert(stripBarcodeSuffix(upperCopy(bc)));
            }
            while (std::getline(in, line)) {
                const size_t commaPos = line.find(',');
                const string barcode = upperCopy(stripBarcodeSuffix(
                    commaPos == string::npos ? line : line.substr(0, commaPos)));
                if (subsetNativeSet.find(barcode) != subsetNativeSet.end()) {
                    out << line << "\n";
                }
            }
        }
    }
    {
        const string rawStatsPath = assignOut + "/stats.txt";
        std::ifstream in(rawStatsPath.c_str());
        std::ofstream out((filteredDir + "/stats.txt").c_str());
        if (!out.is_open()) {
            throw runtime_error("Failed to write filtered stats.txt: " + filteredDir);
        }
        if (in.is_open()) {
            string line;
            while (std::getline(in, line)) {
                if (line.find("Total deduped feature counts ") == 0) {
                    out << "Total deduped feature counts " << dedupedCount << "\n";
                } else if (line.find("Total barcodes ") == 0) {
                    out << "Total barcodes " << subsetBarcodes.size() << "\n";
                } else if (line.find("Total excluded barcodes ") == 0) {
                    const uint64_t excluded =
                        rawData.barcodes.size() >= subsetBarcodes.size()
                        ? static_cast<uint64_t>(rawData.barcodes.size() - subsetBarcodes.size())
                        : 0;
                    out << "Total excluded barcodes " << excluded << "\n";
                } else {
                    out << line << "\n";
                }
            }
        } else {
            out << "Total deduped feature counts " << dedupedCount << "\n";
            out << "Total barcodes " << subsetBarcodes.size() << "\n";
            out << "Total excluded barcodes "
                << (rawData.barcodes.size() >= subsetBarcodes.size()
                        ? static_cast<uint64_t>(rawData.barcodes.size() - subsetBarcodes.size())
                        : 0)
                << "\n";
        }
    }

    logStream << "NOTICE: wrote deferred filtered assign output for " << assignOut
              << " (barcodes=" << subsetBarcodes.size()
              << ", deduped_counts=" << dedupedCount << ")\n";
}

std::shared_ptr<PfMultiPreloadHandle> startPfMultiConfigPreload(const Parameters& P) {
    if (P.runMode != "alignReads") {
        return nullptr;
    }
    if (isUnsetToken(P.pfMulti.pfMultiConfig)) {
        return nullptr;
    }
    if (P.dynamicThreadInterface != 1) {
        return nullptr;
    }

    auto preload = std::make_shared<PfMultiPreloadHandle>();
    preload->started = true;
    const PfMultiPreloadInput input = makePfMultiPreloadInput(P);
    preload->worker = std::thread([preload, input]() {
        const auto t0 = std::chrono::steady_clock::now();
        try {
            preload->context = buildPfMultiPreparedContext(input);
            preload->success = true;
        } catch (const std::exception& e) {
            preload->error = e.what();
            preload->success = false;
        } catch (...) {
            preload->error = "unknown exception";
            preload->success = false;
        }
        const auto t1 = std::chrono::steady_clock::now();
        preload->preloadSec =
            std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
        preload->finished = true;
    });
    return preload;
}

std::shared_ptr<PfMultiAssignPhaseResult> runPfMultiAssignPhase(
    Parameters& P,
    const std::shared_ptr<PfMultiPreloadHandle>& preload) {
    if (isUnsetToken(P.pfMulti.pfMultiConfig)) {
        return nullptr;
    }
    
    P.inOut->logMain << timeMonthDayTime() << " ..... started pf-multi feature assignment\n";
    
    auto phase = std::make_shared<PfMultiAssignPhaseResult>();
    try {
        phase->prepared = resolvePfMultiPreparedContext(P, preload);
        const PfMultiPreparedContext& prepared = phase->prepared;

        const string& effectiveChem = prepared.effectiveChem;
        const string& outputChem = prepared.outputChem;
        const AnchoredReadEstimate mapEstimate = prepared.mapEstimate;

        PfMultiAssign::AssignOptions assignOpts;
        assignOpts.maxHammingDistance = P.pfMulti.crAssignMaxHamming;
        assignOpts.featureConstantOffset = P.pfMulti.crAssignFeatureOffset;
        assignOpts.limitSearch = P.pfMulti.crAssignLimitSearch;
        assignOpts.minCounts = P.pfMulti.crAssignMinCounts;
        assignOpts.maxBarcodeMismatches = P.pfMulti.crAssignMaxBarcodeMismatches;
        assignOpts.featureN = P.pfMulti.crAssignFeatureN;
        assignOpts.barcodeN = P.pfMulti.crAssignBarcodeN;
        assignOpts.consumerThreadsPerSet = P.pfMulti.crAssignConsumerThreads;
        assignOpts.searchThreads = P.pfMulti.crAssignSearchThreads;
        assignOpts.readBufferLines = P.pfMulti.crAssignReadBufferLines;
        assignOpts.cbqMode = lowerCopy(P.pfMulti.crAssignCbqMode);
        if (assignOpts.cbqMode != "auto" &&
            assignOpts.cbqMode != "stream" &&
            assignOpts.cbqMode != "range") {
            throw runtime_error("Invalid crAssignCbqMode '" + P.pfMulti.crAssignCbqMode +
                                "'; expected auto, stream, or range");
        }
        assignOpts.minPosterior = P.pfMulti.crAssignMinPosterior;
        assignOpts.maxReads = (P.readMapNumber > 0) ? static_cast<long long>(P.readMapNumber) : -1;
        assignOpts.legacyCbRescue = (P.pfMulti.crAssignLegacyCbRescue != 0);
        assignOpts.skipQcOutputs = (P.pfMulti.crAssignSkipQcOutputs != 0);
        assignOpts.allowUnionWhitelist = (P.pfMulti.crAssignAllowUnionWhitelist != 0);
        assignOpts.useFeatureAnchorSearch = true;
        assignOpts.requireFeatureAnchorMatch = true;
        assignOpts.featureModeBootstrapReads =
            (assignOpts.featureConstantOffset >= 0) ? 0 : 100000;
        assignOpts.skipHeatmaps = true;
        if (const char* env = std::getenv("STAR_PF_USE_FEATURE_ANCHOR_SEARCH")) {
            assignOpts.useFeatureAnchorSearch = (std::atoi(env) != 0);
        }
        if (const char* env = std::getenv("STAR_PF_REQUIRE_FEATURE_ANCHOR_MATCH")) {
            assignOpts.requireFeatureAnchorMatch = (std::atoi(env) != 0);
        }
        if (const char* env = std::getenv("STAR_PF_FEATURE_BOOTSTRAP_READS")) {
            assignOpts.featureModeBootstrapReads = std::atoi(env);
        }
        if (const char* env = std::getenv("STAR_PF_USE_HOT_HASH")) {
            assignOpts.useHotHash = (std::atoi(env) != 0);
        }
        if (const char* env = std::getenv("STAR_PF_SKIP_HEATMAPS")) {
            assignOpts.skipHeatmaps = (std::atoi(env) != 0);
        }
        assignOpts.enableStarDynamicPermitHooks = (P.dynamicThreadInterface == 1);
        const PfPermitControllerMode pfControllerMode = parsePfPermitControllerMode(P.dynamicThreadPfControllerMode);
        const bool pfControllerEnabled = (pfControllerMode != PfPermitControllerMode::Off);
        const bool pfControllerSequenceMode =
            (pfControllerMode == PfPermitControllerMode::Shadow ||
             pfControllerMode == PfPermitControllerMode::Active);
        const bool pfControllerEtaMode = (pfControllerMode == PfPermitControllerMode::Eta);
        const bool pfControllerChunkedMode = (pfControllerMode == PfPermitControllerMode::Chunked);
        const bool pfControllerAppliesUpdates =
            (pfControllerMode == PfPermitControllerMode::Active ||
             pfControllerMode == PfPermitControllerMode::Eta ||
             pfControllerMode == PfPermitControllerMode::Chunked);
        const int pfControllerIntervalMs = P.dynamicThreadPfControllerIntervalMs;
        const vector<int> pfControllerSequence = P.dynamicThreadPfControllerSequence;
        const int pfControllerMaxUpdates = P.dynamicThreadPfControllerMaxUpdates;
        const string pfControllerSequenceLabel = joinIntVector(pfControllerSequence);
        const int pfSearchThreadsForController =
            (assignOpts.searchThreads > 0) ? assignOpts.searchThreads : 4;
        const int pfSearchThreadQuantum = normalizePfSearchThreadQuantum(pfSearchThreadsForController);
        const int pfChunkPermitStep = pfSearchThreadQuantum + 1;
        const bool pfConsumerThreadsExplicit = (assignOpts.consumerThreadsPerSet > 0);

        const bool hasExplicitAssignFilteredBarcodes =
            (!P.pfMulti.crAssignFilteredBarcodes.empty() && P.pfMulti.crAssignFilteredBarcodes != "-");
        string explicitAssignFilteredBarcodesPath;
        phase->usedExplicitAssignFilteredBarcodes = hasExplicitAssignFilteredBarcodes;
        if (hasExplicitAssignFilteredBarcodes) {
            struct stat stFiltered;
            if (stat(P.pfMulti.crAssignFilteredBarcodes.c_str(), &stFiltered) != 0) {
                throw runtime_error("crAssignFilteredBarcodes file not found: " + P.pfMulti.crAssignFilteredBarcodes);
            }
            explicitAssignFilteredBarcodesPath = P.pfMulti.crAssignFilteredBarcodes;
            P.inOut->logMain << "NOTICE: Using explicit filtered barcodes for assignBarcodes: "
                             << explicitAssignFilteredBarcodesPath << "\n";
        } else {
            P.inOut->logMain << "NOTICE: deferring assign filtered-barcode gating until post-Solo merge\n";
        }
        vector<PfMultiFeatureRun>& featureRuns = phase->featureRuns;

        // Phase 4: Library-aware PF scheduler.
        // Compute per-library thread budgets proportional to remaining work,
        // guaranteeing at least 1 producer per library. Libraries still run
        // sequentially (pf_api is not concurrent-safe), but each library's
        // allocation reflects its share of the total work.
        const size_t numFeatureLibs = prepared.featureLibraries.size();
        struct LibrarySchedule {
            uint64_t estimatedWork = 0;
            int fileCount = 0;
            int threadBudget = 0;
        };
        vector<LibrarySchedule> librarySchedules(numFeatureLibs);
        uint64_t totalEstimatedWork = 0;
        int totalFileCount = 0;
        for (size_t li = 0; li < numFeatureLibs; ++li) {
            const auto& est = prepared.featureLibraries[li].featureEstimate;
            uint64_t work = est.valid ? est.estimatedReads : 0;
            if (work == 0) work = 1;
            librarySchedules[li].estimatedWork = work;
            librarySchedules[li].fileCount = std::max(1, static_cast<int>(est.fileCount));
            totalEstimatedWork += work;
            totalFileCount += librarySchedules[li].fileCount;
        }

        // Libraries run sequentially (pf_api holds a mutex), so each gets
        // the full thread budget.  The permit controller manages contention
        // with concurrent MAP work dynamically.
        for (size_t li = 0; li < numFeatureLibs; ++li) {
            librarySchedules[li].threadBudget = P.runThreadN;
        }
        if (numFeatureLibs > 1) {
            P.inOut->logMain << "\npf-multi library scheduler (Phase 4):\n";
            P.inOut->logMain << "  total_threads=" << P.runThreadN
                             << ", libraries=" << numFeatureLibs
                             << ", total_est_reads=" << totalEstimatedWork
                             << ", policy=full_budget_sequential\n";
            for (size_t li = 0; li < numFeatureLibs; ++li) {
                P.inOut->logMain << "  [" << li << "] library_id="
                                 << prepared.featureLibraries[li].libraryId
                                 << ", est_reads=" << librarySchedules[li].estimatedWork
                                 << ", files=" << librarySchedules[li].fileCount
                                 << ", thread_budget=" << librarySchedules[li].threadBudget
                                 << "\n";
            }
            P.inOut->logMain << "\n";
        }

        for (size_t libIdx = 0; libIdx < numFeatureLibs; ++libIdx) {
            const auto& preparedLib = prepared.featureLibraries[libIdx];
            const LibrarySchedule& sched = librarySchedules[libIdx];
            const string& resolvedFastq = preparedLib.resolvedFastq;
            const string& sampleName = preparedLib.sampleName;
            const string& assignOut = preparedLib.assignOut;
            const string& refPath = preparedLib.featureRefPath;
            const string& whitelist = preparedLib.whitelistPath;
            const string& libraryType = preparedLib.libraryType;
            const string& featureRefType = preparedLib.featureRefType;
            const AnchoredReadEstimate featureEstimate = preparedLib.featureEstimate;

            const int libThreadBudget = std::max(1, sched.threadBudget);
            const int pfProducerSlotsAvailable =
                std::max(1, std::min(std::max(1, libThreadBudget - 1),
                                     std::max(1, sched.fileCount)));
            const bool pfSingleProducerOnly = (pfProducerSlotsAvailable <= 1);
            int pfReaderThreadsReserved = 1;
            int pfConsumerBudgetThreads = std::max(1, libThreadBudget - pfReaderThreadsReserved);
            int pfConsumerThreadsForRun = std::max(1, assignOpts.consumerThreadsPerSet);
            int pfConsumerSoftMax = std::max(1, libThreadBudget);
            bool pfConsumerThreadsAuto = false;
            if (!pfConsumerThreadsExplicit && assignOpts.enableStarDynamicPermitHooks) {
                pfConsumerThreadsForRun = pfConsumerBudgetThreads;
                pfConsumerThreadsAuto = true;
            } else if (pfConsumerThreadsExplicit) {
                pfConsumerBudgetThreads = std::max(pfConsumerThreadsForRun, pfConsumerBudgetThreads);
            }
            PfMultiAssign::AssignOptions runAssignOpts = assignOpts;
            runAssignOpts.adtMexOutput = preparedLib.adtMexOutput;
            if (preparedLib.adtMexOutput) {
                P.inOut->logMain << "NOTICE: protein/ADT library_id=" << preparedLib.libraryId
                                 << " will emit assignBarcodes ADT MEX output (output_mode=adt_mex)\n";
            }
            if (preparedLib.starMaxHamming >= 0) {
                runAssignOpts.maxHammingDistance = preparedLib.starMaxHamming;
            }
            runAssignOpts.sampleName = sampleName;
            runAssignOpts.consumerThreadsPerSet = pfConsumerThreadsForRun;
            runAssignOpts.filteredBarcodesPath.clear();

            const bool splitReadLayout = preparedLib.useSplitReadLayout;
            if (splitReadLayout) {
                runAssignOpts.useSplitReadLayout = true;
                runAssignOpts.splitReadLayout = preparedLib.splitReadLayout;
                runAssignOpts.autodetectChemistry = false;
                runAssignOpts.translateNxt = false;
                const string searchMode = lowerCopy(
                    preparedLib.starFeatureSearchMode.empty()
                        ? "free"
                        : preparedLib.starFeatureSearchMode);
                if (searchMode == "free") {
                    runAssignOpts.useFeatureAnchorSearch = false;
                    runAssignOpts.requireFeatureAnchorMatch = false;
                    runAssignOpts.featureModeBootstrapReads = 0;
                    runAssignOpts.limitSearch = -1;
                }
                P.inOut->logMain << "NOTICE: split-read guide layout enabled for library_id="
                                 << preparedLib.libraryId
                                 << " feature_search_mode=" << searchMode << "\n";
            }

            // --- Phase 1: Probe original whitelist to determine read namespace ---
            PfMultiAssign::WhitelistNormalizationResult wlInfo =
                PfMultiAssign::normalizeWhitelistForAssign(whitelist, assignOut);
            string originalWhitelistNamespace = upperCopy(wlInfo.assignmentNamespace);
            if (splitReadLayout) {
                originalWhitelistNamespace = "ATAC";
            } else if (!isKnownNamespace(originalWhitelistNamespace) && isKnownNamespace(effectiveChem)) {
                originalWhitelistNamespace = effectiveChem;
            }
            const bool originalNamespaceKnown =
                splitReadLayout || isKnownNamespace(originalWhitelistNamespace);

            const bool useAutodetect = !splitReadLayout
                                      && !preparedLib.explicitChem
                                      && (preparedLib.resolvedChemRequest == "auto");
            string effectiveReadNamespace = upperCopy(preparedLib.effectiveChem);
            string detectedMatchMode = "UNKNOWN";
            if (preparedLib.explicitChem) {
                P.inOut->logMain << "NOTICE: " << featureRefType
                                 << " star_chemistry=" << preparedLib.resolvedChemRequest
                                 << " → effectiveChem=" << preparedLib.effectiveChem
                                 << " (auto-detect skipped)\n";
            }
            if (useAutodetect) {
                PfMultiAssign::AssignOptions detectOpts = runAssignOpts;
                detectOpts.autodetectChemistry = true;
                detectOpts.probeOnly = true;
                detectOpts.enableStarDynamicPermitHooks = false;
                detectOpts.filteredBarcodesPath.clear();
                if (detectOpts.maxReads <= 0 ||
                    detectOpts.maxReads > detectOpts.autodetectChemistryReads) {
                    detectOpts.maxReads = detectOpts.autodetectChemistryReads;
                }
                const string detectOut = assignOut + "/.autodetect_probe";
                PfMultiAssign::AssignResult detectResult =
                    PfMultiAssign::runAssignBarcodes(
                        whitelist, refPath, resolvedFastq, detectOut, detectOpts);
                detectedMatchMode = detectResult.detectedMatchMode;

                if (detectedMatchMode != "RAW_MATCH" && detectedMatchMode != "TRANSLATED_MATCH") {
                    throw runtime_error(
                        "Ambiguous or low-confidence chemistry auto-detect for library_id=" + preparedLib.libraryId +
                        " (" + featureRefType + ", match_mode=" + detectedMatchMode +
                        "). Provide explicit star_chemistry in multi config or set --crChemistry NXT/TRU.");
                }
                if (originalNamespaceKnown) {
                    if (detectedMatchMode == "RAW_MATCH") {
                        effectiveReadNamespace = originalWhitelistNamespace;
                    } else {
                        effectiveReadNamespace = oppositeNamespace(originalWhitelistNamespace);
                    }
                    P.inOut->logMain << "NOTICE: auto-detect " << detectedMatchMode
                                     << " for " << featureRefType
                                     << " → effectiveReadNamespace=" << effectiveReadNamespace
                                     << ", originalWhitelistNamespace=" << originalWhitelistNamespace
                                     << "\n";
                } else {
                    P.inOut->logMain << "WARNING: assignment whitelist namespace unresolved for "
                                     << featureRefType
                                     << " (whitelist=" << wlInfo.sourcePath
                                     << "). Cannot normalize whitelist to read namespace.\n";
                }

                const string cleanupCmd = "rm -rf \"" + detectOut + "\"";
                if (system(cleanupCmd.c_str()) != 0) {
                    P.inOut->logMain << "WARNING: failed to clean auto-detect probe dir: "
                                     << detectOut << "\n";
                }
            } else if (!splitReadLayout && !useAutodetect && !preparedLib.explicitChem) {
                if (originalNamespaceKnown && isKnownNamespace(effectiveReadNamespace)) {
                    // Both known, no probe needed.
                } else {
                    std::ostringstream oss;
                    oss << "Namespace relation unresolved for library_id=" << preparedLib.libraryId
                        << " (effectiveReadNamespace=" << effectiveReadNamespace
                        << ", originalWhitelistNamespace=" << originalWhitelistNamespace
                        << "). Provide explicit star_chemistry and/or a whitelist with resolvable NXT/TRU namespace.";
                    throw runtime_error(oss.str());
                }
            }
            if (splitReadLayout) {
                effectiveReadNamespace = "ATAC";
                detectedMatchMode = "RAW_MATCH";
            }
            runAssignOpts.autodetectChemistry = false;

            // --- Phase 2: Normalize whitelist to match read namespace ---
            // After probe, effectiveReadNamespace is determined. Translate
            // whitelist so assignBarcodes matches in read namespace directly.
            if (!splitReadLayout &&
                isKnownNamespace(effectiveReadNamespace) && originalNamespaceKnown &&
                effectiveReadNamespace != originalWhitelistNamespace) {
                PfMultiAssign::WhitelistNormalizationResult nsWlInfo =
                    PfMultiAssign::normalizeWhitelistToNamespace(
                        whitelist, assignOut, effectiveReadNamespace);
                P.inOut->logMain << "NOTICE: translated whitelist from "
                                 << originalWhitelistNamespace << " to "
                                 << effectiveReadNamespace << " for " << featureRefType
                                 << " (" << nsWlInfo.normalizedRowCount << " barcodes"
                                 << " → " << nsWlInfo.normalizedPath << ")\n";
                wlInfo = nsWlInfo;
            }
            // Assignment whitelist is now in read namespace; no runtime
            // translate_NXT needed in the matcher.
            string assignmentWhitelistNamespace = upperCopy(wlInfo.assignmentNamespace);
            const bool assignmentNamespaceKnown = isKnownNamespace(assignmentWhitelistNamespace);
            runAssignOpts.translateNxt = false;

            PfLibraryNamespaceContext nsCtx;
            nsCtx.libraryId = preparedLib.libraryId;
            nsCtx.assignmentWhitelistNamespace = assignmentWhitelistNamespace;
            nsCtx.outputNamespace = outputChem;
            nsCtx.autoDetectMatchMode = detectedMatchMode;
            nsCtx.isChemistryExplicit = preparedLib.explicitChem;
            nsCtx.isNamespaceConfident = assignmentNamespaceKnown && wlInfo.namespaceConfidence;
            nsCtx.effectiveReadNamespace = effectiveReadNamespace;
            nsCtx.translateNxtForAssign = false;
            nsCtx.allowUnionWhitelist = runAssignOpts.allowUnionWhitelist;

            // --- Phase 3: Filtered barcode normalization ---
            // Filtered barcodes (from Solo/GEX) are in TRU namespace.
            // Assignment now happens in effectiveReadNamespace, so filtered
            // barcodes must be normalized to that namespace.
            const bool needAssignFilteredNormalization = hasExplicitAssignFilteredBarcodes;
            if (needAssignFilteredNormalization) {
                std::unordered_set<string> whitelistSet;
                uint64_t wlRowsRead = 0;
                uint64_t wlInvalidRows = 0;
                if (!loadWhitelistBarcodeSet(wlInfo.normalizedPath, whitelistSet, wlRowsRead, wlInvalidRows)) {
                    throw runtime_error("Failed to load normalized assignment whitelist for filtered-barcode normalization: "
                        + wlInfo.normalizedPath);
                }

                vector<string> sourceFiltered;
                const string sourceLabel = "explicit";
                uint64_t explicitRowsRead = 0;
                uint64_t explicitInvalidRows = 0;
                if (!loadBarcodeListFromFile(explicitAssignFilteredBarcodesPath,
                                             sourceFiltered,
                                             explicitRowsRead,
                                             explicitInvalidRows)) {
                    throw runtime_error("Failed to load explicit assign filtered barcode list: "
                        + explicitAssignFilteredBarcodesPath);
                }
                if (explicitInvalidRows > 0) {
                    P.inOut->logMain
                        << "WARNING: explicit assign filtered barcode file had "
                        << explicitInvalidRows
                        << " invalid rows ignored: "
                        << explicitAssignFilteredBarcodesPath << "\n";
                }

                FilteredBarcodeNormalizationStats normStats;
                if (!assignmentNamespaceKnown) {
                    throw runtime_error(
                        "Cannot normalize filtered barcodes for library_id="
                        + preparedLib.libraryId + " (source=" + sourceLabel
                        + "): assignment whitelist namespace is unresolved.  "
                        "Provide a whitelist with resolvable NXT/TRU namespace.");
                }

                vector<string> normalizedFiltered;
                if (runAssignOpts.allowUnionWhitelist && sourceLabel == "explicit") {
                    normStats = FilteredBarcodeNormalizationStats{};
                    normStats.outputCount = sourceFiltered.size();
                    normalizedFiltered.reserve(sourceFiltered.size());
                    for (const auto& rawBc : sourceFiltered) {
                        normalizedFiltered.push_back(upperCopy(rawBc));
                    }
                } else {
                    // Filtered barcodes from Solo/GEX are in Solo's output
                    // namespace (typically TRU for CR-compat). Normalize from
                    // that source into the assignment (read) namespace.
                    const string filteredSourceNs = upperCopy(prepared.soloOutputNamespace);
                    normalizedFiltered =
                        normalizeFilteredBarcodesForAssignNamespace(
                            sourceFiltered, whitelistSet, normStats,
                            filteredSourceNs, assignmentWhitelistNamespace);
                }
                if (normalizedFiltered.empty() && !sourceFiltered.empty()) {
                    std::ostringstream oss;
                    oss << "Filtered barcodes became empty after assignment-namespace normalization for library_id="
                        << preparedLib.libraryId
                        << " (source=" << sourceLabel
                        << ", source_count=" << sourceFiltered.size()
                        << ", whitelist_rows=" << wlRowsRead
                        << ", whitelist_invalid_rows=" << wlInvalidRows
                        << ", unmatched=" << normStats.unmatched << ")";
                    throw runtime_error(oss.str());
                }

                const string filteredPath = assignOut + "/filtered_gex_barcodes.assign_namespace.txt";
                std::ofstream fbc(filteredPath.c_str());
                if (!fbc.is_open()) {
                    throw runtime_error("Failed to write normalized filtered barcode list for assignBarcodes: " + filteredPath);
                }
                for (const auto& bc : normalizedFiltered) {
                    fbc << bc << "\n";
                }
                fbc.close();
                runAssignOpts.filteredBarcodesPath = filteredPath;
                nsCtx.normalizationStats = normStats;

                P.inOut->logMain
                    << "NOTICE: assign filtered-barcode normalization for library_id=" << preparedLib.libraryId
                    << " source=" << sourceLabel
                    << " in_set=" << normStats.inSet
                    << " translated_to_set=" << normStats.translatedToSet
                    << " unmatched=" << normStats.unmatched
                    << " dedup_dropped=" << normStats.dedupDropped
                    << " output=" << normStats.outputCount
                    << " whitelist_rows=" << wlRowsRead
                    << " whitelist_invalid_rows=" << wlInvalidRows
                    << "\n";
            }

            if (!runAssignOpts.filteredBarcodesPath.empty()) {
                enforceUniqueBarcodeFileOrThrow(
                    runAssignOpts.filteredBarcodesPath,
                    "filtered barcode input (library_id=" + preparedLib.libraryId + ")");
            }

            // Assignment whitelist and filtered barcodes are both in
            // effectiveReadNamespace now; source/target for pf_api ingress
            // normalization should both be that same namespace (no-op).
            runAssignOpts.sourceNamespace = assignmentWhitelistNamespace;
            runAssignOpts.targetNamespace = assignmentWhitelistNamespace;

            const bool pfSingleConsumer = (pfConsumerThreadsForRun <= 1);
            const int pfPermitCeilingDuringPf =
                std::max(1, libThreadBudget - pfReaderThreadsReserved);
            const bool pfControllerCpuAware = (P.dynamicThreadPfControllerCpuAware == 1);
            const int pfControllerCpuSampleMs = P.dynamicThreadPfControllerCpuSampleMs;
            const double pfControllerCpuEmaAlpha = P.dynamicThreadPfControllerCpuEmaAlpha;
            const double pfControllerCpuIdleThreshold = P.dynamicThreadPfControllerCpuIdleThreshold;
            const double pfControllerCpuBusyThreshold = P.dynamicThreadPfControllerCpuBusyThreshold;

            std::atomic<bool> stopPfController(false);
            std::atomic<uint64_t> pfControllerTicks(0);
            std::atomic<uint64_t> pfControllerApplied(0);
            std::atomic<int> pfControllerLastTarget(0);
            std::atomic<uint64_t> pfControllerMapEstimateInitial(0);
            std::atomic<uint64_t> pfControllerFeatureEstimateInitial(0);
            std::atomic<uint64_t> pfControllerMapEstimateFinal(0);
            std::atomic<uint64_t> pfControllerFeatureEstimateFinal(0);
            std::thread pfControllerThread;

            if (pfControllerEnabled) {
                const ThreadControl::MapPermitSnapshot controllerBaseline = g_threadChunks.mapPermitSnapshot();
                uint64_t mapEstimateTotal = std::max<uint64_t>(
                    controllerBaseline.mapDomain.workUnitsTotal, mapEstimate.estimatedReads);
                uint64_t featureEstimateTotal = controllerBaseline.featureDomain.workUnitsTotal;
                if (featureEstimate.valid) {
                    featureEstimateTotal += featureEstimate.estimatedReads;
                }
                mapEstimateTotal = inflateEstimateUntilNonNegativeRemaining(
                    mapEstimateTotal, controllerBaseline.mapDomain.workUnitsTotal);
                featureEstimateTotal = inflateEstimateUntilNonNegativeRemaining(
                    featureEstimateTotal, controllerBaseline.featureDomain.workUnitsTotal);
                pfControllerMapEstimateInitial.store(mapEstimateTotal, std::memory_order_relaxed);
                pfControllerFeatureEstimateInitial.store(featureEstimateTotal, std::memory_order_relaxed);
                pfControllerMapEstimateFinal.store(mapEstimateTotal, std::memory_order_relaxed);
                pfControllerFeatureEstimateFinal.store(featureEstimateTotal, std::memory_order_relaxed);

                P.inOut->logMain << "pf-dynamic-controller: start mode="
                                 << pfPermitControllerModeName(pfControllerMode)
                                 << ", intervalMs=" << pfControllerIntervalMs
                                 << ", sequence=" << pfControllerSequenceLabel
                                 << ", maxUpdates=" << pfControllerMaxUpdates
                                 << ", chunkPermitStep=" << pfChunkPermitStep
                                 << ", searchThreadQuantum=" << pfSearchThreadQuantum
                                 << ", searchThreadsForController=" << pfSearchThreadsForController
                                 << ", consumerThreadsForRun=" << pfConsumerThreadsForRun
                                 << ", consumerSoftMax=" << pfConsumerSoftMax
                                 << ", producerSlotsAvailable=" << pfProducerSlotsAvailable
                                 << ", singleProducerOnly=" << (pfSingleProducerOnly ? "yes" : "no")
                                 << ", singleConsumerMode=" << (pfSingleConsumer ? "yes" : "no")
                                 << ", readerThreadsReserved=" << pfReaderThreadsReserved
                                 << ", permitCeilingDuringPf=" << pfPermitCeilingDuringPf
                                 << ", cpuAware=" << (pfControllerCpuAware ? "yes" : "no")
                                 << ", cpuSampleMs=" << pfControllerCpuSampleMs
                                 << ", cpuEmaAlpha=" << pfControllerCpuEmaAlpha
                                 << ", cpuIdleThreshold=" << pfControllerCpuIdleThreshold
                                 << ", cpuBusyThreshold=" << pfControllerCpuBusyThreshold
                                 << ", mapEstimate=" << mapEstimateTotal
                                 << ", featureEstimate=" << featureEstimateTotal
                                 << ", libraryType=" << libraryType
                                 << ", sample=" << sampleName
                                 << "\n";
                if (pfReaderThreadsReserved >= P.runThreadN) {
                    P.inOut->logMain << "pf-dynamic-controller: WARNING readerThreadsReserved="
                                     << pfReaderThreadsReserved
                                     << " >= runThreadN=" << P.runThreadN
                                     << "; permit ceiling clamped to "
                                     << pfPermitCeilingDuringPf
                                     << " (oversubscription risk remains)\n";
                }

                pfControllerThread = std::thread([&]() {
                    size_t seqIndex = 0;
                    uint64_t mapEstimateLoop = mapEstimateTotal;
                    uint64_t featureEstimateLoop = featureEstimateTotal;
                    uint64_t prevMapDone = controllerBaseline.mapDomain.workUnitsTotal;
                    uint64_t prevFeatureDone = controllerBaseline.featureDomain.workUnitsTotal;
                    auto prevTick = std::chrono::steady_clock::now();
                    double mapRateEwma = 0.0;
                    double featureRateEwma = 0.0;
                    constexpr double kEtaAlpha = 0.30;
                    constexpr double kEtaGapThreshold = 0.10;
                    constexpr double kEtaRateEps = 1.0e-9;
                    constexpr long double kChunkGapThreshold = 0.10L;
                    constexpr int kCpuStableTicksRequired = 2;
                    int cpuIdleStableTicks = 0;
                    int cpuBusyStableTicks = 0;

                    while (!stopPfController.load(std::memory_order_relaxed)) {
                        if (pfControllerCpuAware) {
                            g_threadChunks.mapPermitCpuMaybeSample();
                        }
                        const ThreadControl::MapPermitSnapshot snapshot = g_threadChunks.mapPermitSnapshot();
                        int nextTarget = snapshot.targetPermits;
                        uint64_t mapRemainingForCpu = 0;
                        uint64_t featureRemainingForCpu = 0;
                        int permitCeilingForCpu = P.runThreadN;
                        if (pfControllerSequenceMode) {
                            nextTarget = pfControllerSequence[seqIndex % pfControllerSequence.size()];
                            ++seqIndex;
                            permitCeilingForCpu = pfPermitCeilingDuringPf;
                        } else if (pfControllerEtaMode) {
                            const uint64_t mapDone = snapshot.mapDomain.workUnitsTotal;
                            const uint64_t featureDone = snapshot.featureDomain.workUnitsTotal;
                            mapEstimateLoop = inflateEstimateUntilNonNegativeRemaining(mapEstimateLoop, mapDone);
                            featureEstimateLoop = inflateEstimateUntilNonNegativeRemaining(featureEstimateLoop, featureDone);
                            pfControllerMapEstimateFinal.store(mapEstimateLoop, std::memory_order_relaxed);
                            pfControllerFeatureEstimateFinal.store(featureEstimateLoop, std::memory_order_relaxed);

                            uint64_t mapRemaining = (mapEstimateLoop > mapDone) ? (mapEstimateLoop - mapDone) : 0;
                            uint64_t featureRemaining =
                                (featureEstimateLoop > featureDone) ? (featureEstimateLoop - featureDone) : 0;
                            mapRemainingForCpu = mapRemaining;
                            featureRemainingForCpu = featureRemaining;

                            const auto nowTick = std::chrono::steady_clock::now();
                            const double dt = std::chrono::duration_cast<std::chrono::duration<double>>(nowTick - prevTick).count();
                            if (dt > 0.0) {
                                const double mapInstRate = (mapDone >= prevMapDone)
                                    ? static_cast<double>(mapDone - prevMapDone) / dt
                                    : 0.0;
                                const double featureInstRate = (featureDone >= prevFeatureDone)
                                    ? static_cast<double>(featureDone - prevFeatureDone) / dt
                                    : 0.0;
                                mapRateEwma = (mapRateEwma <= 0.0)
                                    ? mapInstRate
                                    : (kEtaAlpha * mapInstRate + (1.0 - kEtaAlpha) * mapRateEwma);
                                featureRateEwma = (featureRateEwma <= 0.0)
                                    ? featureInstRate
                                    : (kEtaAlpha * featureInstRate + (1.0 - kEtaAlpha) * featureRateEwma);
                            }
                            prevTick = nowTick;
                            prevMapDone = mapDone;
                            prevFeatureDone = featureDone;

                            double mapEta = std::numeric_limits<double>::infinity();
                            double featureEta = std::numeric_limits<double>::infinity();
                            if (mapRemaining == 0) {
                                mapEta = 0.0;
                            } else if (mapRateEwma > kEtaRateEps) {
                                mapEta = static_cast<double>(mapRemaining) / mapRateEwma;
                            }
                            if (featureRemaining == 0) {
                                featureEta = 0.0;
                            } else if (featureRateEwma > kEtaRateEps) {
                                featureEta = static_cast<double>(featureRemaining) / featureRateEwma;
                            }

                            nextTarget = snapshot.targetPermits;
                            if (mapRemaining == 0 && featureRemaining > 0) {
                                nextTarget = pfPermitCeilingDuringPf;
                            } else if (featureRemaining == 0 && mapRemaining > 0) {
                                nextTarget = P.runThreadN;
                            } else if (mapRemaining > 0 && featureRemaining > 0 &&
                                       std::isfinite(mapEta) && std::isfinite(featureEta)) {
                                const double maxEta = std::max(mapEta, featureEta);
                                if (maxEta > 0.0) {
                                    const double etaGap = std::fabs(mapEta - featureEta) / maxEta;
                                    if (etaGap > kEtaGapThreshold) {
                                        if (featureEta > mapEta) {
                                            nextTarget = std::min(P.runThreadN, snapshot.targetPermits + 1);
                                        } else {
                                            nextTarget = std::max(1, snapshot.targetPermits - 1);
                                        }
                                    }
                                }
                            } else if (featureRemaining > 0 && !std::isfinite(featureEta)) {
                                nextTarget = std::min(P.runThreadN, snapshot.targetPermits + 1);
                            }
                            const int etaPermitCeiling =
                                (featureRemaining > 0) ? pfPermitCeilingDuringPf : P.runThreadN;
                            nextTarget = std::min(nextTarget, etaPermitCeiling);
                            permitCeilingForCpu = etaPermitCeiling;
                        } else if (pfControllerChunkedMode) {
                            const uint64_t mapDone = snapshot.mapDomain.workUnitsTotal;
                            const uint64_t featureDone = snapshot.featureDomain.workUnitsTotal;
                            mapEstimateLoop = inflateEstimateUntilNonNegativeRemaining(mapEstimateLoop, mapDone);
                            featureEstimateLoop = inflateEstimateUntilNonNegativeRemaining(featureEstimateLoop, featureDone);
                            pfControllerMapEstimateFinal.store(mapEstimateLoop, std::memory_order_relaxed);
                            pfControllerFeatureEstimateFinal.store(featureEstimateLoop, std::memory_order_relaxed);

                            const uint64_t mapRemaining =
                                (mapEstimateLoop > mapDone) ? (mapEstimateLoop - mapDone) : 0;
                            const uint64_t featureRemaining =
                                (featureEstimateLoop > featureDone) ? (featureEstimateLoop - featureDone) : 0;
                            mapRemainingForCpu = mapRemaining;
                            featureRemainingForCpu = featureRemaining;

                            nextTarget = snapshot.targetPermits;
                            if (mapRemaining == 0 && featureRemaining == 0) {
                                // nothing to retune
                            } else if (mapRemaining == 0 || featureRemaining == 0) {
                                if (featureRemaining > 0) {
                                    nextTarget = pfPermitCeilingDuringPf;
                                } else {
                                    // PF side drained: release back to map side
                                    nextTarget = P.runThreadN;
                                }
                            } else {
                                const uint64_t dominantRemaining = std::max(mapRemaining, featureRemaining);
                                const uint64_t remainingGap =
                                    (mapRemaining > featureRemaining)
                                    ? (mapRemaining - featureRemaining)
                                    : (featureRemaining - mapRemaining);
                                const long double gapFrac = (dominantRemaining == 0)
                                    ? 0.0L
                                    : (static_cast<long double>(remainingGap) /
                                       static_cast<long double>(dominantRemaining));
                                if (gapFrac > kChunkGapThreshold) {
                                    if (featureRemaining > mapRemaining) {
                                        nextTarget = std::max(1, snapshot.targetPermits - pfChunkPermitStep);
                                    } else {
                                        nextTarget = std::min(P.runThreadN, snapshot.targetPermits + pfChunkPermitStep);
                                    }
                                }
                            }
                            const int chunkPermitCeiling =
                                (featureRemaining > 0) ? pfPermitCeilingDuringPf : P.runThreadN;
                            nextTarget = std::min(nextTarget, chunkPermitCeiling);
                            permitCeilingForCpu = chunkPermitCeiling;
                        }

                        if (pfControllerSequenceMode) {
                            // Sequence mode has no remaining-work signal; always reserve producer/read threads.
                            nextTarget = std::min(nextTarget, pfPermitCeilingDuringPf);
                            permitCeilingForCpu = pfPermitCeilingDuringPf;
                        }

                        if (pfControllerCpuAware && snapshot.cpuAwareEnabled && snapshot.cpuInitialized) {
                            const bool canNudgeUp =
                                featureRemainingForCpu > 0 &&
                                nextTarget < permitCeilingForCpu &&
                                snapshot.cpuIdleEma >= pfControllerCpuIdleThreshold;
                            const bool canNudgeDown =
                                (mapRemainingForCpu > 0 || featureRemainingForCpu > 0) &&
                                nextTarget > 1 &&
                                snapshot.cpuBusyEma >= pfControllerCpuBusyThreshold;

                            if (canNudgeUp && !canNudgeDown) {
                                ++cpuIdleStableTicks;
                                cpuBusyStableTicks = 0;
                                if (cpuIdleStableTicks >= kCpuStableTicksRequired) {
                                    nextTarget = std::min(permitCeilingForCpu, nextTarget + 1);
                                    cpuIdleStableTicks = 0;
                                }
                            } else if (canNudgeDown) {
                                ++cpuBusyStableTicks;
                                cpuIdleStableTicks = 0;
                                if (cpuBusyStableTicks >= kCpuStableTicksRequired) {
                                    nextTarget = std::max(1, nextTarget - 1);
                                    cpuBusyStableTicks = 0;
                                }
                            } else {
                                cpuIdleStableTicks = 0;
                                cpuBusyStableTicks = 0;
                            }
                        } else {
                            cpuIdleStableTicks = 0;
                            cpuBusyStableTicks = 0;
                        }
                        if (nextTarget < 1) {
                            nextTarget = 1;
                        }
                        pfControllerLastTarget.store(nextTarget, std::memory_order_relaxed);
                        const uint64_t tick =
                            pfControllerTicks.fetch_add(1, std::memory_order_relaxed) + 1;
                        if (pfControllerAppliesUpdates) {
                            g_threadChunks.mapPermitSetTargetPermits(nextTarget);
                            pfControllerApplied.fetch_add(1, std::memory_order_relaxed);
                        }
                        if (pfControllerMaxUpdates > 0 &&
                            tick >= static_cast<uint64_t>(pfControllerMaxUpdates)) {
                            break;
                        }
                        std::this_thread::sleep_for(std::chrono::milliseconds(pfControllerIntervalMs));
                    }
                });
            }

            PfMultiAssign::AssignResult assignResult;
            try {
                assignResult = PfMultiAssign::runAssignBarcodes(wlInfo.normalizedPath, refPath, resolvedFastq, assignOut, runAssignOpts);
            } catch (...) {
                stopPfController.store(true, std::memory_order_relaxed);
                if (pfControllerThread.joinable()) {
                    pfControllerThread.join();
                }
                if (pfControllerAppliesUpdates) {
                    g_threadChunks.mapPermitSetTargetPermits(P.runThreadN);
                }
                throw;
            }

            stopPfController.store(true, std::memory_order_relaxed);
            if (pfControllerThread.joinable()) {
                pfControllerThread.join();
            }
            if (pfControllerAppliesUpdates) {
                // PF stage finished: release full permit budget back to map side.
                g_threadChunks.mapPermitSetTargetPermits(P.runThreadN);
            }

            if (pfControllerEnabled) {
                const string apiRunPath = assignOut + "/assignBarcodes.api_run.txt";
                std::ofstream apiRunOut(apiRunPath.c_str(), std::ios::app);
                if (apiRunOut.is_open()) {
                    apiRunOut << "pfController.mode=" << pfPermitControllerModeName(pfControllerMode) << "\n";
                    apiRunOut << "pfController.intervalMs=" << pfControllerIntervalMs << "\n";
                    apiRunOut << "pfController.sequence=" << pfControllerSequenceLabel << "\n";
                    apiRunOut << "pfController.maxUpdates=" << pfControllerMaxUpdates << "\n";
                    apiRunOut << "pfController.ticks=" << pfControllerTicks.load(std::memory_order_relaxed) << "\n";
                    apiRunOut << "pfController.applied=" << pfControllerApplied.load(std::memory_order_relaxed) << "\n";
                    apiRunOut << "pfController.lastTarget=" << pfControllerLastTarget.load(std::memory_order_relaxed) << "\n";
                    apiRunOut << "pfController.mapEstimateInitial=" << pfControllerMapEstimateInitial.load(std::memory_order_relaxed) << "\n";
                    apiRunOut << "pfController.featureEstimateInitial=" << pfControllerFeatureEstimateInitial.load(std::memory_order_relaxed) << "\n";
                    apiRunOut << "pfController.mapEstimateFinal=" << pfControllerMapEstimateFinal.load(std::memory_order_relaxed) << "\n";
                    apiRunOut << "pfController.featureEstimateFinal=" << pfControllerFeatureEstimateFinal.load(std::memory_order_relaxed) << "\n";
                    apiRunOut << "pfController.estimatePolicy=smallest_anchor_bytescale\n";
                    apiRunOut << "pfController.estimateNegativeRemainingInflationPct=10\n";
                    apiRunOut << "pfController.producerReservePolicy="
                              << "expand_producers_by_consumer_cap"
                              << "\n";
                    apiRunOut << "pfController.searchThreadsForController=" << pfSearchThreadsForController << "\n";
                    apiRunOut << "pfController.searchThreadQuantum=" << pfSearchThreadQuantum << "\n";
                    apiRunOut << "pfController.chunkPermitStep=" << pfChunkPermitStep << "\n";
                    apiRunOut << "pfController.consumerThreadsForRun=" << pfConsumerThreadsForRun << "\n";
                    apiRunOut << "pfController.consumerSoftMax=" << pfConsumerSoftMax << "\n";
                    apiRunOut << "pfController.consumerThreadsExplicit=" << (pfConsumerThreadsExplicit ? 1 : 0) << "\n";
                    apiRunOut << "pfController.consumerThreadsAuto=" << (pfConsumerThreadsAuto ? 1 : 0) << "\n";
                    apiRunOut << "pfController.consumerBudgetThreads=" << pfConsumerBudgetThreads << "\n";
                    apiRunOut << "pfController.producerSlotsAvailable=" << pfProducerSlotsAvailable << "\n";
                    apiRunOut << "pfController.singleProducerOnly=" << (pfSingleProducerOnly ? 1 : 0) << "\n";
                    apiRunOut << "pfController.singleConsumerMode=" << (pfSingleConsumer ? 1 : 0) << "\n";
                    apiRunOut << "pfController.readerThreadsReserved=" << pfReaderThreadsReserved << "\n";
                    apiRunOut << "pfController.permitCeilingDuringPf=" << pfPermitCeilingDuringPf << "\n";
                    apiRunOut << "pfController.cpuAware=" << (pfControllerCpuAware ? 1 : 0) << "\n";
                    apiRunOut << "pfController.cpuSampleMs=" << pfControllerCpuSampleMs << "\n";
                    apiRunOut << "pfController.cpuEmaAlpha=" << pfControllerCpuEmaAlpha << "\n";
                    apiRunOut << "pfController.cpuIdleThreshold=" << pfControllerCpuIdleThreshold << "\n";
                    apiRunOut << "pfController.cpuBusyThreshold=" << pfControllerCpuBusyThreshold << "\n";
                    if (pfControllerEtaMode) {
                        apiRunOut << "pfController.etaGapThresholdPct=10\n";
                    }
                    if (pfControllerChunkedMode) {
                        apiRunOut << "pfController.chunkGapThresholdPct=10\n";
                        apiRunOut << "pfController.chunkPriorityPolicy=remaining_work_coarse_quantum\n";
                    }
                    apiRunOut << "pfController.mapAnchorFile=" << mapEstimate.anchorFile << "\n";
                    apiRunOut << "pfController.mapAnchorReads=" << mapEstimate.anchorReads << "\n";
                    apiRunOut << "pfController.mapAnchorBytes=" << mapEstimate.anchorBytes << "\n";
                    apiRunOut << "pfController.mapTotalBytes=" << mapEstimate.totalBytes << "\n";
                    apiRunOut << "pfController.featureAnchorFile=" << featureEstimate.anchorFile << "\n";
                    apiRunOut << "pfController.featureAnchorReads=" << featureEstimate.anchorReads << "\n";
                    apiRunOut << "pfController.featureAnchorBytes=" << featureEstimate.anchorBytes << "\n";
                    apiRunOut << "pfController.featureTotalBytes=" << featureEstimate.totalBytes << "\n";
                    const ThreadControl::MapPermitSnapshot finalSnapshot = g_threadChunks.mapPermitSnapshot();
                    apiRunOut << "pfController.cpuInitializedFinal=" << (finalSnapshot.cpuInitialized ? 1 : 0) << "\n";
                    apiRunOut << "pfController.cpuSampleCountFinal=" << finalSnapshot.cpuSampleCount << "\n";
                    apiRunOut << "pfController.cpuBusyInstantFinal=" << finalSnapshot.cpuBusyInstant << "\n";
                    apiRunOut << "pfController.cpuBusyEmaFinal=" << finalSnapshot.cpuBusyEma << "\n";
                    apiRunOut << "pfController.cpuIdleEmaFinal=" << finalSnapshot.cpuIdleEma << "\n";
                } else {
                    P.inOut->logMain << "WARNING: failed to append pf controller summary to "
                                     << apiRunPath << "\n";
                }
            }

            if (pfControllerEnabled) {
                const ThreadControl::MapPermitSnapshot controllerStopSnapshot = g_threadChunks.mapPermitSnapshot();
                P.inOut->logMain << "pf-dynamic-controller: stop mode="
                                 << pfPermitControllerModeName(pfControllerMode)
                                 << ", ticks=" << pfControllerTicks.load(std::memory_order_relaxed)
                                 << ", applied=" << pfControllerApplied.load(std::memory_order_relaxed)
                                 << ", lastTarget=" << pfControllerLastTarget.load(std::memory_order_relaxed)
                                 << ", mapEstimateFinal=" << pfControllerMapEstimateFinal.load(std::memory_order_relaxed)
                                 << ", featureEstimateFinal=" << pfControllerFeatureEstimateFinal.load(std::memory_order_relaxed)
                                 << ", cpuInitialized=" << (controllerStopSnapshot.cpuInitialized ? "yes" : "no")
                                 << ", cpuSampleCount=" << controllerStopSnapshot.cpuSampleCount
                                 << ", cpuBusyEma=" << controllerStopSnapshot.cpuBusyEma
                                 << ", cpuIdleEma=" << controllerStopSnapshot.cpuIdleEma
                                 << ", libraryType=" << libraryType
                                 << ", sample=" << sampleName
                                 << "\n";
            }

            if (assignResult.returnCode != 0) {
                throw runtime_error("Failed to process feature library: type=" + libraryType
                    + ", libraryId=" + preparedLib.libraryId
                    + ", sample=" + sampleName
                    + ", featureRef=" + refPath
                    + ", fastq=" + resolvedFastq);
            }
            P.inOut->logMain << "NOTICE: assignBarcodes completed for library_id="
                             << preparedLib.libraryId
                             << " input_format=" << assignResult.inputFormat
                             << " cbq_mode_requested=" << assignResult.cbqModeRequested
                             << " cbq_mode_effective=" << assignResult.cbqModeEffective;
            if (!assignResult.cbqModeFallbackReason.empty()) {
                P.inOut->logMain << " cbq_mode_fallback_reason="
                                 << assignResult.cbqModeFallbackReason;
            }
            P.inOut->logMain << "\n";

            appendAssignNormalizationStats(assignOut, nsCtx, wlInfo);

            PfMultiFeatureRun run;
            run.featureType = featureRefType;
            run.assignOut = assignOut;
            run.featureRefPath = refPath;
            run.whitelistPath = wlInfo.normalizedPath;
            run.barcodeOutputMapPath =
                !preparedLib.starBarcodeOutputMap.empty()
                    ? preparedLib.starBarcodeOutputMap
                    : ((wlInfo.hasTwoColumnSource && assignmentWhitelistNamespace == "NXT")
                        ? wlInfo.sourcePath
                        : "");
            run.effectiveChem = effectiveReadNamespace;
            run.effectiveReadNamespace = nsCtx.effectiveReadNamespace;
            run.assignmentWhitelistNamespace = nsCtx.assignmentWhitelistNamespace;
            run.translateNxtForAssign = false;
            // Whitelist is now in read namespace; translateNxt is false so
            // assignBarcodes outputs barcodes in the assignment/read namespace.
            // A two-column 10x translation whitelist maps NXT assignment
            // barcodes to TRU MEX barcodes at stub generation. Keep the MEX
            // namespace aligned with the barcode TSV that downstream merge reads.
            run.featureMexOutputNamespace =
                !run.barcodeOutputMapPath.empty()
                    ? (splitReadLayout ? "GEX" : "TRU")
                    : nsCtx.assignmentWhitelistNamespace;
            run.outputNamespace = nsCtx.outputNamespace;
            run.namespaceConfident = nsCtx.isNamespaceConfident;
            run.detectedMatchMode = detectedMatchMode;
            run.whitelistNormalization = wlInfo;
            run.normalizationStats = nsCtx.normalizationStats;
            run.libraryId = preparedLib.libraryId;
            run.sampleName = sampleName;
            run.resolvedFastq = resolvedFastq;
            run.resolvedChemRequest = preparedLib.resolvedChemRequest;
            run.explicitChem = preparedLib.explicitChem;
            run.adtMexOutput = preparedLib.adtMexOutput;
            run.returnCode = assignResult.returnCode;
            featureRuns.push_back(run);
        }

        if (featureRuns.empty()) {
            throw runtime_error("No feature libraries found in multi config");
        }

        // Clear stale provenance artifacts before new validation/write cycle.
        // A previous successful run may have left status=OK manifests that would
        // be misleading if this run fails partway through.
        for (auto& run : featureRuns) {
            string realOut = findAssignOutputDir(run.assignOut);
            run.assignOut = realOut;
            for (const char* suffix : {"pf_library_provenance.tsv", "pf_library_provenance.tsv.tmp"}) {
                string path = run.assignOut + "/" + suffix;
                if (remove(path.c_str()) != 0 && errno != ENOENT) {
                    throw runtime_error("Failed to remove stale provenance artifact: " + path
                        + " (" + strerror(errno) + ")");
                }
            }
        }

        for (auto& run : featureRuns) {
            int stubRet = PfMultiMexStub::processAssignOutput(
                run.assignOut, run.featureRefPath, run.featureType,
                true,
                run.barcodeOutputMapPath.empty() ? run.whitelistPath : run.barcodeOutputMapPath,
                run.featureType);
            if (stubRet != 0) {
                throw runtime_error("Failed to generate MEX stub for library: type="
                    + run.featureType + ", library_id=" + run.libraryId
                    + ", assign_out=" + run.assignOut);
            }
            // With the namespace refactor, the whitelist is pre-normalized to
            // the read namespace before assignment. barcodes.txt is already in
            // the assignment (read) namespace. copyBarcodesTsv may attempt a
            // COL1->COL2 remap via the original whitelist, but since
            // barcodes.txt barcodes won't match COL1 keys (different namespace),
            // they pass through unchanged. featureMexOutputNamespace was set
            // correctly in the assignment loop and should not be overridden here.
        }

        // Read feature MEX to validate before writing provenance (fail-fast on read errors)
        vector<PfMultiFeatureMexEntry> featureMexEntries;
        for (const auto& run : featureRuns) {
            try {
                PfMultiFeatureMexEntry entry;
                entry.data = PfMultiMerge::readMex(run.assignOut);
                entry.effectiveChem = run.effectiveChem;
                entry.featureMexOutputNamespace = run.featureMexOutputNamespace;
                entry.featureType = run.featureType;
                entry.libraryId = run.libraryId;
                featureMexEntries.push_back(std::move(entry));
            } catch (const exception& e) {
                throw runtime_error("Failed to read feature MEX for library: type="
                    + run.featureType + ", library_id=" + run.libraryId
                    + ", assign_out=" + run.assignOut + ": " + e.what());
            }
        }

        // All libraries passed assign + stub + MEX read.
        // Write provenance manifests atomically: tmp file + rename.
        // On failure, clean up any .tmp files already written.
        vector<string> manifestTmpPaths;
        vector<string> manifestFinalPaths;
        try {
            for (const auto& run : featureRuns) {
                string finalPath = run.assignOut + "/pf_library_provenance.tsv";
                string tmpPath = finalPath + ".tmp";
                manifestTmpPaths.push_back(tmpPath);
                manifestFinalPaths.push_back(finalPath);

                std::ofstream manifest(tmpPath.c_str());
                if (!manifest.is_open()) {
                    throw runtime_error("Failed to create provenance manifest: " + tmpPath);
                }
                manifest << "key\tvalue\n";
                manifest << "library_id\t" << run.libraryId << "\n";
                manifest << "sample\t" << run.sampleName << "\n";
                manifest << "feature_type\t" << run.featureType << "\n";
                manifest << "adt_mex_output\t" << (run.adtMexOutput ? "yes" : "no") << "\n";
                manifest << "fastq_dir\t" << run.resolvedFastq << "\n";
                manifest << "feature_ref\t" << run.featureRefPath << "\n";
                manifest << "whitelist\t" << run.whitelistPath << "\n";
                manifest << "barcode_output_map\t" << run.barcodeOutputMapPath << "\n";
                manifest << "chemistry_request\t" << run.resolvedChemRequest << "\n";
                manifest << "chemistry_explicit\t" << (run.explicitChem ? "yes" : "no") << "\n";
                manifest << "effective_chemistry\t" << run.effectiveChem << "\n";
                manifest << "detected_match_mode\t" << run.detectedMatchMode << "\n";
                manifest << "effective_read_namespace\t" << run.effectiveReadNamespace << "\n";
                manifest << "assignment_whitelist_namespace\t" << run.assignmentWhitelistNamespace << "\n";
                manifest << "translate_nxt_for_assign\t" << (run.translateNxtForAssign ? "yes" : "no") << "\n";
                manifest << "output_namespace\t" << run.outputNamespace << "\n";
                manifest << "namespace_confidence\t" << (run.namespaceConfident ? "yes" : "no") << "\n";
                manifest << "whitelist_source\t" << run.whitelistNormalization.sourcePath << "\n";
                manifest << "whitelist_normalized\t" << run.whitelistNormalization.normalizedPath << "\n";
                manifest << "whitelist_has_two_columns\t"
                         << (run.whitelistNormalization.hasTwoColumnSource ? "yes" : "no") << "\n";
                manifest << "whitelist_normalized_rows\t" << run.whitelistNormalization.normalizedRowCount << "\n";
                manifest << "barcode_norm_in_set\t" << run.normalizationStats.inSet << "\n";
                manifest << "barcode_norm_translated_to_set\t" << run.normalizationStats.translatedToSet << "\n";
                manifest << "barcode_norm_unmatched\t" << run.normalizationStats.unmatched << "\n";
                manifest << "barcode_norm_dedup_dropped\t" << run.normalizationStats.dedupDropped << "\n";
                manifest << "barcode_norm_output_count\t" << run.normalizationStats.outputCount << "\n";
                manifest << "return_code\t" << run.returnCode << "\n";
                manifest << "status\tOK\n";
                manifest << "assign_output_dir\t" << run.assignOut << "\n";
                manifest.close();
                if (manifest.fail()) {
                    throw runtime_error("Failed to flush provenance manifest: " + tmpPath);
                }
            }
            // All tmp files written successfully; atomically promote them
            for (size_t i = 0; i < manifestTmpPaths.size(); ++i) {
                if (rename(manifestTmpPaths[i].c_str(), manifestFinalPaths[i].c_str()) != 0) {
                    throw runtime_error("Failed to rename provenance manifest: "
                        + manifestTmpPaths[i] + " -> " + manifestFinalPaths[i]);
                }
            }
        } catch (...) {
            for (const auto& tmp : manifestTmpPaths) {
                remove(tmp.c_str());
            }
            throw;
        }

        // Log top-level multi-library summary table.
        // All libraries are OK here (fail-fast throws before reaching this point).
        P.inOut->logMain << "\npf-multi library summary:\n";
        P.inOut->logMain << "  library_id\tfeature_type\tsample\teffective_chem\t"
                         << "match_mode\tstatus\tassign_out\n";
        for (const auto& run : featureRuns) {
            P.inOut->logMain << "  " << run.libraryId
                             << "\t" << run.featureType
                             << "\t" << run.sampleName
                             << "\t" << run.effectiveChem
                             << "\t" << run.detectedMatchMode
                             << "\tOK"
                             << "\t" << run.assignOut
                             << "\n";
        }
        P.inOut->logMain << "\n";
        
        phase->featureMexEntries = std::move(featureMexEntries);
        P.inOut->logMain << timeMonthDayTime() << " ..... finished pf-multi feature assignment\n";
        return phase;
    } catch (const exception& e) {
        P.inOut->logMain << "ERROR in pf-multi feature assignment: " << e.what() << "\n";
        return nullptr;
    }
}

std::shared_ptr<PfFeatureAssignHandle> startPfFeatureAssignment(
    Parameters& P,
    const std::shared_ptr<PfMultiPreloadHandle>& preload) {
    if (P.runMode != "alignReads") {
        return nullptr;
    }
    if (isUnsetToken(P.pfMulti.pfMultiConfig)) {
        return nullptr;
    }

    auto handle = std::make_shared<PfFeatureAssignHandle>();
    handle->started = true;
    handle->worker = std::thread([handle, &P, preload]() {
        try {
            handle->result = runPfMultiAssignPhase(P, preload);
            handle->success = (handle->result != nullptr);
            if (!handle->success) {
                handle->error = "raw pf-multi feature-assignment phase returned no result";
            }
        } catch (const std::exception& e) {
            handle->success = false;
            handle->error = e.what();
        } catch (...) {
            handle->success = false;
            handle->error = "unknown exception";
        }
        handle->finished = true;
    });
    return handle;
}

int finalizePfMultiConfig(Parameters& P,
                          const Solo* solo,
                          const std::shared_ptr<PfMultiAssignPhaseResult>& assignPhase) {
    if (!assignPhase) {
        return 0;
    }
    try {
        const PfMultiPreparedContext& prepared = assignPhase->prepared;
        const PfMultiConfig::Config& config = prepared.config;
        const string& soloWhitelist = prepared.soloWhitelist;
        const string gexWhitelistNamespace = upperCopy(prepared.soloEffectiveChem);
        const string gexNormalizationChem = upperCopy(prepared.soloOutputNamespace);
        const string& outputChem = prepared.outputChem;
        const string& outPrefix = prepared.outPrefix;
        const vector<PfMultiFeatureRun>& featureRuns = assignPhase->featureRuns;
        vector<PfMultiFeatureMexEntry> featureMexEntries = assignPhase->featureMexEntries;

        P.inOut->logMain << timeMonthDayTime() << " ..... started pf-multi merge/finalize\n";

        string soloOut = outPrefix + "/Solo.out";
        string geneOut = soloOut + "/Gene";
        string geneFullOut = soloOut + "/GeneFull";
        string filteredOut = geneOut + "/filtered";
        string rawOut = geneOut + "/raw";
        string geneFullFiltered = geneFullOut + "/filtered";
        string geneFullRaw = geneFullOut + "/raw";

        if (P.pSolo.crGexFeature == ParametersSolo::CrGexGeneFull) {
            rawOut = geneFullRaw;
            filteredOut = geneFullFiltered;
            P.inOut->logMain << "NOTICE: --soloCrGexFeature=genefull (using GeneFull MEX for CR-compat merge)\n";
        } else if (P.pSolo.crGexFeature == ParametersSolo::CrGexGene) {
            rawOut = geneOut + "/raw";
            filteredOut = geneOut + "/filtered";
            P.inOut->logMain << "NOTICE: --soloCrGexFeature=gene (using Gene MEX for CR-compat merge)\n";
        }

        bool hasRaw = (P.pSolo.type != 0 && hasMexFiles(rawOut));
        bool hasFiltered = (P.pSolo.type != 0 && hasMexFiles(filteredOut));
        if (!hasRaw && !hasFiltered) {
            if (P.pSolo.crGexFeature == ParametersSolo::CrGexGeneFull) {
                P.inOut->logMain << "ERROR: GeneFull MEX directory not found for CR-compat merge\n";
                return 1;
            }
            if (P.pSolo.crGexFeature == ParametersSolo::CrGexGene) {
                P.inOut->logMain << "ERROR: Gene MEX directory not found for CR-compat merge\n";
                return 1;
            }
            bool hasGeneFullRaw = (P.pSolo.type != 0 && hasMexFiles(geneFullRaw));
            bool hasGeneFullFiltered = (P.pSolo.type != 0 && hasMexFiles(geneFullFiltered));
            if (hasGeneFullRaw || hasGeneFullFiltered) {
                rawOut = geneFullRaw;
                filteredOut = geneFullFiltered;
                hasRaw = hasGeneFullRaw;
                hasFiltered = hasGeneFullFiltered;
                P.inOut->logMain << "NOTICE: Using GeneFull MEX for CR-compat merge (Gene missing)\n";
            } else if (P.pSolo.type != 0 && hasMexFiles(geneOut)) {
                hasRaw = true;
                rawOut = geneOut;
            } else {
                P.inOut->logMain << "WARNING: GEX MEX directory not found, skipping merge\n";
                return 0;
            }
        }

        vector<string> filteredGexBarcodes;
        bool useFilteredGex = false;
        std::unordered_set<string> gexWhitelistSet;
        uint64_t gexWhitelistRows = 0;
        uint64_t gexWhitelistInvalidRows = 0;
        const bool haveGexWhitelistSet =
            loadWhitelistBarcodeSet(soloWhitelist, gexWhitelistSet, gexWhitelistRows, gexWhitelistInvalidRows);
        if (!haveGexWhitelistSet) {
            P.inOut->logMain << "WARNING: failed to load whitelist set for GEX barcode normalization: "
                             << soloWhitelist << "\n";
        }
        if (solo && getFilteredBarcodesFromSolo(solo, P, filteredGexBarcodes, true)) {
            if (haveGexWhitelistSet && isKnownNamespace(gexWhitelistNamespace)) {
                FilteredBarcodeNormalizationStats gexNormStats;
                filteredGexBarcodes =
                    normalizeFilteredBarcodesForAssignNamespace(
                        filteredGexBarcodes, gexWhitelistSet, gexNormStats,
                        gexNormalizationChem, gexWhitelistNamespace);
                P.inOut->logMain
                    << "NOTICE: gex filtered-barcode normalization from Solo"
                    << " sourceNS=" << gexNormalizationChem
                    << " whitelistNS=" << gexWhitelistNamespace
                    << " in_set=" << gexNormStats.inSet
                    << " translated_to_set=" << gexNormStats.translatedToSet
                    << " unmatched=" << gexNormStats.unmatched
                    << " dedup_dropped=" << gexNormStats.dedupDropped
                    << " output=" << gexNormStats.outputCount
                    << "\n";
                normalizeBarcodeVecToTru(filteredGexBarcodes, gexWhitelistNamespace);
            } else {
                normalizeBarcodeVecToTru(filteredGexBarcodes, gexNormalizationChem);
            }
            useFilteredGex = true;
            P.inOut->logMain << "NOTICE: Using GEX filtered barcodes from Solo "
                             << "(normalized to TRU, in-memory)\n";
        }

        PfMultiMerge::MexData gexRawData;
        try {
            if (hasRaw) {
                gexRawData = PfMultiMerge::readMex(rawOut);
                bool hasMultipleTypes = false;
                for (const auto& type : gexRawData.featureTypes) {
                    if (type != "Gene Expression" && !type.empty()) {
                        hasMultipleTypes = true;
                        break;
                    }
                }
                if (hasMultipleTypes) {
                    gexRawData = PfMultiMerge::filterByFeatureType(gexRawData, "Gene Expression");
                }
            }
        } catch (const exception& e) {
            P.inOut->logMain << "WARNING: Failed to read raw GEX MEX: " << e.what() << "\n";
            if (!hasFiltered) {
                return 0;
            }
        }

        PfMultiMerge::MexData gexFilteredData;
        bool loadedFilteredMex = false;
        bool needFilteredMexForCounts = (!hasRaw || gexRawData.features.empty());
        bool needFilteredMexForBarcodes = (!useFilteredGex && hasFiltered);
        if (hasFiltered && (needFilteredMexForCounts || needFilteredMexForBarcodes)) {
            try {
                gexFilteredData = PfMultiMerge::readMex(filteredOut);
                bool hasMultipleTypes = false;
                for (const auto& type : gexFilteredData.featureTypes) {
                    if (type != "Gene Expression" && !type.empty()) {
                        hasMultipleTypes = true;
                        break;
                    }
                }
                if (hasMultipleTypes) {
                    gexFilteredData = PfMultiMerge::filterByFeatureType(gexFilteredData, "Gene Expression");
                }
                loadedFilteredMex = true;
            } catch (const exception& e) {
                P.inOut->logMain << "WARNING: Failed to read filtered GEX MEX: " << e.what()
                                 << ", will use observed raw GEX barcodes for filtered output\n";
            }
        }

        PfMultiMerge::MexData gexData;
        if (hasRaw && !gexRawData.features.empty()) {
            gexData = gexRawData;
        } else if (hasFiltered && loadedFilteredMex && !gexFilteredData.features.empty()) {
            gexData = gexFilteredData;
            P.inOut->logMain << "WARNING: Raw GEX MEX not available, using filtered GEX MEX as merge base\n";
        } else {
            P.inOut->logMain << "ERROR: No valid GEX MEX data available for merging\n";
            return 1;
        }

        normalizeMexBarcodesToTru(gexData, gexNormalizationChem);
        normalizeMexBarcodesToTru(gexRawData, gexNormalizationChem);
        if (loadedFilteredMex) {
            normalizeMexBarcodesToTru(gexFilteredData, gexNormalizationChem);
        }
        for (auto& entry : featureMexEntries) {
            normalizeMexBarcodesToTru(entry.data, entry.featureMexOutputNamespace);
        }
        P.inOut->logMain << "Barcode namespace: all inputs normalized to TRU before merge"
                         << " (GEX soloOutputNS=" << gexNormalizationChem
                         << ", whitelistNS=" << gexWhitelistNamespace
                         << ", outputChem=" << outputChem << ")\n";

        vector<PfMultiMerge::MexData> featureDataVec;
        featureDataVec.reserve(featureMexEntries.size());
        for (auto& entry : featureMexEntries) {
            featureDataVec.push_back(std::move(entry.data));
        }
        PfMultiMerge::MexData mergedData = PfMultiMerge::mergeMex(gexData, featureDataVec);
        enforceUniqueBarcodesOrThrow(mergedData.barcodes, "merged MEX barcode list");

        {
            map<string, size_t> typeCountMap;
            for (const auto& ft : mergedData.featureTypes) {
                typeCountMap[ft]++;
            }
            P.inOut->logMain << "pf-multi merged feature breakdown:\n";
            for (const auto& pair : typeCountMap) {
                P.inOut->logMain << "  " << pair.first << ": " << pair.second << " features\n";
            }
            P.inOut->logMain << "  total barcodes: " << mergedData.barcodes.size()
                             << ", total triplets: " << mergedData.triplets.size() << "\n";
        }

        string gemWell = "1";
        vector<PfMultiConfig::LibraryEntry> gexLibs = config.getGexLibraries();
        if (!gexLibs.empty()) {
            gemWell = gexLibs[0].gem_well;
            for (size_t i = 1; i < gexLibs.size(); ++i) {
                if (gexLibs[i].gem_well != gemWell) {
                    P.inOut->logMain << "WARNING: Multiple GEX libraries have different gem_well values ("
                                     << gemWell << " vs " << gexLibs[i].gem_well
                                     << "), using first: " << gemWell << "\n";
                    break;
                }
            }
        } else if (!config.libraries.empty()) {
            gemWell = config.libraries[0].gem_well;
        }

        vector<string> observedRawGexBarcodes;
        if (hasRaw && !gexRawData.features.empty()) {
            observedRawGexBarcodes = PfMultiMerge::computeObservedGexBarcodes(gexRawData);
        } else if (hasFiltered && !gexFilteredData.features.empty()) {
            observedRawGexBarcodes = PfMultiMerge::computeObservedGexBarcodes(gexFilteredData);
            P.inOut->logMain << "WARNING: Raw GEX MEX not available, using filtered GEX barcodes for raw output\n";
        } else {
            P.inOut->logMain << "ERROR: Cannot compute observed GEX barcodes - no valid GEX data\n";
            return 1;
        }

        if (!useFilteredGex && hasFiltered && loadedFilteredMex && !gexFilteredData.features.empty()) {
            filteredGexBarcodes = PfMultiMerge::computeObservedGexBarcodes(gexFilteredData);
            useFilteredGex = true;
        }
        if (!useFilteredGex) {
            filteredGexBarcodes = observedRawGexBarcodes;
            P.inOut->logMain << "WARNING: Filtered GEX barcodes not available, using observed raw GEX barcodes for filtered output\n";
        }
        enforceUniqueBarcodesOrThrow(observedRawGexBarcodes, "observed raw GEX barcode list");
        enforceUniqueBarcodesOrThrow(filteredGexBarcodes, "filtered GEX barcode list");

        if (!assignPhase->usedExplicitAssignFilteredBarcodes) {
            for (const auto& run : featureRuns) {
                writeDeferredFilteredAssignOutput(
                    run.assignOut, filteredGexBarcodes, run.featureMexOutputNamespace,
                    run.barcodeOutputMapPath.empty() ? run.whitelistPath : run.barcodeOutputMapPath,
                    P.inOut->logMain);
            }
        }

        string rawOutDir = outPrefix + "/outs/raw_feature_bc_matrix";
        int ret = PfMultiMerge::writeCombinedMex(rawOutDir,
                                                 mergedData,
                                                 gemWell,
                                                 P.inOut->logMain,
                                                 observedRawGexBarcodes,
                                                 "TRU",
                                                 outputChem);
        if (ret != 0) {
            throw runtime_error("Failed to write raw combined MEX");
        }
        P.inOut->logMain << "Raw MEX written to: " << rawOutDir << "\n";

        string filteredOutDir = outPrefix + "/outs/filtered_feature_bc_matrix";
        ret = PfMultiMerge::writeCombinedMex(filteredOutDir,
                                             mergedData,
                                             gemWell,
                                             P.inOut->logMain,
                                             filteredGexBarcodes,
                                             "TRU",
                                             outputChem);
        if (ret != 0) {
            throw runtime_error("Failed to write filtered combined MEX");
        }
        P.inOut->logMain << "Filtered MEX written to: " << filteredOutDir << "\n";

        bool hasCrisprFeatures = false;
        for (const auto& run : featureRuns) {
            if (run.featureType == "CRISPR Guide Capture") {
                hasCrisprFeatures = true;
                break;
            }
        }
        if (hasCrisprFeatures) {
            string outsDir = outPrefix + "/outs";
            ret = runCrisprFeatureCalling(filteredOutDir, outsDir, P.pfMulti.crMinUmi, P.inOut->logMain);
            if (ret != 0) {
                P.inOut->logMain << "WARNING: CRISPR feature calling failed, continuing without crispr_analysis/\n";
            }
        }

        P.inOut->logMain << timeMonthDayTime() << " ..... finished pf-multi config processing\n";
    } catch (const exception& e) {
        P.inOut->logMain << "ERROR in pf-multi config processing: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

int processPfMultiConfig(Parameters& P,
                         const Solo* solo,
                         const std::shared_ptr<PfMultiPreloadHandle>& preload,
                         const std::shared_ptr<PfFeatureAssignHandle>& assignHandle) {
    std::shared_ptr<PfMultiAssignPhaseResult> phase;
    if (isUnsetToken(P.pfMulti.pfMultiConfig)) {
        return 0;
    }

    if (assignHandle && assignHandle->started) {
        assignHandle->join();
        if (!assignHandle->isSuccess()) {
            P.inOut->logMain << "ERROR: pf-feature-assign Phase A failed: "
                             << (assignHandle->getError().empty() ? string("unknown") : assignHandle->getError())
                             << "\n";
            return 1;
        }
        phase = assignHandle->getResult();
    } else {
        phase = runPfMultiAssignPhase(P, preload);
    }

    if (!phase) {
        return 1;
    }
    return finalizePfMultiConfig(P, solo, phase);
}

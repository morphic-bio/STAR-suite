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
#include <sys/stat.h>
#include <dirent.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cctype>
#include <cmath>
#include <chrono>
#include <thread>
#include <atomic>
#include <limits>
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
    string assignOut;
    string featureRefPath;
    string effectiveChem;
    string resolvedChemRequest;   // "auto", "nxt", or "tru" — after column > global resolution
    bool explicitChem = false;    // true when star_chemistry column provided NXT/TRU
    bool usedFilteredRef = false;
    AnchoredReadEstimate featureEstimate;
};

struct PfMultiPreparedContext {
    PfMultiConfig::Config config;
    string featureRef;
    string whitelist;
    string requestedChem;
    string requestedOutputChem;
    string inferredChem;
    string inferredReason;
    bool inferredChemConfident = false;
    string effectiveChem;
    string outputChem;
    string outPrefix;
    string crAssignRoot;
    AnchoredReadEstimate mapEstimate;
    vector<PfPreparedFeatureLibrary> featureLibraries;
    string prepLog;
};

class PfMultiPreloadHandle {
public:
    ~PfMultiPreloadHandle() {
        if (worker.joinable()) {
            worker.join();
        }
    }

private:
    friend std::shared_ptr<PfMultiPreloadHandle> startPfMultiConfigPreload(const Parameters& P);
    friend int processPfMultiConfig(Parameters& P,
                                    const Solo* solo,
                                    const std::shared_ptr<PfMultiPreloadHandle>& preload);

    std::thread worker;
    bool started = false;
    bool finished = false;
    bool success = false;
    double preloadSec = 0.0;
    string error;
    PfMultiPreparedContext context;
};

namespace {

struct FeatureSpec {
    string libraryType;
    string featureRefType;
};

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

static void normalizeMexBarcodesToTru(PfMultiMerge::MexData& data, const string& fromChem) {
    string norm = upperCopy(fromChem);
    if (norm != "NXT") {
        return;
    }
    for (auto& bc : data.barcodes) {
        if (bc.size() >= 9) {
            bc[7] = complementBase(bc[7]);
            bc[8] = complementBase(bc[8]);
        }
    }
}

static void normalizeBarcodeVecToTru(vector<string>& barcodes, const string& fromChem) {
    string norm = upperCopy(fromChem);
    if (norm != "NXT") {
        return;
    }
    for (auto& bc : barcodes) {
        if (bc.size() >= 9) {
            bc[7] = complementBase(bc[7]);
            bc[8] = complementBase(bc[8]);
        }
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
        return "unknown";
    }

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
                reason = msg.str();
                confident = true;
                return "NXT";
            }
            msg << " — below 80% threshold, falling through to filename/default";
            reason = msg.str();
        }
    }

    // Fallback to filename hints.
    string base = lowerCopy(basenameOf(whitelistPath));
    if (base.find("nxt") != string::npos) {
        reason = "whitelist filename contains 'nxt'";
        confident = true;
        return "NXT";
    }
    if (base.find("tru") != string::npos) {
        reason = "whitelist filename contains 'tru'";
        confident = true;
        return "TRU";
    }

    reason = "no chemistry marker found in whitelist content or filename; defaulting to TRU";
    confident = false;
    return "TRU";
}

static string oppositeNamespace(const string& ns) {
    if (ns == "NXT") return "TRU";
    if (ns == "TRU") return "NXT";
    return ns;
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

static vector<FeatureSpec> defaultPfFeatureSpecs() {
    return {
        {"CRISPR Guide Capture", "CRISPR Guide Capture"},
        {"Antibody Capture", "Antibody Capture"},
        {"CellPlex (CMO)", "Multiplexing Capture"},
        {"Multiplexing Capture", "Multiplexing Capture"}
    };
}

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
    if (isUnsetToken(context.featureRef)) {
        throw runtime_error("Feature reference not provided (use --crFeatureRef or set in config)");
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
            << " output_requested=" << context.requestedOutputChem
            << " output_effective=" << context.outputChem
            << "\n";

    // Per-library GEX chemistry: star_chemistry column on GEX row overrides the global.
    // Must update both effectiveChem AND inferredChem so that auto-detect composition
    // for feature libraries uses the correct wlNamespace anchor.
    {
        vector<PfMultiConfig::LibraryEntry> gexLibs = context.config.getGexLibraries();
        for (const auto& gex : gexLibs) {
            if (!gex.starChemistry.empty() && gex.starChemistry != "auto") {
                string prevEffective = context.effectiveChem;
                string prevInferred = context.inferredChem;
                if (gex.starChemistry == "nxt") {
                    context.effectiveChem = "NXT";
                    context.inferredChem = "NXT";
                } else if (gex.starChemistry == "tru") {
                    context.effectiveChem = "TRU";
                    context.inferredChem = "TRU";
                }
                context.inferredChemConfident = true;
                prepLog << "  GEX star_chemistry=" << gex.starChemistry
                        << " overrides namespace: effectiveChem " << prevEffective
                        << " → " << context.effectiveChem
                        << ", inferredChem " << prevInferred
                        << " → " << context.inferredChem
                        << " (confident=yes, user-specified)\n";
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

    const vector<FeatureSpec> featureSpecs = defaultPfFeatureSpecs();
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
            prepared.resolvedFastq =
                PfMultiConfig::resolveFastqDir(lib.fastqs, input.crFastqRoot, fastqMap);
            prepared.sampleName = lib.sample.empty()
                ? basenameOf(prepared.resolvedFastq)
                : lib.sample;
            prepared.assignOut =
                assignBase + "/" + sanitizeDirName(prepared.sampleName);

            const string mkAssignOutCmd = "mkdir -p \"" + prepared.assignOut + "\"";
            if (system(mkAssignOutCmd.c_str()) != 0) {
                throw runtime_error("Failed to create assign output directory: " + prepared.assignOut);
            }

            const string filteredRef = prepared.assignOut + "/feature_reference.filtered.csv";
            prepared.usedFilteredRef =
                filterFeatureRefCsv(context.featureRef, spec.featureRefType, filteredRef);
            prepared.featureRefPath = prepared.usedFilteredRef ? filteredRef : context.featureRef;

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
            prepLog << "  library " << prepared.sampleName << "/" << spec.libraryType
                    << ": star_chemistry=" << (lib.starChemistry.empty() ? "(empty)" : lib.starChemistry)
                    << ", global=" << context.requestedChem
                    << ", resolved=" << libChemRequest
                    << ", effectiveChem=" << prepared.effectiveChem
                    << (prepared.explicitChem ? " (explicit)" : " (auto-detect eligible)")
                    << "\n";

            if (!prepared.usedFilteredRef) {
                prepLog << "WARNING: feature reference not filtered for "
                        << spec.libraryType << "; using full reference\n";
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

int processPfMultiConfig(Parameters& P,
                         const Solo* solo,
                         const std::shared_ptr<PfMultiPreloadHandle>& preload) {
    if (isUnsetToken(P.pfMulti.pfMultiConfig)) {
        return 0; // Not enabled
    }
    
    P.inOut->logMain << timeMonthDayTime() << " ..... started pf-multi config processing\n";
    
    try {
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

        const PfMultiConfig::Config& config = prepared.config;
        const string& whitelist = prepared.whitelist;
        const string& effectiveChem = prepared.effectiveChem;
        const string& outputChem = prepared.outputChem;
        const string& outPrefix = prepared.outPrefix;
        const string& crAssignRoot = prepared.crAssignRoot;
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
        assignOpts.minPosterior = P.pfMulti.crAssignMinPosterior;
        assignOpts.maxReads = (P.readMapNumber > 0) ? static_cast<long long>(P.readMapNumber) : -1;
        assignOpts.legacyCbRescue = (P.pfMulti.crAssignLegacyCbRescue != 0);
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
        const int pfConsumerCapPerProducer = 8;

        if (!P.pfMulti.crAssignFilteredBarcodes.empty() && P.pfMulti.crAssignFilteredBarcodes != "-") {
            struct stat stFiltered;
            if (stat(P.pfMulti.crAssignFilteredBarcodes.c_str(), &stFiltered) != 0) {
                throw runtime_error("crAssignFilteredBarcodes file not found: " + P.pfMulti.crAssignFilteredBarcodes);
            }
            assignOpts.filteredBarcodesPath = P.pfMulti.crAssignFilteredBarcodes;
            P.inOut->logMain << "NOTICE: Using explicit filtered barcodes for assignBarcodes: "
                             << assignOpts.filteredBarcodesPath << "\n";
        } else {
            vector<string> filteredAssignBarcodes;
            if (getFilteredBarcodesFromSolo(solo, P, filteredAssignBarcodes, false)) {
                string filteredPath = crAssignRoot + "/filtered_gex_barcodes.txt";
                std::ofstream fbc(filteredPath.c_str());
                if (!fbc.is_open()) {
                    throw runtime_error("Failed to write filtered barcode list for assignBarcodes: " + filteredPath);
                }
                for (const auto& bc : filteredAssignBarcodes) {
                    fbc << bc << "\n";
                }
                fbc.close();
                assignOpts.filteredBarcodesPath = filteredPath;
                P.inOut->logMain << "NOTICE: Using Solo filtered barcodes for assignBarcodes "
                                 << "(matching namespace, " 
                                 << filteredAssignBarcodes.size() << " barcodes)\n";
            }
        }
        
        struct FeatureRun {
            string featureType;
            string assignOut;
            string featureRefPath;
            string effectiveChem;
            string detectedMatchMode;
        };
        vector<FeatureRun> featureRuns;

        for (const auto& preparedLib : prepared.featureLibraries) {
            const string& resolvedFastq = preparedLib.resolvedFastq;
            const string& sampleName = preparedLib.sampleName;
            const string& assignOut = preparedLib.assignOut;
            const string& refPath = preparedLib.featureRefPath;
            const string& libraryType = preparedLib.libraryType;
            const string& featureRefType = preparedLib.featureRefType;
            const AnchoredReadEstimate featureEstimate = preparedLib.featureEstimate;
            const int pfProducerSlotsAvailable =
                std::max(1, std::min(std::max(1, P.runThreadN - 1),
                                     std::max(1, static_cast<int>(featureEstimate.fileCount))));
            int pfReaderThreadsReserved = 1;
            int pfConsumerBudgetThreads = std::max(1, P.runThreadN - pfReaderThreadsReserved);
            int pfConsumerThreadsForRun = std::max(1, assignOpts.consumerThreadsPerSet);
            int pfConsumerSoftMax = pfConsumerCapPerProducer;
            bool pfConsumerThreadsAuto = false;
            if (!pfConsumerThreadsExplicit && assignOpts.enableStarDynamicPermitHooks) {
                // Choose producer reservation and consumer pool jointly:
                // maximize consumers, cap consumers per producer, then increase producers as needed.
                int bestProducerReserve = 1;
                int bestConsumerThreads = 1;
                for (int producerReserve = 1; producerReserve <= pfProducerSlotsAvailable; ++producerReserve) {
                    const int consumerBudget =
                        std::max(1, P.runThreadN - producerReserve);
                    const int consumerByBudget =
                        std::max(1, consumerBudget / std::max(1, pfSearchThreadQuantum));
                    const int consumerByProducerCap =
                        std::max(1, producerReserve * pfConsumerCapPerProducer);
                    const int candidateConsumers =
                        std::max(1, std::min(consumerByBudget, consumerByProducerCap));
                    if (candidateConsumers > bestConsumerThreads ||
                        (candidateConsumers == bestConsumerThreads &&
                         producerReserve < bestProducerReserve)) {
                        bestConsumerThreads = candidateConsumers;
                        bestProducerReserve = producerReserve;
                    }
                }
                pfReaderThreadsReserved = bestProducerReserve;
                pfConsumerThreadsForRun = bestConsumerThreads;
                pfConsumerBudgetThreads = std::max(1, P.runThreadN - pfReaderThreadsReserved);
                pfConsumerSoftMax =
                    std::max(1, pfReaderThreadsReserved * pfConsumerCapPerProducer);
                pfConsumerThreadsAuto = true;
            } else if (pfConsumerThreadsExplicit) {
                // Respect explicit consumer count, but scale producer reserve to maintain per-producer cap.
                const int neededProducerReserve = std::max(
                    1, (pfConsumerThreadsForRun + pfConsumerCapPerProducer - 1) / pfConsumerCapPerProducer);
                pfReaderThreadsReserved = std::min(pfProducerSlotsAvailable, neededProducerReserve);
                pfConsumerBudgetThreads = std::max(1, P.runThreadN - pfReaderThreadsReserved);
                pfConsumerSoftMax =
                    std::max(1, pfReaderThreadsReserved * pfConsumerCapPerProducer);
            }
            PfMultiAssign::AssignOptions runAssignOpts = assignOpts;
            runAssignOpts.consumerThreadsPerSet = pfConsumerThreadsForRun;

            const bool useAutodetect = !preparedLib.explicitChem
                                      && (preparedLib.resolvedChemRequest == "auto");
            const string wlNamespace = prepared.inferredChem;
            const bool wlNamespaceConfident = prepared.inferredChemConfident;
            if (useAutodetect) {
                runAssignOpts.autodetectChemistry = true;
            }
            if (preparedLib.explicitChem) {
                P.inOut->logMain << "NOTICE: " << featureRefType
                                 << " star_chemistry=" << preparedLib.resolvedChemRequest
                                 << " → effectiveChem=" << preparedLib.effectiveChem
                                 << " (auto-detect skipped)\n";
            } else if (useAutodetect && !wlNamespaceConfident) {
                P.inOut->logMain << "WARNING: Whitelist namespace uncertain ('"
                                 << wlNamespace << "', inferred by default), "
                                 << "auto-detect will determine effective chemistry for "
                                 << featureRefType << "\n";
            }
            const bool pfSingleConsumer = (pfConsumerThreadsForRun <= 1);
            const int pfPermitCeilingDuringPf =
                std::max(1, P.runThreadN - pfReaderThreadsReserved);

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
                                 << ", consumerCapPerProducer=" << pfConsumerCapPerProducer
                                 << ", producerSlotsAvailable=" << pfProducerSlotsAvailable
                                 << ", singleConsumerMode=" << (pfSingleConsumer ? "yes" : "no")
                                 << ", readerThreadsReserved=" << pfReaderThreadsReserved
                                 << ", permitCeilingDuringPf=" << pfPermitCeilingDuringPf
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

                    while (!stopPfController.load(std::memory_order_relaxed)) {
                        int nextTarget = g_threadChunks.mapPermitSnapshot().targetPermits;
                        if (pfControllerSequenceMode) {
                            nextTarget = pfControllerSequence[seqIndex % pfControllerSequence.size()];
                            ++seqIndex;
                        } else if (pfControllerEtaMode) {
                            const ThreadControl::MapPermitSnapshot snapshot = g_threadChunks.mapPermitSnapshot();
                            const uint64_t mapDone = snapshot.mapDomain.workUnitsTotal;
                            const uint64_t featureDone = snapshot.featureDomain.workUnitsTotal;
                            mapEstimateLoop = inflateEstimateUntilNonNegativeRemaining(mapEstimateLoop, mapDone);
                            featureEstimateLoop = inflateEstimateUntilNonNegativeRemaining(featureEstimateLoop, featureDone);
                            pfControllerMapEstimateFinal.store(mapEstimateLoop, std::memory_order_relaxed);
                            pfControllerFeatureEstimateFinal.store(featureEstimateLoop, std::memory_order_relaxed);

                            uint64_t mapRemaining = (mapEstimateLoop > mapDone) ? (mapEstimateLoop - mapDone) : 0;
                            uint64_t featureRemaining =
                                (featureEstimateLoop > featureDone) ? (featureEstimateLoop - featureDone) : 0;

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
                        } else if (pfControllerChunkedMode) {
                            const ThreadControl::MapPermitSnapshot snapshot = g_threadChunks.mapPermitSnapshot();
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
                        }

                        if (pfControllerSequenceMode) {
                            // Sequence mode has no remaining-work signal; always reserve producer/read threads.
                            nextTarget = std::min(nextTarget, pfPermitCeilingDuringPf);
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
                assignResult = PfMultiAssign::runAssignBarcodes(whitelist, refPath, resolvedFastq, assignOut, runAssignOpts);
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
                    apiRunOut << "pfController.consumerCapPerProducer=" << pfConsumerCapPerProducer << "\n";
                    apiRunOut << "pfController.consumerThreadsExplicit=" << (pfConsumerThreadsExplicit ? 1 : 0) << "\n";
                    apiRunOut << "pfController.consumerThreadsAuto=" << (pfConsumerThreadsAuto ? 1 : 0) << "\n";
                    apiRunOut << "pfController.consumerBudgetThreads=" << pfConsumerBudgetThreads << "\n";
                    apiRunOut << "pfController.producerSlotsAvailable=" << pfProducerSlotsAvailable << "\n";
                    apiRunOut << "pfController.singleConsumerMode=" << (pfSingleConsumer ? 1 : 0) << "\n";
                    apiRunOut << "pfController.readerThreadsReserved=" << pfReaderThreadsReserved << "\n";
                    apiRunOut << "pfController.permitCeilingDuringPf=" << pfPermitCeilingDuringPf << "\n";
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
                } else {
                    P.inOut->logMain << "WARNING: failed to append pf controller summary to "
                                     << apiRunPath << "\n";
                }
            }

            if (pfControllerEnabled) {
                P.inOut->logMain << "pf-dynamic-controller: stop mode="
                                 << pfPermitControllerModeName(pfControllerMode)
                                 << ", ticks=" << pfControllerTicks.load(std::memory_order_relaxed)
                                 << ", applied=" << pfControllerApplied.load(std::memory_order_relaxed)
                                 << ", lastTarget=" << pfControllerLastTarget.load(std::memory_order_relaxed)
                                 << ", mapEstimateFinal=" << pfControllerMapEstimateFinal.load(std::memory_order_relaxed)
                                 << ", featureEstimateFinal=" << pfControllerFeatureEstimateFinal.load(std::memory_order_relaxed)
                                 << ", libraryType=" << libraryType
                                 << ", sample=" << sampleName
                                 << "\n";
            }

            if (assignResult.returnCode != 0) {
                throw runtime_error("Failed to process feature library: " + libraryType);
            }

            FeatureRun run;
            run.featureType = featureRefType;
            run.assignOut = assignOut;
            run.featureRefPath = refPath;
            run.effectiveChem = preparedLib.effectiveChem;
            run.detectedMatchMode = assignResult.detectedMatchMode;

            if (useAutodetect && wlNamespaceConfident) {
                if (assignResult.detectedMatchMode == "RAW_MATCH") {
                    run.effectiveChem = wlNamespace;
                    P.inOut->logMain << "NOTICE: auto-detect RAW_MATCH for "
                                     << featureRefType
                                     << " → effectiveChem=" << run.effectiveChem << "\n";
                } else if (assignResult.detectedMatchMode == "TRANSLATED_MATCH") {
                    run.effectiveChem = oppositeNamespace(wlNamespace);
                    P.inOut->logMain << "NOTICE: auto-detect TRANSLATED_MATCH for "
                                     << featureRefType
                                     << " → effectiveChem=" << run.effectiveChem << "\n";
                } else {
                    P.inOut->logMain << "WARNING: auto-detect "
                                     << assignResult.detectedMatchMode << " for "
                                     << featureRefType << ", keeping inferred: "
                                     << run.effectiveChem << "\n";
                }
            } else if (useAutodetect && !wlNamespaceConfident) {
                P.inOut->logMain << "WARNING: auto-detect " << assignResult.detectedMatchMode
                                 << " for " << featureRefType
                                 << " but whitelist namespace is uncertain ('"
                                 << wlNamespace << "' by default). "
                                 << "effectiveChem kept as " << run.effectiveChem
                                 << " to stay consistent with GEX. "
                                 << "Use star_chemistry column, --crChemistry NXT/TRU, "
                                 << "or rename whitelist to enable absolute namespace resolution.\n";
            }
            featureRuns.push_back(run);
        }

        if (featureRuns.empty()) {
            throw runtime_error("No feature libraries found in multi config");
        }

        for (auto& run : featureRuns) {
            string realOut = findAssignOutputDir(run.assignOut);
            run.assignOut = realOut;
            int ret = PfMultiMexStub::processAssignOutput(run.assignOut, run.featureRefPath, run.featureType, false, whitelist);
            if (ret != 0) {
                P.inOut->logMain << "WARNING: Failed to generate MEX stub for " << run.featureType << "\n";
            }
        }
        
        // Read GEX MEX from STARsolo output (both raw and filtered if available)
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
            // Fallback to GeneFull if Gene is absent (auto mode)
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
        
        // Prefer in-memory filtered barcodes from Solo if available.
        // Fetch read-space barcodes (useOutputNamespace=false) so that the
        // namespace is unambiguous regardless of 1-col vs 2-col whitelist,
        // then normalize to TRU to match the merge-boundary convention.
        vector<string> filteredGexBarcodes;
        bool useFilteredGex = false;
        if (getFilteredBarcodesFromSolo(solo, P, filteredGexBarcodes, false)) {
            normalizeBarcodeVecToTru(filteredGexBarcodes, effectiveChem);
            useFilteredGex = true;
            P.inOut->logMain << "NOTICE: Using GEX filtered barcodes from Solo "
                             << "(read-space, normalized to TRU, in-memory)\n";
        }

        // Read raw GEX MEX (required for raw_feature_bc_matrix)
        PfMultiMerge::MexData gexRawData;
        try {
            if (hasRaw) {
                gexRawData = PfMultiMerge::readMex(rawOut);
                // Filter to Gene Expression only if needed
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
            // If raw failed but filtered exists, we'll use filtered as fallback below
            if (!hasFiltered) {
                return 0;
            }
        }
        
        // Read filtered GEX MEX if needed (for counts fallback or barcode list)
        PfMultiMerge::MexData gexFilteredData;
        bool loadedFilteredMex = false;
        bool needFilteredMexForCounts = (!hasRaw || gexRawData.features.empty());
        bool needFilteredMexForBarcodes = (!useFilteredGex && hasFiltered);
        if (hasFiltered && (needFilteredMexForCounts || needFilteredMexForBarcodes)) {
            try {
                gexFilteredData = PfMultiMerge::readMex(filteredOut);
                // Filter to Gene Expression only if needed
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
        
        // Determine primary GEX data for merging (use raw if available, otherwise filtered)
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
        
        // Read feature MEX files, pairing each with its chemistry so that
        // a failed read does not cause index misalignment.
        struct FeatureMexEntry {
            PfMultiMerge::MexData data;
            string effectiveChem;
            string featureType;
        };
        vector<FeatureMexEntry> featureMexEntries;
        for (const auto& run : featureRuns) {
            try {
                FeatureMexEntry entry;
                entry.data = PfMultiMerge::readMex(run.assignOut);
                entry.effectiveChem = run.effectiveChem;
                entry.featureType = run.featureType;
                featureMexEntries.push_back(std::move(entry));
            } catch (const exception& e) {
                P.inOut->logMain << "WARNING: Failed to read feature MEX for " << run.featureType
                                << ": " << e.what() << "\n";
            }
        }
        
        if (featureMexEntries.empty()) {
            P.inOut->logMain << "WARNING: No feature MEX files found, skipping merge\n";
            return 0;
        }
        
        // Normalize all barcode namespaces to TRU at the integration boundary.
        // Each source may be NXT or TRU depending on its whitelist; normalizing
        // here ensures cross-library barcode joins and output are in a single
        // canonical namespace (TRU) regardless of per-library chemistry.
        normalizeMexBarcodesToTru(gexData, effectiveChem);
        normalizeMexBarcodesToTru(gexRawData, effectiveChem);
        if (loadedFilteredMex) {
            normalizeMexBarcodesToTru(gexFilteredData, effectiveChem);
        }
        for (auto& entry : featureMexEntries) {
            normalizeMexBarcodesToTru(entry.data, entry.effectiveChem);
        }
        P.inOut->logMain << "Barcode namespace: all inputs normalized to TRU before merge"
                         << " (GEX effectiveChem=" << effectiveChem
                         << ", outputChem=" << outputChem << ")\n";

        // Merge MEX files (all inputs now in TRU namespace)
        vector<PfMultiMerge::MexData> featureDataVec;
        featureDataVec.reserve(featureMexEntries.size());
        for (auto& entry : featureMexEntries) {
            featureDataVec.push_back(std::move(entry.data));
        }
        PfMultiMerge::MexData mergedData = PfMultiMerge::mergeMex(gexData, featureDataVec);
        
        // Extract GEM well from GEX library (with fallback logic)
        string gemWell = "1"; // Default
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
        
        // Compute observed GEX barcodes from raw GEX triplets (already in TRU)
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
        
        // Fallback filtered barcodes from MEX data (already TRU-normalized).
        if (!useFilteredGex && hasFiltered && loadedFilteredMex && !gexFilteredData.features.empty()) {
            filteredGexBarcodes = PfMultiMerge::computeObservedGexBarcodes(gexFilteredData);
            useFilteredGex = true;
        }
        if (!useFilteredGex) {
            filteredGexBarcodes = observedRawGexBarcodes;
            P.inOut->logMain << "WARNING: Filtered GEX barcodes not available, using observed raw GEX barcodes for filtered output\n";
        }
        
        // Write raw_feature_bc_matrix — merged barcodes are in TRU;
        // outputChem controls the final namespace (default TRU, --crOutputChemistry NXT reverses).
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
        
        // Write filtered_feature_bc_matrix
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
        
        // Run CRISPR feature calling if CRISPR Guide Capture features were processed
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

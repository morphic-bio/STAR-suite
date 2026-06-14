#include "PfMultiTableImport.h"
#ifndef PF_TABLE_IMPORT_NO_PERMITS
#include "PfMultiAssign.h"
#endif
#include "PfMultiMexStub.h"
#include "MexWriter.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <cstdint>
#include <stdexcept>

namespace PfMultiTableImport {

namespace {

static constexpr size_t kPermitRowChunk = 1024;

static void trimInPlace(string& s) {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        s.clear();
        return;
    }
    size_t last = s.find_last_not_of(" \t\r\n");
    s = s.substr(first, last - first + 1);
}

static string upperCopy(string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return value;
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

static bool hasDigitDashSuffix(const string& barcode) {
    if (barcode.size() < 3 || barcode.back() < '0' || barcode.back() > '9') {
        return false;
    }
    size_t pos = barcode.size() - 1;
    while (pos > 0 && barcode[pos] >= '0' && barcode[pos] <= '9') {
        --pos;
    }
    return pos > 0 && barcode[pos] == '-';
}

static string stripGemWellSuffix(const string& barcode) {
    if (!hasDigitDashSuffix(barcode)) {
        return barcode;
    }
    size_t pos = barcode.size() - 1;
    while (pos > 0 && barcode[pos] >= '0' && barcode[pos] <= '9') {
        --pos;
    }
    if (pos > 0 && barcode[pos] == '-') {
        return barcode.substr(0, pos);
    }
    return barcode;
}

static string namespaceLabel(bool sawSuffixed, bool sawUnsuffixed) {
    if (sawSuffixed && sawUnsuffixed) {
        return "MIXED";
    }
    if (sawSuffixed) {
        return "SUFFIXED";
    }
    if (sawUnsuffixed) {
        return "UNSUFFIXED";
    }
    return "UNKNOWN";
}

static char detectDelimiterFromHeader(const string& headerLine) {
    size_t tabs = 0;
    size_t commas = 0;
    for (char c : headerLine) {
        if (c == '\t') {
            ++tabs;
        } else if (c == ',') {
            ++commas;
        }
    }
    if (tabs > 0 && tabs >= commas) {
        return '\t';
    }
    if (commas > 0) {
        return ',';
    }
    return '\t';
}

static char detectDelimiterFromPath(const string& path) {
    if (path.size() >= 4) {
        string ext = path.substr(path.size() - 4);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".csv") {
            return ',';
        }
    }
    return '\t';
}

static vector<string> splitFields(const string& line, char delim) {
    vector<string> fields;
    string field;
    bool inQuotes = false;
    for (char c : line) {
        if (c == '"') {
            inQuotes = !inQuotes;
        } else if (c == delim && !inQuotes) {
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(field);
    for (auto& value : fields) {
        trimInPlace(value);
        if (!value.empty() && value.front() == '"' && value.back() == '"') {
            value = value.substr(1, value.length() - 2);
            trimInPlace(value);
        }
    }
    return fields;
}

static bool parseNonNegativeInt(const string& text, uint64_t& out) {
    if (text.empty()) {
        return false;
    }
    uint64_t value = 0;
    for (char c : text) {
        if (c < '0' || c > '9') {
            return false;
        }
        value = value * 10 + static_cast<uint64_t>(c - '0');
    }
    out = value;
    return true;
}

static bool loadWhitelistBarcodes(const string& path,
                                  vector<string>& orderedBarcodes,
                                  std::unordered_map<string, uint32_t>& indexByBarcode) {
    orderedBarcodes.clear();
    indexByBarcode.clear();
    std::ifstream in(path.c_str());
    if (!in.is_open()) {
        return false;
    }
    string line;
    while (std::getline(in, line)) {
        trimInPlace(line);
        if (line.empty()) {
            continue;
        }
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t end = line.find_first_of("\t, \r\n", first);
        string token = (end == string::npos) ? line.substr(first) : line.substr(first, end - first);
        token = upperCopy(token);
        if (!isValidBarcodeSeq(token)) {
            continue;
        }
        if (indexByBarcode.find(token) != indexByBarcode.end()) {
            continue;
        }
        indexByBarcode[token] = static_cast<uint32_t>(orderedBarcodes.size());
        orderedBarcodes.push_back(token);
    }
    return !orderedBarcodes.empty();
}

static bool loadFilteredBarcodeSet(const string& path, std::unordered_set<string>& out) {
    out.clear();
    if (path.empty()) {
        return true;
    }
    std::ifstream in(path.c_str());
    if (!in.is_open()) {
        return false;
    }
    string line;
    while (std::getline(in, line)) {
        trimInPlace(line);
        if (line.empty()) {
            continue;
        }
        out.insert(upperCopy(line));
    }
    return true;
}

static int resolveBarcodeIndex(const string& rawBarcode,
                               const std::unordered_map<string, uint32_t>& indexByBarcode,
                               bool whitelistUsesSuffix,
                               bool& suffixNormalized) {
    suffixNormalized = false;
    string bc = upperCopy(rawBarcode);
    auto it = indexByBarcode.find(bc);
    if (it != indexByBarcode.end()) {
        return static_cast<int>(it->second);
    }

    if (hasDigitDashSuffix(bc)) {
        string stripped = stripGemWellSuffix(bc);
        it = indexByBarcode.find(stripped);
        if (it != indexByBarcode.end()) {
            suffixNormalized = true;
            return static_cast<int>(it->second);
        }
    } else if (!whitelistUsesSuffix) {
        string withSuffix = bc + "-1";
        it = indexByBarcode.find(withSuffix);
        if (it != indexByBarcode.end()) {
            suffixNormalized = true;
            return static_cast<int>(it->second);
        }
    }
    return -1;
}

static bool acquirePermitChunk(bool enabled, uint64_t& waitNsOut) {
#ifndef PF_TABLE_IMPORT_NO_PERMITS
    return PfMultiAssign::acquireFeaturePermitChunk(enabled, waitNsOut);
#else
    (void)enabled;
    waitNsOut = 0;
    return false;
#endif
}

static void releasePermitChunk(bool enabled,
                               uint64_t waitNs,
                               uint64_t workUnits,
                               uint64_t workBytes) {
#ifndef PF_TABLE_IMPORT_NO_PERMITS
    PfMultiAssign::releaseFeaturePermitChunk(enabled, waitNs, workUnits, workBytes);
#else
    (void)enabled;
    (void)waitNs;
    (void)workUnits;
    (void)workBytes;
#endif
}

static void writeApiRunSummary(const string& assignOut,
                               const TableImportStats& stats,
                               const TableImportOptions& options,
                               const string& whitelistPath,
                               const string& featureRefPath
#ifndef PF_TABLE_IMPORT_NO_PERMITS
                               ,
                               const ThreadControl::MapPermitSnapshot* permitBefore,
                               const ThreadControl::MapPermitSnapshot* permitAfter
#endif
                               ) {
    std::ofstream out((assignOut + "/table_feature_import.api_run.txt").c_str());
    if (!out.is_open()) {
        return;
    }
    out << "mode=table_feature_import\n";
    out << "whitelist=" << whitelistPath << "\n";
    out << "feature_ref=" << featureRefPath << "\n";
    out << "table_input_path=" << stats.inputPath << "\n";
    out << "table_delimiter=" << (stats.delimiter == "\t" ? "tab" : "comma") << "\n";
    out << "sample=" << options.sampleName << "\n";
    out << "feature_type=" << options.featureTypeLabel << "\n";
    out << "assignment_whitelist_namespace=" << options.assignmentWhitelistNamespace << "\n";
    out << "filtered_barcodes_path=" << options.filteredBarcodesPath << "\n";
    out << "enableStarDynamicPermitHooks=" << (options.enableStarDynamicPermitHooks ? 1 : 0) << "\n";
    out << "table_rows_read=" << stats.rowsRead << "\n";
    out << "table_rows_retained=" << stats.rowsRetained << "\n";
    out << "table_rows_rejected_barcode=" << stats.rowsRejectedBarcode << "\n";
    out << "table_rows_rejected_feature=" << stats.rowsRejectedFeature << "\n";
    out << "table_rows_rejected_count=" << stats.rowsRejectedCount << "\n";
    out << "table_duplicate_pairs_collapsed=" << stats.duplicatePairsCollapsed << "\n";
    out << "table_rows_zero_skipped=" << stats.rowsZeroSkipped << "\n";
    out << "table_rows_suffix_normalized=" << stats.rowsSuffixNormalized << "\n";
    out << "table_permit_chunks_processed=" << stats.permitChunksProcessed << "\n";
    out << "table_feature_permit_acquires=" << stats.featurePermitAcquires << "\n";
    out << "barcode_namespace_input=" << stats.barcodeNamespaceInput << "\n";
    out << "barcode_namespace_output=" << stats.barcodeNamespaceOutput << "\n";

#ifndef PF_TABLE_IMPORT_NO_PERMITS
    if (permitBefore != nullptr && permitAfter != nullptr) {
        auto delta = [](uint64_t after, uint64_t before) -> uint64_t {
            return after >= before ? (after - before) : 0;
        };
        out << "dynamicPermitDelta.feature.acquires="
            << delta(permitAfter->featureDomain.acquireCalls, permitBefore->featureDomain.acquireCalls) << "\n";
        out << "dynamicPermitDelta.feature.workUnits="
            << delta(permitAfter->featureDomain.workUnitsTotal, permitBefore->featureDomain.workUnitsTotal) << "\n";
        out << "dynamicPermitDelta.feature.workBytes="
            << delta(permitAfter->featureDomain.workBytesTotal, permitBefore->featureDomain.workBytesTotal) << "\n";
    }
#endif
}

static bool copyFeatureRefSnapshot(const string& featureRefPath, const string& assignOut) {
    std::ifstream in(featureRefPath.c_str());
    if (!in.is_open()) {
        return false;
    }
    std::ofstream out((assignOut + "/feature_reference.csv").c_str());
    if (!out.is_open()) {
        return false;
    }
    out << in.rdbuf();
    return true;
}

} // namespace

TableImportResult runTableFeatureImport(const string& whitelistNormalizedPath,
                                        const string& featureRefPath,
                                        const string& tablePath,
                                        const string& assignOut,
                                        const TableImportOptions& options) {
    TableImportResult result;
    result.stats.inputPath = tablePath;

#ifndef PF_TABLE_IMPORT_NO_PERMITS
    if (options.enableStarDynamicPermitHooks) {
        PfMultiAssign::waitForFeaturePermitInterface(true);
    }
    ThreadControl::MapPermitSnapshot permitBefore{};
    ThreadControl::MapPermitSnapshot permitAfter{};
    const bool capturePermitDelta =
        PfMultiAssign::featurePermitTelemetryEnabled(options.enableStarDynamicPermitHooks);
    if (capturePermitDelta) {
        permitBefore = PfMultiAssign::featurePermitSnapshot();
    }
#else
    const bool capturePermitDelta = false;
#endif

    vector<PfMultiMexStub::FeatureRow> featureRows;
    try {
        featureRows = PfMultiMexStub::loadFeatureCsv(featureRefPath);
    } catch (const std::exception& e) {
        throw std::runtime_error(string("table import: failed to load feature ref: ") + e.what());
    }
    if (featureRows.empty()) {
        throw std::runtime_error("table import: feature reference is empty: " + featureRefPath);
    }

    std::unordered_map<string, uint32_t> featureIndexById;
    vector<MexWriter::Feature> mexFeatures;
    mexFeatures.reserve(featureRows.size());
    for (size_t i = 0; i < featureRows.size(); ++i) {
        const string id = featureRows[i].id.empty() ? featureRows[i].name : featureRows[i].id;
        if (id.empty()) {
            continue;
        }
        featureIndexById[id] = static_cast<uint32_t>(mexFeatures.size());
        string ftype = options.featureTypeLabel;
        if (ftype.empty() && !featureRows[i].featureType.empty()) {
            ftype = featureRows[i].featureType;
        }
        if (ftype.empty()) {
            ftype = "Custom";
        }
        string name = featureRows[i].name.empty() ? id : featureRows[i].name;
        mexFeatures.emplace_back(id, name, ftype);
    }

    vector<string> orderedBarcodes;
    std::unordered_map<string, uint32_t> barcodeIndexBySeq;
    if (!loadWhitelistBarcodes(whitelistNormalizedPath, orderedBarcodes, barcodeIndexBySeq)) {
        throw std::runtime_error("table import: failed to load whitelist: " + whitelistNormalizedPath);
    }

    bool whitelistUsesSuffix = false;
    for (const auto& bc : orderedBarcodes) {
        if (hasDigitDashSuffix(bc)) {
            whitelistUsesSuffix = true;
            break;
        }
    }
    result.stats.barcodeNamespaceInput = "UNKNOWN";
    result.stats.barcodeNamespaceOutput = whitelistUsesSuffix ? "SUFFIXED" : "UNSUFFIXED";

    std::unordered_set<string> filteredSet;
    if (!options.filteredBarcodesPath.empty()) {
        if (!loadFilteredBarcodeSet(options.filteredBarcodesPath, filteredSet)) {
            throw std::runtime_error("table import: failed to load filtered barcodes: "
                + options.filteredBarcodesPath);
        }
    }

    std::ifstream tableIn(tablePath.c_str());
    if (!tableIn.is_open()) {
        throw std::runtime_error("table import: failed to open table: " + tablePath);
    }

    string headerLine;
    if (!std::getline(tableIn, headerLine)) {
        throw std::runtime_error("table import: table is empty: " + tablePath);
    }
    trimInPlace(headerLine);
    char delim = detectDelimiterFromHeader(headerLine);
    if (delim == '\t' && headerLine.find('\t') == string::npos && headerLine.find(',') != string::npos) {
        delim = ',';
    }
    if (headerLine.find('\t') == string::npos && headerLine.find(',') == string::npos) {
        delim = detectDelimiterFromPath(tablePath);
    }
    result.stats.delimiter = (delim == '\t') ? string("\t") : string(",");

    vector<string> headers = splitFields(headerLine, delim);
    int barcodeCol = -1;
    int featureCol = -1;
    int countCol = -1;
    for (size_t i = 0; i < headers.size(); ++i) {
        string h = headers[i];
        std::transform(h.begin(), h.end(), h.begin(), ::tolower);
        if (h == "barcode") {
            barcodeCol = static_cast<int>(i);
        } else if (h == "feature_id" || h == "featureid") {
            featureCol = static_cast<int>(i);
        } else if (h == "count") {
            countCol = static_cast<int>(i);
        }
    }
    if (barcodeCol < 0 || featureCol < 0 || countCol < 0) {
        throw std::runtime_error("table import: required columns barcode, feature_id, count missing in "
            + tablePath);
    }

    std::unordered_map<uint64_t, uint64_t> aggregated;
    aggregated.reserve(4096);
    auto pairKey = [](uint32_t bcIdx, uint32_t featIdx) -> uint64_t {
        return (static_cast<uint64_t>(bcIdx) << 32) | static_cast<uint64_t>(featIdx);
    };

    string line;
    size_t rowsSincePermit = 0;
    uint64_t bytesSincePermit = 0;
    uint64_t permitWaitNs = 0;
    bool permitChunkActive = false;
    bool sawSuffixedTableBarcode = false;
    bool sawUnsuffixedTableBarcode = false;
    while (std::getline(tableIn, line)) {
        bytesSincePermit += line.size() + 1;
        if (line.empty()) {
            continue;
        }
        trimInPlace(line);
        if (line.empty()) {
            continue;
        }

        if (rowsSincePermit == 0) {
            permitChunkActive =
                acquirePermitChunk(options.enableStarDynamicPermitHooks, permitWaitNs);
            if (permitChunkActive) {
                result.stats.featurePermitAcquires++;
            }
        }
        ++rowsSincePermit;
        result.stats.rowsRead++;

        vector<string> fields = splitFields(line, delim);
        if (static_cast<int>(fields.size()) <= std::max({barcodeCol, featureCol, countCol})) {
            result.stats.rowsRejectedCount++;
            continue;
        }

        const string tableBarcode = upperCopy(fields[barcodeCol]);
        if (!tableBarcode.empty()) {
            if (hasDigitDashSuffix(tableBarcode)) {
                sawSuffixedTableBarcode = true;
            } else {
                sawUnsuffixedTableBarcode = true;
            }
        }

        uint64_t count = 0;
        if (!parseNonNegativeInt(fields[countCol], count)) {
            result.stats.rowsRejectedCount++;
            continue;
        }
        if (count == 0) {
            result.stats.rowsZeroSkipped++;
            continue;
        }

        const string& featureId = fields[featureCol];
        auto featIt = featureIndexById.find(featureId);
        if (featIt == featureIndexById.end()) {
            result.stats.rowsRejectedFeature++;
            continue;
        }

        bool suffixNormalized = false;
        const int bcIdx = resolveBarcodeIndex(
            tableBarcode, barcodeIndexBySeq, whitelistUsesSuffix, suffixNormalized);
        if (bcIdx < 0) {
            result.stats.rowsRejectedBarcode++;
            continue;
        }
        if (suffixNormalized) {
            result.stats.rowsSuffixNormalized++;
        }

        if (!filteredSet.empty()) {
            const string& canonBc = orderedBarcodes[static_cast<size_t>(bcIdx)];
            if (filteredSet.find(canonBc) == filteredSet.end()) {
                result.stats.rowsRejectedBarcode++;
                continue;
            }
        }

        const uint64_t key = pairKey(static_cast<uint32_t>(bcIdx), featIt->second);
        auto aggIt = aggregated.find(key);
        if (aggIt == aggregated.end()) {
            aggregated[key] = count;
        } else {
            aggIt->second += count;
            result.stats.duplicatePairsCollapsed++;
        }

        if (rowsSincePermit >= kPermitRowChunk || bytesSincePermit >= 65536) {
            if (permitChunkActive) {
                releasePermitChunk(options.enableStarDynamicPermitHooks,
                                   permitWaitNs, rowsSincePermit, bytesSincePermit);
                result.stats.permitChunksProcessed++;
            }
            rowsSincePermit = 0;
            bytesSincePermit = 0;
            permitWaitNs = 0;
            permitChunkActive = false;
        }
    }
    if (rowsSincePermit > 0 && permitChunkActive) {
        releasePermitChunk(options.enableStarDynamicPermitHooks,
                           permitWaitNs, rowsSincePermit, bytesSincePermit);
        result.stats.permitChunksProcessed++;
    }

    result.stats.barcodeNamespaceInput =
        namespaceLabel(sawSuffixedTableBarcode, sawUnsuffixedTableBarcode);
    result.stats.rowsRetained = aggregated.size();

    vector<MexWriter::Triplet> triplets;
    triplets.reserve(aggregated.size());
    for (const auto& entry : aggregated) {
        const uint32_t bcIdx = static_cast<uint32_t>(entry.first >> 32);
        const uint32_t featIdx = static_cast<uint32_t>(entry.first & 0xffffffffu);
        if (entry.second == 0) {
            continue;
        }
        if (entry.second > std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("table import: count overflow for barcode/feature pair");
        }
        triplets.push_back({bcIdx, featIdx,
                            static_cast<uint32_t>(entry.second)});
    }

    vector<uint32_t> usedBarcodeIdx;
    usedBarcodeIdx.reserve(orderedBarcodes.size());
    std::unordered_set<uint32_t> seenBc;
    for (const auto& trip : triplets) {
        if (seenBc.insert(trip.cell_idx).second) {
            usedBarcodeIdx.push_back(trip.cell_idx);
        }
    }
    std::sort(usedBarcodeIdx.begin(), usedBarcodeIdx.end());

    vector<string> outputBarcodes;
    vector<uint32_t> oldToNew(orderedBarcodes.size(), UINT32_MAX);
    outputBarcodes.reserve(usedBarcodeIdx.size());
    for (size_t i = 0; i < usedBarcodeIdx.size(); ++i) {
        const uint32_t oldIdx = usedBarcodeIdx[i];
        oldToNew[oldIdx] = static_cast<uint32_t>(i);
        outputBarcodes.push_back(orderedBarcodes[oldIdx]);
    }

    vector<MexWriter::Triplet> remapped;
    remapped.reserve(triplets.size());
    for (const auto& trip : triplets) {
        const uint32_t newIdx = oldToNew[trip.cell_idx];
        if (newIdx == UINT32_MAX) {
            continue;
        }
        remapped.push_back({newIdx, trip.gene_idx, trip.count});
    }

    if (outputBarcodes.empty() || mexFeatures.empty()) {
        throw std::runtime_error("table import: no retained barcode/feature pairs in " + tablePath);
    }

    const string mexPrefix = assignOut + "/";
    if (MexWriter::writeMex(mexPrefix, outputBarcodes, mexFeatures, remapped, -1) != 0) {
        throw std::runtime_error("table import: failed to write MEX under " + assignOut);
    }
    copyFeatureRefSnapshot(featureRefPath, assignOut);

#ifndef PF_TABLE_IMPORT_NO_PERMITS
    if (capturePermitDelta) {
        permitAfter = PfMultiAssign::featurePermitSnapshot();
    }
    writeApiRunSummary(assignOut, result.stats, options, whitelistNormalizedPath, featureRefPath,
                       capturePermitDelta ? &permitBefore : nullptr,
                       capturePermitDelta ? &permitAfter : nullptr);
#else
    writeApiRunSummary(assignOut, result.stats, options, whitelistNormalizedPath, featureRefPath);
#endif

    result.returnCode = 0;
    return result;
}

} // namespace PfMultiTableImport

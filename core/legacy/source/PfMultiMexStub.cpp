#include "PfMultiMexStub.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <stdexcept>
#include <dirent.h>
#include <unordered_map>
#include <unordered_set>
#include <cstdio>
using std::cerr;
using std::endl;

namespace PfMultiMexStub {

static void trimInPlace(string& s) {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        s.clear();
        return;
    }
    size_t last = s.find_last_not_of(" \t\r\n");
    s = s.substr(first, last - first + 1);
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

static bool parseWhitelistOutputMap(const string& whitelistPath, std::unordered_map<string, string>& outMap) {
    outMap.clear();
    if (whitelistPath.empty()) {
        return false;
    }

    ifstream in(whitelistPath.c_str());
    if (!in.is_open()) {
        return false;
    }

    string line;
    bool sawSecondColumn = false;
    while (getline(in, line)) {
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t end1 = line.find_first_of(" \t,\r\n", first);
        string col1 = (end1 == string::npos) ? line.substr(first) : line.substr(first, end1 - first);
        size_t second = (end1 == string::npos) ? string::npos : line.find_first_not_of(" \t,\r\n", end1);
        if (second == string::npos) {
            continue;
        }
        size_t end2 = line.find_first_of(" \t,\r\n", second);
        string col2 = (end2 == string::npos) ? line.substr(second) : line.substr(second, end2 - second);
        sawSecondColumn = true;

        std::transform(col1.begin(), col1.end(), col1.begin(), ::toupper);
        std::transform(col2.begin(), col2.end(), col2.begin(), ::toupper);
        if (!isValidBarcodeSeq(col1) || !isValidBarcodeSeq(col2)) {
            continue;
        }
        if (col1.size() != col2.size()) {
            continue;
        }
        outMap[col1] = col2;
    }

    return sawSecondColumn && !outMap.empty();
}

vector<FeatureRow> loadFeatureCsv(const string& csvPath) {
    vector<FeatureRow> result;
    ifstream file(csvPath);
    if (!file.is_open()) {
        ostringstream err;
        err << "Failed to open feature CSV: " << csvPath;
        throw runtime_error(err.str());
    }
    
    // Read header
    string headerLine;
    if (!getline(file, headerLine)) {
        throw runtime_error("Feature CSV is empty");
    }
    
    // Parse header
    vector<string> headers;
    istringstream headerStream(headerLine);
    string headerField;
    while (getline(headerStream, headerField, ',')) {
        trimInPlace(headerField);
        std::transform(headerField.begin(), headerField.end(), headerField.begin(), ::tolower);
        headers.push_back(headerField);
    }
    
    // Find column indices
    int nameIdx = -1, idIdx = -1, ftypeIdx = -1;
    for (size_t i = 0; i < headers.size(); ++i) {
        if (headers[i] == "name") {
            nameIdx = i;
        } else if (headers[i] == "id") {
            idIdx = i;
        } else if (headers[i] == "feature_type" || headers[i] == "type") {
            ftypeIdx = i;
        }
    }
    
    if (nameIdx == -1 && idIdx == -1) {
        throw runtime_error("Feature CSV header must include 'name' or 'id'");
    }
    
    // Read data rows (first-definition-wins on duplicate feature names; matches io.c)
    string line;
    std::unordered_set<string> seenNames;
    int lineNumber = 1;
    while (getline(file, line)) {
        ++lineNumber;
        if (line.empty()) continue;
        
        vector<string> fields;
        string field;
        bool inQuotes = false;
        
        // Simple CSV parsing
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
        if (!field.empty() || line.back() == ',') {
            fields.push_back(field);
        }
        
        FeatureRow row;
        if (nameIdx >= 0 && (size_t)nameIdx < fields.size()) {
            row.name = fields[nameIdx];
            // Remove quotes if present
            if (!row.name.empty() && row.name.front() == '"' && row.name.back() == '"') {
                row.name = row.name.substr(1, row.name.length() - 2);
            }
            trimInPlace(row.name);
        }
        if (idIdx >= 0 && (size_t)idIdx < fields.size()) {
            row.id = fields[idIdx];
            if (!row.id.empty() && row.id.front() == '"' && row.id.back() == '"') {
                row.id = row.id.substr(1, row.id.length() - 2);
            }
            trimInPlace(row.id);
        }
        if (ftypeIdx >= 0 && (size_t)ftypeIdx < fields.size()) {
            row.featureType = fields[ftypeIdx];
            if (!row.featureType.empty() && row.featureType.front() == '"' && row.featureType.back() == '"') {
                row.featureType = row.featureType.substr(1, row.featureType.length() - 2);
            }
            trimInPlace(row.featureType);
        }
        
        // Use name as id if id is missing, and vice versa
        if (row.id.empty() && !row.name.empty()) {
            row.id = row.name;
        }
        if (row.name.empty() && !row.id.empty()) {
            row.name = row.id;
        }
        
        const string dedupeKey = !row.name.empty() ? row.name : row.id;
        if (dedupeKey.empty()) {
            continue;
        }
        if (seenNames.count(dedupeKey)) {
            fprintf(stderr,
                    "Warning: duplicate feature name '%s' in %s at line %d; "
                    "ignoring later definition and keeping the first\n",
                    dedupeKey.c_str(), csvPath.c_str(), lineNumber);
            continue;
        }
        seenNames.insert(dedupeKey);
        result.push_back(row);
    }
    
    if (result.empty()) {
        throw runtime_error("Feature CSV has no data rows");
    }
    
    return result;
}

vector<string> readFeaturesTxt(const string& txtPath) {
    vector<string> result;
    ifstream file(txtPath);
    if (!file.is_open()) {
        return result; // File doesn't exist, return empty
    }
    
    string line;
    while (getline(file, line)) {
        trimInPlace(line);
        if (!line.empty()) {
            result.push_back(line);
        }
    }
    return result;
}

string compareFeatureNames(const vector<FeatureRow>& featureRows, 
                          const vector<string>& featuresTxt) {
    if (featureRows.size() != featuresTxt.size()) {
        ostringstream err;
        err << "feature count mismatch: csv=" << featureRows.size() 
            << " features.txt=" << featuresTxt.size();
        return err.str();
    }
    
    for (size_t i = 0; i < featureRows.size(); ++i) {
        if (featureRows[i].name != featuresTxt[i]) {
            ostringstream err;
            err << "name mismatch at row " << (i + 1) << ": csv='" << featureRows[i].name 
                << "' features.txt='" << featuresTxt[i] << "'";
            return err.str();
        }
    }
    
    return "";
}

bool writeFeaturesTsv(const string& outPath, const vector<FeatureRow>& featureRows,
                     const string& defaultType, bool force,
                     const string& featureTypeOverride) {
    struct stat st;
    if (stat(outPath.c_str(), &st) == 0 && !force) {
        return false;
    }
    
    ofstream out(outPath);
    if (!out.is_open()) {
        ostringstream err;
        err << "Failed to write features.tsv: " << outPath;
        throw runtime_error(err.str());
    }
    
    for (const auto& row : featureRows) {
        string id = row.id.empty() ? row.name : row.id;
        string name = row.name.empty() ? row.id : row.name;
        string ftype;
        if (!featureTypeOverride.empty()) {
            ftype = featureTypeOverride;
        } else if (!row.featureType.empty()) {
            ftype = row.featureType;
        } else {
            ftype = defaultType;
        }
        out << id << "\t" << name << "\t" << ftype << "\n";
    }
    
    return true;
}

bool copyBarcodesTsv(const string& barcodesTxt, const string& barcodesTsv, bool force,
                     const string& whitelistPath) {
    struct stat st;
    if (stat(barcodesTxt.c_str(), &st) != 0) {
        return false; // Source file doesn't exist
    }
    
    if (stat(barcodesTsv.c_str(), &st) == 0 && !force) {
        return false; // Destination exists and not forcing
    }
    
    ifstream src(barcodesTxt);
    if (!src.is_open()) {
        return false;
    }
    
    ofstream dst(barcodesTsv);
    if (!dst.is_open()) {
        return false;
    }

    std::unordered_map<string, string> outputMap;
    bool useOutputMap = parseWhitelistOutputMap(whitelistPath, outputMap);

    if (!useOutputMap) {
        dst << src.rdbuf();
        return true;
    }

    string line;
    while (getline(src, line)) {
        trimInPlace(line);
        if (line.empty()) {
            continue;
        }
        auto it = outputMap.find(line);
        if (it != outputMap.end()) {
            dst << it->second << "\n";
        } else {
            dst << line << "\n";
        }
    }
    return true;
}

static bool remapFeaturePerCellCsv(const string& featurePerCellCsv,
                                   const string& whitelistPath,
                                   bool force) {
    struct stat st;
    if (stat(featurePerCellCsv.c_str(), &st) != 0) {
        return false;
    }

    std::unordered_map<string, string> outputMap;
    bool useOutputMap = parseWhitelistOutputMap(whitelistPath, outputMap);
    if (!useOutputMap) {
        return false;
    }

    std::unordered_set<string> outputValues;
    outputValues.reserve(outputMap.size());
    for (const auto& item : outputMap) {
        outputValues.insert(item.second);
    }

    ifstream in(featurePerCellCsv.c_str());
    if (!in.is_open()) {
        return false;
    }

    vector<string> lines;
    lines.reserve(1024);
    string line;
    size_t assignmentHits = 0;
    size_t outputHits = 0;
    while (getline(in, line)) {
        if (!lines.empty() || line.rfind("barcode,", 0) != 0) {
            string barcode = line.substr(0, line.find(','));
            trimInPlace(barcode);
            std::transform(barcode.begin(), barcode.end(), barcode.begin(), ::toupper);
            if (outputMap.find(barcode) != outputMap.end()) {
                assignmentHits++;
            }
            if (outputValues.find(barcode) != outputValues.end()) {
                outputHits++;
            }
        }
        lines.push_back(line);
    }
    in.close();

    if (assignmentHits == 0 || outputHits >= assignmentHits) {
        return false;
    }
    if (!force && outputHits > 0) {
        return false;
    }

    const string tmpPath = featurePerCellCsv + ".tmp";
    ofstream out(tmpPath.c_str());
    if (!out.is_open()) {
        return false;
    }

    for (const auto& row : lines) {
        if (row.rfind("barcode,", 0) == 0 || row.empty()) {
            out << row << "\n";
            continue;
        }
        size_t commaPos = row.find(',');
        string barcode = commaPos == string::npos ? row : row.substr(0, commaPos);
        string suffix = commaPos == string::npos ? "" : row.substr(commaPos);
        trimInPlace(barcode);
        std::transform(barcode.begin(), barcode.end(), barcode.begin(), ::toupper);
        auto it = outputMap.find(barcode);
        if (it != outputMap.end()) {
            out << it->second << suffix << "\n";
        } else {
            out << barcode << suffix << "\n";
        }
    }
    out.close();

    if (std::rename(tmpPath.c_str(), featurePerCellCsv.c_str()) != 0) {
        std::remove(tmpPath.c_str());
        return false;
    }
    return true;
}

int processAssignOutput(const string& assignOutDir, const string& featureCsvPath,
                       const string& defaultFeatureType, bool force,
                       const string& whitelistPath,
                       const string& featureTypeOverride) {
    vector<string> outDirs;
    outDirs.push_back(assignOutDir);
    
    string filteredDir = assignOutDir + "/filtered";
    struct stat st;
    if (stat(filteredDir.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        outDirs.push_back(filteredDir);
    }
    
    vector<FeatureRow> featureRows;
    try {
        featureRows = loadFeatureCsv(featureCsvPath);
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
    
    vector<string> warnings;
    bool wroteAny = false;
    
    for (const auto& outDir : outDirs) {
        string featuresTxt = outDir + "/features.txt";
        string barcodesTxt = outDir + "/barcodes.txt";
        string barcodesTsv = outDir + "/barcodes.tsv";
        string featuresTsv = outDir + "/features.tsv";
        
        vector<string> featuresTxtRows = readFeaturesTxt(featuresTxt);
        if (!featuresTxtRows.empty()) {
            string warn = compareFeatureNames(featureRows, featuresTxtRows);
            if (!warn.empty()) {
                warnings.push_back(outDir + ": " + warn);
            }
        } else {
            warnings.push_back(outDir + ": features.txt not found");
        }
        
        try {
            if (writeFeaturesTsv(featuresTsv, featureRows, defaultFeatureType, force,
                                 featureTypeOverride)) {
                wroteAny = true;
            }
            if (copyBarcodesTsv(barcodesTxt, barcodesTsv, force, whitelistPath)) {
                wroteAny = true;
            }
            if (remapFeaturePerCellCsv(outDir + "/feature_per_cell.csv", whitelistPath, force)) {
                wroteAny = true;
            }
        } catch (const exception& e) {
            cerr << "ERROR: " << e.what() << endl;
            return 1;
        }
    }
    
    for (const auto& warn : warnings) {
        cerr << "WARNING: " << warn << endl;
    }
    
    if (!wroteAny) {
        if (!force) {
            bool allExist = true;
            for (const auto& outDir : outDirs) {
                struct stat fst;
                if (stat((outDir + "/features.tsv").c_str(), &fst) != 0 ||
                    stat((outDir + "/barcodes.tsv").c_str(), &fst) != 0) {
                    allExist = false;
                    break;
                }
            }
            if (allExist) {
                return 0;
            }
        }
        cerr << "No outputs written and expected files missing." << endl;
        return 1;
    }
    
    return 0;
}

} // namespace PfMultiMexStub

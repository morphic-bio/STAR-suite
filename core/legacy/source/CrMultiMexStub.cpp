#include "CrMultiMexStub.h"
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
using std::cerr;
using std::endl;

namespace CrMultiMexStub {

static void trimInPlace(string& s) {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        s.clear();
        return;
    }
    size_t last = s.find_last_not_of(" \t\r\n");
    s = s.substr(first, last - first + 1);
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
    
    // Read data rows
    string line;
    while (getline(file, line)) {
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
        
        if (!row.id.empty() || !row.name.empty()) {
            result.push_back(row);
        }
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
                     const string& defaultType, bool force) {
    struct stat st;
    if (stat(outPath.c_str(), &st) == 0 && !force) {
        return false; // File exists and not forcing
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
        string ftype = row.featureType.empty() ? defaultType : row.featureType;
        out << id << "\t" << name << "\t" << ftype << "\n";
    }
    
    return true;
}

bool copyBarcodesTsv(const string& barcodesTxt, const string& barcodesTsv, bool force) {
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
    
    dst << src.rdbuf();
    return true;
}

int processAssignOutput(const string& assignOutDir, const string& featureCsvPath,
                       const string& defaultFeatureType, bool force) {
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
            if (writeFeaturesTsv(featuresTsv, featureRows, defaultFeatureType, force)) {
                wroteAny = true;
            }
            if (copyBarcodesTsv(barcodesTxt, barcodesTsv, force)) {
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
        cerr << "No outputs written (files may already exist)." << endl;
        return 1;
    }
    
    return 0;
}

} // namespace CrMultiMexStub

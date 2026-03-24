#include "CrMultiConfig.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <stdexcept>

namespace CrMultiConfig {

static void trimInPlace(string& s) {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        s.clear();
        return;
    }
    size_t last = s.find_last_not_of(" \t\r\n");
    s = s.substr(first, last - first + 1);
}

string LibraryEntry::normalizedFeatureType() const {
    string normalized;
    normalized.reserve(feature_types.size());
    for (unsigned char c : feature_types) {
        if (std::isalnum(c)) {
            normalized.push_back(static_cast<char>(std::tolower(c)));
        }
    }
    return normalized;
}

vector<LibraryEntry> Config::getGexLibraries() const {
    vector<LibraryEntry> result;
    for (const auto& lib : libraries) {
        string norm = lib.normalizedFeatureType();
        if (norm.find("geneexpression") != string::npos || norm.find("gex") != string::npos) {
            result.push_back(lib);
        }
    }
    return result;
}

vector<LibraryEntry> Config::getFeatureLibraries(const string& featureType) const {
    vector<LibraryEntry> result;
    string targetNorm = featureType;
    std::transform(targetNorm.begin(), targetNorm.end(), targetNorm.begin(), 
                  [](unsigned char c) { return std::tolower(c); });
    targetNorm.erase(std::remove_if(targetNorm.begin(), targetNorm.end(), ::isspace), targetNorm.end());
    
    for (const auto& lib : libraries) {
        string norm = lib.normalizedFeatureType();
        if (norm.find(targetNorm) != string::npos) {
            result.push_back(lib);
        }
    }
    return result;
}

Config parseConfig(const string& configPath) {
    Config config;
    ifstream file(configPath);
    if (!file.is_open()) {
        ostringstream err;
        err << "Failed to open multi config file: " << configPath;
        throw runtime_error(err.str());
    }
    
    string currentSection;
    vector<string> currentLines;
    bool inLibraries = false;
    vector<string> librariesHeader;
    
    string line;
    while (getline(file, line)) {
        // Remove comments and whitespace
        size_t commentPos = line.find_first_of("#;");
        if (commentPos != string::npos) {
            line = line.substr(0, commentPos);
        }
        
        // Trim whitespace safely
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t last = line.find_last_not_of(" \t\r\n");
        line = line.substr(first, last - first + 1);
        
        if (line.empty()) {
            continue;
        }
        
        // Check for section header
        if (line.front() == '[' && line.back() == ']') {
            // Process previous section
            if (currentSection == "libraries" && inLibraries) {
                // Parse CSV
                if (librariesHeader.empty()) {
                    ostringstream err;
                    err << "libraries section has no header";
                    throw runtime_error(err.str());
                }
                
                for (const auto& row : currentLines) {
                    if (row.empty()) continue;
                    istringstream ss(row);
                    vector<string> fields;
                    string field;
                    bool inQuotes = false;
                    
                    // Simple CSV parsing (handles quoted fields)
                    for (char c : row) {
                        if (c == '"') {
                            inQuotes = !inQuotes;
                        } else if (c == ',' && !inQuotes) {
                            fields.push_back(field);
                            field.clear();
                        } else {
                            field += c;
                        }
                    }
                    if (!field.empty() || row.back() == ',') {
                        fields.push_back(field);
                    }
                    
                    if (fields.size() < librariesHeader.size()) {
                        continue; // Skip incomplete rows
                    }
                    
                    LibraryEntry entry;
                    entry.gem_well = "1"; // Default to "1" if not specified
                    for (size_t i = 0; i < librariesHeader.size() && i < fields.size(); ++i) {
                        string header = librariesHeader[i];
                        std::transform(header.begin(), header.end(), header.begin(), ::tolower);
                        string value = fields[i];
                        // Remove quotes if present
                        if (!value.empty() && value.front() == '"' && value.back() == '"') {
                            value = value.substr(1, value.length() - 2);
                        }
                        
                        if (header == "fastqs") {
                            entry.fastqs = value;
                        } else if (header == "feature_types" || header == "feature_type") {
                            entry.feature_types = value;
                        } else if (header == "sample") {
                            entry.sample = value;
                        } else if (header == "library_type" || header == "librarytype") {
                            entry.library_type = value;
                        } else if (header == "gem_well") {
                            entry.gem_well = value.empty() ? "1" : value;
                        }
                    }
                    if (!entry.fastqs.empty()) {
                        config.libraries.push_back(entry);
                    }
                }
            } else if (currentSection == "feature") {
                // Parse key=value pairs
                for (const auto& kvLine : currentLines) {
                    size_t commaPos = kvLine.find(',');
                    if (commaPos == string::npos) continue;
                    string key = kvLine.substr(0, commaPos);
                    string value = kvLine.substr(commaPos + 1);
                    // Trim
                    trimInPlace(key);
                    trimInPlace(value);
                    
                    string keyLower = key;
                    std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
                    if (keyLower == "ref" || keyLower == "reference") {
                        config.featureRef = value;
                    }
                }
            } else if (currentSection == "reference") {
                // Parse key=value pairs
                for (const auto& kvLine : currentLines) {
                    size_t commaPos = kvLine.find(',');
                    if (commaPos == string::npos) continue;
                    string key = kvLine.substr(0, commaPos);
                    string value = kvLine.substr(commaPos + 1);
                    // Trim
                    trimInPlace(key);
                    trimInPlace(value);
                    
                    string keyLower = key;
                    std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
                    if (keyLower == "path" || keyLower == "genome") {
                        config.referencePath = value;
                    }
                }
            }
            
            // Start new section
            currentSection = line.substr(1, line.length() - 2);
            std::transform(currentSection.begin(), currentSection.end(), currentSection.begin(), ::tolower);
            currentLines.clear();
            inLibraries = false;
            
            if (currentSection == "libraries") {
                inLibraries = true;
                // Read header
                if (!getline(file, line)) {
                    break;
                }
                // Remove comments
                size_t commentPos = line.find_first_of("#;");
                if (commentPos != string::npos) {
                    line = line.substr(0, commentPos);
                }
                // Parse header
                istringstream headerStream(line);
                string headerField;
                while (getline(headerStream, headerField, ',')) {
                    headerField.erase(0, headerField.find_first_not_of(" \t"));
                    headerField.erase(headerField.find_last_not_of(" \t") + 1);
                    std::transform(headerField.begin(), headerField.end(), headerField.begin(), ::tolower);
                    librariesHeader.push_back(headerField);
                }
            }
        } else {
            currentLines.push_back(line);
        }
    }
    
    // Process final section
    if (currentSection == "libraries" && inLibraries) {
        // Same parsing as above
        if (!librariesHeader.empty()) {
            for (const auto& row : currentLines) {
                if (row.empty()) continue;
                istringstream ss(row);
                vector<string> fields;
                string field;
                bool inQuotes = false;
                
                for (char c : row) {
                    if (c == '"') {
                        inQuotes = !inQuotes;
                    } else if (c == ',' && !inQuotes) {
                        fields.push_back(field);
                        field.clear();
                    } else {
                        field += c;
                    }
                }
                if (!field.empty() || row.back() == ',') {
                    fields.push_back(field);
                }
                
                if (fields.size() < librariesHeader.size()) {
                    continue;
                }
                
                LibraryEntry entry;
                entry.gem_well = "1"; // Default to "1" if not specified
                for (size_t i = 0; i < librariesHeader.size() && i < fields.size(); ++i) {
                    string header = librariesHeader[i];
                    std::transform(header.begin(), header.end(), header.begin(), ::tolower);
                    string value = fields[i];
                    if (!value.empty() && value.front() == '"' && value.back() == '"') {
                        value = value.substr(1, value.length() - 2);
                    }
                    
                    if (header == "fastqs") {
                        entry.fastqs = value;
                    } else if (header == "feature_types" || header == "feature_type") {
                        entry.feature_types = value;
                    } else if (header == "sample") {
                        entry.sample = value;
                    } else if (header == "library_type" || header == "librarytype") {
                        entry.library_type = value;
                    } else if (header == "gem_well") {
                        entry.gem_well = value.empty() ? "1" : value;
                    }
                }
                if (!entry.fastqs.empty()) {
                    config.libraries.push_back(entry);
                }
            }
        }
    } else if (currentSection == "feature") {
        for (const auto& kvLine : currentLines) {
            size_t commaPos = kvLine.find(',');
            if (commaPos == string::npos) continue;
            string key = kvLine.substr(0, commaPos);
            string value = kvLine.substr(commaPos + 1);
            trimInPlace(key);
            trimInPlace(value);
            
            string keyLower = key;
            std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
            if (keyLower == "ref" || keyLower == "reference") {
                config.featureRef = value;
            }
        }
    } else if (currentSection == "reference") {
        for (const auto& kvLine : currentLines) {
            size_t commaPos = kvLine.find(',');
            if (commaPos == string::npos) continue;
            string key = kvLine.substr(0, commaPos);
            string value = kvLine.substr(commaPos + 1);
            trimInPlace(key);
            trimInPlace(value);
            
            string keyLower = key;
            std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
            if (keyLower == "path" || keyLower == "genome") {
                config.referencePath = value;
            }
        }
    }
    
    return config;
}

string resolveFastqDir(const string& configPath, const string& fastqRoot, 
                       const map<string, string>& fastqMap) {
    // Check fastq_map first
    auto it = fastqMap.find(configPath);
    if (it != fastqMap.end()) {
        return it->second;
    }
    
    // Check if path exists as-is
    struct stat st;
    if (stat(configPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        return configPath;
    }
    
    // Try fastq_root + basename
    if (!fastqRoot.empty()) {
        size_t lastSlash = configPath.find_last_of("/\\");
        string basename = (lastSlash == string::npos) ? configPath : configPath.substr(lastSlash + 1);
        // Remove trailing slash if present
        while (!basename.empty() && (basename.back() == '/' || basename.back() == '\\')) {
            basename.pop_back();
        }
        string candidate = fastqRoot;
        if (candidate.back() != '/' && candidate.back() != '\\') {
            candidate += "/";
        }
        candidate += basename;
        if (stat(candidate.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            return candidate;
        }
    }
    
    // Return original path (will fail later if invalid)
    return configPath;
}

map<string, string> parseFastqMap(const vector<string>& fastqMapVec) {
    map<string, string> result;
    for (const auto& entry : fastqMapVec) {
        size_t eqPos = entry.find('=');
        if (eqPos == string::npos) {
            continue; // Skip invalid entries
        }
        string key = entry.substr(0, eqPos);
        string value = entry.substr(eqPos + 1);
        // Trim whitespace
        trimInPlace(key);
        trimInPlace(value);
        result[key] = value;
    }
    return result;
}

} // namespace CrMultiConfig

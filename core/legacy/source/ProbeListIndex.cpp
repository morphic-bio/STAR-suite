#include "ProbeListIndex.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

bool ProbeListIndex::containsIgnoreCase(const std::string& haystack, const std::string& needle) {
    std::string h = haystack;
    std::string n = needle;
    std::transform(h.begin(), h.end(), h.begin(), [](unsigned char c) { return std::toupper(c); });
    std::transform(n.begin(), n.end(), n.begin(), [](unsigned char c) { return std::toupper(c); });
    return h.find(n) != std::string::npos;
}

bool ProbeListIndex::load(const std::string &path, bool removeDeprecated, uint32_t *deprecatedCount) {
    geneIdToIndex_.clear();
    if (deprecatedCount) *deprecatedCount = 0;
    
    if (path.empty() || path == "-") return true; // treat empty as not provided
    std::ifstream in(path.c_str());
    if (!in.is_open()) return false;
    std::string line;
    uint32_t lineNo = 0;
    uint32_t deprecatedRemoved = 0;
    
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') continue;
        // trim
        size_t beg = line.find_first_not_of(" \t\r\n");
        size_t end = line.find_last_not_of(" \t\r\n");
        if (beg == std::string::npos) continue;
        std::string geneId = line.substr(beg, end - beg + 1);
        
        // Skip deprecated entries if removeDeprecated is enabled
        if (removeDeprecated && containsIgnoreCase(geneId, "DEPRECATED")) {
            deprecatedRemoved++;
            continue;
        }
        
        lineNo++;
        // 1-based index stored; guard against overflow but defer to geneIndex15()
        geneIdToIndex_[geneId] = lineNo;
    }
    
    if (deprecatedCount) *deprecatedCount = deprecatedRemoved;
    return true;
}



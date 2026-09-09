#ifndef STAR_SUITE_MITOCHONDRIAL_RANK_MASK_H
#define STAR_SUITE_MITOCHONDRIAL_RANK_MASK_H

#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Exact IDs from the active feature annotation, never display-name heuristics.
// An empty path disables masking; a supplied but unusable list is an error.
inline std::vector<uint8_t> loadMitochondrialRankMask(
    const std::vector<std::string>& features, const std::string& path)
{
    if (path.empty() || path == "-" || path == "None") return {};
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open mitochondrial gene list: " + path);
    std::set<std::string> ids, matched;
    std::string line, id;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        if (row >> id && id[0] != '#') ids.insert(id);
    }
    if (ids.empty()) throw std::runtime_error("Mitochondrial gene list is empty: " + path);
    std::vector<uint8_t> mask(features.size(), 0);
    size_t rows = 0;
    for (size_t i = 0; i < features.size(); ++i) {
        if (ids.count(features[i])) { mask[i] = 1; matched.insert(features[i]); ++rows; }
    }
    if (matched.empty()) throw std::runtime_error("No mitochondrial gene IDs match matrix feature IDs");
    std::cerr << "[OrdMag MT mask] " << matched.size() << "/" << ids.size()
              << " gene IDs matched; " << rows << " feature rows; original counts retained\n";
    return mask;
}

#endif

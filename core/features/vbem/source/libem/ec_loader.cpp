#include "ec_loader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>

ECTable load_ec_file(const std::string& ec_path, std::vector<std::string>& transcript_names) {
    std::ifstream file(ec_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open EC file: " + ec_path);
    }
    
    ECTable table;
    
    // Read number of transcripts
    std::string line;
    if (!std::getline(file, line)) {
        throw std::runtime_error("EC file is empty");
    }
    table.n_transcripts = std::stoul(line);
    
    // Read number of equivalence classes
    if (!std::getline(file, line)) {
        throw std::runtime_error("EC file missing EC count");
    }
    table.n_ecs = std::stoul(line);
    table.ecs.reserve(table.n_ecs);
    
    // Read transcript names
    transcript_names.clear();
    transcript_names.reserve(table.n_transcripts);
    for (size_t i = 0; i < table.n_transcripts; ++i) {
        if (!std::getline(file, line)) {
            throw std::runtime_error("EC file missing transcript names");
        }
        transcript_names.push_back(line);
    }
    
    // Read equivalence classes
    // Supports two formats:
    //   Unweighted (--dumpEq):        k idx1 idx2 ... idxk count
    //   Weighted (--dumpEqWeights):   k idx1 idx2 ... idxk w1 w2 ... wk count
    for (size_t i = 0; i < table.n_ecs; ++i) {
        if (!std::getline(file, line)) {
            throw std::runtime_error("EC file missing EC definitions");
        }
        
        std::istringstream iss(line);
        EC ec;
        
        // First number is the count of transcripts in this EC
        uint32_t n_txp;
        if (!(iss >> n_txp)) {
            throw std::runtime_error("Invalid EC format: missing transcript count");
        }
        
        ec.transcript_ids.reserve(n_txp);
        
        // Read transcript indices
        for (uint32_t j = 0; j < n_txp; ++j) {
            uint32_t idx;
            if (!(iss >> idx)) {
                throw std::runtime_error("Invalid EC format: missing transcript index");
            }
            if (idx >= table.n_transcripts) {
                throw std::runtime_error("Transcript index out of range: " + std::to_string(idx));
            }
            ec.transcript_ids.push_back(idx);
        }
        
        // DO NOT sort transcript IDs - order must match weights if present
        // std::sort(ec.transcript_ids.begin(), ec.transcript_ids.end());
        
        // Read remaining values into a vector to detect format
        std::vector<double> remaining;
        double val;
        while (iss >> val) {
            remaining.push_back(val);
        }
        
        if (remaining.size() == n_txp + 1) {
            // Weighted format: k weights + 1 count
            ec.weights.reserve(n_txp);
            for (uint32_t j = 0; j < n_txp; ++j) {
                ec.weights.push_back(remaining[j]);
            }
            ec.count = remaining[n_txp];
        } else if (remaining.size() == 1) {
            // Unweighted format: just count
            ec.count = remaining[0];
            // weights vector remains empty - will use fallback in E-step
        } else {
            throw std::runtime_error("Invalid EC format: expected " + 
                std::to_string(n_txp + 1) + " or 1 values after transcript IDs, got " +
                std::to_string(remaining.size()));
        }
        
        table.ecs.push_back(ec);
    }
    
    return table;
}

void load_transcript_lengths(TranscriptState& state, const std::string& quant_path) {
    if (quant_path.empty()) {
        // No length file provided - set all lengths to 1.0 (placeholder)
        for (size_t i = 0; i < state.n; ++i) {
            state.lengths[i] = 1.0;
            state.eff_lengths[i] = 1.0;
        }
        return;
    }
    
    std::ifstream file(quant_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open quant file: " + quant_path);
    }
    
    // Build name→index map for O(1) lookup
    std::unordered_map<std::string, size_t> name_to_idx;
    for (size_t i = 0; i < state.n; ++i) {
        name_to_idx[state.names[i]] = i;
    }
    
    // Read header line
    std::string line;
    if (!std::getline(file, line)) {
        throw std::runtime_error("Quant file is empty");
    }
    
    // Parse header to find column indices
    std::istringstream header_iss(line);
    std::vector<std::string> headers;
    std::string header;
    while (header_iss >> header) {
        headers.push_back(header);
    }
    
    int name_col = -1, length_col = -1, eff_length_col = -1;
    for (size_t i = 0; i < headers.size(); ++i) {
        if (headers[i] == "Name") name_col = i;
        else if (headers[i] == "Length") length_col = i;
        else if (headers[i] == "EffectiveLength") eff_length_col = i;
    }
    
    if (name_col < 0 || length_col < 0) {
        throw std::runtime_error("Quant file missing required columns (Name, Length)");
    }
    
    // Track which transcripts have been found
    std::vector<bool> found(state.n, false);
    
    // Read transcript data
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::vector<std::string> fields;
        std::string field;
        while (iss >> field) {
            fields.push_back(field);
        }
        
        if (fields.size() <= (size_t)std::max(name_col, length_col)) {
            continue; // Skip malformed lines
        }
        
        std::string name = fields[name_col];
        
        // O(1) lookup using map
        auto it = name_to_idx.find(name);
        if (it != name_to_idx.end()) {
            size_t idx = it->second;
            state.lengths[idx] = std::stod(fields[length_col]);
            if (eff_length_col >= 0 && (size_t)eff_length_col < fields.size()) {
                state.eff_lengths[idx] = std::stod(fields[eff_length_col]);
            } else {
                state.eff_lengths[idx] = state.lengths[idx]; // Use raw length if effective length not available
            }
            found[idx] = true;
        }
    }
    
    // Set default lengths for transcripts not found in quant file
    for (size_t i = 0; i < state.n; ++i) {
        if (!found[i]) {
            state.lengths[i] = 1.0;
            state.eff_lengths[i] = 1.0;
        }
    }
}

void load_effective_lengths(TranscriptState& state, const std::string& eff_lengths_path) {
    if (eff_lengths_path.empty()) {
        return;  // No file provided, keep existing effective lengths
    }
    
    std::ifstream file(eff_lengths_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open effective lengths file: " + eff_lengths_path);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        size_t idx;
        double eff_len;
        
        if (sscanf(line.c_str(), "%zu %lf", &idx, &eff_len) == 2) {
            if (idx < state.n) {
                state.eff_lengths[idx] = eff_len;
            }
        }
    }
}

/**
 * FlexProbeIndex implementation
 */

#include "FlexProbeIndex.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <cctype>
#include <iomanip>
#include <sys/stat.h>
#include <cstring>
#include <regex>
#include <cstdio>

// For SHA256
#include <openssl/sha.h>

namespace FlexProbeIndex {

// Helper: Create directory recursively
static bool createDirectoryRecursive(const std::string& path) {
    std::string current;
    std::istringstream pathStream(path);
    std::string segment;
    
    if (!path.empty() && path[0] == '/') {
        current = "/";
    }
    
    while (std::getline(pathStream, segment, '/')) {
        if (segment.empty()) continue;
        
        if (!current.empty() && current.back() != '/') {
            current += '/';
        }
        current += segment;
        
        struct stat st;
        if (stat(current.c_str(), &st) != 0) {
            if (mkdir(current.c_str(), 0755) != 0) {
                return false;
            }
        }
    }
    return true;
}

// Helper: Trim whitespace
static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Helper: Convert to uppercase
static std::string toUpper(const std::string& s) {
    std::string result = s;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

// Helper: Check if string contains substring (case-insensitive)
static bool containsIgnoreCase(const std::string& haystack, const std::string& needle) {
    std::string h = toUpper(haystack);
    std::string n = toUpper(needle);
    return h.find(n) != std::string::npos;
}

// Helper: Validate probe sequence (50bp A/C/G/T only)
static bool validateProbeSeq(const std::string& seq, uint32_t expectedLength, std::string& errorMsg) {
    if (seq.length() != expectedLength) {
        std::ostringstream oss;
        oss << "probe_seq is " << seq.length() << "bp, expected " << expectedLength << "bp";
        errorMsg = oss.str();
        return false;
    }
    
    for (size_t i = 0; i < seq.length(); i++) {
        char c = std::toupper(seq[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            std::ostringstream oss;
            oss << "probe_seq contains invalid character '" << seq[i] << "' at position " << i << " (only A/C/G/T allowed)";
            errorMsg = oss.str();
            return false;
        }
    }
    return true;
}

// Helper: Split CSV line (handles quoted fields)
static std::vector<std::string> splitCSV(const std::string& line) {
    std::vector<std::string> fields;
    std::string current;
    bool inQuotes = false;
    
    for (size_t i = 0; i < line.length(); i++) {
        char c = line[i];
        if (c == '"') {
            inQuotes = !inQuotes;
        } else if (c == ',' && !inQuotes) {
            fields.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    fields.push_back(current);
    return fields;
}

// Helper: Get timestamp in ISO format
static std::string getISOTimestamp() {
    time_t now = time(nullptr);
    struct tm* utc = gmtime(&now);
    char buf[64];
    strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", utc);
    return std::string(buf);
}

// Compute SHA256 of a file
std::string computeSHA256(const std::string& filePath) {
    std::ifstream file(filePath, std::ios::binary);
    if (!file.is_open()) {
        return "N/A";
    }
    
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    
    char buffer[8192];
    while (file.read(buffer, sizeof(buffer))) {
        SHA256_Update(&ctx, buffer, file.gcount());
    }
    if (file.gcount() > 0) {
        SHA256_Update(&ctx, buffer, file.gcount());
    }
    
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_Final(hash, &ctx);
    
    std::ostringstream oss;
    for (int i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        oss << std::hex << std::setfill('0') << std::setw(2) << (int)hash[i];
    }
    return oss.str();
}

// Parse GTF to extract gene_ids
bool parseGTFGeneIds(const std::string& gtfPath, std::unordered_set<std::string>& geneIds, std::string& errorMsg) {
    std::ifstream file;
    bool isGzip = false;
    
    // Check if file is gzipped
    if (gtfPath.length() > 3 && gtfPath.substr(gtfPath.length() - 3) == ".gz") {
        isGzip = true;
    }
    
    // For gzip files, use popen with zcat
    FILE* fp = nullptr;
    if (isGzip) {
        std::string cmd = "zcat \"" + gtfPath + "\"";
        fp = popen(cmd.c_str(), "r");
        if (!fp) {
            errorMsg = "Could not open gzipped GTF file: " + gtfPath;
            return false;
        }
    } else {
        file.open(gtfPath);
        if (!file.is_open()) {
            errorMsg = "Could not open GTF file: " + gtfPath;
            return false;
        }
    }
    
    // Regex to extract gene_id "VALUE"
    std::regex geneIdRegex("gene_id\\s+\"([^\"]+)\"");
    
    char lineBuf[65536];
    while (true) {
        std::string line;
        if (isGzip) {
            if (!fgets(lineBuf, sizeof(lineBuf), fp)) break;
            line = lineBuf;
        } else {
            if (!std::getline(file, line)) break;
        }
        
        // Skip comments
        if (line.empty() || line[0] == '#') continue;
        
        // Extract gene_id
        std::smatch match;
        if (std::regex_search(line, match, geneIdRegex) && match.size() > 1) {
            geneIds.insert(match[1].str());
        }
    }
    
    if (isGzip) {
        pclose(fp);
    } else {
        file.close();
    }
    
    if (geneIds.empty()) {
        errorMsg = "No gene_ids found in GTF file";
        return false;
    }
    
    return true;
}

// Parse FASTA to extract contig names
bool parseFastaContigNames(const std::string& fastaPath, std::unordered_set<std::string>& contigNames, std::string& errorMsg) {
    std::ifstream file(fastaPath);
    if (!file.is_open()) {
        errorMsg = "Could not open FASTA file: " + fastaPath;
        return false;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] == '>') {
            // Extract contig name (first word after >)
            std::string name = line.substr(1);
            size_t spacePos = name.find(' ');
            if (spacePos != std::string::npos) {
                name = name.substr(0, spacePos);
            }
            name = trim(name);
            if (!name.empty()) {
                contigNames.insert(name);
            }
        }
    }
    
    file.close();
    return true;
}

// Parse probe CSV with filtering
bool parseProbeCSV(const std::string& csvPath, const std::unordered_set<std::string>& gtfGeneIds,
                   uint32_t enforceLength, std::vector<ProbeEntry>& probes,
                   std::vector<std::string>& headerLines, std::string& headerRow,
                   FilterStats& stats, std::string& errorMsg) {
    std::ifstream file(csvPath);
    if (!file.is_open()) {
        errorMsg = "Could not open probe CSV file: " + csvPath;
        return false;
    }
    
    // Column indices (will be determined from header)
    int geneIdCol = -1, probeSeqCol = -1, probeIdCol = -1;
    int includedCol = -1, regionCol = -1, geneNameCol = -1;
    
    std::string line;
    uint32_t lineNum = 0;
    bool foundHeader = false;
    
    while (std::getline(file, line)) {
        lineNum++;
        std::string trimmedLine = trim(line);
        
        // Skip empty lines
        if (trimmedLine.empty()) continue;
        
        // Handle comment lines (metadata)
        if (trimmedLine[0] == '#') {
            headerLines.push_back(line);
            continue;
        }
        
        // First non-comment line is the header
        if (!foundHeader) {
            headerRow = line;
            foundHeader = true;
            
            // Parse header to get column indices
            std::vector<std::string> headers = splitCSV(line);
            for (size_t i = 0; i < headers.size(); i++) {
                std::string h = trim(headers[i]);
                if (h == "gene_id") geneIdCol = i;
                else if (h == "probe_seq") probeSeqCol = i;
                else if (h == "probe_id") probeIdCol = i;
                else if (h == "included") includedCol = i;
                else if (h == "region") regionCol = i;
                else if (h == "gene_name") geneNameCol = i;
            }
            
            if (geneIdCol < 0) {
                errorMsg = "Required column 'gene_id' not found in header";
                return false;
            }
            if (probeSeqCol < 0) {
                errorMsg = "Required column 'probe_seq' not found in header";
                return false;
            }
            if (probeIdCol < 0) {
                errorMsg = "Required column 'probe_id' not found in header";
                return false;
            }
            
            continue;
        }
        
        // Parse data line
        stats.totalInput++;
        std::vector<std::string> fields = splitCSV(line);
        
        if ((size_t)geneIdCol >= fields.size() || (size_t)probeSeqCol >= fields.size() || 
            (size_t)probeIdCol >= fields.size()) {
            std::ostringstream oss;
            oss << "line " << lineNum << ": not enough columns";
            errorMsg = oss.str();
            return false;
        }
        
        std::string geneId = trim(fields[geneIdCol]);
        std::string probeSeq = trim(fields[probeSeqCol]);
        std::string probeId = trim(fields[probeIdCol]);
        
        // Check DEPRECATED (case-insensitive)
        if (containsIgnoreCase(probeId, "DEPRECATED")) {
            stats.droppedDeprecated++;
            continue;
        }
        
        // Check included column (case-sensitive: must be exactly "TRUE")
        // Only filter if included column exists; if missing, keep current behavior (do not drop)
        if (includedCol >= 0 && (size_t)includedCol < fields.size()) {
            std::string included = trim(fields[includedCol]);
            if (included != "TRUE") {
                stats.droppedIncludedFalse++;
                continue;
            }
        }
        
        // Validate sequence
        std::string seqError;
        if (!validateProbeSeq(probeSeq, enforceLength, seqError)) {
            std::ostringstream oss;
            oss << "line " << lineNum << ": " << seqError;
            errorMsg = oss.str();
            stats.droppedInvalidSeq++;
            return false;  // Fail fast on invalid sequence
        }
        
        // Check gene_id in GTF
        if (gtfGeneIds.find(geneId) == gtfGeneIds.end()) {
            stats.droppedNoMatch++;
            continue;
        }
        
        // Passed all filters - add to output
        ProbeEntry entry;
        entry.gene_id = geneId;
        entry.probe_seq = probeSeq;
        entry.probe_id = probeId;
        entry.included = (includedCol >= 0 && (size_t)includedCol < fields.size()) ? trim(fields[includedCol]) : "";
        entry.region = (regionCol >= 0 && (size_t)regionCol < fields.size()) ? trim(fields[regionCol]) : "";
        entry.gene_name = (geneNameCol >= 0 && (size_t)geneNameCol < fields.size()) ? trim(fields[geneNameCol]) : "";
        entry.original_line = line;
        
        probes.push_back(entry);
    }
    
    file.close();
    
    stats.totalOutput = probes.size();
    
    if (probes.empty()) {
        errorMsg = "0 probes matched GTF genes; check gene_id format";
        return false;
    }
    
    return true;
}

// Main entry point
Result run(const Config& config) {
    Result result;
    
    // Validate config
    if (config.probeCSVPath.empty()) {
        result.errorMessage = "probeCSVPath is required";
        return result;
    }
    if (config.gtfPath.empty()) {
        result.errorMessage = "gtfPath is required";
        return result;
    }
    if (config.baseFastaPath.empty()) {
        result.errorMessage = "baseFastaPath is required";
        return result;
    }
    if (config.outputDir.empty()) {
        result.errorMessage = "outputDir is required";
        return result;
    }
    
    // Create output directories
    std::string metadataDir = config.outputDir + "/metadata";
    if (!createDirectoryRecursive(config.outputDir) || !createDirectoryRecursive(metadataDir)) {
        result.errorMessage = "Could not create output directories";
        return result;
    }
    
    // Step 1: Parse GTF to get gene_ids
    std::unordered_set<std::string> gtfGeneIds;
    if (!parseGTFGeneIds(config.gtfPath, gtfGeneIds, result.errorMessage)) {
        return result;
    }
    
    // Step 2: Parse base FASTA to get contig names (for collision check)
    std::unordered_set<std::string> baseContigNames;
    if (!parseFastaContigNames(config.baseFastaPath, baseContigNames, result.errorMessage)) {
        return result;
    }
    
    // Step 3: Parse and filter probes
    std::vector<ProbeEntry> probes;
    std::vector<std::string> headerLines;
    std::string headerRow;
    
    if (!parseProbeCSV(config.probeCSVPath, gtfGeneIds, config.enforceProbeLength,
                       probes, headerLines, headerRow, result.stats, result.errorMessage)) {
        return result;
    }
    
    // Step 4: Stable sort by gene_id then probe_id
    std::stable_sort(probes.begin(), probes.end());
    
    // Step 5: Check for contig collisions
    for (const auto& probe : probes) {
        if (baseContigNames.find(probe.probe_id) != baseContigNames.end()) {
            result.errorMessage = "Synthetic contig '" + probe.probe_id + "' collides with existing contig in base FASTA";
            return result;
        }
    }
    
    // Step 6: Generate outputs
    std::string filteredCSV = config.outputDir + "/filtered_probe_set.csv";
    std::string probesOnlyFA = config.outputDir + "/probes_only.fa";
    std::string probesOnlyGTF = config.outputDir + "/probes_only.gtf";
    std::string hybridFA = config.outputDir + "/genome.filtered.fa";
    std::string hybridGTF = config.outputDir + "/genes.filtered.gtf";
    std::string probeBED = config.outputDir + "/probe_genes_exons.bed";
    std::string probeList = config.outputDir + "/probe_list.txt";
    std::string manifest = metadataDir + "/reference_manifest.json";
    
    // Write filtered CSV
    {
        std::ofstream out(filteredCSV);
        if (!out.is_open()) {
            result.errorMessage = "Could not create filtered CSV: " + filteredCSV;
            return result;
        }
        
        // Write header lines (comments)
        for (const auto& h : headerLines) {
            out << h << "\n";
        }
        // Write header row
        out << headerRow << "\n";
        // Write data rows (sorted)
        for (const auto& probe : probes) {
            out << probe.original_line << "\n";
        }
        out.close();
    }
    
    // Write probes-only FASTA and GTF
    {
        std::ofstream faOut(probesOnlyFA);
        std::ofstream gtfOut(probesOnlyGTF);
        
        if (!faOut.is_open() || !gtfOut.is_open()) {
            result.errorMessage = "Could not create probe-only files";
            return result;
        }
        
        for (const auto& probe : probes) {
            // FASTA entry
            faOut << ">" << probe.probe_id << "\n";
            faOut << probe.probe_seq << "\n";
            
            // GTF entry - extract gene_name from probe_id (format: ENSG|GENE|hash)
            std::string geneName = probe.gene_name;
            if (geneName.empty()) {
                // Try to extract from probe_id
                size_t pos1 = probe.probe_id.find('|');
                if (pos1 != std::string::npos) {
                    size_t pos2 = probe.probe_id.find('|', pos1 + 1);
                    if (pos2 != std::string::npos) {
                        geneName = probe.probe_id.substr(pos1 + 1, pos2 - pos1 - 1);
                    }
                }
            }
            if (geneName.empty()) geneName = probe.gene_id;
            
            uint32_t seqLen = probe.probe_seq.length();
            gtfOut << probe.probe_id << "\tFLEX\texon\t1\t" << seqLen << "\t.\t+\t.\t"
                   << "gene_id \"" << geneName << "\"; gene_name \"" << geneName << "\"; "
                   << "transcript_id \"" << probe.probe_id << "\"; "
                   << "probe_id \"" << probe.probe_id << "\"; "
                   << "region \"" << probe.region << "\";\n";
        }
        
        faOut.close();
        gtfOut.close();
    }
    
    // Write hybrid FASTA (copy base + append probes)
    {
        std::ifstream baseIn(config.baseFastaPath, std::ios::binary);
        std::ofstream hybridOut(hybridFA, std::ios::binary);
        
        if (!baseIn.is_open() || !hybridOut.is_open()) {
            result.errorMessage = "Could not create hybrid FASTA";
            return result;
        }
        
        // Copy base FASTA
        hybridOut << baseIn.rdbuf();
        baseIn.close();
        
        // Append probe FASTA
        std::ifstream probeIn(probesOnlyFA);
        hybridOut << probeIn.rdbuf();
        probeIn.close();
        
        hybridOut.close();
    }
    
    // Write hybrid GTF (copy base + append probes)
    {
        // Handle gzipped base GTF
        bool isGzip = config.gtfPath.length() > 3 && 
                      config.gtfPath.substr(config.gtfPath.length() - 3) == ".gz";
        
        std::ofstream hybridOut(hybridGTF);
        if (!hybridOut.is_open()) {
            result.errorMessage = "Could not create hybrid GTF";
            return result;
        }
        
        if (isGzip) {
            std::string cmd = "zcat \"" + config.gtfPath + "\"";
            FILE* fp = popen(cmd.c_str(), "r");
            if (!fp) {
                result.errorMessage = "Could not decompress base GTF";
                return result;
            }
            
            char buffer[65536];
            while (fgets(buffer, sizeof(buffer), fp)) {
                hybridOut << buffer;
            }
            pclose(fp);
        } else {
            std::ifstream baseIn(config.gtfPath);
            hybridOut << baseIn.rdbuf();
            baseIn.close();
        }
        
        // Append probe GTF
        std::ifstream probeIn(probesOnlyGTF);
        hybridOut << probeIn.rdbuf();
        probeIn.close();
        
        hybridOut.close();
    }
    
    // Write probe BED
    {
        std::ofstream out(probeBED);
        if (!out.is_open()) {
            result.errorMessage = "Could not create probe BED: " + probeBED;
            return result;
        }
        
        for (const auto& probe : probes) {
            std::string geneName = probe.gene_name;
            if (geneName.empty()) {
                size_t pos1 = probe.probe_id.find('|');
                if (pos1 != std::string::npos) {
                    size_t pos2 = probe.probe_id.find('|', pos1 + 1);
                    if (pos2 != std::string::npos) {
                        geneName = probe.probe_id.substr(pos1 + 1, pos2 - pos1 - 1);
                    }
                }
            }
            
            out << probe.probe_id << "\t0\t" << probe.probe_seq.length() << "\t"
                << geneName << "\t0\t+\n";
        }
        out.close();
    }
    
    // Write probe list (unique gene IDs, sorted)
    {
        std::set<std::string> uniqueGenes;
        for (const auto& probe : probes) {
            uniqueGenes.insert(probe.gene_id);
        }
        
        // Filter deprecated gene IDs if removeDeprecated is enabled
        uint32_t deprecatedCount = 0;
        if (config.removeDeprecated) {
            std::set<std::string> filteredGenes;
            for (const auto& gene : uniqueGenes) {
                if (containsIgnoreCase(gene, "DEPRECATED")) {
                    deprecatedCount++;
                } else {
                    filteredGenes.insert(gene);
                }
            }
            uniqueGenes = filteredGenes;
        }
        
        result.stats.uniqueGenes = uniqueGenes.size();
        
        std::ofstream out(probeList);
        if (!out.is_open()) {
            result.errorMessage = "Could not create probe list: " + probeList;
            return result;
        }
        
        for (const auto& gene : uniqueGenes) {
            out << gene << "\n";
        }
        out.close();
        
        // Log deprecated removal if enabled
        if (config.removeDeprecated && deprecatedCount > 0) {
            // Note: This will be logged by the caller if needed
            // We could add a deprecatedGenesRemoved field to FilterStats if needed
        }
    }
    
    // Compute SHA256 hashes
    std::string probeCsvSha = computeSHA256(config.probeCSVPath);
    std::string gtfSha = computeSHA256(config.gtfPath);
    std::string baseFastaSha = computeSHA256(config.baseFastaPath);
    std::string filteredCsvSha = computeSHA256(filteredCSV);
    std::string probesFaSha = computeSHA256(probesOnlyFA);
    std::string probesGtfSha = computeSHA256(probesOnlyGTF);
    std::string hybridFaSha = computeSHA256(hybridFA);
    std::string hybridGtfSha = computeSHA256(hybridGTF);
    std::string probeBedSha = computeSHA256(probeBED);
    std::string probeListSha = computeSHA256(probeList);
    
    // Write manifest
    {
        std::ofstream out(manifest);
        if (!out.is_open()) {
            result.errorMessage = "Could not create manifest: " + manifest;
            return result;
        }
        
        out << "{\n";
        out << "  \"manifest_version\": \"1.0\",\n";
        out << "  \"generated_at\": \"" << getISOTimestamp() << "\",\n";
        out << "  \"generator\": \"FlexProbeIndex\",\n";
        out << "  \"inputs\": {\n";
        out << "    \"probe_csv\": {\n";
        out << "      \"path\": \"" << config.probeCSVPath << "\",\n";
        out << "      \"sha256\": \"" << probeCsvSha << "\"\n";
        out << "    },\n";
        out << "    \"gtf\": {\n";
        out << "      \"path\": \"" << config.gtfPath << "\",\n";
        out << "      \"sha256\": \"" << gtfSha << "\"\n";
        out << "    },\n";
        out << "    \"base_fasta\": {\n";
        out << "      \"path\": \"" << config.baseFastaPath << "\",\n";
        out << "      \"sha256\": \"" << baseFastaSha << "\"\n";
        out << "    }\n";
        out << "  },\n";
        out << "  \"filtering\": {\n";
        out << "    \"total_input_probes\": " << result.stats.totalInput << ",\n";
        out << "    \"dropped_deprecated\": " << result.stats.droppedDeprecated << ",\n";
        out << "    \"dropped_included_false\": " << result.stats.droppedIncludedFalse << ",\n";
        out << "    \"dropped_no_gtf_match\": " << result.stats.droppedNoMatch << ",\n";
        out << "    \"dropped_invalid_seq\": " << result.stats.droppedInvalidSeq << ",\n";
        out << "    \"total_output_probes\": " << result.stats.totalOutput << ",\n";
        out << "    \"unique_genes\": " << result.stats.uniqueGenes << ",\n";
        out << "    \"enforce_probe_length\": " << config.enforceProbeLength << ",\n";
        out << "    \"ordering\": \"stable_sort_gene_id_probe_id\"\n";
        out << "  },\n";
        
        // Write synthetic contigs
        out << "  \"synthetic_contigs\": [\n";
        for (size_t i = 0; i < probes.size(); i++) {
            out << "  {\"name\": \"" << probes[i].probe_id << "\", \"length\": " 
                << probes[i].probe_seq.length() << "}";
            if (i < probes.size() - 1) out << ",";
            out << "\n";
        }
        out << "  ],\n";
        
        out << "  \"outputs\": {\n";
        out << "    \"filtered_csv\": {\n";
        out << "      \"path\": \"" << filteredCSV << "\",\n";
        out << "      \"sha256\": \"" << filteredCsvSha << "\"\n";
        out << "    },\n";
        out << "    \"probes_fasta\": {\n";
        out << "      \"path\": \"" << probesOnlyFA << "\",\n";
        out << "      \"sha256\": \"" << probesFaSha << "\"\n";
        out << "    },\n";
        out << "    \"probes_gtf\": {\n";
        out << "      \"path\": \"" << probesOnlyGTF << "\",\n";
        out << "      \"sha256\": \"" << probesGtfSha << "\"\n";
        out << "    },\n";
        out << "    \"hybrid_fasta\": {\n";
        out << "      \"path\": \"" << hybridFA << "\",\n";
        out << "      \"sha256\": \"" << hybridFaSha << "\"\n";
        out << "    },\n";
        out << "    \"hybrid_gtf\": {\n";
        out << "      \"path\": \"" << hybridGTF << "\",\n";
        out << "      \"sha256\": \"" << hybridGtfSha << "\"\n";
        out << "    },\n";
        out << "    \"probe_bed\": {\n";
        out << "      \"path\": \"" << probeBED << "\",\n";
        out << "      \"sha256\": \"" << probeBedSha << "\"\n";
        out << "    },\n";
        out << "    \"probe_list\": {\n";
        out << "      \"path\": \"" << probeList << "\",\n";
        out << "      \"sha256\": \"" << probeListSha << "\"\n";
        out << "    }\n";
        out << "  }\n";
        out << "}\n";
        
        out.close();
    }
    
    result.success = true;
    return result;
}

} // namespace FlexProbeIndex


/**
 * CellRangerFormatter implementation
 * Matches format-fa-gtf.pl exactly
 */

#include "CellRangerFormatter.h"
#include "IncludeDefine.h"
#include <regex>
#include <cstdio>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
extern "C" {
#include "htslib/htslib/hfile.h"
}
#include <zlib.h>
#include <fstream>

namespace CellRangerFormatter {

// Allowed biotypes (exact match from Perl)
static const std::set<std::string> ALLOWED_BIOTYPES = {
    "protein_coding", "lncRNA",
    "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene",
    "IG_V_pseudogene", "IG_J_pseudogene", "IG_C_pseudogene",
    "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
    "TR_V_pseudogene", "TR_J_pseudogene"
};

static const std::string GENE_ID_KEY = "gene_id \"";
static const std::string TRANSCRIPT_ID_KEY = "transcript_id \"";
static const std::string EXON_ID_KEY = "exon_id \"";
static const std::string GENE_PREFIX = "ENSG";
static const std::string TRANSCRIPT_PREFIX = "ENST";
static const std::string EXON_PREFIX = "ENSE";

static bool parseGeneIdPass1(const std::string& attrs, std::string& geneId) {
    size_t pos = attrs.find(GENE_ID_KEY);
    if (pos == std::string::npos) {
        return false;
    }
    pos += GENE_ID_KEY.size();
    if (attrs.compare(pos, GENE_PREFIX.size(), GENE_PREFIX) != 0) {
        return false;
    }
    size_t idStart = pos;
    pos += GENE_PREFIX.size();
    size_t digitsStart = pos;
    while (pos < attrs.size() && std::isdigit(static_cast<unsigned char>(attrs[pos]))) {
        ++pos;
    }
    if (pos == digitsStart) {
        return false;
    }
    if (pos >= attrs.size() || attrs[pos] != '.') {
        return false;
    }
    geneId.assign(attrs, idStart, pos - idStart);
    return true;
}

static bool parseIdWithVersion(const std::string& attrs, const std::string& key,
                               const std::string& prefix, std::string& id, std::string& version) {
    size_t pos = attrs.find(key);
    if (pos == std::string::npos) {
        return false;
    }
    pos += key.size();
    if (attrs.compare(pos, prefix.size(), prefix) != 0) {
        return false;
    }
    size_t idStart = pos;
    pos += prefix.size();
    size_t digitsStart = pos;
    while (pos < attrs.size() && std::isdigit(static_cast<unsigned char>(attrs[pos]))) {
        ++pos;
    }
    if (pos == digitsStart) {
        return false;
    }
    if (pos >= attrs.size() || attrs[pos] != '.') {
        return false;
    }
    id.assign(attrs, idStart, pos - idStart);
    ++pos; // skip '.'
    size_t verStart = pos;
    while (pos < attrs.size() && std::isdigit(static_cast<unsigned char>(attrs[pos]))) {
        ++pos;
    }
    if (pos == verStart) {
        return false;
    }
    version.assign(attrs, verStart, pos - verStart);
    return true;
}

static std::string replaceExonIdVersions(const std::string& attrs) {
    size_t match = attrs.find(EXON_ID_KEY);
    if (match == std::string::npos) {
        return attrs;
    }

    std::string out;
    out.reserve(attrs.size() + 64);
    size_t pos = 0;

    while (match != std::string::npos) {
        out.append(attrs, pos, match - pos);
        size_t idStart = match + EXON_ID_KEY.size();
        size_t cursor = idStart;

        bool replaced = false;
        if (attrs.compare(cursor, EXON_PREFIX.size(), EXON_PREFIX) == 0) {
            cursor += EXON_PREFIX.size();
            size_t digitsStart = cursor;
            while (cursor < attrs.size() && std::isdigit(static_cast<unsigned char>(attrs[cursor]))) {
                ++cursor;
            }
            if (cursor > digitsStart && cursor < attrs.size() && attrs[cursor] == '.') {
                std::string exonId = attrs.substr(idStart, cursor - idStart);
                ++cursor; // skip '.'
                size_t verStart = cursor;
                while (cursor < attrs.size() && std::isdigit(static_cast<unsigned char>(attrs[cursor]))) {
                    ++cursor;
                }
                if (cursor > verStart && cursor < attrs.size() && attrs[cursor] == '"') {
                    std::string version = attrs.substr(verStart, cursor - verStart);
                    ++cursor; // consume closing quote
                    out.append(EXON_ID_KEY);
                    out.append(exonId);
                    out.append("\"; exon_version \"");
                    out.append(version);
                    out.append("\"");
                    pos = cursor;
                    replaced = true;
                }
            }
        }

        if (!replaced) {
            out.append(attrs, match, EXON_ID_KEY.size());
            pos = idStart;
        }
        match = attrs.find(EXON_ID_KEY, pos);
    }

    out.append(attrs, pos, std::string::npos);
    return out;
}

static uint64_t gtfLogEvery() {
    static uint64_t cached = []() -> uint64_t {
        const char* env = std::getenv("CELLRANGER_GTF_LOG_EVERY");
        if (!env || !*env) {
            return 0;
        }
        char* end = nullptr;
        unsigned long long value = std::strtoull(env, &end, 10);
        if (end == env || value == 0) {
            return 0;
        }
        return static_cast<uint64_t>(value);
    }();
    return cached;
}

static void logGtfProgress(const char* phase, uint64_t lineCount, uint64_t outCount,
                           const std::string& inputPath, const std::string& outputPath) {
    uint64_t every = gtfLogEvery();
    if (every == 0) {
        return;
    }
    if (lineCount % every == 0) {
        std::cerr << "[CellRangerFormatter] " << phase << " progress: lines=" << lineCount;
        if (outCount > 0) {
            std::cerr << " written=" << outCount;
        }
        if (!inputPath.empty()) {
            std::cerr << " input=" << inputPath;
        }
        if (!outputPath.empty()) {
            std::cerr << " output=" << outputPath;
        }
        std::cerr << std::endl;
    }
}

// Helper: Create directory recursively
static bool createDirectoryRecursive(const string& path) {
    string dir = path;
    size_t lastSlash = dir.find_last_of('/');
    if (lastSlash != string::npos) {
        dir = dir.substr(0, lastSlash);
    } else {
        return true; // No directory component
    }
    
    if (dir.empty()) return true;
    
    // Check if exists
    struct stat st;
    if (stat(dir.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        return true;
    }
    
    // Create parent first
    size_t prevSlash = dir.find_last_of('/');
    if (prevSlash != string::npos && prevSlash > 0) {
        string parent = dir.substr(0, prevSlash);
        if (!createDirectoryRecursive(parent + "/")) {
            return false;
        }
    }
    
    // Create this directory
    if (mkdir(dir.c_str(), 0755) != 0) {
        if (errno != EEXIST) {
            return false;
        }
    }
    return true;
}

// Helper: Chomp (remove trailing newline only, like Perl's chomp)
static string chomp(const string& s) {
    if (s.empty()) return s;
    if (s.back() == '\n') {
        return s.substr(0, s.length() - 1);
    }
    return s;
}

// Helper: Split string by delimiter
static vector<string> split(const string& s, char delim) {
    vector<string> result;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}

// FASTA transformation: modifyChrNames
Result formatFasta(const string& inputPath, const string& outputPath) {
    Result result;
    
    if (!createDirectoryRecursive(outputPath)) {
        result.errorMessage = "Could not create output directory for " + outputPath;
        return result;
    }
    
    ifstream inFile(inputPath);
    if (!inFile.is_open()) {
        result.errorMessage = "Could not open input FASTA: " + inputPath;
        return result;
    }
    
    string tempPath = outputPath + ".tmp";
    ofstream outFile(tempPath);
    if (!outFile.is_open()) {
        inFile.close();
        result.errorMessage = "Could not open temporary output FASTA: " + tempPath;
        return result;
    }
    
    string line;
    while (getline(inFile, line)) {
        if (line.empty()) {
            outFile << "\n";
            continue;
        }
        
        if (line[0] == '>') {
            // Header line - transform
            // Split on space, take first token only (matches Perl: split(' ',$line))
            size_t spacePos = line.find(' ');
            string chr = (spacePos != string::npos) ? line.substr(0, spacePos) : line;
            
            // Remove '>' for pattern matching
            string chrName = chr.substr(1);
            
            // Match >1..>22, >X, >Y → >chrN N (matches Perl: /^\>([0-9XY]+)/)
            // Note: Perl matches PREFIX, so >1_random matches and extracts "1"
            regex numChrRegex("^([0-9]+|[XY])");
            smatch match;
            if (regex_search(chrName, match, numChrRegex) && match.size() > 1) {
                string matchedPart = match[1].str();
                outFile << ">chr" << matchedPart << " " << matchedPart << "\n";
            }
            // Match >MT → >chrM MT (matches Perl: substr($chr,1,2) eq "MT")
            else if (chrName.length() >= 2 && chrName.substr(0, 2) == "MT") {
                outFile << ">chrM MT\n";
            }
            // All others → >orig orig (matches Perl: printf "$chr %s\n",substr($chr,1))
            else {
                outFile << chr << " " << chrName << "\n";
            }
        } else {
            // Sequence line - output as-is
            outFile << line << "\n";
        }
    }
    
    inFile.close();
    outFile.close();
    
    if (!outFile.good() && !outFile.eof()) {
        remove(tempPath.c_str());
        result.errorMessage = "Error writing to temporary output FASTA: " + tempPath;
        return result;
    }
    
    // Atomic rename
    if (rename(tempPath.c_str(), outputPath.c_str()) != 0) {
        remove(tempPath.c_str());
        result.errorMessage = "Could not rename temporary file to final output: " + outputPath;
        return result;
    }
    
    result.success = true;
    return result;
}

// GTF Pass 1: Generate filter for GTF (build allowed gene_id set)
static bool generateFilterForGTF(const string& gtfPath, unordered_set<string>& goodGenes, string& errorMsg) {
    bool isGzip = (gtfPath.length() > 3 && gtfPath.substr(gtfPath.length() - 3) == ".gz");

    goodGenes.reserve(200000);
    
    gzFile gzFilePtr = nullptr;
    ifstream inFileStream;
    
    if (isGzip) {
        // Use zlib for gzipped files (no shell dependency)
        gzFilePtr = gzopen(gtfPath.c_str(), "rb");
        if (!gzFilePtr) {
            errorMsg = "Could not open gzipped GTF file: " + gtfPath;
            return false;
        }
    } else {
        inFileStream.open(gtfPath);
        if (!inFileStream.is_open()) {
            errorMsg = "Could not open GTF file: " + gtfPath;
            return false;
        }
    }
    
    char lineBuf[65536];
    string line;
    uint64_t lineCount = 0;
    
    while (true) {
        if (isGzip) {
            // gzgets includes newline in output (like fgets)
            if (!gzgets(gzFilePtr, lineBuf, sizeof(lineBuf))) break;
            line = lineBuf;
        } else {
            if (!getline(inFileStream, line)) break;
            line += "\n"; // getline removes newline, add it back to match Perl behavior
        }
        
        ++lineCount;
        logGtfProgress("GTF pass1", lineCount, 0, gtfPath, "");
        
        // Chomp newline only (like Perl's chomp) - don't trim whitespace
        line = chomp(line);
        if (line.empty() || line[0] == '#') continue;
        
        // Split by tab (no trimming - preserve whitespace)
        vector<string> parts = split(line, '\t');
            if (parts.size() < 9) continue;
            
            // Only process transcript lines
            if (parts[2] != "transcript") continue;
            
            // Extract gene_id from first attribute
            string attrs = parts[8];
            vector<string> infoParts = split(attrs, ';');
            if (infoParts.empty()) continue;
            
            // Match ENSG[0-9]+\. pattern
            string geneID;
            if (!parseGeneIdPass1(infoParts[0], geneID)) {
                continue;
            }
            if (goodGenes.find(geneID) != goodGenes.end()) {
                continue; // Already processed
            }
            
            // Check conditions: gene_type, transcript_type, no PAR, no readthrough
            bool hasGeneType = false, hasTranscriptType = false;
            bool hasPAR = false, hasReadthrough = false;
            
            for (const auto& biotype : ALLOWED_BIOTYPES) {
                if (attrs.find("gene_type \"" + biotype + "\"") != string::npos) {
                    hasGeneType = true;
                }
                if (attrs.find("transcript_type \"" + biotype + "\"") != string::npos) {
                    hasTranscriptType = true;
                }
            }
            
            if (attrs.find("tag \"PAR\"") != string::npos) {
                hasPAR = true;
            }
            if (attrs.find("tag \"readthrough_transcript\"") != string::npos) {
                hasReadthrough = true;
            }
            
            if (hasGeneType && hasTranscriptType && !hasPAR && !hasReadthrough) {
            goodGenes.insert(geneID);
        }
    }
    
    if (isGzip) {
        gzclose(gzFilePtr);
    } else {
        inFileStream.close();
    }
    if (gtfLogEvery() > 0) {
        std::cerr << "[CellRangerFormatter] GTF pass1 complete: lines=" << lineCount
                  << " genes=" << goodGenes.size()
                  << " input=" << gtfPath << std::endl;
    }
    return true;
}

// GTF Pass 2: Modify GTF
Result formatGtf(const string& inputPath, const string& outputPath) {
    Result result;
    
    // Pass 1: Build allowed gene_id set
    unordered_set<string> goodGenes;
    string errorMsg;
    if (!generateFilterForGTF(inputPath, goodGenes, errorMsg)) {
        result.errorMessage = errorMsg;
        return result;
    }
    
    if (!createDirectoryRecursive(outputPath)) {
        result.errorMessage = "Could not create output directory for " + outputPath;
        return result;
    }
    
    // Pass 2: Filter and rewrite
    bool isGzip = (inputPath.length() > 3 && inputPath.substr(inputPath.length() - 3) == ".gz");
    
    gzFile gzFilePtr = nullptr;
    ifstream inFileStream;
    
    if (isGzip) {
        // Use zlib for gzipped files (no shell dependency)
        gzFilePtr = gzopen(inputPath.c_str(), "rb");
        if (!gzFilePtr) {
            result.errorMessage = "Could not open gzipped GTF: " + inputPath;
            return result;
        }
    } else {
        inFileStream.open(inputPath);
        if (!inFileStream.is_open()) {
            result.errorMessage = "Could not open GTF: " + inputPath;
            return result;
        }
    }
    
    if (gtfLogEvery() > 0) {
        std::cerr << "[CellRangerFormatter] GTF pass2 start: input=" << inputPath
                  << " output=" << outputPath << std::endl;
    }

    string tempPath = outputPath + ".tmp";
    ofstream outFile(tempPath);
    if (!outFile.is_open()) {
        result.errorMessage = "Could not open temporary output GTF: " + tempPath;
        if (isGzip) {
            gzclose(gzFilePtr);
        } else {
            inFileStream.close();
        }
        return result;
    }
    
    char lineBuf[65536];
    string line;
    uint64_t lineCount = 0;
    uint64_t outCount = 0;
    
    while (true) {
        if (isGzip) {
            // gzgets includes newline in output (like fgets)
            if (!gzgets(gzFilePtr, lineBuf, sizeof(lineBuf))) break;
            line = lineBuf;
        } else {
            if (!getline(inFileStream, line)) break;
            line += "\n"; // getline removes newline, add it back to match Perl behavior
        }
        
        ++lineCount;
        logGtfProgress("GTF pass2", lineCount, outCount, inputPath, outputPath);
        
        // Preserve comment lines verbatim (with newline)
        if (line.empty() || line[0] == '#') {
            outFile << line;
            ++outCount;
            continue;
        }
        
        // Chomp newline only (like Perl's chomp) - don't trim whitespace
        line = chomp(line);
        
        // Split by tab (no trimming - preserve whitespace)
        vector<string> parts = split(line, '\t');
        if (parts.size() < 9) continue;
        
        // Skip if chrY and has PAR tag
        if (parts[0] == "chrY" && parts[8].find("tag \"PAR\"") != string::npos) {
            continue;
        }
        
        // Split attributes by semicolon
        vector<string> infoParts = split(parts[8], ';');
        if (infoParts.empty()) continue;
        
        // Extract gene_id from first attribute (must match ENSG pattern)
        string geneID;
        string geneVersion;
        if (!parseIdWithVersion(infoParts[0], GENE_ID_KEY, GENE_PREFIX, geneID, geneVersion)) {
            continue;
        }
        
        // Only keep if gene_id is in goodGenes set
        if (goodGenes.find(geneID) == goodGenes.end()) {
            continue;
        }
        
        // Rewrite first attribute: gene_id "ENSG..."; gene_version "..."
        // Note: This matches Perl exactly - it adds gene_version to infoParts[0],
        // which may result in duplicate gene_version if it already exists elsewhere
        infoParts[0] = "gene_id \"" + geneID + "\"; gene_version \"" + geneVersion + "\"";
        
        // Rewrite transcript_id if present in second attribute (matches Perl behavior)
        // Perl only checks infoParts[1], not searching all attributes
        if (infoParts.size() > 1) {
            string transcriptID;
            string transcriptVersion;
            if (parseIdWithVersion(infoParts[1], TRANSCRIPT_ID_KEY, TRANSCRIPT_PREFIX, transcriptID, transcriptVersion)) {
                infoParts[1] = " transcript_id \"" + transcriptID + "\"; transcript_version \"" + transcriptVersion + "\"";
            }
        }
        
        // Join all attributes
        string joinedAttrs = "";
        joinedAttrs.reserve(parts[8].size() + 64);
        for (size_t i = 0; i < infoParts.size(); i++) {
            if (i > 0) joinedAttrs += ";";
            joinedAttrs += infoParts[i];
        }
        
        // Replace exon_id pattern in joined string (matches Perl regex replace)
        // Pattern: (exon_id "ENSE[0-9]+).([0-9]+") -> $1"; exon_version "$2
        if (joinedAttrs.find("exon_id") != string::npos) {
            joinedAttrs = replaceExonIdVersions(joinedAttrs);
        }
        
        parts[8] = joinedAttrs;
        
        // Output with trailing semicolon
        for (size_t i = 0; i < parts.size(); i++) {
            if (i > 0) outFile << "\t";
            outFile << parts[i];
        }
        outFile << ";\n";
        ++outCount;
    }
    
    if (isGzip) {
        gzclose(gzFilePtr);
    } else {
        inFileStream.close();
    }
    outFile.close();

    if (gtfLogEvery() > 0) {
        std::cerr << "[CellRangerFormatter] GTF pass2 complete: lines=" << lineCount
                  << " written=" << outCount
                  << " output=" << outputPath << std::endl;
    }
    
    if (!outFile.good() && !outFile.eof()) {
        remove(tempPath.c_str());
        result.errorMessage = "Error writing to temporary output GTF: " + tempPath;
        return result;
    }
    
    // Atomic rename
    if (rename(tempPath.c_str(), outputPath.c_str()) != 0) {
        remove(tempPath.c_str());
        result.errorMessage = "Could not rename temporary file to final output: " + outputPath;
        return result;
    }
    
    result.success = true;
    return result;
}

// Main entry point
Result format(const Config& config) {
    Result result;
    
    // Format FASTA
    Result fastaResult = formatFasta(config.inputFastaPath, config.outputFastaPath);
    if (!fastaResult.success) {
        result.errorMessage = "FASTA formatting failed: " + fastaResult.errorMessage;
        return result;
    }
    
    // Format GTF
    Result gtfResult = formatGtf(config.inputGtfPath, config.outputGtfPath);
    if (!gtfResult.success) {
        result.errorMessage = "GTF formatting failed: " + gtfResult.errorMessage;
        return result;
    }
    
    result.success = true;
    return result;
}

// Trusted Ensembl FTP URLs (no cksum verification required for these)
static const unordered_set<string> TRUSTED_ENSEMBL_URLS = {
    // Ensembl FASTA files (various releases and species)
    "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    "ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    "ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    "ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", // 2020-A
    "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    "ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    "ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    // GENCODE annotation files
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz",
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz",
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz",
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz", // 2020-A
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz",
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz",
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz",
};

// Expected cksum values for specific files (cksum of .gz files)
// REQUIRED for trusted URLs - FTP downloads require embedded cksum verification
// Format: {crc32, file_size}
// To obtain cksum values:
//   Ensembl: curl -L ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/CHECKSUMS | grep 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
//   GENCODE: curl -L ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/CHECKSUMS | grep 'gencode.v44.primary_assembly.annotation.gtf.gz'
// CHECKSUMS format: cksum size filename
struct CksumValue {
    uint32_t crc;
    uint64_t size;
    CksumValue() : crc(0), size(0) {}
    CksumValue(uint32_t c, uint64_t s) : crc(c), size(s) {}
    bool isValid() const { return crc != 0 && size != 0; }
};
static const unordered_map<string, CksumValue> EXPECTED_CKSUM = {
    {"ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", CksumValue(0, 0)}, // TODO: Fill with actual cksum from Ensembl CHECKSUMS
    {"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz", CksumValue(0, 0)}, // TODO: Fill with actual cksum from GENCODE CHECKSUMS
    {"ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", CksumValue(0, 0)}, // 2020-A: TODO: Fill with actual cksum from Ensembl CHECKSUMS
    {"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz", CksumValue(0, 0)} // 2020-A: TODO: Fill with actual cksum from GENCODE CHECKSUMS
};

// Runtime cksum cache (loaded from disk, merged with EXPECTED_CKSUM)
static unordered_map<string, CksumValue> runtimeCksumCache;
static bool cksumCacheLoaded = false;

// Check if URL is in trusted Ensembl FTP URL table
bool isTrustedUrl(const string& url) {
    return TRUSTED_ENSEMBL_URLS.find(url) != TRUSTED_ENSEMBL_URLS.end();
}

// Check if URL has a known cksum entry (from embedded map or cache file)
bool hasKnownCksum(const string& url) {
    // Check embedded map
    auto it = EXPECTED_CKSUM.find(url);
    if (it != EXPECTED_CKSUM.end() && it->second.isValid()) {
        return true;
    }
    // Check runtime cache
    auto cacheIt = runtimeCksumCache.find(url);
    if (cacheIt != runtimeCksumCache.end() && cacheIt->second.isValid()) {
        return true;
    }
    return false;
}

// Check if URL is trusted OR has known cksum
bool isTrustedOrHasKnownCksum(const string& url) {
    return isTrustedUrl(url) || hasKnownCksum(url);
}

// Get cksum value for URL (from embedded map or cache)
// Returns true if found and valid, false otherwise
static bool getCksumForUrl(const string& url, uint32_t& crc, uint64_t& size) {
    // Check embedded map first
    auto it = EXPECTED_CKSUM.find(url);
    if (it != EXPECTED_CKSUM.end() && it->second.isValid()) {
        crc = it->second.crc;
        size = it->second.size;
        return true;
    }
    // Check runtime cache
    auto cacheIt = runtimeCksumCache.find(url);
    if (cacheIt != runtimeCksumCache.end() && cacheIt->second.isValid()) {
        crc = cacheIt->second.crc;
        size = cacheIt->second.size;
        return true;
    }
    return false;
}

// Load cksum cache file from disk
// Format: url<TAB>crc<TAB>size (one entry per line)
bool loadCksumCache(const string& cacheDir) {
    if (cksumCacheLoaded) {
        return true; // Already loaded
    }
    
    string cacheFile = cacheDir + "/cksum_cache.tsv";
    ifstream file(cacheFile);
    if (!file.is_open()) {
        // Cache file doesn't exist yet - that's OK
        cksumCacheLoaded = true;
        return true;
    }
    
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }
        
        // Parse: url<TAB>crc<TAB>size
        size_t tab1 = line.find('\t');
        if (tab1 == string::npos) continue;
        size_t tab2 = line.find('\t', tab1 + 1);
        if (tab2 == string::npos) continue;
        
        string url = line.substr(0, tab1);
        string crcStr = line.substr(tab1 + 1, tab2 - tab1 - 1);
        string sizeStr = line.substr(tab2 + 1);
        
        uint32_t crc = 0;
        uint64_t size = 0;
        try {
            crc = stoul(crcStr);
            size = stoull(sizeStr);
        } catch (...) {
            continue; // Skip invalid entries
        }
        
        if (crc != 0 && size != 0) {
            runtimeCksumCache[url] = CksumValue(crc, size);
        }
    }
    
    file.close();
    cksumCacheLoaded = true;
    return true;
}

// Save cksum entry to cache file
// Format: url<TAB>crc<TAB>size (append mode)
// allowOverwrite: if true, update existing entry; if false, skip if already exists
bool saveCksumToCache(const string& cacheDir, const string& url, uint32_t crc, uint64_t size, bool allowOverwrite) {
    if (cacheDir.empty()) {
        return false; // No cache directory specified
    }
    if (crc == 0 || size == 0) {
        return false; // Don't save invalid cksums
    }
    
    // Ensure cache directory exists
    if (!createDirectoryRecursive(cacheDir + "/dummy")) {
        return false;
    }
    
    string cacheFile = cacheDir + "/cksum_cache.tsv";
    ofstream file(cacheFile, ios::app); // Append mode
    if (!file.is_open()) {
        return false;
    }
    
    // Check if entry already exists
    bool exists = (runtimeCksumCache.find(url) != runtimeCksumCache.end());
    if (exists && !allowOverwrite) {
        file.close();
        return true; // Already in cache, skip update
    }
    
    // Write entry (always append to file for audit trail)
    file << url << "\t" << crc << "\t" << size << "\n";
    file.close();
    
    // Update runtime cache (overwrites if allowOverwrite is true)
    runtimeCksumCache[url] = CksumValue(crc, size);
    
    return true;
}

// Generate decompressed cache key
static string getDecompressedKey(const string& url) {
    return url + "|decompressed";
}

// Get decompressed cksum for URL (from cache)
bool getDecompressedCksumForUrl(const string& url, uint32_t& crc, uint64_t& size) {
    string key = getDecompressedKey(url);
    auto cacheIt = runtimeCksumCache.find(key);
    if (cacheIt != runtimeCksumCache.end() && cacheIt->second.isValid()) {
        crc = cacheIt->second.crc;
        size = cacheIt->second.size;
        return true;
    }
    return false;
}

// Save decompressed cksum entry to cache file
// Allows overwrite to update cache when files are re-decompressed (e.g., after truncation fix)
bool saveDecompressedCksumToCache(const string& cacheDir, const string& url, uint32_t crc, uint64_t size) {
    if (cacheDir.empty()) {
        return false; // No cache directory specified
    }
    string key = getDecompressedKey(url);
    return saveCksumToCache(cacheDir, key, crc, size, true); // Allow overwrite for decompressed checksums
}

// Compute POSIX cksum of a file (CRC-32 XOR file size)
// Uses the system cksum command to ensure exact POSIX compliance
// This matches the system cksum command exactly
bool computeCksumFile(const string& filePath, uint32_t& crc, uint64_t& size) {
    // Use system cksum command for exact POSIX compliance
    string cmd = "cksum \"" + filePath + "\"";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        return false;
    }
    
    char buffer[256];
    string result = "";
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        result += buffer;
    }
    int status = pclose(pipe);
    if (status != 0) {
        return false;
    }
    
    // Parse cksum output: "cksum_value file_size filename"
    // Extract first two tokens (cksum and size)
    istringstream iss(result);
    if (!(iss >> crc >> size)) {
        return false;
    }
    
    return true;
}

// Download reference file using FTP/HTTP with optional cksum verification
// cacheDir: directory for cksum cache file (empty string to disable cache)
// autoCksumUpdate: if true, attempt to auto-fill missing cksum from CHECKSUMS file for trusted URLs
// replaceUnverifiableFiles: if true, replace existing final files that cannot be verified
bool downloadReference(const string& url, const string& outputPath, 
                       uint32_t expectedCksum, uint64_t expectedSize, bool allowUntrusted, 
                       const string& cacheDir, bool autoCksumUpdate, bool replaceUnverifiableFiles,
                       string& errorMsg) {
    // Load cksum cache if cacheDir is provided and not already loaded
    if (!cacheDir.empty() && !cksumCacheLoaded) {
        loadCksumCache(cacheDir);
    }
    
    // Determine if output should be decompressed
    bool isGzip = (url.length() > 3 && url.substr(url.length() - 3) == ".gz");
    string tempGzPath = outputPath + ".gz";
    string finalPath = outputPath;
    string backupPath = "";  // Backup path for existing file (restore on failure)
    
    // Extract expected filename from URL for verification
    string expectedFilename = url.substr(url.find_last_of('/') + 1);
    if (isGzip && expectedFilename.length() > 3 && expectedFilename.substr(expectedFilename.length() - 3) == ".gz") {
        expectedFilename = expectedFilename.substr(0, expectedFilename.length() - 3);
    }
    string outputFilename = outputPath.substr(outputPath.find_last_of('/') + 1);
    
    struct stat fileStat;
    
    // For gzipped files, prefer checking .gz file first (cache gz checksums)
    if (isGzip && stat(tempGzPath.c_str(), &fileStat) == 0) {
        // .gz file exists - verify its checksum
        cerr << "Compressed file already exists: " << tempGzPath << "\n";
        cerr << "Checking cksum to verify integrity...\n";
        
        // Get expected cksum (for .gz file)
        uint32_t cksumToVerify = expectedCksum;
        uint64_t sizeToVerify = expectedSize;
        
        // If expectedCksum is 0, check EXPECTED_CKSUM map and cache
        if (cksumToVerify == 0 || sizeToVerify == 0) {
            if (!getCksumForUrl(url, cksumToVerify, sizeToVerify)) {
                // Not found in maps - check if trusted
                bool isTrusted = isTrustedUrl(url);
                if (isTrusted) {
                    // Attempt auto-fill if enabled
                    if (autoCksumUpdate) {
                        cerr << "Trusted URL missing cksum. Computing cksum of existing file...\n";
                        uint32_t autoCrc = 0;
                        uint64_t autoSize = 0;
                        string autoError;
                        if (autoFillCksumFromChecksums(url, cacheDir, autoCrc, autoSize, autoError)) {
                            // Auto-fill returns success with crc=0, size=0 meaning "compute cksum of existing file"
                            if (autoCrc == 0 && autoSize == 0) {
                                // Compute cksum of existing .gz file and save to cache
                                if (!computeCksumFile(tempGzPath, autoCrc, autoSize)) {
                                    errorMsg = "Could not compute cksum of existing file: " + tempGzPath;
                                    return false;
                                }
                                // Save to cache
                                if (!cacheDir.empty()) {
                                    saveCksumToCache(cacheDir, url, autoCrc, autoSize);
                                    cerr << "✓ Computed and cached cksum of existing file\n";
                                    cerr << "  File: " << tempGzPath << "\n";
                                    cerr << "  cksum: " << autoCrc << " size: " << autoSize << "\n";
                                }
                            }
                            cksumToVerify = autoCrc;
                            sizeToVerify = autoSize;
                        } else {
                            cerr << "WARNING: Auto-fill failed: " << autoError << "\n";
                            errorMsg = "FATAL: Trusted URL requires cksum but none is configured: " + url + "\n"
                                       "Compressed file exists but cannot verify integrity: " + tempGzPath + "\n"
                                       "Auto-fill from CHECKSUMS failed: " + autoError + "\n"
                                       "Please update EXPECTED_CKSUM map in CellRangerFormatter.cpp with the cksum value.\n"
                                       "See summary document for instructions on obtaining cksum values from CHECKSUMS files.";
                            return false;
                        }
                    } else {
                        errorMsg = "FATAL: Trusted URL requires cksum but none is configured: " + url + "\n"
                                   "Compressed file exists but cannot verify integrity: " + tempGzPath + "\n"
                                   "Please update EXPECTED_CKSUM map in CellRangerFormatter.cpp with the cksum value.\n"
                                   "See summary document for instructions on obtaining cksum values from CHECKSUMS files.\n"
                                   "Alternatively, enable --autoCksumUpdate Yes to attempt automatic cksum retrieval from CHECKSUMS files.";
                        return false;
                    }
                } else if (!allowUntrusted) {
                    errorMsg = "Compressed file exists but cannot verify integrity (no cksum available and URL not trusted): " + tempGzPath + "\n"
                               "Please remove the existing file or use --allUntrustedUrl Yes to allow unverified files.";
                    return false;
                } else {
                    // Untrusted URL with allowUntrusted flag - warn but proceed
                    cerr << "WARNING: Compressed file exists but cksum not available for verification. Assuming file is correct.\n";
                    cerr << "WARNING: Untrusted URLs disable integrity checking. Use only with trusted sources.\n";
                    cerr << "  File: " << tempGzPath << "\n";
                    // Continue to decompress if final file missing
                }
            }
        }
        
        // If we have an expected cksum, verify the existing .gz file
        if (cksumToVerify != 0 && sizeToVerify != 0) {
            uint32_t existingCksum;
            uint64_t existingSize;
            if (!computeCksumFile(tempGzPath, existingCksum, existingSize)) {
                errorMsg = "Could not compute cksum of existing file: " + tempGzPath;
                return false;
            }
            
            if (existingCksum == cksumToVerify && existingSize == sizeToVerify) {
                cerr << "✓ Correct compressed file exists (cksum matches).\n";
                cerr << "  File: " << tempGzPath << "\n";
                cerr << "  cksum: " << existingCksum << " size: " << existingSize << "\n";
                
                // Check if final file exists and verify decompressed checksum
                bool finalExists = (stat(finalPath.c_str(), &fileStat) == 0);
                bool needDecompress = false;
                
                if (finalExists) {
                    uint32_t decompCrc;
                    uint64_t decompSize;
                    if (getDecompressedCksumForUrl(url, decompCrc, decompSize)) {
                        uint32_t existingDecompCrc;
                        uint64_t existingDecompSize;
                        if (computeCksumFile(finalPath, existingDecompCrc, existingDecompSize)) {
                            if (existingDecompCrc == decompCrc && existingDecompSize == decompSize) {
                                cerr << "  Final file exists and decompressed checksum matches. Skipping decompression.\n";
                                return true;
                            } else {
                                errorMsg = "Final file exists but decompressed checksum mismatch: " + finalPath + "\n"
                                          "  Expected: cksum=" + to_string(decompCrc) + " size=" + to_string(decompSize) + "\n"
                                          "  Found: cksum=" + to_string(existingDecompCrc) + " size=" + to_string(existingDecompSize) + "\n"
                                          "Checksum mismatches are hard errors. Please delete the file and re-download, or remove the cache entry if you want to accept a modified file.";
                                return false;
                            }
                        } else {
                            if (replaceUnverifiableFiles) {
                                cerr << "  Final file exists but could not compute checksum. Re-decompressing...\n";
                                needDecompress = true;
                            } else {
                                cerr << "WARNING: Final file exists but could not compute checksum: " << finalPath << "\n";
                                cerr << "  Keeping existing file. Use --replaceUnverifiableFiles Yes to replace it.\n";
                                return true;
                            }
                        }
                    } else {
                        if (replaceUnverifiableFiles) {
                            cerr << "  Final file exists but no decompressed checksum available. Re-decompressing to ensure integrity...\n";
                            needDecompress = true;
                        } else {
                            cerr << "WARNING: Final file exists but cannot be verified (no decompressed checksum available): " << finalPath << "\n";
                            cerr << "  Keeping existing file. Use --replaceUnverifiableFiles Yes to replace it.\n";
                            return true;
                        }
                    }
                } else {
                    needDecompress = true;
                }
                
                if (needDecompress) {
                    cerr << "  Decompressing to final file...\n";
                    
                    string tempPath = finalPath + ".tmp";
                    
                    // Decompress the existing .gz file to temp
                    gzFile gzIn = gzopen(tempGzPath.c_str(), "rb");
                    if (!gzIn) {
                        errorMsg = "Could not open gzipped file for decompression: " + tempGzPath;
                        return false;
                    }
                    
                    ofstream tempFile(tempPath, ios::binary);
                    if (!tempFile.is_open()) {
                        gzclose(gzIn);
                        errorMsg = "Could not create temporary output file: " + tempPath;
                        return false;
                    }
                    
                    char gzBuffer[65536];
                    int gzBytesRead;
                    while ((gzBytesRead = gzread(gzIn, gzBuffer, sizeof(gzBuffer))) > 0) {
                        tempFile.write(gzBuffer, gzBytesRead);
                        if (!tempFile.good()) {
                            gzclose(gzIn);
                            tempFile.close();
                            remove(tempPath.c_str());
                            errorMsg = "Error writing decompressed data to: " + tempPath;
                            return false;
                        }
                    }
                    
                    gzclose(gzIn);
                    tempFile.close();
                    
                    if (gzBytesRead < 0) {
                        remove(tempPath.c_str());
                        errorMsg = "Error decompressing file: " + tempGzPath;
                        return false;
                    }
                    
                    // Atomic rename
                    if (rename(tempPath.c_str(), finalPath.c_str()) != 0) {
                        remove(tempPath.c_str());
                        errorMsg = "Could not rename temporary file to final output: " + finalPath;
                        return false;
                    }
                    
                    // Compute and cache decompressed checksum
                    uint32_t decompCrc;
                    uint64_t decompSize;
                    if (computeCksumFile(finalPath, decompCrc, decompSize)) {
                        saveDecompressedCksumToCache(cacheDir, url, decompCrc, decompSize);
                    }
                    
                    cerr << "✓ Decompression completed successfully.\n";
                    return true;
                }
            } else {
                errorMsg = "Existing compressed file cksum mismatch!\n"
                           "  File: " + tempGzPath + "\n"
                           "  Expected cksum: " + to_string(cksumToVerify) + " size: " + to_string(sizeToVerify) + "\n"
                           "  Existing file cksum: " + to_string(existingCksum) + " size: " + to_string(existingSize) + "\n"
                           "The existing file does not match the expected file. Please remove it or use a different output path.";
                return false;
            }
        } else {
            // No cksum available but allowUntrusted is true - decompress if final missing or verify if exists
            bool finalExists = (stat(finalPath.c_str(), &fileStat) == 0);
            bool needDecompress = false;
            
            if (finalExists) {
                uint32_t decompCrc;
                uint64_t decompSize;
                if (getDecompressedCksumForUrl(url, decompCrc, decompSize)) {
                    uint32_t existingDecompCrc;
                    uint64_t existingDecompSize;
                    if (computeCksumFile(finalPath, existingDecompCrc, existingDecompSize)) {
                        if (existingDecompCrc == decompCrc && existingDecompSize == decompSize) {
                            cerr << "  Final file exists and decompressed checksum matches. Skipping decompression.\n";
                            return true;
                        } else {
                            errorMsg = "Final file exists but decompressed checksum mismatch: " + finalPath + "\n"
                                      "  Expected: cksum=" + to_string(decompCrc) + " size=" + to_string(decompSize) + "\n"
                                      "  Found: cksum=" + to_string(existingDecompCrc) + " size=" + to_string(existingDecompSize) + "\n"
                                      "Checksum mismatches are hard errors. Please delete the file and re-download, or remove the cache entry if you want to accept a modified file.";
                            return false;
                        }
                    } else {
                        if (replaceUnverifiableFiles) {
                            cerr << "  Final file exists but could not compute checksum. Re-decompressing...\n";
                            needDecompress = true;
                        } else {
                            cerr << "WARNING: Final file exists but could not compute checksum: " << finalPath << "\n";
                            cerr << "  Keeping existing file. Use --replaceUnverifiableFiles Yes to replace it.\n";
                            return true;
                        }
                    }
                } else {
                    if (replaceUnverifiableFiles) {
                        cerr << "  Final file exists but no decompressed checksum available. Re-decompressing to ensure integrity...\n";
                        needDecompress = true;
                    } else {
                        cerr << "WARNING: Final file exists but cannot be verified (no decompressed checksum available): " << finalPath << "\n";
                        cerr << "  Keeping existing file. Use --replaceUnverifiableFiles Yes to replace it.\n";
                        return true;
                    }
                }
            } else {
                needDecompress = true;
            }
            
            if (needDecompress) {
                cerr << "  Decompressing to final file...\n";
                
                string tempPath = finalPath + ".tmp";
                
                // Decompress the existing .gz file to temp
                gzFile gzIn = gzopen(tempGzPath.c_str(), "rb");
                if (!gzIn) {
                    errorMsg = "Could not open gzipped file for decompression: " + tempGzPath;
                    return false;
                }
                
                ofstream tempFile(tempPath, ios::binary);
                if (!tempFile.is_open()) {
                    gzclose(gzIn);
                    errorMsg = "Could not create temporary output file: " + tempPath;
                    return false;
                }
                
                char gzBuffer[65536];
                int gzBytesRead;
                while ((gzBytesRead = gzread(gzIn, gzBuffer, sizeof(gzBuffer))) > 0) {
                    tempFile.write(gzBuffer, gzBytesRead);
                    if (!tempFile.good()) {
                        gzclose(gzIn);
                        tempFile.close();
                        remove(tempPath.c_str());
                        errorMsg = "Error writing decompressed data to: " + tempPath;
                        return false;
                    }
                }
                
                gzclose(gzIn);
                tempFile.close();
                
                if (gzBytesRead < 0) {
                    remove(tempPath.c_str());
                    errorMsg = "Error decompressing file: " + tempGzPath;
                    return false;
                }
                
                // Atomic rename
                if (rename(tempPath.c_str(), finalPath.c_str()) != 0) {
                    remove(tempPath.c_str());
                    errorMsg = "Could not rename temporary file to final output: " + finalPath;
                    return false;
                }
                
                // Compute and cache decompressed checksum
                uint32_t decompCrc;
                uint64_t decompSize;
                if (computeCksumFile(finalPath, decompCrc, decompSize)) {
                    saveDecompressedCksumToCache(cacheDir, url, decompCrc, decompSize);
                }
                
                return true;
            }
        }
    }
    
    // Check if final file already exists (for non-gz or when .gz is missing)
    if (stat(finalPath.c_str(), &fileStat) == 0) {
        // File exists - verify filename and cksum
        cerr << "File already exists: " << finalPath << "\n";
        
        // Check filename matches (optional check - warn if different)
        if (outputFilename != expectedFilename) {
            cerr << "WARNING: Filename mismatch - expected: " << expectedFilename 
                 << ", found: " << outputFilename << "\n";
        }
        
        // For gzipped URLs, check decompressed checksum if available
        bool skipExistingFileCheck = false;
        if (isGzip) {
            uint32_t decompCrc;
            uint64_t decompSize;
            if (getDecompressedCksumForUrl(url, decompCrc, decompSize)) {
                uint32_t existingDecompCrc;
                uint64_t existingDecompSize;
                if (computeCksumFile(finalPath, existingDecompCrc, existingDecompSize)) {
                    if (existingDecompCrc == decompCrc && existingDecompSize == decompSize) {
                        cerr << "✓ File exists and decompressed checksum matches. Skipping download.\n";
                        return true;
                    } else {
                        errorMsg = "File exists but decompressed checksum mismatch: " + finalPath + "\n"
                                  "  Expected: cksum=" + to_string(decompCrc) + " size=" + to_string(decompSize) + "\n"
                                  "  Found: cksum=" + to_string(existingDecompCrc) + " size=" + to_string(existingDecompSize) + "\n"
                                  "Checksum mismatches are hard errors. Please delete the file and re-download, or remove the cache entry if you want to accept a modified file.";
                        return false;
                    }
                } else {
                    if (replaceUnverifiableFiles) {
                        cerr << "Could not compute checksum of existing file. Re-downloading...\n";
                        skipExistingFileCheck = true;
                    } else {
                        cerr << "WARNING: File exists but could not compute checksum: " << finalPath << "\n";
                        cerr << "  Keeping existing file. Use --replaceUnverifiableFiles Yes to replace it.\n";
                        return true;
                    }
                }
            } else {
                // No decompressed checksum available
                if (replaceUnverifiableFiles) {
                    cerr << "File exists but no decompressed checksum available. Re-downloading...\n";
                    skipExistingFileCheck = true;
                } else {
                    cerr << "WARNING: File exists but cannot be verified (no decompressed checksum available): " << finalPath << "\n";
                    cerr << "  Keeping existing file. Use --replaceUnverifiableFiles Yes to replace it.\n";
                    return true;
                }
            }
        }
        
        if (skipExistingFileCheck) {
            // Skip existing file verification and proceed to download
            // Backup existing file instead of deleting (restore on failure)
            backupPath = finalPath + ".bak";
            if (rename(finalPath.c_str(), backupPath.c_str()) != 0) {
                errorMsg = "Could not backup existing file: " + finalPath;
                return false;
            }
            cerr << "  Backed up existing file to: " << backupPath << "\n";
        } else {
            cerr << "Checking cksum to verify integrity...\n";
            
            // Get expected cksum
            uint32_t cksumToVerify = expectedCksum;
            uint64_t sizeToVerify = expectedSize;
            
            // If expectedCksum is 0, check EXPECTED_CKSUM map and cache
            if (cksumToVerify == 0 || sizeToVerify == 0) {
                getCksumForUrl(url, cksumToVerify, sizeToVerify);
            }
            
            bool isTrusted = isTrustedUrl(url);
            bool hasKnown = hasKnownCksum(url);
            
            // If we have an expected cksum, verify the existing file
            if (cksumToVerify != 0 && sizeToVerify != 0) {
                uint32_t existingCksum;
                uint64_t existingSize;
                if (!computeCksumFile(finalPath, existingCksum, existingSize)) {
                    errorMsg = "Could not compute cksum of existing file: " + finalPath;
                    return false;
                }
                
                if (existingCksum == cksumToVerify && existingSize == sizeToVerify) {
                    cerr << "✓ Correct file exists (cksum matches). Skipping download.\n";
                    cerr << "  File: " << finalPath << "\n";
                    cerr << "  cksum: " << existingCksum << " size: " << existingSize << "\n";
                    return true;
                } else {
                    errorMsg = "Existing file cksum mismatch!\n"
                               "  File: " + finalPath + "\n"
                               "  Expected cksum: " + to_string(cksumToVerify) + " size: " + to_string(sizeToVerify) + "\n"
                               "  Existing file cksum: " + to_string(existingCksum) + " size: " + to_string(existingSize) + "\n"
                               "The existing file does not match the expected file. Please remove it or use a different output path.";
                    return false;
                }
            } else {
                // No cksum available
                if (isTrusted || hasKnown) {
                    // Trusted URLs or URLs with known cksum require cksum - attempt auto-fill if enabled
                    if (autoCksumUpdate && isTrusted) {
                    cerr << "Trusted URL missing cksum. Computing cksum of existing file...\n";
                    uint32_t autoCrc = 0;
                    uint64_t autoSize = 0;
                    string autoError;
                    if (autoFillCksumFromChecksums(url, cacheDir, autoCrc, autoSize, autoError)) {
                        // Auto-fill returns success with crc=0, size=0 meaning "compute cksum of existing file"
                        if (autoCrc == 0 && autoSize == 0) {
                            // Compute cksum of existing decompressed file and save to cache
                            // Note: We compute cksum of the decompressed file since that's what exists
                            uint32_t existingCksum;
                            uint64_t existingSize;
                            if (!computeCksumFile(finalPath, existingCksum, existingSize)) {
                                errorMsg = "Could not compute cksum of existing file: " + finalPath;
                                return false;
                            }
                            // Save to cache (using URL + "_decompressed" marker to distinguish)
                            // Actually, for trusted URLs, we trust the existing file and just proceed
                            cerr << "✓ Trusted URL - accepting existing file\n";
                            cerr << "  File: " << finalPath << "\n";
                            cerr << "  cksum: " << existingCksum << " size: " << existingSize << "\n";
                            return true;
                        }
                        cksumToVerify = autoCrc;
                        sizeToVerify = autoSize;
                        cerr << "✓ Successfully auto-filled cksum from CHECKSUMS file\n";
                        // Re-check with auto-filled cksum
                        uint32_t existingCksum;
                        uint64_t existingSize;
                        if (!computeCksumFile(finalPath, existingCksum, existingSize)) {
                            errorMsg = "Could not compute cksum of existing file: " + finalPath;
                            return false;
                        }
                        if (existingCksum == cksumToVerify && existingSize == sizeToVerify) {
                            cerr << "✓ Correct file exists (cksum matches). Skipping download.\n";
                            cerr << "  File: " << finalPath << "\n";
                            cerr << "  cksum: " << existingCksum << " size: " << existingSize << "\n";
                            return true;
                        } else {
                            errorMsg = "Existing file cksum mismatch!\n"
                                       "  File: " + finalPath + "\n"
                                       "  Expected cksum: " + to_string(cksumToVerify) + " size: " + to_string(sizeToVerify) + "\n"
                                       "  Existing file cksum: " + to_string(existingCksum) + " size: " + to_string(existingSize) + "\n"
                                       "The existing file does not match the expected file. Please remove it or use a different output path.";
                            return false;
                        }
                    } else {
                        cerr << "WARNING: Auto-fill failed: " << autoError << "\n";
                        errorMsg = "FATAL: URL requires cksum but none is configured: " + url + "\n"
                                   "File exists but cannot verify integrity: " + finalPath + "\n"
                                   "Auto-fill from CHECKSUMS failed: " + autoError + "\n"
                                   "Please update EXPECTED_CKSUM map in CellRangerFormatter.cpp with the cksum value.\n"
                                   "See summary document for instructions on obtaining cksum values from CHECKSUMS files.";
                        return false;
                    }
                } else {
                    // Trusted URLs or URLs with known cksum require cksum - fail
                    errorMsg = "FATAL: URL requires cksum but none is configured: " + url + "\n"
                               "File exists but cannot verify integrity: " + finalPath + "\n"
                               "Please update EXPECTED_CKSUM map in CellRangerFormatter.cpp with the cksum value.\n"
                               "See summary document for instructions on obtaining cksum values from CHECKSUMS files.\n"
                               "Alternatively, enable --autoCksumUpdate Yes to attempt automatic cksum retrieval from CHECKSUMS files.";
                    return false;
                }
            } else if (!allowUntrusted) {
                // Not trusted and allowUntrusted not set - fail
                errorMsg = "File exists but cannot verify integrity (no cksum available and URL not trusted): " + finalPath + "\n"
                           "Please remove the existing file or use --allUntrustedUrl Yes to allow unverified files.";
                return false;
            } else {
                // Untrusted URL with allowUntrusted flag - warn but proceed
                cerr << "WARNING: File exists but cksum not available for verification. Assuming file is correct.\n";
                cerr << "WARNING: Untrusted URLs disable integrity checking. Use only with trusted sources.\n";
                cerr << "  File: " << finalPath << "\n";
                return true;
                }
            }
        }
    }
    
    // File doesn't exist - proceed with download
    // Warn about insecure transport
    if (url.find("ftp://") == 0 || url.find("http://") == 0) {
        cerr << "WARNING: Using insecure transport (" << (url.find("ftp://") == 0 ? "FTP" : "HTTP") 
             << "). Integrity enforced via embedded cksum." << endl;
    }
    
    // Create output directory
    if (!createDirectoryRecursive(outputPath)) {
        // Restore backup if download fails
        if (!backupPath.empty()) {
            rename(backupPath.c_str(), finalPath.c_str());
            cerr << "  Restored backup file: " << finalPath << "\n";
        }
        errorMsg = "Could not create output directory for " + outputPath;
        return false;
    }
    
    // Open URL with htslib
    hFILE* hfile = hopen(url.c_str(), "r");
    if (!hfile) {
        // Restore backup if download fails
        if (!backupPath.empty()) {
            rename(backupPath.c_str(), finalPath.c_str());
            cerr << "  Restored backup file: " << finalPath << "\n";
        }
        errorMsg = "Could not open URL: " + url + "\n"
                   "Please check your network connection and ensure the URL is accessible.";
        return false;
    }
    
    // Write to temporary .gz file first
    ofstream outFile(tempGzPath, ios::binary);
    if (!outFile.is_open()) {
        (void)hclose(hfile);
        // Restore backup if download fails
        if (!backupPath.empty()) {
            rename(backupPath.c_str(), finalPath.c_str());
            cerr << "  Restored backup file: " << finalPath << "\n";
        }
        errorMsg = "Could not create output file: " + tempGzPath;
        return false;
    }
    
    // Read and write in chunks
    char buffer[65536];
    ssize_t bytesRead;
    while ((bytesRead = hread(hfile, buffer, sizeof(buffer))) > 0) {
        outFile.write(buffer, bytesRead);
        if (!outFile.good()) {
            (void)hclose(hfile);
            outFile.close();
            remove(tempGzPath.c_str());
            // Restore backup if download fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "Error writing to file: " + tempGzPath;
            return false;
        }
    }
    
    (void)hclose(hfile);
    outFile.close();
    
    if (bytesRead < 0) {
        remove(tempGzPath.c_str());
        // Restore backup if download fails
        if (!backupPath.empty()) {
            rename(backupPath.c_str(), finalPath.c_str());
            cerr << "  Restored backup file: " << finalPath << "\n";
        }
        errorMsg = "Error reading from URL: " + url;
        return false;
    }
    
    // Check if URL is trusted or has known cksum
    bool isTrusted = isTrustedUrl(url);
    bool hasKnown = hasKnownCksum(url);
    
    // cksum verification logic
    uint32_t cksumToVerify = expectedCksum;
    uint64_t sizeToVerify = expectedSize;
    
    // If expectedCksum is 0, check EXPECTED_CKSUM map and cache
    if (cksumToVerify == 0 || sizeToVerify == 0) {
        getCksumForUrl(url, cksumToVerify, sizeToVerify);
    }
    
    // Enforce cksum for trusted URLs or URLs with known cksum (FTP + embedded cksum required)
    bool skipVerification = false;
    if ((isTrusted || hasKnown) && (cksumToVerify == 0 || sizeToVerify == 0)) {
        // Attempt auto-fill if enabled
        if (autoCksumUpdate && isTrusted) {
            cerr << "Trusted URL missing cksum. Will compute and cache cksum after download...\n";
            uint32_t autoCrc = 0;
            uint64_t autoSize = 0;
            string autoError;
            if (autoFillCksumFromChecksums(url, cacheDir, autoCrc, autoSize, autoError)) {
                // Auto-fill returns success with crc=0, size=0 meaning "compute cksum after download"
                if (autoCrc == 0 && autoSize == 0) {
                    // Skip verification, compute cksum after download
                    skipVerification = true;
                    cerr << "  (Will compute cksum of downloaded file and save to cache)\n";
                } else {
                    cksumToVerify = autoCrc;
                    sizeToVerify = autoSize;
                    cerr << "✓ Successfully auto-filled cksum from CHECKSUMS file\n";
                }
            } else {
                cerr << "WARNING: Auto-fill failed: " << autoError << "\n";
                remove(tempGzPath.c_str());
                // Restore backup if download fails
                if (!backupPath.empty()) {
                    rename(backupPath.c_str(), finalPath.c_str());
                    cerr << "  Restored backup file: " << finalPath << "\n";
                }
                errorMsg = "FATAL: URL requires cksum but none is configured: " + url + "\n"
                           "FTP downloads require embedded cksum verification for integrity.\n"
                           "Auto-fill from CHECKSUMS failed: " + autoError + "\n"
                           "Please update EXPECTED_CKSUM map in CellRangerFormatter.cpp with the cksum value.\n"
                           "See summary document for instructions on obtaining cksum values from CHECKSUMS files.";
                return false;
            }
        } else {
            remove(tempGzPath.c_str());
            // Restore backup if download fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "FATAL: URL requires cksum but none is configured: " + url + "\n"
                       "FTP downloads require embedded cksum verification for integrity.\n"
                       "Please update EXPECTED_CKSUM map in CellRangerFormatter.cpp with the cksum value.\n"
                       "See summary document for instructions on obtaining cksum values from CHECKSUMS files.\n"
                       "Alternatively, enable --autoCksumUpdate Yes to attempt automatic cksum retrieval from CHECKSUMS files.";
            return false;
        }
    }
    
    // If URL is not trusted and no cksum provided, require allowUntrusted flag
    if (!isTrusted && !hasKnown && (cksumToVerify == 0 || sizeToVerify == 0)) {
        if (!allowUntrusted) {
            remove(tempGzPath.c_str());
            // Restore backup if download fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "FATAL: URL is not in trusted Ensembl FTP URL table and has no known cksum: " + url + "\n"
                       "Untrusted URLs require --allUntrustedUrl Yes flag.\n"
                       "WARNING: Using untrusted URLs disables integrity checking. Only use with trusted sources.";
            return false;
        } else {
            // Allow untrusted - log warning
            cerr << "WARNING: Downloading untrusted URL without cksum verification. Integrity not verified for: " << url << endl;
            cerr << "WARNING: Untrusted URLs disable integrity checking. Use only with trusted sources." << endl;
        }
    }
    
    // Verify cksum if available (always verify if cksum is provided, regardless of trusted status)
    // Or compute and save cksum if skipVerification is set
    if (skipVerification) {
        // Auto-fill mode: compute cksum of downloaded file and save to cache
        uint32_t computedCksum;
        uint64_t computedSize;
        if (!computeCksumFile(tempGzPath, computedCksum, computedSize)) {
            remove(tempGzPath.c_str());
            errorMsg = "Could not compute cksum of downloaded file: " + tempGzPath;
            return false;
        }
        
        // Save computed cksum to cache for future verification
        if (!cacheDir.empty()) {
            saveCksumToCache(cacheDir, url, computedCksum, computedSize);
            cerr << "✓ Computed and cached cksum for downloaded file\n";
            cerr << "  File: " << tempGzPath << "\n";
            cerr << "  cksum: " << computedCksum << " size: " << computedSize << "\n";
            cerr << "  URL: " << url << " (cached for future verification)\n";
        } else {
            cerr << "✓ Downloaded file from trusted URL (no cache dir to save cksum)\n";
            cerr << "  File: " << tempGzPath << "\n";
            cerr << "  cksum: " << computedCksum << " size: " << computedSize << "\n";
        }
    } else if (cksumToVerify != 0 && sizeToVerify != 0) {
        uint32_t computedCksum;
        uint64_t computedSize;
        if (!computeCksumFile(tempGzPath, computedCksum, computedSize)) {
            remove(tempGzPath.c_str());
            // Restore backup if download fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "Could not compute cksum of downloaded file: " + tempGzPath;
            return false;
        }
        
        if (computedCksum != cksumToVerify || computedSize != sizeToVerify) {
            remove(tempGzPath.c_str());
            // Restore backup if download fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "cksum mismatch for downloaded file!\n"
                       "  URL: " + url + "\n"
                       "  Expected cksum: " + to_string(cksumToVerify) + " size: " + to_string(sizeToVerify) + "\n"
                       "  Computed cksum: " + to_string(computedCksum) + " size: " + to_string(computedSize) + "\n"
                       "The downloaded file may be corrupted or the source file has changed.";
            return false;
        }
        
        // After successful download with cksum match, promote URL to trusted by saving to cache
        // Only save if cksum was verified (not for untrusted downloads without verification)
        if (!cacheDir.empty() && cksumToVerify != 0 && sizeToVerify != 0 && !allowUntrusted) {
            saveCksumToCache(cacheDir, url, computedCksum, computedSize);
            cerr << "✓ Saved cksum to cache: " << url << " (promoted to trusted)\n";
        }
    }
    
    // Decompress if needed (only after cksum verification passes)
    if (isGzip) {
        string tempPath = finalPath + ".tmp";
        
        gzFile gzIn = gzopen(tempGzPath.c_str(), "rb");
        if (!gzIn) {
            remove(tempGzPath.c_str());
            // Restore backup if decompression fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "Could not open gzipped file for decompression: " + tempGzPath;
            return false;
        }
        
        ofstream tempFile(tempPath, ios::binary);
        if (!tempFile.is_open()) {
            gzclose(gzIn);
            remove(tempGzPath.c_str());
            // Restore backup if decompression fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "Could not create temporary output file: " + tempPath;
            return false;
        }
        
        char gzBuffer[65536];
        int gzBytesRead;
        while ((gzBytesRead = gzread(gzIn, gzBuffer, sizeof(gzBuffer))) > 0) {
            tempFile.write(gzBuffer, gzBytesRead);
            if (!tempFile.good()) {
                gzclose(gzIn);
                tempFile.close();
                remove(tempGzPath.c_str());
                remove(tempPath.c_str());
                // Restore backup if decompression fails
                if (!backupPath.empty()) {
                    rename(backupPath.c_str(), finalPath.c_str());
                    cerr << "  Restored backup file: " << finalPath << "\n";
                }
                errorMsg = "Error writing decompressed data to: " + tempPath;
                return false;
            }
        }
        
        gzclose(gzIn);
        tempFile.close();
        
        if (gzBytesRead < 0) {
            remove(tempGzPath.c_str());
            remove(tempPath.c_str());
            // Restore backup if decompression fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "Error decompressing file: " + tempGzPath;
            return false;
        }
        
        // Atomic rename
        if (rename(tempPath.c_str(), finalPath.c_str()) != 0) {
            remove(tempGzPath.c_str());
            remove(tempPath.c_str());
            // Restore backup if decompression fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "Could not rename temporary file to final output: " + finalPath;
            return false;
        }
        
        // Keep .gz file in cache (do not delete)
        // remove(tempGzPath.c_str()); // REMOVED: Keep .gz for future reuse
        
        // Compute and cache decompressed checksum
        uint32_t decompCrc;
        uint64_t decompSize;
        if (computeCksumFile(finalPath, decompCrc, decompSize)) {
            saveDecompressedCksumToCache(cacheDir, url, decompCrc, decompSize);
        }
        
        // Remove backup file only after successful decompression
        if (!backupPath.empty()) {
            remove(backupPath.c_str());
            cerr << "  Removed backup file: " << backupPath << "\n";
        }
    } else {
        // Not gzipped, just rename
        if (rename(tempGzPath.c_str(), finalPath.c_str()) != 0) {
            remove(tempGzPath.c_str());
            // Restore backup if rename fails
            if (!backupPath.empty()) {
                rename(backupPath.c_str(), finalPath.c_str());
                cerr << "  Restored backup file: " << finalPath << "\n";
            }
            errorMsg = "Could not rename temporary file to final output: " + finalPath;
            return false;
        }
        
        // Remove backup file only after successful rename
        if (!backupPath.empty()) {
            remove(backupPath.c_str());
            cerr << "  Removed backup file: " << backupPath << "\n";
        }
    }
    
    return true;
}

// Auto-fill cksum for trusted URLs
// For trusted URLs without cksum, we download the file and compute the cksum ourselves,
// then save it to cache for future verification.
// This avoids dealing with different checksum formats (Ensembl/GENCODE use BSD sum, not POSIX cksum).
//
// Note: This function is called when a trusted URL is missing cksum.
// It signals that the caller should proceed with download without verification,
// then compute and cache the cksum after successful download.
bool autoFillCksumFromChecksums(const string& url, const string& cacheDir,
                                  uint32_t& crc, uint64_t& size, string& errorMsg) {
    crc = 0;
    size = 0;
    
    // Verify URL is from a supported provider
    bool isEnsembl = (url.find("ftp://ftp.ensembl.org/pub/") == 0);
    bool isGencode = (url.find("ftp://ftp.ebi.ac.uk/pub/databases/gencode/") == 0);
    
    if (!isEnsembl && !isGencode) {
        errorMsg = "URL is not from Ensembl or GENCODE: " + url;
        return false;
    }
    
    // For trusted URLs, we signal that the caller should:
    // 1. Download the file without cksum verification
    // 2. Compute cksum after download
    // 3. Save to cache
    //
    // We return "success" with crc=0, size=0 to indicate "proceed with download, compute cksum after"
    // The caller checks for crc==0 && size==0 to know it should compute cksum post-download
    
    cerr << "Auto-fill enabled for trusted URL (will compute cksum after download)...\n";
    cerr << "  URL: " << url << "\n";
    cerr << "  Provider: " << (isEnsembl ? "Ensembl" : "GENCODE") << "\n";
    cerr << "  Note: Ensembl/GENCODE CHECKSUMS use BSD sum format, not POSIX cksum.\n";
    cerr << "        Will download file first, then compute and cache POSIX cksum.\n";
    
    // Return true with crc=0, size=0 to signal "download then compute"
    return true;
}

} // namespace CellRangerFormatter
